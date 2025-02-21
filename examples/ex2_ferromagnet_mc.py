# Copyright 2025 Google LLC
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#     https://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Monte Carlo based simulation of a classical ferromagnetic spin chain."""

import collections
import os

import numpy
from matplotlib import pyplot
import scipy.integrate

import matplotlib

import ex2_ferromagnet_analytic

# matplotlib.rcParams.update({'font.size': 14})
matplotlib.rcParams.update({'font.size': 18})

# Parameters

NUM_TRIALS = 100
CHAIN_LENGTH = 50000
NUM_TIME_STEPS = 4000
SITES_PER_PAIR = 250
NUM_TRIALS_PER_TIME_STEP = CHAIN_LENGTH // 100
beta = 1.0  
J = 1.0
h = -0.25  # Pushing magnetizations towards "down" state.
t_max = 40
t_steps = 4000

DATA_FILE = 'ferromagnet_mc_chain_counts.npz'

def simulate(current_chain,
             num_time_steps,
             num_trials_per_time_step=1000,
             # Energy conventions in alignment with the common Ising model
             # definitions; the sigma_i are {-1, +1}.
             # E_total = -J sigma_i sigma_{i+1} - h sigma_i
             J=1,
             h=0,
             beta=1,
             rng=None):
  """Simulates temporal evolution of a classical 'spin-chain'.

  At every time step, `num_trials_per_time_step` sites are picked at
  random (with repetitions) and flipped if a randomly rolled number is
  below the corresponding energy-dependent rejection-sampling
  threshold.

  Args:
    current_chain: 1-axis numpy ArrayLike with 0=down, 1=up entries
      representing the initial chain.
    num_time_steps: Number of time steps to perform.
    num_trials_per_time_step: How many sites to look at for one time-step.
    J: nearest-neighbor coupling strength.
    h: external field strength.
    beta: Inverse temperature 1/(k_B T).
    rng: Optional numpy.random.RandomState instance to use
      for generating random numbers reproducibly. If `None`,
      a new one is allocated without explicit seeding.
  """
  if rng is None:
    rng = numpy.random.RandomState()
  beta_J = beta * J
  beta_h = beta * h
  current_chain = numpy.asarray(current_chain, dtype=numpy.int8)
  chain_length = current_chain.size
  offsets = numpy.array((-1, 0, 1))[numpy.newaxis, :]
  # We have the following relevant transitions:
  # udu => uuu,  delta_E = -4J - 2h
  # dud => ddd,  delta_E = -4J + 2h
  # uud => udd,  delta_E = +2h
  # ddu => duu,  delta_E = -2h
  # (and the symmetric forms)
  #
  # A flip at time t will adjust the state at that time and into all future.
  # So, the array has to start out with copies of the initial state for
  # all times, since this is the evolution we would get if zero flips
  # were successful.
  result = numpy.pad(current_chain[numpy.newaxis, :],
                     [(0, num_time_steps-1), (0, 0)])
  for nt in range(1, num_time_steps):
    result[nt, :] = result[nt - 1, :]  # Copy state.
    indices = rng.randint(0, chain_length, size=num_trials_per_time_step)
    random01s = rng.uniform(0, 1, size=num_trials_per_time_step)
    strides = (indices[:, numpy.newaxis] + offsets) % chain_length
    # We need to process strides one-by-one, since first stride-adjustment
    # may impact every subsequent stride-adjustment.
    for random01, stride3 in zip(random01s, strides):
      ijk = result[nt-1, stride3]
      # We consider flipping the middle element:
      # uuu <-> udu must come with an energy-change of -2J <-> 2J
      # in (-1, 1)-encoding. Total energy change hence is -4, 0, or +4.
      e_neigbors_after_minus_before = 2 * (
        (int(ijk[0] == ijk[1]) + int(ijk[1] == ijk[2])) -
        (int(ijk[0] != ijk[1]) + int(ijk[1] != ijk[2])))
      # Biggest energy change is ddd <-> dud: 
      e_J_relative_rate_factor = numpy.exp(
        -beta_J * (e_neigbors_after_minus_before + 4))
      # Case h>0: If the flipping cell was "up", flipping to the energetically
      # disfavored "down" state gets punished with rate-reduction.
      e_h_relative_rate_factor = (
        numpy.exp(-2 * beta_h * ijk[1]) if h > 0 else
        numpy.exp(+2 * beta_h * (1 - ijk[1])))
      relative_rate_factor = e_J_relative_rate_factor * e_h_relative_rate_factor
      assert relative_rate_factor < 1.0001
      if random01 < relative_rate_factor:
          result[nt, stride3[1]] ^= 1
  return result


def energy(chains, J, h):
  """Computes the total energy of all chains."""
  chains_pm1 = chains.astype(numpy.float32) * 2 - 1
  E_J = -J * ((chains_pm1[..., 1:] * chains_pm1[..., :-1]).sum(axis=-1) +
              (chains_pm1[..., 0] * chains_pm1[..., -1]))
  E_h = - h * chains_pm1.sum(axis=-1)
  return E_J + E_h


def island_length_stats(chain, is_up=True):
  """Counts island lengths (with wraparound).

  Args:
    chain: 1-axis numpy ArrayLike with 0=down, 1=up entries
      representing the full chain.

  Returns:
    A new `{island_length: count}` dictionary with island counts.
  """
  chain = numpy.asarray(chain).astype(numpy.int8)
  eff_chain = chain if is_up else 1 - chain
  wrapping_chain_length_prefix = eff_chain.argmin()
  wrapping_chain_length_suffix = eff_chain[::-1].argmin()
  wrapping_chain_length = (wrapping_chain_length_prefix +
                           wrapping_chain_length_suffix)
  eff_chain_reduced = eff_chain[wrapping_chain_length_prefix:
                                chain.size - wrapping_chain_length_suffix]
  stats = {wrapping_chain_length: int(wrapping_chain_length > 0)}
  if eff_chain_reduced.size == 0:
    return stats
  swap_positions = eff_chain_reduced[:-1] ^ eff_chain_reduced[1:]
  swap_indices = numpy.arange(eff_chain_reduced.size-1)[swap_positions == 1]
  # Invariant: We started at 0 and ended at 0,
  # so swapped an even number of times.
  assert len(swap_indices) % 2 == 0
  for low, high in swap_indices.reshape(-1, 2):
    chain_length = high - low
    stats[chain_length] = 1 + stats.get(chain_length, 0)
  return stats

###

if not os.access(DATA_FILE, os.R_OK):
  # 
  seed_offset = 1000
  # We are only interested in chain lengths 1-5.
  # Counts indexing: [num_trial, time_step, chain_length]    
  chain_counts = numpy.zeros([NUM_TRIALS, NUM_TIME_STEPS, 6])
  for n_trial in range(NUM_TRIALS):
    print('Doing trial:', n_trial)
    rng = numpy.random.RandomState(seed=n_trial + seed_offset)
    # Every site has to have the given probability to be
    # the left element of a pair.
    pair_positions = rng.uniform(0, 1, size=CHAIN_LENGTH) < 1/SITES_PER_PAIR
    chain0 = (pair_positions | numpy.roll(pair_positions, 1)).astype(numpy.int8)
    history = simulate(
      chain0,
      NUM_TIME_STEPS,
      num_trials_per_time_step=NUM_TRIALS_PER_TIME_STEP,
      J=J, h=h, beta=beta,
      rng=rng)
    for n_time, chain in enumerate(history):
      stats = island_length_stats(chain)
      for c_len in range(1, 6):
        chain_counts[n_trial, n_time, c_len] = stats.get(c_len, 0)
  numpy.savez_compressed(DATA_FILE,
                         chain_counts=chain_counts)

  
if os.access(DATA_FILE, os.R_OK):
  chain_counts = numpy.load(DATA_FILE)['chain_counts']
  p10 = numpy.percentile(chain_counts, 10, axis=0) / CHAIN_LENGTH
  p50 = numpy.percentile(chain_counts, 50, axis=0) / CHAIN_LENGTH
  p90 = numpy.percentile(chain_counts, 90, axis=0) / CHAIN_LENGTH
  p_history_analytic = ex2_ferromagnet_analytic.get_p_history(
    beta=beta, J=J, h=h, t_max=t_max, t_steps=t_steps,
    p0_pair_start=1/SITES_PER_PAIR)
  fig = pyplot.figure(figsize=(16, 12))
  ax = fig.gca()
  ts = numpy.linspace(0, t_max, t_steps)
  ax.grid()
  for n, style in enumerate(('-k', '-r', '-b', '-g'), 1):
    scaling, scaling_text = (0.25, '*0.25') if n == 2 else (1, '')
    ax.plot(ts, p10[:, n] * scaling, style, label=f'p(L={n}){scaling_text}, MC')
    ax.plot(ts, p50[:, n] * scaling, style)
    ax.plot(ts, p90[:, n] * scaling, style)
    ax.plot(ts, scaling * p_history_analytic[:, n-1],
            style, label=f'p(L={n}){scaling_text}, AA', linewidth=3)
  ax.set_title('p(L), Monte-Carlo 10th/50th/90th percentile '
               'vs. analytic approximation')
  ax.set_xlabel('Time')
  ax.set_ylabel('p(L)')
  ax.legend(loc='best')
  fig.savefig('ferromagnet_mc_plot.pdf')
  fig.show()
