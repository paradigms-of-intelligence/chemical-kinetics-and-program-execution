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

"""Ferromagnet example.

Execute as e.g.: PYTHONPATH=../framework python3 ex2_ferromagnet_tape.py


This script produces the graphs for validating the Markov process
dynamics calculation against the analytic approximation, and exploring
the impact of Markov process context length size.
"""

import itertools
import markov_tapes

import ex2_ferromagnet_analytic

import numpy

from matplotlib import pyplot
import matplotlib

matplotlib.rcParams.update({'font.size': 18})


_DEBUG = False

T_MAX=60


def get_p0(cl_k, p_pair=0.01):
  """Produces the initial probability distribution."""
  p0 = numpy.zeros([2]*cl_k, dtype=numpy.float64)
  p0_flat = p0.ravel()  # New view on same data.
  p0_flat[0] = 1.0 - p_pair * (cl_k + 1)  
  for k in range(cl_k - 1):
    p0_flat[0b11 << k] = p_pair
  p0_flat[1] = p_pair
  p0_flat[1 << (cl_k-1)] = p_pair
  return p0


def get_p0_v1(cl_k, p_pair=0.01):
  """Produces the initial probability distribution (fixed version)."""
  p0 = numpy.zeros([2]*cl_k, dtype=numpy.float64)
  p0_flat = p0.ravel()  # New view on same data.
  for k in range(cl_k - 1):
    p0_flat[0b11 << k] = p_pair
  p0_flat[1] = p_pair
  p0_flat[1 << (cl_k-1)] = p_pair
  p0_flat[(1 << (cl_k-1)) | 1] = p_pair**2  # This is the essential correction!
  p0_flat[0] = 1.0 - p0_flat.sum()
  return p0


# ode_ts = numpy.linspace(0, 40, 1001)
ode_ts = numpy.linspace(0, T_MAX, 1001)

p_history_by_cl_k = {}


for cl_k in range(3, 8):
  print(f'Doing {cl_k=}...')
  p_at_t0 = get_p0(cl_k, p_pair=1/250)
  ode_ys = markov_tapes.ode_integrate(
    tag='ex2-ferromagnetic-chain',
    size_a=2, cl_k=cl_k,
    p0=p_at_t0,
    ts=ode_ts,
    debug=_DEBUG,
    odeint_kwargs=dict(rtol=1e-9, atol=1e-9))
  p_history_by_cl_k[cl_k] = ode_ys.reshape(ode_ts.shape + (2,) * cl_k)


if 'plot-cl_k-comparison':
  # Plotting a comparison for different cl_k.
  fig = pyplot.figure(figsize=(16, 12))
  ax = fig.gca()
  for length, color in ((1, 'k'), (2, 'r'), (3, 'b'), (4, 'g'), (5, 'm')):
    for cl_k, style in ((7, '-'), (5, '--'), (4, '-.'), (3, ':')):
      probs = markov_tapes.seq_prob(p_history_by_cl_k[cl_k],
                                    (0, *((1,)*length), 0),
                                    num_prefix_indices=1)[0][1:]
      to_plot = probs
      to_plot = numpy.log(numpy.clip(probs, 1e-30, numpy.inf)) / numpy.log(10)
      ax.plot(ode_ts[1:], to_plot,
              style + color,
              label=f'L={length}, cl_k={cl_k}',
              linewidth=(3 if style == '--' else 1)
              )
  ax.legend(loc='best')
  ax.set_title('Impact of Context Length')
  ax.set_xlabel('Time')
  ax.set_ylabel('log10(p)')
  ax.grid()
  fig.savefig('ferromagnet_mpd_cl_k_comparison.pdf')
  fig.show()


if 'plot-analytic_comparison':
  cl_k = 7
  p_history_analytic = ex2_ferromagnet_analytic.get_p_history(t_max=T_MAX)
  #
  fig = pyplot.figure(figsize=(16, 12))
  ax = fig.gca()
  ax.grid()
  for length, color in ((1, 'k'), (2, 'r'), (3, 'b'), (4, 'g'), (5, 'm')):
    scaling, scaling_text = (0.25, '*0.25') if length == 2 else (1, '')
    probs = markov_tapes.seq_prob(p_history_by_cl_k[cl_k],
                                  (0, *((1,)*length), 0),
                                  num_prefix_indices=1)[0][1:]
    ax.plot(ode_ts[1:], scaling * probs, f'-{color}', label=f'p(L={length}){scaling_text}, MPD')
    ax.plot(ode_ts[1:], scaling * p_history_analytic[1:, length-1],
            f'--{color}',
            label=f'p(L={length}){scaling_text}, AA',
            linewidth=3)
  ax.set_ylabel(r'p')
  ax.set_xlabel('Time')
  ax.legend(loc='best')
  ax.set_title('Comparison of Analytic Approximation (AA) and Markov Process '
               f'Dynamics (MPD)')  
  fig.savefig('ferromagnet_mpd_aa_comparison.pdf')
  fig.show()
