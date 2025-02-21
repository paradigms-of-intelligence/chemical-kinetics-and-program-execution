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

"""Chemical Turing Machine example

Alphabet is: #(A B C D I O P X S)

Execute as e.g.: PYTHONPATH=../framework python3 ex4_chemical_turing.py
"""


import itertools
import markov_tapes

import numpy

import matplotlib
from matplotlib import pyplot

matplotlib.rcParams.update({'font.size': 18})


# _DEBUG = False
_DEBUG = True
SIZE_A, CL_K = 9, 5
TAG = 'ex4-chemical-turing'

fn_dy_dt = markov_tapes.get_dy_dt(
  tag=TAG,
  size_a=SIZE_A, cl_k=CL_K, debug=_DEBUG)


def get_p0(cl_k=CL_K,
           tape_fraction=0.25,
           # We need to supply only a limited amount of "powered" stuff
           # due to the competition with fixation.
           cursor_fraction=0.01,
           # With powered_fraction == cursor_fraction, we would have
           # just enough oomph to finish every calculation, except for
           # entropy, due to tape_fraction=1/4.  We need some excess,
           # but not too much, to not stabilize the A-state too much.
           powered_fraction=0.05,
           random01=False
           ):
  SYM_A, SYM_B, SYM_C, SYM_D, SYM_I, SYM_O, SYM_P, SYM_X, SYM_S = range(9)
  p0 = numpy.zeros([SIZE_A]*cl_k, dtype=numpy.float64)
  p0r = p0.ravel()
  g0_solvent = numpy.full([cl_k], fill_value=SYM_S)
  for n, group in enumerate(itertools.product(range(SIZE_A), repeat=cl_k)):
    sg = numpy.array(sorted(group))
    if (sg[1:] == g0_solvent[1:]).all():
      if sg[0] == SYM_P:  # Energized molecule in solvent.
        p0r[n] = (1-tape_fraction) * powered_fraction
      elif sg[0] == SYM_S:
        # Just solvent.
        p0r[n] = (1-tape_fraction) * (1 - powered_fraction * cl_k)
    elif (sg <= SYM_O).all():  # On-tape.
      if random01:
        if sg[0] == SYM_A and (sg[1:] >= SYM_I).all():
          # Cursor on tape
          p0r[n] = tape_fraction * cursor_fraction * 0.5**(cl_k-1)
        elif (sg >= SYM_I).all():
          # Just tape.
          p0r[n] = tape_fraction * (1 - cursor_fraction * cl_k) * 0.5**cl_k
      else:  # Only zeros and A-cursors on tape.
        if sg[0] == SYM_A and (sg[1:] == SYM_O).all():
          # Cursor on tape
          p0r[n] = tape_fraction * cursor_fraction
        elif (sg == SYM_O).all():
          # Just tape.
          p0r[n] = tape_fraction * (1 - cursor_fraction * cl_k)
  return p0


p0_a = get_p0(powered_fraction=0.04)

if 'do_check_p0' and False:
  # Disabled by default. We can enable this to verify that the
  # initial-state is compatible with describing a Markov process. For
  # our parameters, the check succeeds, but takes some time (and a
  # nontrivial amount of RAM).
  delta, eigenspace = markov_tapes.get_ctm_eigenvalue1_eigenspace(p0_a)
  if delta > 1e-10:
    raise ValueError('Impossible p0.')


p0_b = get_p0(powered_fraction=0.01)  # "Starved"


ode_ts = numpy.linspace(0, 2000.0, 2001)
# We also should verify that nothing changes in the very long run.
# (This is indeed the case).
# ode_ts = numpy.linspace(0, 50000.0, 2001)

if True:
  for (p0, filename) in ((p0_a, 'ex4_chemical_turing_a.pdf'),
                         (p0_b, 'ex4_chemical_turing_b.pdf')):
    ode_ys = markov_tapes.ode_integrate_ivp(
      tag=TAG,
      size_a=SIZE_A, cl_k=CL_K,
      p0=p0,
      ts=ode_ts,
      debug=_DEBUG,
      ivp_kwargs=dict(rtol=1e-13, atol=1e-13, method='DOP853')
    )
    ####
    fig = pyplot.figure(figsize=(16, 8))
    ax = fig.gca()
    ax.grid()
    def get_seq_prob(ode_ys, seq):
      return [markov_tapes.seq_prob(spd.reshape(*[SIZE_A]*CL_K), seq)[0]
              for spd in ode_ys]
    def plot_seq_prob(ode_ys, ts, seq, style, label, **extra):
      ys = get_seq_prob(ode_ys, seq)
      ax.plot(ts[1:], (1e-100 + numpy.log(ys[1:])) / numpy.log(10),
              style, label=label, **extra)
      print(f'{label}: p_final={ys[-1]}')
    plot_seq_prob(ode_ys, ode_ts, [5, 0, 5, 5, 5], '-k', label='p(OAOOO)')
    plot_seq_prob(ode_ys, ode_ts, [5, 4, 1, 5, 5], '-b', label='p(OIBOO)')
    plot_seq_prob(ode_ys, ode_ts, [5, 4, 1, 4, 5], '--b', label='p(OIBIO)')    
    plot_seq_prob(ode_ys, ode_ts, [5, 4, 5, 2, 5], '-m', label='p(OIOCO)')
    plot_seq_prob(ode_ys, ode_ts, [5, 4, 5, 2, 4], '--m', label='p(OIOCI)')    
    plot_seq_prob(ode_ys, ode_ts, [5, 4, 5, 4, 3], '-r', label='p(OIOID)')
    plot_seq_prob(ode_ys, ode_ts, [6], ':k',
                  label='p(P)')
    plot_seq_prob(ode_ys, ode_ts, [7], ':b',
                  label='p(X)')
    # Total cursor
    c_tot = sum(numpy.asarray(get_seq_prob(ode_ys, seq))
                for seq in [[0], [1], [2], [3]])
    ax.plot(ode_ts[1:], numpy.log(1e-100 + c_tot[1:]) / numpy.log(10),
            '--r', label='{total cursor}')
    ax.set_ylabel(r'$\log_{10}(p)$')
    ax.set_xlabel('time')
    ax.legend(loc='best')
    fig.savefig(filename)
    fig.show()

## The p0_a scenario gives us at t=t_final:
#
# p(OAOOO): p_final=1.069972289390935e-08
# p(OIBOO): p_final=6.515573824924313e-07
# p(OIBIO): p_final=6.515311604360241e-07
# p(OIOCO): p_final=3.968674272397802e-05
# p(OIOCI): p_final=3.968643987041947e-05
# p(OIOID): p_final=0.00241751541540069
# p(P): p_final=0.02258485544510012
# p(X): p_final=0.007415144554899872
#
## For the p0_b scenario, we find at t=t_final:
#
# p(OAOOO): p_final=0.00012550563638350954
# p(OIBOO): p_final=0.00031502540335240174
# p(OIBIO): p_final=5.084130198577003e-05
# p(OIOCO): p_final=0.0005186964734668385
# p(OIOCI): p_final=9.96749791258151e-05
# p(OIOID): p_final=0.0013280547249873754
# p(P): p_final=0.0019018941966848447
# p(X): p_final=0.005598105803315155
