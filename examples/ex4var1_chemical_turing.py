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

Execute as e.g.: PYTHONPATH=../framework python3 ex4var1_chemical_turing.py
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
TAG = 'ex4var1-chemical-turing'

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


p0 = get_p0(cursor_fraction=0.001, powered_fraction=0.05, random01=True)

if 'do_check_p0' and False:
  # Disabled by default. We can enable this to verify that the
  # initial-state is compatible with describing a Markov process. For
  # our parameters, the check succeeds, but takes some time (and a
  # nontrivial amounts of RAM).
  delta, eigenspace = markov_tapes.get_ctm_eigenvalue1_eigenspace(p0_a)
  if delta > 1e-10:
    raise ValueError('Impossible p0.')


ode_ts = numpy.linspace(0, 2000.0, 2001)


if True:
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
  plot_seq_prob(ode_ys, ode_ts, [0], '-k', label='p(A)')
  plot_seq_prob(ode_ys, ode_ts, [1], '-b', label='p(B)')
  plot_seq_prob(ode_ys, ode_ts, [2], '-m', label='p(C)')    
  plot_seq_prob(ode_ys, ode_ts, [3], '-r', label='p(D)')
  plot_seq_prob(ode_ys, ode_ts, [4, 5, 4, 3], '--r', label='p(IOID)', linewidth=3)
  plot_seq_prob(ode_ys, ode_ts, [6], ':k',
                label='p(P)')
  plot_seq_prob(ode_ys, ode_ts, [7], ':b',
                label='p(X)')
  # Total cursor
  c_tot = sum(numpy.asarray(get_seq_prob(ode_ys, seq))
              for seq in [[0], [1], [2], [3]])
  ax.plot(ode_ts[1:], numpy.log(1e-100 + c_tot[1:]) / numpy.log(10),
          ':m', label='{total cursor}')
  ax.set_ylabel(r'$\log_{10}(p)$')
  ax.set_xlabel('time')
  ax.legend(loc='best')
  fig.savefig('ex4var1_chemical_turing.pdf')
  fig.show()
  # Entropy
  entropies = [markov_tapes.markov_entropy(ys.reshape([SIZE_A]*CL_K))
               for ys in ode_ys]
  print('S_initial: {entropies[0]}, S_final: '
        f'{entropies[-1]}, '
        f'S_delta: {entropies[-1]-entropies[0]}')
  fig_S = pyplot.figure(figsize=(16, 8))
  ax_S = fig_S.gca()
  ax_S.plot(ode_ts, entropies, '-k')
  ax_S.set_ylabel('Markov entropy')
  ax_S.set_xlabel('time')
  ax_S.grid()
  fig_S.savefig('ex4var1_chemical_turing_s.pdf')
  fig_S.show()

  
