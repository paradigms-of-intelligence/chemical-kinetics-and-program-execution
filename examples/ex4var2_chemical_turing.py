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

Alphabet is: #(A B C D I O P X S E)

Execute as e.g.: PYTHONPATH=../framework python3 ex4var2_chemical_turing.py
"""


import itertools
import os
import time

import markov_tapes

import numpy

import matplotlib
from matplotlib import pyplot

matplotlib.rcParams.update({'font.size': 18})


# _DEBUG = False
_DEBUG = True
SIZE_A, CL_K = 10, 5
# SIZE_A, CL_K = 10, 4
TAG = 'ex4var2-chemical-turing'
DATA_FILENAME = 'ex4var2.npz'

fn_dy_dt = markov_tapes.get_dy_dt(
  tag=TAG,
  size_a=SIZE_A, cl_k=CL_K, debug=_DEBUG)


def get_p0(cl_k=CL_K,
           tape_fraction=0.25,
           # We need to supply only a limited amount of "powered" stuff
           # due to the competition with fixation.
           cursor_fraction=0.03,
           powered_fraction=0.06,
           random01=False
           ):
  SYM_A, SYM_B, SYM_C, SYM_D, SYM_I, SYM_O, SYM_P, SYM_X, SYM_S, SYM_E = range(10)
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


# DDD Initially, evaluator-molecules are only in solution.
def get_p0e(cl_k=CL_K,
            tape_fraction=0.25,
            # We need to supply only a limited amount of "powered" stuff
            # due to the competition with fixation.
            cursor_fraction=0.04,
            powered_fraction=0.1,
            random01=False
            ):
  SYM_A, SYM_B, SYM_C, SYM_D, SYM_I, SYM_O, SYM_P, SYM_X, SYM_S, SYM_E = range(10)
  p0 = numpy.zeros([SIZE_A]*cl_k, dtype=numpy.float64)
  p0r = p0.ravel()
  g0_solvent = numpy.full([cl_k], fill_value=SYM_S)
  for n, group in enumerate(itertools.product(range(SIZE_A), repeat=cl_k)):
    sg = numpy.array(sorted(group))
    if (sg == SYM_S).all():
      p0r[n] = (1-tape_fraction) * (1 - powered_fraction * cl_k - cursor_fraction * cl_k)
    elif (sg[1:] == SYM_S).all() and sg[0] == SYM_P:
        p0r[n] = (1-tape_fraction) * powered_fraction
    elif (sg[:-1] == SYM_S).all() and sg[-1] == SYM_E:
        p0r[n] = (1-tape_fraction) * cursor_fraction
    elif (sg <= SYM_O).all():  # On-tape.
      if random01:
        if (sg >= SYM_I).all():
          # Just tape.
          p0r[n] = tape_fraction * 0.5**cl_k
      elif (sg == SYM_O).all():  # Only zeros on tape.
          p0r[n] = tape_fraction
  return p0


p0 = get_p0e()


if 'do_check_p0' and False:
  # Disabled by default. We can enable this to verify that the
  # initial-state is compatible with describing a Markov process. For
  # our parameters, the check succeeds, but takes some time (and a
  # nontrivial amount of RAM).
  delta, eigenspace = markov_tapes.get_ctm_eigenvalue1_eigenspace(p0)
  if delta > 1e-6:
    raise ValueError('Impossible p0.')


ode_ts = numpy.linspace(0, 10000.0, 5001)

if not os.access(DATA_FILENAME, os.R_OK):
  # Making the costly computation idempotent.
  t0 = time.monotonic()
  ode_ys = markov_tapes.ode_integrate_ivp(
    tag=TAG,
    size_a=SIZE_A, cl_k=CL_K,
    p0=p0,
    ts=ode_ts,
    debug=_DEBUG,
    # f1:
    # mv ex4var2.npz ex4var2f1.npz
    # ivp_kwargs=dict(rtol=1e-11, atol=1e-11, method='DOP853')
    # f2:
    # ivp_kwargs=dict(rtol=1.5e-11, atol=1.5e-11, method='DOP853')
    # f3:
    # ivp_kwargs=dict(rtol=2e-11, atol=2e-11, method='DOP853')
    ivp_kwargs=dict(rtol=1e-11, atol=1e-11, method='DOP853')
  )
  numpy.savez_compressed(DATA_FILENAME, ode_ys=ode_ys)

if True:
  ode_ys = numpy.load(DATA_FILENAME)['ode_ys']
  fig = pyplot.figure(figsize=(16, 8))
  ax = fig.gca()
  ax.grid()
  def get_seq_prob(ode_ys, seq):
    return [markov_tapes.seq_prob(spd.reshape(*[SIZE_A]*CL_K), seq)[0]
            for spd in ode_ys]
  def plot_seq_prob(ode_ys, ts, seq, style, label, offset=0, **extra):
    ys = get_seq_prob(ode_ys, seq)
    ax.plot(numpy.log(ts[1:]) / numpy.log(10),
            offset + (1e-100 + numpy.log(ys[1:])) / numpy.log(10),
            style, label=label, **extra)
    print(f'{label}: p_final={ys[-1]}')
  plot_seq_prob(ode_ys, ode_ts, [0], '-k', label='p(A)')
  plot_seq_prob(ode_ys, ode_ts, [1], '-b', label='p(B)')
  plot_seq_prob(ode_ys, ode_ts, [2], '-m', label='p(C)')    
  plot_seq_prob(ode_ys, ode_ts, [3], '-r', label='p(D)')
  plot_seq_prob(ode_ys, ode_ts, [4, 5, 4, 3], '--r', label='p(IOID)', linewidth=3)
  plot_seq_prob(ode_ys, ode_ts, [5, 4, 5, 4, 5], '--k', label='p(OIOIO)', linewidth=3)
  # plot_seq_prob(ode_ys, ode_ts, [4, 5, 4], '--b', label='p(IOI)', linewidth=3)
  # plot_seq_prob(ode_ys, ode_ts, [4, 4, 4], '-g', label='p(III)', linewidth=3)
  plot_seq_prob(ode_ys, ode_ts, [4, 4, 4, 4], '--g', label='p(IIII)', linewidth=3)    
  plot_seq_prob(ode_ys, ode_ts, [6], '-.k', label='p(P)')
  plot_seq_prob(ode_ys, ode_ts, [7], '-.b', label='p(X)')
  plot_seq_prob(ode_ys, ode_ts, [9], '-.g', label='p(E)')
  # Total cursor
  c_tot = sum(numpy.asarray(get_seq_prob(ode_ys, seq))
              for seq in [[0], [1], [2], [3]])
  ax.plot(numpy.log(ode_ts[1:]) / numpy.log(10),
          numpy.log(1e-100 + c_tot[1:]) / numpy.log(10),
          '-.m', label='{total cursor}')
  ax.set_ylabel(r'$\log_{10}(p)$')
  ax.set_xlabel(r'$\log_{10}$(time)')
  ax.legend(loc='best')
  fig.savefig('ex4var2_chemical_turing.pdf')
  fig.show()
