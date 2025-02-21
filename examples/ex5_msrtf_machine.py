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

"""Simple Machine Language Example

Alphabet is: #(M S R T F)

Execute as e.g.: PYTHONPATH=../framework python3 ex5_msrtf_machine.py
"""


import itertools
import os
import markov_tapes

import numpy

import matplotlib
from matplotlib import pyplot

matplotlib.rcParams.update({'font.size': 18})


_DEBUG = True
SIZE_A, CL_K = 5, 5
TAG = 'ex5-msrtf-machine'
DATA_FILENAME = 'ex5_msrtf.npz'

fn_dy_dt = markov_tapes.get_dy_dt(
  tag=TAG,
  size_a=SIZE_A, cl_k=CL_K, debug=_DEBUG)


def get_p0(cl_k=CL_K):
  p0 = numpy.zeros([SIZE_A]*cl_k, dtype=numpy.float64)
  # First three symbols have equal probability, others have zero probability.
  p0[(slice(0, 3),)*CL_K] = 3.0**(-cl_k)
  return p0


p0 = get_p0()


delta, eigenspace = markov_tapes.get_ctm_eigenvalue1_eigenspace(p0)
if delta > 1e-10:
    raise ValueError('Impossible p0.')


ode_ts = numpy.linspace(0, 500.0, 4001)


if not os.access(DATA_FILENAME, os.R_OK):
  ode_ys = markov_tapes.ode_integrate_ivp(
    tag=TAG,
    size_a=SIZE_A, cl_k=CL_K,
    p0=p0,
    ts=ode_ts,
    debug=_DEBUG,
    ivp_kwargs=dict(rtol=1e-13, atol=1e-13, method='DOP853')
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
  def plot_seq_prob(ode_ys, ts, seq, style, label, **extra):
    ys = get_seq_prob(ode_ys, seq)
    ax.plot(ts[1:], ys[1:], style, label=label, **extra)
    print(f'{label}: p_final={ys[-1]}')
  plot_seq_prob(ode_ys, ode_ts, [1, 4, 3, 0], '-k', label='p(SFTM)')
  plot_seq_prob(ode_ys, ode_ts, [1, 3, 0, 1], '-b', label='p(STMS)')
  plot_seq_prob(ode_ys, ode_ts, [0, 0, 0, 0], '-r', label='p(MMMM)')    
  plot_seq_prob(ode_ys, ode_ts, [2, 2, 2, 2], '--r', label='p(RRRR)')
  plot_seq_prob(ode_ys, ode_ts, [0, 2, 0, 0], '-m', label='p(MRMM)')
  plot_seq_prob(ode_ys, ode_ts, [0, 1, 2, 3], '--m', label='p(MSRT)')
  plot_seq_prob(ode_ys/50, ode_ts, [0], ':k', label='p(M)/50')
  plot_seq_prob(ode_ys/50, ode_ts, [1], ':b', label='p(S)/50')
  plot_seq_prob(ode_ys/50, ode_ts, [2], ':c', label='p(R)/50')
  plot_seq_prob(ode_ys/50, ode_ts, [3], ':r', label='p(T)/50')
  plot_seq_prob(ode_ys/50, ode_ts, [4], ':m', label='p(F)/50')
  ax.set_ylabel('probability')
  ax.set_xlabel('time')
  ax.legend(loc='best')
  fig.savefig('ex5_msrtf.pdf')
  fig.show()
