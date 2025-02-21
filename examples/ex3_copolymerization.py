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

"""Nylon Copolymerization example.

Execute as e.g.: PYTHONPATH=../framework python3 ex3_copolymerization.py
"""

import itertools
import os

import numpy
import matplotlib
from matplotlib import pyplot

import markov_tapes


matplotlib.rcParams.update({'font.size': 18})


_DEBUG = False

CL_K = 6
DATA_FILE = 'ex3var2.npz'

def get_p0(cl_k=CL_K, p_a=0.02):
  p0 = numpy.zeros([4]*cl_k, dtype=numpy.float64)
  sym_o, sym_a, sym_m, sym_n = range(4)
  set_o = {0}
  # Inefficiency of this piece of code is irrelevant here.
  # (We could do this better with some more advanced Python idioms.)
  for xs in itertools.product(range(4), repeat=cl_k):
    # Need to have three O-s in the four-symbol sequence.
    if set(sorted(xs)[:cl_k-1]) != set_o: continue
    if sym_a in xs:
        p0[xs] = p_a
    elif sym_m in xs or sym_n in xs:
        p0[xs] = 0.5 * p_a
    else:  # all-O sequence
        p0[xs] = 1 - cl_k * p_a * 2
  return p0


p0 = get_p0()


delta, eigenspace = markov_tapes.get_ctm_eigenvalue1_eigenspace(p0)
if delta > 1e-10:
    raise ValueError('Impossible p0.')


odeint_kwargs = dict(rtol=1e-9, atol=1e-9)


ode_ts = numpy.linspace(0, 1000, 1001)
ode_ts2 = numpy.linspace(0, 200, 1001)

ode_ys = markov_tapes.ode_integrate(
    tag='ex3-copolymerization',
    size_a=4, cl_k=CL_K,
    p0=p0,
    ts=ode_ts,
    odeint_kwargs=odeint_kwargs,
    debug=_DEBUG)

ode_ys_var1 = markov_tapes.ode_integrate(
    tag='ex3var1-copolymerization',
    size_a=4, cl_k=CL_K,
    p0=p0,
    ts=ode_ts,
    odeint_kwargs=odeint_kwargs,
    debug=_DEBUG)


if not os.access(DATA_FILENAME, os.R_OK):
    ode_ys_var2 = markov_tapes.ode_integrate(
        tag='ex3var2-copolymerization',
        size_a=4, cl_k=CL_K,
        p0=p0,
        ts=ode_ts2,
        odeint_kwargs=odeint_kwargs,        
        debug=_DEBUG)
    numpy.savez_compressed(DATA_FILE, ode_ys_var2=ode_ys_var2)    


ode_ys_var2 = numpy.load(DATA_FILE)['ode_ys_var2']


for suffix, history, ts in (('', ode_ys, ode_ts),
                        ('_var1', ode_ys_var1, ode_ts),
                        ('_var2', ode_ys_var2, ode_ts2)):
  fig = pyplot.figure(figsize=(16, 8))
  ax = fig.gca()
  ax.grid()
  def plot_seq_prob(ode_ys, ts, seq, style, label, **extra):
    ys = [markov_tapes.seq_prob(spd.reshape(*[4]*CL_K), seq)[0]
          for spd in ode_ys[1:, :]]
    ax.plot(ts[1:], (1e-100 + numpy.log(ys)) / numpy.log(10),
            style, label=label, **extra)
  plot_seq_prob(history, ts, [0, 1, 0], '-k', label='p(OAO)')
  plot_seq_prob(history, ts, [0, 2, 0], '-g', label='p(OMO)')
  plot_seq_prob(history, ts, [0, 1, 2, 0], '-m', label='p(OAMO)')
  plot_seq_prob(history, ts, [0, 2, 1, 3, 0], '-c', label='p(OMANO)')
  plot_seq_prob(history, ts, [0, 2, 1, 2, 0], '-r', label='p(OMAMO)')
  plot_seq_prob(history, ts, [1, 3, 1, 2], '-b', label='p(ANAM)')
  plot_seq_prob(history, ts, [1, 3, 1, 3], '--b', label='p(ANAN)', linewidth=3)
  ax.set_ylabel(r'$\log_{10}(p)$')
  ax.set_xlabel('time')
  ax.legend(loc='best')
  fig.savefig(f'ex3_copolymerization{suffix}.pdf')
  fig.show()

