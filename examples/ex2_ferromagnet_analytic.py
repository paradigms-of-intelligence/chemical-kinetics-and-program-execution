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

"""Ferromagnetic Chain Example: Approximative analytic treatment.

This is the benchmark against which we want to compare both the
Monte-Carlo as well as the Markov-Parameter-Dynamics results.

"""

import numpy
import scipy.integrate


def get_dy_dt_matrix(param_a, param_b, num_lengths_to_track):
  result = numpy.zeros([num_lengths_to_track, num_lengths_to_track])
  result[0, 0] = -1  # Melting of a length-1 chain.
  for k in range(1, num_lengths_to_track):
    # k -> k-1 "melting"
    result[k-1, k] += 2 * param_a
    result[k, k]   -= 2 * param_a
    # Growing.
    result[k, k-1] += 2 * param_a * param_b
    result[k, k] -= 2 * param_a * param_b
  return result


def get_p_history(*,
                  beta=1.0, J=1.0, h=-0.25,
                  num_lengths_to_track=20,
                  t_max=40, t_steps=1001,
                  p0_pair_start=1/250,
                  rtol=1e-10, atol=1e-10
                  ):
  # Fastest rate is dud -> ddd: delta-E = -4J+2h
  # Second fastest rate is uud -> udd: delta-E = 2h
  # Third fastest rate is udd -> uud: delta-E = -2h
  mm = get_dy_dt_matrix(param_a=numpy.exp(-beta*4*J),
                        param_b=numpy.exp(beta*2*h),  # Note h<0 here!
                        num_lengths_to_track=num_lengths_to_track)
  ts = numpy.linspace(0, t_max, t_steps)
  # We do have to take spontaneous creation of 1-chains into account.
  r_formation = numpy.array([numpy.exp(-8*beta*J + 2*beta*h)] +
                            [0.0] * (num_lengths_to_track-1))
  p_history = scipy.integrate.odeint(
      lambda y, t: mm.dot(y) + r_formation,
      numpy.array([p0_pair_start if k==1 else 0.0
                   for k in range(num_lengths_to_track)]),
      ts, rtol=rtol, atol=atol)
  return numpy.clip(p_history, 0, numpy.inf)
