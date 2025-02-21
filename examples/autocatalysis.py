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


# JAX requires a hack to enable float64 support.
if 'HACK':
  import os    
  if not os.environ.get('JAX_ENABLE_X64', ''):
    # If this has been set, but explicitly to '0', we go with that
    # user-override.
    import warnings
    warnings.warn(
      'Setting JAX_ENABLE_X64=1 environment variable to enable '
      'jax.numpy.float64/.complex128 arrays')
    os.environ['JAX_ENABLE_X64'] = '1'
  #
  if not os.environ.get('JAX_TRACEBACK_FILTERING', ''):
    # If this has been set, but explicitly to '0', we go with that
    # user-override.
    import warnings
    warnings.warn(
      'Setting JAX_TRACEBACK_FILTERING=off environment variable to enable '
      'meaningful tracebacks.')
    os.environ['JAX_TRACEBACK_FILTERING'] = 'off'


import jax
from jax import numpy as jnp
import numpy
import scipy.integrate
import scipy.optimize
import matplotlib
from matplotlib import pyplot

matplotlib.rcParams.update({'font.size': 18})


# ===
PARAMS10 = jnp.array([
  0.0, 0.0, 1.0,
  0.001, 20.0, 10.0,
  0.001, 50.0, 20.0,
  0.0, 0.0], dtype=jnp.float64)
PARAMS11 = jnp.array([
  0.2, 0.1, 0.4,
  0.001, 20.0, 10.0,
  0.001, 50.0, 20.0,
  0.0, 0.0], dtype=jnp.float64)
PARAMS12 = jnp.array([
  0.0, 0.0, 1.0,
  0.001, 20.0, 10.0,
  0.001, 80.0, 20.0,
  0.0, 0.0], dtype=jnp.float64)
PARAMS13 = jnp.array([
  0.0, 0.0, 1.0,
  0.001, 50.0, 10.0,
  0.001, 20.0, 20.0,
  0.0, 0.0], dtype=jnp.float64)
PARAM_SET1 = ((0, '-', PARAMS10), (1, '--', PARAMS11),
              (2, '-.', PARAMS12), (3, ':', PARAMS13))

PARAMS20 = jnp.array([
  0.0, 0.0, 1.0,  
  0.001, 20.0, 10.0,
  0.001, 50.0, 20.0,
  0.0, 0.0], dtype=jnp.float64)
PARAMS21 = jnp.array([
  0.0, 0.0, 1.0,
  0.01, 20.0, 10.0,
  0.01, 50.0, 20.0,
  0.1, 0.1], dtype=jnp.float64)
PARAMS22 = jnp.array([
  0.0, 0.0, 1.0,  
  0.01, 20.0, 10.0,
  0.01, 50.0, 20.0,
  0.5, 0.5], dtype=jnp.float64)
PARAMS23 = jnp.array([
  0.0, 0.0, 1.0,
  0.01, 20.0, 10.0,
  0.01, 50.0, 20.0,
  10, 10], dtype=jnp.float64)
PARAM_SET2 = ((0, '-', PARAMS20), (1, '--', PARAMS21),
              (2, '-.', PARAMS22), (3, ':', PARAMS23))

PARAMS30 = jnp.array([
  0.0, 0.0, 1.0,
  0.05, 20.0, 10.0,
  0.05, 25.0, 10.0,
  0.1, 0.1], dtype=jnp.float64)
PARAMS31 = jnp.array([
  0.0, 0.0, 1.0,
  0.05, 20.0, 10.0,
  0.05, 25.0, 10.0,
  1.0, 1.0], dtype=jnp.float64)
PARAMS32 = jnp.array([
  0.0, 0.0, 1.0,
  0.05, 20.0, 10.0,
  0.05, 25.0, 10.0,
  5.0, 5.0], dtype=jnp.float64)
PARAMS33 = jnp.array([
    0.0, 0.0, 1.0,
  0.05, 20.0, 10.0,
  0.05, 25.0, 10.0,
  30.0, 30.0], dtype=jnp.float64)
# High flow rate, no autocatalysis.
# PARAMS33 = jnp.array([
#     0.0, 0.0, 1.0,
#   0.05, 0.0, 10.0,
#   0.05, 0.0, 10.0,
#   100.0, 100.0], dtype=jnp.float64)
PARAM_SET3 = ((0, '-', PARAMS30), (1, '--', PARAMS31),
              (2, '-.', PARAMS32), (3, ':', PARAMS33))


@jax.jit
def fn_dy_dt(y, params):
  (c_form_a, c_auto_a, c_stab_a, c_form_b, c_auto_b, c_stab_b, c_add, c_remove
   ) = params
  # Spontaneuos-pathway and autocatalytic-pathway dissociation
  c_sdiss_a = c_form_a / c_stab_a
  c_adiss_a = c_auto_a / c_stab_a
  c_sdiss_b = c_form_b / c_stab_b
  c_adiss_b = c_auto_b / c_stab_b
  #
  ca, cb, cm = y  # A-dimer, B-dimer, Monomer.
  return jnp.array([
      (+ c_form_a * cm * cm + c_auto_a * ca * cm * cm
       - c_sdiss_a * ca - c_adiss_a * ca * ca
       - c_remove * ca
       ),
      (+ c_form_b * cm * cm + c_auto_b * cb * cm * cm
       - c_sdiss_b * cb - c_adiss_b * cb * cb
       - c_remove * cb
       ),
      (+ 2 * (c_sdiss_a * ca + c_sdiss_b * cb)
       + 2 * (c_adiss_a * ca * ca + c_adiss_b * cb * cb)
       - 2 * (c_form_a * cm * cm + c_form_b * cm * cm)
       - 2 * (c_auto_a * ca * cm * cm + c_auto_b * cb * cm * cm)
       - c_remove * cm + c_add
       )], dtype=jnp.float64)


ts = numpy.linspace(0, 100, 10001)


for filename, param_set in (
    ('autocatalysis1.pdf', PARAM_SET1),
    ('autocatalysis2.pdf', PARAM_SET2),
    ('autocatalysis3.pdf', PARAM_SET3),
):
  fig = pyplot.figure(figsize=(12, 8))
  ax = fig.gca()
  ax.grid()
  def aplot(ts, ys, *args, **kwargs):
    ax.plot(numpy.log(ts) / numpy.log(10), ys, *args, **kwargs)
  for n, style, y0_and_params in param_set:
    y0 = y0_and_params[:3]
    params = y0_and_params[3:]
    ode_ys = scipy.integrate.odeint(lambda y, t: fn_dy_dt(y, params),
                                    y0, ts)
    aplot(ts[1:], ode_ys[1:, 0], style + 'b', label=('A' if n==0 else None))
    aplot(ts[1:], ode_ys[1:, 1], style + 'g', label=('B' if n==0 else None))
    aplot(ts[1:], ode_ys[1:, 2], style + 'r', label=('M' if n==0 else None))
    aplot(ts[1:], ode_ys[1:, 0] * 2 + ode_ys[1:, 1] * 2 + ode_ys[1:, 2],
            '-k', label=('M(total)' if n==0 else None))
  ax.set_ylabel('Concentration')
  ax.set_xlabel(r'$\log_{10}$(time)')
  ax.legend(loc='upper right')
  fig.savefig(filename)
  fig.show()


## Tools for finding the equilibrium. We generally want to start from
## an approximate solution, such as obtained as the endpoint of
## finite-time ODE-integration. In general, there will be multiple
## equilibria, and any given one might not be sensible (such as:
## having negative concentrations), or perhaps unstable under
## perturbations.

def get_equilibrium_fn(fn_dy_dt):
  @jax.jit
  def fn_f(y, params):
    dy_dt = fn_dy_dt(y, params)
    return dy_dt @ dy_dt
  fn_fprime = jax.grad(fn_f)
  #
  def fn_opt(y0, params):
    y0 = jnp.asarray(y0, dtype=jnp.float64)
    opt = scipy.optimize.fmin_bfgs(
      fn_f, y0, fprime=fn_fprime, gtol=1e-10,
      args=(params,),
      disp=0)
    return opt, float(fn_f(opt, params))
  #
  return fn_opt


fn_eq = get_equilibrium_fn(fn_dy_dt)

