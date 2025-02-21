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

"""Interfacing compiled Gambit-C code."""


# Ideally, a Python/C-interface is built such that even providing bad
# arguments from the Python side cannot crash the process.
#
# This generally requires dynamic PyObject type checking in C code.
#
# Here, we instead go with "users must not interfere with the contents
# of this module, since doing so may endanger memory-correctness",
# since this simplifies the code.

import atexit
import ctypes
import itertools
import os
import sys
import types

import numpy
import scipy.integrate


IS_DEBUG = bool(int(os.getenv('MARKOV_TAPES_DEBUG', '0')))

u_lib = ctypes.CDLL(os.path.join(os.path.dirname(__file__),
                                 'tapes_py_interface.so'))

u_lib.setup_gambit.restype = ctypes.c_void_p
u_lib.setup_gambit.argtypes = []

u_lib.cleanup_gambit.restype = None
u_lib.cleanup_gambit.argtypes = [ctypes.c_void_p]

u_lib.c_register_problems.restype = ctypes.c_int64
u_lib.c_register_problems.argtypes = [ctypes.c_int64]


u_lib.c_compute_dy_dt.restype = None
u_lib.c_compute_dy_dt.argtypes = [ctypes.c_void_p,
                                  ctypes.c_int64, ctypes.c_int64,
                                  ctypes.c_void_p, ctypes.c_void_p]

def init_gambit():
  """Initializes the Gambit-C runtime."""
  # TODO: Depending on how Gambit-C is exited (especially after an error that
  # opens the Gambit REPL), the process may segfault. See if this can be fixed.
  if IS_DEBUG:
    print('DEBUG: initializing gambit.', file=sys.stderr, flush=True)
  gambit = u_lib.setup_gambit()
  def _cleanup_gambit():
    if IS_DEBUG:
      print('DEBUG: cleanup gambit.', file=sys.stderr, flush=True)
    u_lib.cleanup_gambit(gambit)
  atexit.register(_cleanup_gambit)
  if IS_DEBUG:
    print('DEBUG: registering problems.', file=sys.stderr, flush=True)
  if 124 != u_lib.c_register_problems(123):
    # The `c_register_problems` function has "canary" behavior and returns
    # one more than its input. If this does not work, something is off
    # with the Scheme<->C interface, such as misaligned data-marshalling.
    raise ValueError('Registering problems failed.')


### Helpers

def mpp_from_spd(spd, eps=None):
  """Computes Markov Process Parameters from Sequence Probability Distribution.
 
  Args:
    spd: numpy ArrayLike of shape `B+(N,)*k`, where `N` is th alphabet size,
      `B` is a shape-prefix (such as for batch indices), and `spd[b + indices]`
      is the probability to find the sequence denoted by the `k`-index-tuple
      `indices` when probing a tape at a random position and reading `k` symbols
      in total, while `b` is an index-prefix. 
    eps: Optional offset to add to subsequence probabilities to handle cases
      with an impossible prefix, i.e. p({prefix} + (s,)) = 0 for all s.
      In this case, it is better to have each follow-up symbol have the same
      probability than encountering division-by-zero. If `None`, a default is
      used that is a positive offset far smaller than float64 ulp(1).
  
  Returns:
    numpy.ndarray `r` of same shape as `spd` where `r[prefix + (s,)]` with a
    (k-1)-index-tuple `prefix` is the probability for the sequence denoted by
    `prefix` to be followed by `s`.
  """
  if eps is None: eps = 1e-100
  spd_normalized = numpy.clip(
    numpy.asarray(spd).astype(numpy.float64), eps, 1)
  return spd_normalized / spd_normalized.sum(axis=-1, keepdims=True)


def ctm_from_mpp(num_alphabet, num_context, mpp):
  """Computes a Context Transfer Matrix from Markov Process Parameters.

  Args:
    num_alphabet: Number of symbols in the alphabet.
    num_context: length of the context-prefix subsequence from which
      the next symbol's probability is determined.
    mpp: [num_alphabet]*(num_context+1) numpy.ndarray such that
      `mpp[*prefix, i]` is the probability for the symbol-index
      sequence `prefix` to be followed by the symbol with index `i`.

  Returns:
    [num_alphabet**num_context, num_alphabet**num_context]-ndarray
    with prefix-index-sequence transition probabilities.
  """
  result = numpy.zeros([num_alphabet ** num_context] * 2)
  result_stepwise = result.reshape([num_alphabet] * (2 * num_context))
  mp_stepwise = mpp.reshape([num_alphabet] * (1 + num_context))
  # There may be more elegant ways to express this multiindex operation,
  # but this is likely clearest:
  for indices in itertools.product(range(num_alphabet), repeat=num_context + 1):
    prob = mp_stepwise[indices]
    result_stepwise[indices[1:] + indices[:-1]] += prob
  return result


def get_ctm_eigenvalue1_eigenspace(spd,
                                   eps_mpp=None,
                                   eps=1e-7):
  """Computes the eigenvalue=1 transfer matrix eigenspace.

  Args:
    spd: sequence probability distribution, as for `mpp_from_spd()`.
    eps_mpp: `eps`-parameter for `mpp_from_spd()`.
    eps: Maximal magnitude-deviation of the eigenvalue from 1,
      or maximal deviation between left and right marginals of `spd`.

  Returns:
    If left and right marginals are compatible, a pair
    `(deviation, eigenspace)`, where `deviation` is the least-squares
    distance for obtaining the original `spd` as a linear combination of
    eigenvalue-1 eigenvectors, and `eigenspace` is a
    `[num_alphabet**num_context, dim_eigenspace]-numpy.ndarray` where each
    column-vector is a generalized eigenvalue-1-eigenspace basis vector.
    If they are not, returns `(marginals_distance, None)`
  Raises:
    ValueError: If marginal subsequence-probabilities over leading and trailing
      index differ by more than `eps`.
  """
  spd = numpy.asarray(spd, dtype=numpy.float64)
  num_alphabet = spd.shape[0]
  num_context = spd.ndim - 1
  spd_marginal_right = spd.sum(axis=-1)
  spd_marginal_left = spd.sum(axis=0)
  marginal_distance = numpy.linalg.norm(spd_marginal_right.ravel() -
                                        spd_marginal_left.ravel())
  if not marginal_distance <= eps:
    return marginal_distance, None
  mpp = mpp_from_spd(spd, eps=eps_mpp)
  ctm = ctm_from_mpp(num_alphabet, num_context, mpp)
  # Need .eig() since in general ctm != ctm.T  
  eigvals, eigvecs = numpy.linalg.eig(ctm)
  good_eigvals = abs(eigvals - 1.0) <= eps
  eigenspace = eigvecs[:, good_eigvals]
  coeffs, residuals, *_ = numpy.linalg.lstsq(eigenspace,
                                             spd_marginal_left.ravel(),
                                             rcond=None)
  # Residuals are squares-of-coordinate-distances.
  return numpy.linalg.norm(residuals**.5), eigenspace


def markov_entropy(spd):
  """Computes the Markov chain entropy of the subseq p.d."""
  # See: https://en.wikipedia.org/wiki/Entropy_rate
  eps = 1e-280
  spd_normalized = numpy.clip(
    numpy.asarray(spd).astype(numpy.float64), eps, 1)
  spd_reduced = spd_normalized.sum(axis=-1)
  p_conditional = spd_normalized / spd_reduced[..., numpy.newaxis]
  return (-p_conditional * numpy.log(p_conditional)
          ).sum(axis=-1).ravel().dot(spd_reduced.ravel())


def seq_prob(spd, seq, *,
             num_prefix_indices=0, eps=None,
             mpp=None, want_mpp=False):
  """Computes the probability of a sequence given a probability distribution.

  Args:
    spd: sequence probability distribution, as for `mpp_from_spd()`.
    seq: symbol-sequence to compute the probability for.
    num_prefix_indices: Number of leading non-sequence indices on the
      subsequence probability distribution (such as: a single time-step index).
    eps: Epsilon to use for `mpp_from_spd()`.
    mpp: `None` or Markov process parameters. If `None`, these will be
      (re-)computed from `spd`. (Providing this parameter can speed up
      re-evaluation.)
    want_mpp: Whether Markov process parameters should be returned.

  Returns:
    A pair `(probability, mpp)`, where `probability` is the
    subsequence-probability (with `num_prefix_indices` many
    inert prefix-indices), and `mpp` are the Markov process
    parameters - either as provided or as computed via
    `mpp_from_spd()`.
  """
  spd = numpy.asarray(spd, dtype=numpy.float64)
  num_sequence_indices = spd.ndim - num_prefix_indices
  num_excess_sequence_indices = num_sequence_indices - len(seq)
  if num_excess_sequence_indices >= 0:
   # We have sufficiently many indices to satisfy the request.
   return (
     spd[..., *seq].sum(
       axis=tuple(range(num_prefix_indices,
                        num_prefix_indices + num_excess_sequence_indices))),
     mpp_from_spd(spd, eps=eps) if want_mpp else mpp)
  # Otherwise, we will have to compute probabilities via iterative expansion.
  # This requires knowing the Markov Process Parameters.
  if mpp is None:  # Compute if not provided.
    mpp = mpp_from_spd(spd, eps=eps)
  p_current = spd[..., *seq[:num_sequence_indices]]
  tail_todo = seq[1:]
  while len(tail_todo) >= num_sequence_indices:
    context = tail_todo[:num_sequence_indices]
    p_current = mpp[..., *context] * p_current
    tail_todo = tail_todo[1:]
  return p_current, mpp
    
 
def tprint(size_a, cl_k, adata, epsilon=1e-10, nmax=float('inf'), file=None):
  """Debug-prints non-zero entries of a Markov transition table.

  Args:
    size_a: Size of the symbol-alphabet.
    cl_k: Length of the subsequences for which we have probabilities.
    adata: probability-array, must be reshapeable to
      `[size_a] * (2 * (cl_k - 1))`.
    epsilon: Magnitude threshold for printing an entry.
    nmax: Maximal number of entries to print.
    file: `file=` parameter to forward to `print()`.
  """
  num_indices_in = cl_k - 1
  a_multiindex = numpy.asarray(adata).reshape([size_a] * (2 * num_indices_in))
  for n, indices in enumerate(itertools.product(range(size_a),
                                                repeat=2*num_indices_in)):
    if n >= nmax:
      print('... more entries...', file=file)
    val = a_multiindex[indices]
    if not abs(val) < epsilon:
      print(f'{indices[:num_indices_in]} {indices[num_indices_in:]}: {val}')


def get_dy_dt(*, tag, size_a, cl_k, debug=False):
  """Returns the dy/dt-function to be used for ODE-integration.

  Args:
    tag: Name of the problem under which it was registered by Scheme code.
    size_a: Size of the symbol-alphabet.
    cl_k: Length of the subsequences for which we have probabilities.
    debug: whether to turn on debugging output for Scheme.

  Returns:
    The `(probabilities_in, t) -> d_dt_probabilities` ODE right-hand-side function
    for this problem.
  """
  a_tag = numpy.frombuffer(tag.encode() + b'\x00', dtype=numpy.uint8)
  do_debug = 1 if debug else 0
  expected_size = size_a ** cl_k 
  def dy_dt(a_probs_in, t):
    # del t  # Unused, required by ODE-solver interface.
    print(f'DDD {t=:.10g}')
    c_probs_in = numpy.asarray(a_probs_in, dtype=numpy.float64).ravel()
    c_probs_out = numpy.zeros_like(c_probs_in)
    if c_probs_in.size != expected_size:
      raise ValueError(f'probability-array should have size {expected_size}, '
                       f'observed: {c_probs_in.size}')
    u_lib.c_compute_dy_dt(
      a_tag.__array_interface__['data'][0],  # tag-string.
      cl_k, do_debug,
      c_probs_in.__array_interface__['data'][0],
      c_probs_out.__array_interface__['data'][0])
    return c_probs_out
  return dy_dt


def ode_integrate(*, tag, size_a, cl_k, p0, ts,
                  odeint_kwargs=types.MappingProxyType({}),
                  debug=False):
  """Performs ODE-integration via `scipy.integrate.odeint()`.

  Args:
    tag: Name of the problem under which it was registered by Scheme code.
    size_a: Size of the symbol-alphabet.
    cl_k: Length of the subsequences for which we have probabilities.
    p0: Initial subsequence probability-distribution. Array-like, must have
      total size `size_a**cl_k`.
    ts: ArrayLike of requested points in time for ODE-integration.
    odeint_kwargs: Mapping with extra `scipy.integrate.odeint()` keyword
      arguments.
    debug: Whether to get a `dy/dt`-function that has debug-printing enabled.

  Returns:
    The array returned by calling `scipy.integrate.odeint()` for this problem.
  """
  p0 = numpy.asarray(p0, dtype=numpy.float64).ravel()
  if not (p0.size == size_a**cl_k and
          (0 <= p0).all() and (p0 <= 1).all() and
          abs(p0.sum() - 1) < 1e-10):
    raise ValueError(
      'Parameter p0 is not a subsequence probability distribution.')
  dy_dt = get_dy_dt(tag=tag, size_a=size_a, cl_k=cl_k, debug=debug)
  return scipy.integrate.odeint(dy_dt, p0, ts, **odeint_kwargs)


def ode_integrate_ivp(*, tag, size_a, cl_k, p0, ts,
                     ivp_kwargs=types.MappingProxyType({}),
                     debug=False):
  """Performs ODE-integration via `scipy.integrate.solve_ivp()`.

  Args:
    tag: Name of the problem under which it was registered by Scheme code.
    size_a: Size of the symbol-alphabet.
    cl_k: Length of the subsequences for which we have probabilities.
    p0: Initial subsequence probability-distribution. Array-like, must have
      total size `size_a**cl_k`.
    ts: ArrayLike of requested points in time for ODE-integration.
    ivp_kwargs: Mapping with extra `scipy.integrate.solve_ivp()` keyword
      arguments.
    debug: Whether to get a `dy/dt`-function that has debug-printing enabled.

  Returns:
    The data-array returned by calling `scipy.integrate.solve_ivp()` for this
    problem, brought into the shape that `scipy.integrate.odeint()` would
    return.
  """
  p0 = numpy.asarray(p0, dtype=numpy.float64).ravel()
  if not (p0.size == size_a**cl_k and
          (0 <= p0).all() and (p0 <= 1).all() and
          abs(p0.sum() - 1) < 1e-10):
    raise ValueError(
      'Parameter p0 is not a subsequence probability distribution.')
  dy_dt = get_dy_dt(tag=tag, size_a=size_a, cl_k=cl_k, debug=debug)
  return scipy.integrate.solve_ivp(
    lambda t, y: dy_dt(y, t),
    (ts[0], ts[-1]),
    p0,
    t_eval=ts,
    **ivp_kwargs).y.T


def _run_validation():
  fn_dy_dt = get_dy_dt(tag='__canary_problem_radioactive_decay',
                       size_a=2, cl_k=3, debug=False)
  observed = fn_dy_dt(numpy.full([8], fill_value=0.125,
                                 dtype=numpy.float64), 0.0).tolist()
  expected = [0.375, 0.125, 0.125, -0.125, 0.125, -0.125, -0.125, -0.375]
  if expected != observed:
    raise RuntimeError(
      'Load-time validation problem failed to produce the expected result.')
  

######

# We unconditionally initialize Gambit at import time...
init_gambit()

# ...and validate that we can call Scheme code by running a smoke test example.
_run_validation()
