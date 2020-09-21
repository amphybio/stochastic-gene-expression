# vim: fileencoding=utf-8

# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, you can obtain one at https://mozilla.org/MPL/2.0/
#
# Copyright 2020 Alexandre Ferreira Ramos - AMPhyBio Laboratory
#
# Project:  github.com/amphybio/stochastic-gene-expression
# Version:  1.0
# Created:  20-09-2019
# Authors:  Leonardo R. Gama <leonardo.gama@usp.br>

"""
Shannon's entropy for stochastic gene models.

Functions to calculate the Shannon's entropy, entropy conditional to promoter
state and mutual information for constitutive and binary stochastic gene models
using SymPy, Maple or compiled C code.
"""

__all__ = [
        'H_constitutive',
        'H_external', 'H_ON_external', 'H_OFF_external', 'I_external',
        'set_n_processes', 'symbolic_H',
]

import atexit
import logging
import math
import multiprocessing as mp
import os
import pathlib
import signal
import subprocess as sub
from multiprocessing import pool

import mpmath
from sympy import *
from sympy.abc import *
from sympy import E, N as evalf

import utils
from functions import log2
from steady_state import alpha_n_external, beta_n_external, palpha, phi_n_external


def set_n_processes(n):
    """Initialize number of worker processes in pool. Should be called right after imports."""
    global N_PROCESSES
    N_PROCESSES = n
set_n_processes(os.cpu_count())


DOC = {
    'epsilon': "ratio between promotor switching rates and protein degradation rate",
    'palpha': "probability of finding the promotor at the ON state",
    'N': "mean number of proteins of a constitutive gene with the same synthesis/degradation rates",
    'k': "upper bound of summation for the entropy calculation",
    'precision': "number of decimal digits of precision",
    'method': "either 'C', 'maple', 'maple-async', 'sympy' or 'sympy-parallel'",
    'backup_method': "(list of) backup method(s) to try if 'method' fails",
    'func': "function to be computed",
    'subs': "dictionay with parameters to 'func'",
}

symbolic_H = {}


### Constitutive gene ###

""" Shannon entropy for a constitutive gene with mean expression N.

                      oo     n
               ⎛N⎞   \⎺⎺⎺` ⎛N   -N         ⎞
    H = -N⋅log₂⎜─⎟ +  ⟩    ⎜──⋅e  ⋅log₂(n!)⎟
               ⎝e⎠   /⎽⎽⎽, ⎝n!             ⎠
                     n = 0
"""
H_poisson_const = -N*log2(N/E)
H_poisson_sum_term = N**n / factorial(n) * exp(-N) * log2(factorial(n))
H_poisson_sum = Sum(H_poisson_sum_term, (n, 0, oo))
symbolic_H['constitutive'] = H_poisson_const + H_poisson_sum

@utils.memoized
def H_constitutive(N, precision=mpmath.mp.dps):
    """Shannon entropy for the constitutive gene model.

    :N: mean number of proteins
    :precision: {precision}
    :returns: entropy of a gene with parameter 'N'
    """
    return symbolic_H['constitutive'].evalf(precision, subs={'N': N})


### Binary gene ###

# Auxiliary expressions for parallel computation on SymPy.
parallel_H = {}
parallel_n = i*N_PROCESSES + c  # i: index; c: constant


## External Regulation Gene ##

""" Shannon entropy for the binary stochastic gene model.

          ∞
    H = - ∑  φₙ⋅log₂(φₙ)
         n=0
"""

H_external_sum_term = phi_n_external*log2(phi_n_external)
symbolic_H['external'] = -Sum(H_external_sum_term, (n, 0, oo))
parallel_H['external'] = -Sum(H_external_sum_term.subs(n, parallel_n), (i, 0, oo))

def H_external(epsilon, palpha, N, k=oo, precision=mpmath.mp.dps, method='sympy-parallel', backup_method=None):
    """Shannon entropy for the externally regulated gene model.

    :epsilon: {epsilon}
    :palpha: {palpha}
    :N: {N}
    :k: {k}
    :precision: {precision}
    :method: {method}
    :backup_method: {backup_method}
    :returns: shannon entropy of gene with parameters ε, pₐ and N
    """
    subs = {'epsilon': epsilon, 'p_a': palpha, 'N': N}
    return _H_dispatch('external', subs, k, precision, method, backup_method)


"""Shannon entropy of the number of proteins conditional to the promoter state being ON or OFF.

            oo
           \⎺⎺⎺` αₙ     ⎛αₙ⎞
    H   = - ⟩    ──⋅log₂⎜──⎟
     ON    /⎽⎽⎽, pₐ     ⎝pₐ⎠
           n = 0
"""

H_ON_external_sum_term = alpha_n_external/palpha*log2(alpha_n_external/palpha)
symbolic_H['ON_external'] = -Sum(H_ON_external_sum_term, (n, 0, oo))
parallel_H['ON_external'] = -Sum(H_ON_external_sum_term.subs(n, parallel_n), (i, 0, oo))

def H_ON_external(epsilon, palpha, N, k=oo, precision=mpmath.mp.dps, method='sympy-parallel', backup_method=None):
    """Entropy conditional to ON state for the externally regulated gene model.

    :epsilon: {epsilon}
    :palpha: {palpha}
    :N: {N}
    :k: {k}
    :precision: {precision}
    :method: {method}
    :backup_method: {backup_method}
    :returns: shannon entropy of the distribution given the promoter is at ON state
    """
    subs = {'epsilon': epsilon, 'p_a': palpha, 'N': N}
    return _H_dispatch('ON_external', subs, k, precision, method, backup_method)

H_OFF_external_sum_term = beta_n_external/(1 - palpha)*log2(beta_n_external/(1 - palpha))
symbolic_H['OFF_external'] = -Sum(H_OFF_external_sum_term, (n, 0, oo))
parallel_H['OFF_external'] = -Sum(H_OFF_external_sum_term.subs(n, parallel_n), (i, 0, oo))

def H_OFF_external(epsilon, palpha, N, k=oo, precision=mpmath.mp.dps, method='sympy-parallel', backup_method=None):
    """Entropy conditional to OFF state for the externally regulated gene model.

    :epsilon: {epsilon}
    :palpha: {palpha}
    :N: {N}
    :k: {k}
    :precision: {precision}
    :method: {method}
    :backup_method: {backup_method}
    :returns: shannon entropy of the distribution given the promoter is at OFF state
    """
    subs = {'epsilon': epsilon, 'p_a': palpha, 'N': N}
    return _H_dispatch('OFF_external', subs, k, precision, method, backup_method)


"""Mutual information for the externally regulated gene model.

                                  ⎛                   ⎞
    I(X; Y) = H(X) - H(X|Y) = H - ⎜pₐ⋅H  + (1-pₐ)⋅H   ⎟
                                  ⎝    ON          OFF⎠
"""

class IAsyncResult(pool.AsyncResult):
    """Wrapper to multiple AsyncResult's"""
    def __init__(self, palpha, hs):
        self.palpha = palpha
        is_async = (isinstance(res, pool.AsyncResult) for res in hs)
        self.results = tuple(zip(hs, is_async))

    def ready(self):
        return all(not is_async or res.ready() for res, is_async in self.results)

    def successful(self):
        return all(not is_async or res.successful() for res, is_async in self.results)

    def wait(self, timeout=None):
        for res, is_async in self.results:
            if is_async:
                res.wait(timeout)

    def get(self, timeout=None):
        h, h_on, h_off = (res.get(timeout) if is_async else res for res, is_async in self.results)
        return h - self.palpha*h_on - (1 - self.palpha)*h_off

def I_external(epsilon, palpha, N, k=oo, precision=mpmath.mp.dps, method='sympy-parallel', backup_method=None):
    """Mutual information for the externally regulated gene model.

    :epsilon: {epsilon}
    :palpha: {palpha}
    :N: {N}
    :k: {k}
    :precision: {precision}
    :method: {method}
    :backup_method: {backup_method}
    :returns: mutual information of gene with parameters ε, pₐ and N
    """
    h = H_external(epsilon, palpha, N, k, precision, method, backup_method)
    h_on = H_ON_external(epsilon, palpha, N, k, precision, method, backup_method)
    h_off = H_OFF_external(epsilon, palpha, N, k, precision, method, backup_method)
    hs = (h, h_on, h_off)

    if any(res is None for res in hs):
        return None
    elif any(isinstance(res, pool.AsyncResult) for res in hs):
        return IAsyncResult(palpha, hs)
    else:
        return h - palpha*h_on - (1 - palpha)*h_off


### Internals ###

def _H_dispatch(func, subs, k, precision, method, backup_method):
    """Backend calculation dispatcher.

    :func: {func}
    :subs: {subs}
    :k: {k}
    :precision: {precision}
    :method: {method}
    :backup_method: {backup_method}
    :returns: result of 'func' calculation using the required method(s)
    """
    assert method in {'maple', 'maple-async', 'sympy', 'sympy-parallel'}
    if isinstance(backup_method, str):
        backup_method = [backup_method]

    elif method == 'maple':
        res = _H_maple(func, subs, k, precision)
    elif method == 'maple-async':
        # FIXME: avoid 'backup_method' bypass by a None returned later by an pool.AsyncResult.
        # Things should work in a program rerun since (some of) these None async results would be
        # already cached.
        res = global_pool().apply_async(_H_maple, [func, subs, k, precision])
        try:
            res = res.get(timeout=0.1)
        except mp.TimeoutError:
            return res

    else:  # method is 'sympy' or 'sympy-parallel'
        res = _H_sympy(func, subs, k, precision, parallel=method.endswith('parallel'))

    if res is None and backup_method:
        log_msg = "_H_dispatch: method '%s' failed, trying '%s' for '%s' with parameters %s"
        logging.debug(log_msg, method, backup_method[0], func, str(subs))
        return _H_dispatch(func, subs, k, precision, method=backup_method[0], backup_method=backup_method[1:])

    return res

mp_pool = None
def global_pool():
    """Initialize pool of worker processes only when a function requires it."""
    global mp_pool
    if mp_pool is None:
        mp_pool = mp.Pool(processes=N_PROCESSES)
        atexit.register(mp_pool.close)
    return mp_pool

def _map_evalf(arg):
    """Auxiliary function for parallel numeric evaluation.

    :arg: 2-tuple with the variable x and the precision
    :returns: result equivalent to expr.evalf(x, n)
    """
    return evalf(x=arg[0], n=arg[1])

@utils.memoized
def _H_sympy(func, subs, k, precision, parallel):
    """Calculate entropy in SymPy.

    :func: {func}
    :subs: {subs}
    :k: {k}
    :precision: {precision}
    :parallel: wether to run summation in parallel
    :returns: result of 'func' evaluation in SymPy with parameters in 'subs'
    """
    expr = parallel_H[func] if parallel else symbolic_H[func]

    try:
        if not parallel:
            expr = expr.replace(oo, k)
            res = expr.evalf(precision, subs)
            if res == 0:
                raise RuntimeError
        else:
            # Parallel evaluation requires integers or fractions(?).
            expr = expr.subs({key: Rational(str(val)) for key, val in subs.items()})
            expr = expr.replace(oo, k/N_PROCESSES)
            args = [(expr.subs('c', c), precision) for c in range(N_PROCESSES)]
            partial_sums = global_pool().map(_map_evalf, args)
            if any(x == 0 for x in partial_sums):
                raise RuntimeError
            res = sum(partial_sums)
            if logging.getLogger().level >= logging.INFO:
                print(".", end="", flush=True)  # show progress
        return res

    # Note: NaN is stored as None (NULL) in diskcache (SQLite).
    except (RuntimeError, TypeError):
        logging.debug("_H_sympy: invalid result with precision = {}.".format(precision))
        if precision >= 75:
            return None
        return _H_sympy(func, subs, k, precision + 15, parallel)
    except (mpmath.libmp.NoConvergence, ValueError):
        logging.debug("_H_sympy: convergence exception with parameters ε = %(epsilon)f, pₐ = %(p_a)f, N = %(N)d", subs)
        return None

root = pathlib.Path(__file__).parent.resolve()
maple_external = root/'entropy_external.mpl'
@utils.memoized()
def _H_maple(func, subs, k, precision):
    """Calculate entropy using Maple.

    :func: {func}
    :subs: {subs}
    :k: {k}
    :precision: {precision}
    :returns: result of 'func' evaluation in Maple with parameters in 'subs'
    """
    maple_func = 'H_' + func
    args = (maple_func, subs['epsilon'], subs['p_a'], subs['N'], precision)
    args = '-cp:=' + ','.join(str(a) for a in args)
    if not k is oo:
        args += ',' + str(k)
    logging.debug("_H_maple: calling Maple with command: %s %s", maple_external, args)

    with sub.Popen([maple_external, args], stdout=sub.PIPE,
            start_new_session=True, universal_newlines=True) as proc:
        pgid = os.getpgid(proc.pid)

        # Maple subprocesses like to lie around forever... So we KILL it!
        @atexit.register
        def kill_maple(who='atexit'):
            try:
                os.killpg(pgid, signal.SIGKILL)
                logging.debug("%s: killing process group %d", who, pgid)
            except ProcessLookupError:
                pass

        try:
            proc.wait(600)
            res = Float(proc.stdout.readline())
            if math.isnan(res):
                logging.debug("_H_sympy: invalid result with precision = {}.".format(precision))
                if precision >= 75:
                    return None
                return _H_maple(func, subs, k, precision + 15)
        except sub.TimeoutExpired:
            proc.terminate()
            res = None
        finally:
            atexit.unregister(kill_maple)
            kill_maple(who='_H_maple')

    if logging.getLogger().level >= logging.INFO:
        print(".", end="", flush=True)  # show progress
    return res


for func in (H_constitutive, H_external, H_ON_external, H_OFF_external, I_external, _H_dispatch, _H_maple, _H_C):
    func.__doc__ = func.__doc__.format(**DOC)
