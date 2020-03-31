# vim: fileencoding=utf-8

# Copyright 2020 Alexandre Ferreira Ramos
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# authors: Leonardo Gama (29/05/2019)

"""
Expressions for the steady-state solutions from the binary stochastic gene model.
"""

from sympy import *
from sympy.abc import *
import logging
import math
import mpmath as mp

from functions import KummerM, pochhammer


"""Distributions for the externally regulated gene from Ramos et al. (2007 & 2010).

Parameters:
    N = k/ρ
    ε = (f + h)/ρ
    pₐ = f/(f + h)
"""
palpha = Symbol('p_a')

def dist(x, y, N, n):
    """        n
              N  (x)ₙ
    dist(n) = ──⋅────⋅M(x + n, y + n, -N)
              n! (y)ₙ
    """
    return (N**n/factorial(n)) * (pochhammer(x, n)/pochhammer(y, n)) * KummerM(x + n, y + n, -N)

"""Marginal probability of finding n gene products at any moment.
         n
        N  (ε⋅pₐ)ₙ
   φₙ = ──⋅───────⋅M(ε⋅pₐ + n, ε + n, -N)
        n!  (ε)ₙ
"""
phi_n_external = dist(epsilon*palpha, epsilon, N, n)

"""Joint probability of finding n gene products *and* the promoter in the ON state.
             n
            N  (1 + ε⋅pₐ)ₙ
    αₙ = pₐ⋅──⋅───────────⋅M(1 + ε⋅pₐ + n, 1 + ε + n, -N)
            n!  (1 + ε)ₙ
"""
alpha_n_external = palpha * dist(1 + epsilon*palpha, 1 + epsilon, N, n)

"""Joint probability of finding n gene products *and* the promoter in the OFF state.
                   n
                  N  (ε⋅pₐ)ₙ
    βₙ = (1 - pₐ)⋅──⋅────────⋅M(ε⋅pₐ + n, 1 + ε + n, -N)
                  n! (1 + ε)ₙ
"""
beta_n_external = (1 - palpha) * dist(epsilon*palpha, 1 + epsilon, N, n)

"""Mean number of gene products at the ON state.
           ε⋅pₐ + 1
    <nₐ> = ────────⋅N
            1 + ε
"""
mu_alpha_external = (epsilon*palpha + 1)/(epsilon + 1)*N

"""Mean number of gene products at the OFF state.
            ε⋅pₐ
    <nᵦ> = ─────⋅N
           1 + ε
"""
mu_beta_external = (epsilon*palpha)/(epsilon + 1)*N

"""Fano factor of marginal distribution.
            N - μ
    F = 1 + ─────
            1 + ε
"""
fano_external = 1 + N*(1 - palpha)/(1 + epsilon)

"""Variance of marginal distribution.
             μ⋅(N - μ)
    σ² = μ + ─────────
               1 + ε
"""
sigma__2_external = N*palpha*fano_external
