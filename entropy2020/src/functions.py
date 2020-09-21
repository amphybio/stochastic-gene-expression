# vim: fileencoding=utf-8

# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, you can obtain one at https://mozilla.org/MPL/2.0/
#
# Copyright 2020 Alexandre Ferreira Ramos - AMPhyBio Laboratory
#
# Project:  github.com/amphybio/stochastic-gene-expression
# Version:  1.0
# Created:  10-06-2019
# Authors:  Leonardo R. Gama <leonardo.gama@usp.br>

"""
Special mathematical functions for SymPy.
"""

from sympy import Function, hyper, log, prod


def log2(x):
    return log(x, 2)


#TODO: Make evaluation work and remove lambda monkey patch.
class KummerM(hyper):
    """
    KummerM function with pretty printing capabilities, based on SymPy's hypergeometric function.
    """
    def __new__(cls, a, b, z):
        return super().__new__(cls, (a,), (b,), z)

KummerM = lambda a, b, z: hyper((a,), (b,), z)


class pochhammer(Function):
    """
    Pochhammer symbol, raising factorial.

    pochhammer(x, n) = x⋅(x+1)⋅(x+2)⋅⋅⋅(x+n-1)
    """
    @classmethod
    def eval(cls, x, n):
        #NOTE: As of SymPy version 1.4, Float.is_integer == Float.is_zero
        try:
            if float(n).is_integer():
                return prod(x + i for i in range(int(n)))
            else:
                raise NotImplementedError("'n' in pochhammer symbol must be an integer")
        except TypeError:
            pass
