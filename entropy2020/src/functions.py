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

# authors: Leonardo Gama (10/06/2019)

"""
Special mathematical functions.
"""

from sympy import Function, hyper, log, prod


def log2(x):
    return log(x, 2)


#TODO: Make evaluation work.
class KummerM(hyper):
    """
    KummerM function with pretty printing capabilities, based on SymPy's hypergeometric function.
    """
    def __new__(cls, a, b, z):
        return super().__new__(cls, (a,), (b,), z)

#TODO: Remove monkey patch.
KummerM = lambda a, b, z: hyper((a,), (b,), z)


class pochhammer(Function):
    """
    Pochhammer symbol, raising factorial.

    pochhammer(x, n) = x⋅(x+1)⋅(x+2)⋅⋅⋅(x+n-1)
    """
    @classmethod
    def eval(cls, x, n):
        # As of SymPy version 1.4, Float.is_integer == Float.is_zero
        try:
            if float(n).is_integer():
                return prod(x + i for i in range(int(n)))
            else:
                raise NotImplementedError("'n' in pochhammer symbol must be an integer")
        except TypeError:
            pass
