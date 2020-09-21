#!/usr/local/bin/maple -qs
# vim: fileencoding=utf-8

# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, you can obtain one at https://mozilla.org/MPL/2.0/
#
# Copyright 2020 Alexandre Ferreira Ramos - AMPhyBio Laboratory
#
# Project:  github.com/amphybio/stochastic-gene-expression
# Version:  1.0
# Created:  12-12-2019
# Authors:  Leonardo R. Gama <leonardo.gama@usp.br>


# Usage: entropy_external.mpl -cp:=func,epsilon,palpha,N,precision,k

# Convergence limit 'k' estimation for "infinite" sum:
# As the CDF for this distribution is bounded from above by the poissonian distribution at the
# higher tail, it is guarateed to be close to unity with any desired number of decimal places.
# An empirical estimation shows this proximity is of roughly two decimal places for each multiple
# of standard deviations greater than 2: 1 - CDF(N + (2 + i)*σ) < 10^(-2*i)

phi_n_external := N^n*pochhammer(epsilon*palpha, n)*KummerM(epsilon*palpha + n, epsilon + n, -N)/(factorial(n)*pochhammer(epsilon, n)):
alpha_n_external := palpha*N^n*pochhammer(epsilon*palpha + 1, n)*KummerM(epsilon*palpha + n + 1, 1 + epsilon + n, -N)/(factorial(n)*pochhammer(1 + epsilon, n)):
beta_n_external := (1 - palpha)*N^n*pochhammer(epsilon*palpha, n)*KummerM(epsilon*palpha + n, 1 + epsilon + n, -N)/(factorial(n)*pochhammer(1 + epsilon, n)):

H_external := -sum(phi_n_external*log[2](phi_n_external), n=0..k):
H_ON_external := -sum(alpha_n_external*log[2](alpha_n_external/palpha)/palpha, n=0..k):
H_OFF_external := -sum(beta_n_external*log[2](beta_n_external/(1 - palpha))/(1 - palpha), n=0..k):

try
    func := p[1]:
    epsilon := p[2]:
    palpha := p[3]:
    N := p[4]:
    Digits := p[5]:

    if numelems([p]) >= 6 then
        k := p[6]:
    else
        # ~6 digits of precision (see discussion above)
        precision_estimate := 6:
        k := N + ceil((2 + precision_estimate/2)*sqrt(N)):
    end if:

    res := evalf(func):
    if type(res, 'numeric') then
        printf(cat("%.", Digits, "e\n"), res):
    else
        printf("nan\n"):
    end if:

catch:
    stderr := fopen("/dev/stderr", WRITE):
    fprintf(stderr, cat("Maple error: ", StringTools[FormatMessage](lastexception[2], lastexception[3..-1]), "\n")):
    fclose(stderr):
    `quit`(1):
end try:
