#!/usr/bin/env python3
# vim: fileencoding=utf-8

# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, you can obtain one at https://mozilla.org/MPL/2.0/
#
# Copyright 2020 Alexandre Ferreira Ramos - AMPhyBio Laboratory
#
# Project:  github.com/amphybio/stochastic-gene-expression
# Version:  1.0
# Created:  29-07-2019
# Authors:  Leonardo R. Gama <leonardo.gama@usp.br>

import os
import sys
from glob import glob
from runpy import run_path
from subprocess import check_call

old_files = glob('fig_*.eps') + glob('fig_*.png')
if old_files:
    message = ">> These old files were found and may cause errors:"
    prompt = ">> Should they be removed before running? [y/N] "
    print(message, " ".join(old_files), prompt, sep="\n\n", end="")
    if input().strip() == 'y':
        for file in old_files:
            os.remove(file)
        print(">> Files removed.")

for script in sorted(glob('fig_*.py')):
    print("\n>> Running", script)
    run_path(script)

try:
    if sys.argv[1] == '-eps':
        print("\n>> Converting and appending images...")
        check_call('bash -x compile_eps.sh 2>&1 | grep -v "^+ for"', shell=True)
except IndexError:
    pass
