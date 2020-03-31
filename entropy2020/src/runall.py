#!/usr/bin/env python3
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

# authors: Leonardo Gama (29/07/2019)

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
