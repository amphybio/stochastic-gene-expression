# vim: fileencoding=utf-8

# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, you can obtain one at https://mozilla.org/MPL/2.0/
#
# Copyright 2020 Alexandre Ferreira Ramos - AMPhyBio Laboratory
#
# Project:  github.com/amphybio/stochastic-gene-expression
# Version:  1.1
# Created:  14-11-2019
# Authors:  Leonardo R. Gama <leonardo.gama@usp.br>

"""
Common code used by figure generating scripts.

Reads and processes settings (parameters) from JSON file, defines constants and
functions, and does the imports common to all figure scripts.  The usage below
also applies to all figure scripts.

Usage:
    {fig_x}        -- files in EPS format
    {fig_x} -png   -- files in PNG format
    {fig_x} -term  -- plot in "term" window (either 'x', 'qt' or 'wxt')
"""

import json
import logging
import os
import sys
from os import path


with open('config.json') as file:
    config = json.load(file)


# Plot options setup
DEFAULT_OPTIONS = config['default_options']

def default_options():
    """
    Return a 2nd level deep copy of DEFAULT_OPTIONS.
    """
    return {k: v.copy() if 'copy' in dir(v) else v for k, v in DEFAULT_OPTIONS.items()}


# Terminal setup
term_options = config['term_options']

available_terms = {
    'eps': 'postscript eps color {}'.format(term_options['eps']),
    'gp': 'gp',
    'png': 'pngcairo {}'.format(term_options['png']),
    'qt': 'qt {}'.format(term_options['qt']),
    'wxt': 'wxt',
    'x': 'x11'
}

try:
    term = sys.argv[1].lstrip('-')
except IndexError:
    term = 'eps'

if term not in available_terms:
    print(__doc__.format(fig_x=sys.argv[0]).strip(), file=sys.stderr)
    exit(2)

DEFAULT_OPTIONS['terminal'] = available_terms[term]
if term in ('wxt', 'x'):
    DEFAULT_OPTIONS['wait'] = True


def output(curves, options, filename=None, *args):
    """
    Generate output file from curves.

    :curves: list of curves for the gnuplotlib plot function
    :options: options dictionary for the gnuplotlib plot function
    :filename: output file name as a (format) string *without* the extension
    :args: substitutions to be made in 'filename'
    """
    if logger.level >= logging.INFO: print(flush=True)
    if term in ('eps', 'gp', 'png'):
        filename = filename.format(*args)
        options['output'] = filename + '.' + term
        data_file = filename + '.data'

        if logger.level == logging.DEBUG:
            # Save plotted data to file.
            curves = [curve if len(curve) == 3 else curve + ({},) for curve in curves]
            header = "\t".join("x\t" + str(opt.get('legend', 'unlabeled')) for x, y, opt in curves)
            matrix = [row.astype(np.float64) for x, y, opt in curves for row in (x, y)]
            maxlen = max(len(row) for row in matrix)
            matrix = [np.concatenate([row, np.full(maxlen - len(row), np.nan)]) for row in matrix]
            matrix = np.array(matrix).transpose()
            logger.info("Generating data file %s", options['output'])
            np.savetxt(data_file, matrix, fmt="%.3e", delimiter="\t", header=header)

        logger.info("Generating file %s", options['output'])
    gp.gnuplotlib(**options).plot(*curves)


def name(filename):
    return path.splitext(path.basename(filename))[0]


# Text
title = config['title']
label = config['label']
PLOT_LETTER = label['PLOT_LETTER']

# Gnuplot point types
POINT_SQUARE = 5
POINT_BALL = 7
POINT_TRIANGLE = 9
POINT_INVTRIANGLE = 11
POINT_DIAMOND = 13

# Modules used by figures
import math
import multiprocessing
import numpy as np
import sympy as sym
import gnuplotlib as gp
from multiprocessing import pool

import utils

logging.basicConfig(format='[%(levelname)s] %(message)s', level=logging.INFO)
logger = logging.getLogger()

if logger.level == logging.DEBUG:
    DEFAULT_OPTIONS['with'] = 'linespoints '
