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
# Created:  14-11-2019
# Authors:  Leonardo R. Gama <leonardo.gama@usp.br>

"""Graphs of entropy and mutual information vs. mean number of gene products."""

from __init__ import *
from amphybio.stochastic_gene.entropy import *
from amphybio.stochastic_gene.steady_state import mu_alpha_external

options = default_options()
params = config[name(__file__)]

xmax = 50
options['xmin'] = 0
options['xmax'] = xmax
options['ymin'] = 0

# options['set'].append('key')  # modified latter

x = utils.plot_points(0.01, xmax, 25, logspace=True)

const_y = np.array([H_constitutive(N) for N in x])
const_style = options['with'] + 'linecolor black'

plots = [
    {'epsilon': params['epsilon'][0], 'key': 'left top'    },
    {'epsilon': params['epsilon'][1], 'key': 'right bottom'},
    {'epsilon': params['epsilon'][2], 'key': 'right bottom'},
]

for index, plot in enumerate(plots, start=1):
    options['ymax'] = 7.5

    # H vs. <n>
    EPSILON = " ({} = {:.2f})" if plot['epsilon'] < 1 else " ({} = {:.0f})"
    options['title'] = PLOT_LETTER.format('A') + title['H'] + EPSILON.format(label['epsilon'], plot['epsilon'])
    options['xlabel'] = label['mu'] + "_{ }"  # spacing hack
    options['ylabel'] = label['H']
    # options['set'][-1] = 'key ' + plot['key']

    ys = {}
    for palpha in params['palpha']:
        Ns = [mu/palpha for mu in x]  # mu == palpha*N
        ys[palpha] = [H_external(plot['epsilon'], palpha, N, method='maple-async', backup_method='sympy-parallel') for N in Ns]

    curves = []
    for palpha in params['palpha']:
        y = np.array([val.get() if isinstance(val, pool.AsyncResult) else val for val in ys[palpha]])
        curves.append((x, y)) #, {'legend': palpha}))
    curves.append((x, const_y, {'with': const_style}))

    if logger.level >= logging.INFO: print(flush=True)
    output(curves, options, '{}_{}_H', name(__file__), index)

    # Hₒₙ vs. <nₐ>
    options['title'] = PLOT_LETTER.format('B') + title['H_ON']
    options['xlabel'] = label['mu_alpha']
    options['ylabel'] = label['H_ON']
    # options['set'][-1] = 'key right bottom'
    options['set'].append('key opaque font ",18"')
    options['set'].append('key title "{} value"'.format(label['palpha']))
    options['set'].append('key right bottom')
    if logger.level == logging.DEBUG:
        options['title'] += " ({} = {:.1f})".format(label['epsilon'], plot['epsilon'])

    ys = {}
    for palpha in params['palpha']:
        expr = mu_alpha_external.subs({'epsilon': plot['epsilon'], 'p_a': palpha})
        Ns = [sym.solve(sym.Eq(expr, mu_alpha), sym.Symbol('N'))[0] for mu_alpha in x]
        ys[palpha] = [H_ON_external(plot['epsilon'], palpha, N, method='maple-async', backup_method='sympy-parallel') for N in Ns]

    curves = []
    for palpha in params['palpha']:
        y = np.array([val.get() if isinstance(val, pool.AsyncResult) else val for val in ys[palpha]])
        curves.append((x, y, {'legend': palpha}))
    curves.append((x, const_y, {'with': const_style, 'legend': 'const.'}))

    if logger.level >= logging.INFO: print(flush=True)
    output(curves, options, '{}_{}_H_ON', name(__file__), index)
    options['set'][-3:] = []  # remove key

    # I vs. <n>
    options['ymax'] = 1
    options['title'] = PLOT_LETTER.format('C') + title['I']
    options['xlabel'] = label['mu'] + "_{ }"  # spacing hack
    options['ylabel'] = label['I']
    options['ymin'] = 0.001
    options['ymax'] = 1.2
    options['set'].append('ytics ("" 0.005, "0.01" 0.01, "0.05" 0.05, "0.1" 0.1, "0.5" 0.5, "1" 1)')
    options['set'].append('logscale y')
    # options['set'][-1] = 'key right top'
    if logger.level == logging.DEBUG:
        options['title'] += " ({} = {:.1f})".format(label['epsilon'], plot['epsilon'])

    ys = {}
    for palpha in params['palpha']:
        Ns = [mu/palpha for mu in x]  # mu == palpha*N
        ys[palpha] = [I_external(plot['epsilon'], palpha, N, method='maple-async', backup_method='sympy-parallel') for N in Ns]

    curves = []
    for palpha in params['palpha']:
        y = np.array([val.get() if isinstance(val, pool.AsyncResult) else val for val in ys[palpha]])
        curves.append((x, y)) #, {'legend': palpha}))

    if logger.level >= logging.INFO: print(flush=True)
    output(curves, options, '{}_{}_I', name(__file__), index)
    options['set'][-2:] = []  # remove logscale y and ytics
