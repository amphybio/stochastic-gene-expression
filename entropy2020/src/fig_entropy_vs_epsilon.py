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
#           Guiherme Giovanini <ggiovanini@usp.br>

"""Graphs of entropy vs. epsilon."""

from __init__ import *
from entropy import *


options = default_options()
params = config[name(__file__)]

options['xlabel'] = label['epsilon']
options['ylabel'] = label['H']

options['ymin'] = 0
options['ymax'] = 7.5
options['set'].append('logscale x')
options['set'].append('format x ""')

# options['set'].append('key title "{} value"'.format(label['palpha']))
options['set'].append('key left top')


plots = [
    {'mu': params['mu'], 'xmin': 0.01, 'xmax': 100}
]

curves = []
for index, plot in enumerate(plots, start=1):

    # H vs. ε
    options['title'] = PLOT_LETTER.format('A') + title['H_E'] + " ({} = {})".format(label['mu'], plot['mu'])

    partition = 5.
    epsilon1 = np.logspace(math.log10(plot['xmin']), math.log10(partition + 1), 15)
    epsilon2 = np.logspace(math.log10(partition), math.log10(plot['xmax'] + 1), 15)
    epsilon3 = np.logspace(math.log10(100), math.log10(1000), 10)[:-1]
    epsilon = np.concatenate((epsilon1, epsilon2, epsilon3))

    y = {}
    for palpha in params['palpha']:
        y[palpha] = [H_external(e, palpha, plot['mu']/palpha, method='maple-async', backup_method='sympy-parallel') for e in epsilon]
    for palpha in params['palpha']:
        curves.append((
            epsilon,
            np.array([val.get() if isinstance(val, pool.AsyncResult) else val for val in y[palpha]]),
            #{'legend': palpha}
        ))

    # Constitutive reference
    const_style = 'dashtype "-" linecolor "black"'
    const_line = 'arrow from {x0},{y} to {x1},{y} nohead ' + const_style
    options['set'].append(const_line.format(x0=plot['xmin'], x1=plot['xmax'], y=H_constitutive(params['mu'])))
    dummy = np.array([np.nan])
    curves.append((dummy, dummy, {'with': options['with'] + 'dashtype 2 linecolor black', 'legend': 'const.'}))

    # Custom xtics with important epsilon values highlighted
    options['set'].append('xtics add (\
            "{/Helvetica-Bold: 0.01}" 0.01,\
            "0.1" 0.1,\
            "1" 1,\
            "{/Helvetica-Bold: 2}" 2,\
            "{/Helvetica-Bold: 10}" 10,\
            "100" 100,\
            )')

    output(curves, options, '{}_{}_mu_{}', name(__file__), index, plot['mu'])
