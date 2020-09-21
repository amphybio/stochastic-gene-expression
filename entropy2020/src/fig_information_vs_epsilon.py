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
options['ylabel'] = label['I']

options['ymax'] = 1.
options['set'].append('logscale x')
options['set'].append('format x ""') #"10^{%T}"')

options['set'].append('key title "{} value"'.format(label['palpha']))
options['set'].append('key opaque font ",18"')
options['set'].append('key right top')


plots = [

    {'mu': params['mu'], 'xmin': 0.01, 'xmax': 100}

]

curves = []

for index, plot in enumerate(plots, start=1):

    # I vs. ε

    options['title'] = PLOT_LETTER.format('B') + title['I_E']
    if logger.level == logging.DEBUG:
        options['title'] += " ({} = {})".format(label['mu'], plot['mu'])

    partition = 5.
    epsilon1 = np.logspace(math.log10(plot['xmin']), math.log10(partition + 1), 15)
    epsilon2 = np.logspace(math.log10(partition), math.log10(plot['xmax'] + 1), 8)
    # epsilon3 = np.logspace(math.log10(100), math.log10(1000), 10)[:-1]
    epsilon = np.concatenate((epsilon1, epsilon2)) #, epsilon3))

    y = {}
    for palpha in params['palpha']:
        y[palpha] = [I_external(e, palpha, plot['mu']/palpha, method='maple-async', backup_method='sympy-parallel') for e in epsilon]
    for palpha in params['palpha']:
        curves.append((
            epsilon,
            np.array([val.get() if isinstance(val, pool.AsyncResult) else val for val in y[palpha]]),
            {'legend': palpha}
        ))

    marks = config['fig_entropy_vs_mu']['epsilon']
    options['set'].append('xtics add (\
            "{/Helvetica-Bold: 0.01}" 0.01,\
            "0.1" 0.1,\
            "1" 1,\
            "{/Helvetica-Bold: 2}" 2,\
            "{/Helvetica-Bold: 10}" 10,\
            "100" 100,\
            )')

    output(curves, options, '{}_{}_mu_{}', name(__file__), index, plot['mu'])
