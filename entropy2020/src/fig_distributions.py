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

# authors: Leonardo Gama (09-12-2019)

"""Graph of probability distributions."""

import multiprocessing as mp
from __init__ import *
from sympy.abc import lamda, n
from steady_state import alpha_n_external, phi_n_external

poisson = lamda**n*sym.exp(-lamda)/sym.factorial(n)

options = default_options()
params = config[name(__file__)]

options['title'] = title['phi']
options['xlabel'] = label['N'] + "_{ }"  # spacing hack
options['ylabel'] = label['alpha_phi']

options['ymin'] = 1e-5
options['ymax'] = 1e-1
options['set'].append('logscale')
options['set'].append('format y "10^{%T}"')

options['set'].append('label at graph 0.76,graph 0.8 "{} = {}'.format(label['mu'], params['mu']))

# options['set'].append('key title "{} value"'.format(label['palpha']))
# options['set'].append('key opaque font ",18"')
# options['set'].append('key left bottom')

@utils.memoized
def point(args):
    """Parallel computation. Argument for 'func' can be 'phi' or 'alpha'."""
    a, n = args
    func = phi_n_external if a['func'] == 'phi' else alpha_n_external
    return func.subs({'epsilon': a['epsilon'], 'p_a': a['palpha'], 'N': a['N'], 'n': n}).evalf()

plots = [
    {'epsilon': params['epsilon'][0], 'xmax': 1200},
    {'epsilon': params['epsilon'][1], 'xmax': 1200},
    {'epsilon': params['epsilon'][2], 'xmax': 600},
]

with mp.Pool() as pool:
    for index, plot in enumerate(plots, start=1):
        options['title'] = PLOT_LETTER.format('D') + title['phi']
        options['xmax'] = plot['xmax']
        if logger.level == logging.DEBUG:
            options['title'] += " ({} = {:.1f}; {} = {})".format(label['epsilon'], plot['epsilon'], label['mu'], params['mu'])

        x = utils.plot_points(1, plot['xmax'], 200, logspace=True)
        x = np.array(sorted(set(math.ceil(n) for n in x)))

        curves = []
        for color, palpha in enumerate(params['palpha'], start=1):

            # φₙ distributions
            map_args = dict(func='phi', epsilon=plot['epsilon'], palpha=palpha, N=params['mu']/palpha)
            y = np.array(pool.map(point, [(map_args, n) for n in x]))
            curves.append((x, y, {'with': options['with'].replace('2', '1.5') + 'linecolor {}'.format(color)})) #, 'legend': palpha}))

            # αₙ distributions
            map_args['func'] = 'alpha'
            y = np.array(pool.map(point, [(map_args, n) for n in x]))
            curves.append((x, y, {'with': options['with'] + 'dashtype "-" linecolor {}'.format(color)}))

        dummy = np.array([float('nan')])
        # curves.append((dummy, dummy, {'legend': ' ', 'with': 'dots linecolor "white"'}))
        curves.append((dummy, dummy, {'legend': label['phi'], 'with': options['with'] + 'linecolor "gray50"'}))
        curves.append((dummy, dummy, {'legend': label['alpha'], 'with': options['with'] + 'linecolor "gray50" dashtype "-"'}))

        const_y = np.array([poisson.evalf(subs={'lamda': params['mu'], 'n': n}) for n in x])
        curves.append((x, const_y, {'with': options['with'] + 'linecolor "black"'})) #, 'legend': 'const.'}))

        output(curves, options, '{}_{}', name(__file__), index)
