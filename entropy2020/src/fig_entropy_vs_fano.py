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

# authors: Leonardo Gama (10-12-2019)

"""Graphs of entropy and mutual information vs. Fano factor of distribution."""

from __init__ import *
from entropy import *
from steady_state import fano_external
from collections import namedtuple

fano = sym.lambdify(['epsilon', 'p_a', 'N'], fano_external, math)


options = default_options()
params = config[name(__file__)]

options['xlabel'] = label['fano']
options['xmax'] = 4e3
options['set'].append('logscale x')
options['set'].append('format x "10^{%T}"')

epsilon_style = options['with'] + 'linetype {}'
palpha_style = options['with'].replace('2', '1.5') + 'dashtype "-"'

epsilons = utils.plot_points(0.2, 8, 16, logspace=True).tolist()
epsilons = sorted(epsilons + params['epsilon'])
palphas = utils.plot_points(0.9999, 0.015625, 64).tolist()
Ns = [params['mu']/palpha for palpha in palphas]


### Mutual information ###

options['title'] = PLOT_LETTER.format('B') + title['I_vs_fano']
options['ylabel'] = label['I']

options['set'].append('key title "{} value"'.format(label['epsilon']))
options['set'].append('key at graph 0.99,0.98')

curves = []

# pₐ curves
plots = {p: {} for p in params['palpha']}

for palpha, plot in plots.items():
    N = params['mu']/palpha
    plot['x'] = [fano(epsilon, palpha, N) for epsilon in epsilons]
    plot['y'] = [I_external(epsilon, palpha, N, method='maple-async', backup_method='sympy-parallel') for epsilon in epsilons]

for palpha, plot in plots.items():
    x = np.array(plot['x'])
    y = np.array([val.get() if isinstance(val, pool.AsyncResult) else val for val in plot['y']])
    curves.append((x, y, {'with': palpha_style}))

# ε curves
plots = {
    params['epsilon'][0]: {"extra_palphas": []},#0.0175, 0.01, 0.005]},
    params['epsilon'][1]: {"extra_palphas": []},#0.0175, 0.01, 0.005]},
    params['epsilon'][2]: {"extra_palphas": []},#0.0175, 0.01, 0.0075, 0.005, 0.0025]},
    params['epsilon'][3]: {"extra_palphas": []},#0.0175, 0.01, 0.0075, 0.005, 0.0025]},
    params['epsilon'][4]: {"extra_palphas": []},#0.0175, 0.01, 0.0075, 0.005, 0.0025, 0.00175, 0.001, 0.00095, 0.0009, 0.00085]},
}

for epsilon, plot in plots.items():
    extra_Ns = [params['mu']/palpha for palpha in plot['extra_palphas']]
    zip_params = list(zip(palphas + plot['extra_palphas'], Ns + extra_Ns))
    plot['x'] = [fano(epsilon, palpha, N) for palpha, N in zip_params]
    plot['y'] = [I_external(epsilon, palpha, N, method='maple-async', backup_method='sympy-parallel') for palpha, N in zip_params]

for color, (epsilon, plot) in enumerate(plots.items(), start=1):
    x = np.array(plot['x'])
    y = np.array([val.get() if isinstance(val, pool.AsyncResult) else val for val in plot['y']])
    curves.append((x, y, {'with': epsilon_style.format(color), 'legend': "{:.2f}".format(epsilon)}))

if logger.level >= logging.INFO: print(flush=True)
output(curves, options, '{}_I', name(__file__))


### Entropy ###

options['title'] = PLOT_LETTER.format('A') + title['H_vs_fano'] + " ({} = {})".format(label['mu'], params['mu'])
options['ylabel'] = label['H']

options['ymax'] = 7.5

options['set'][-2] = 'key title "{} value"'.format(label['palpha'])
options['set'][-1] = 'key left bottom'

curves = []

# pₐ curves
plots = {p: {} for p in params['palpha']}

for palpha, plot in plots.items():
    N = params['mu']/palpha
    plot['x'] = [fano(epsilon, palpha, N) for epsilon in epsilons]
    plot['y'] = [H_external(epsilon, palpha, N, method='maple-async', backup_method='sympy-parallel') for epsilon in epsilons]

for palpha, plot in plots.items():
    x = np.array(plot['x'])
    y = np.array([val.get() if isinstance(val, pool.AsyncResult) else val for val in plot['y']])
    curves.append((x, y, {'with': palpha_style.format(color), 'legend': '{:.2f}'.format(palpha)}))


# ε curves
plots = {
    params['epsilon'][0]: {"extra_palphas": []},#0.0175, 0.01, 0.005]},
    params['epsilon'][1]: {"extra_palphas": []},#0.0175, 0.01, 0.005]},
    params['epsilon'][2]: {"extra_palphas": []},#0.0175, 0.01, 0.0075, 0.005, 0.0025]},
    params['epsilon'][3]: {"extra_palphas": []},#0.0175, 0.01, 0.0075, 0.005, 0.0025]},
    params['epsilon'][4]: {"extra_palphas": []},#0.0175, 0.01, 0.0075, 0.005, 0.0025, 0.00175, 0.001, 0.00095, 0.0009, 0.00085]},
}

for epsilon, plot in plots.items():
    extra_Ns = [params['mu']/palpha for palpha in plot['extra_palphas']]
    extra_params = list(zip(palphas + plot['extra_palphas'], Ns + extra_Ns))
    plot['x'] = [fano(epsilon, palpha, N) for palpha, N in extra_params]
    plot['y'] = [H_external(epsilon, palpha, N, method='maple-async', backup_method='sympy-parallel') for palpha, N in extra_params]

for color, (epsilon, plot) in enumerate(plots.items(), start=1):
    x = np.array(plot['x'])
    y = np.array([val.get() if isinstance(val, pool.AsyncResult) else val for val in plot['y']])
    curves.append((x, y, {'with': epsilon_style.format(color)}))

## Custom ε key
# key_box = 'object rectangle from graph {x0},{y0} to graph {x1},{y1} fillstyle noborder'
# key_title = 'label at graph {x},graph {y} center "{label} value"'
# key_legend = 'label at graph {x},graph {y} right "{legend:.2f}"'
# key_line = 'arrow from graph {x},graph {y} rto graph 0.1,graph 0 linetype {color} linewidth 2 nohead'

# Point = namedtuple('Point', ['x', 'y'])
# orig = Point(0.98, 0.97)
# options['set'].append(key_box.format(x0=orig.x - 0.24, y0=orig.y - 0.07 - 0.05*len(params['epsilon']), x1=orig.x, y1=orig.y))
# options['set'].append(key_title.format(x=orig.x - 0.12, y=orig.y - 0.035, label=label['epsilon']))
# for i, epsilon in enumerate(params['epsilon'], start=1):
    # options['set'].append(key_legend.format(x=orig.x - 0.14, y=orig.y - 0.04 - i*0.05, legend=epsilon))
    # options['set'].append(key_line.format(x=orig.x - 0.115, y=orig.y - 0.04 - i*0.05, color=i))

# Constitutive gene point
const_style = 'points pointtype {} pointsize 2.5 linecolor "black"'.format(POINT_DIAMOND)
const_line = 'arrow from {x0},{y} to {x1},{y} dashtype "-" linecolor "black" nohead'

const_y = H_constitutive(params['mu'])
curves.append((np.array([1]), np.array([const_y]), {'legend': 'const.', 'with': const_style}))
options['set'].append(const_line.format(x0=1, x1=options['xmax'], y=const_y))

if logger.level >= logging.INFO: print(flush=True)
output(curves, options, "{}_H", name(__file__))
