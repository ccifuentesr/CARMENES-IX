#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: ccifuentesr
"""

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
from matplotlib import rc
from pathlib import Path

fpath = Path(matplotlib.get_data_path(), "/Users/ccifuentesr/Library/Fonts/MinionPro-MediumCapt.otf")
plt.rcParams["font.family"] = "Georgia"
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['axes.linewidth'] = 1.5

# =============================================================================
# DATA
# =============================================================================

input_file = 'cif03.v30'
df = pd.read_csv('Data/'+input_file+'.csv', sep=",", header=0)

filter_single = (df['Component'] == '-')
filter_multiple = (df['Component'] != '-')
filter_crit_1 = (df['ruwe'] >= 2)
filter_crit_2 = (df['rv_error'] >= 10)
filter_crit_3 = ((df['rv_chisq_pvalue'] < 0.01) & (df['rv_renormalised_gof'] > 4) & (df['rv_nb_transits'] >= 10))
filter_crit_4 = (df['rho01'] > 0.1)

#
RV_single = df[filter_single]['rv']
eRV_single = df[filter_single]['rv_error']
G_mag_single = df[filter_single]['G_mag']
#
RV_multiple = df[filter_multiple]['rv']
eRV_multiple = df[filter_multiple]['rv_error']
G_mag_multiple = df[filter_multiple]['G_mag']
#
RV_crit_23 = df[(filter_single | filter_crit_4) & (filter_crit_2 | filter_crit_3)]['rv']
eRV_crit_23 = df[(filter_single | filter_crit_4) & (filter_crit_2 | filter_crit_3)]['rv_error']
G_mag_crit_23 = df[(filter_single | filter_crit_4) & (filter_crit_2 | filter_crit_3)]['G_mag']

x_single = G_mag_single
y_single = eRV_single
x_crit_23 = G_mag_crit_23
y_crit_23 = eRV_crit_23
x_multiple = G_mag_multiple
y_multiple = eRV_multiple

# =============================================================================
# PLOT
# =============================================================================

figsize = (12, 10)
pointsize = 60
linewidth = 2
elinewidth = 2
tickssize = 22
labelsize = 22
legendsize = 18
cblabsize = 18

xlabel = r'$G$ 'f'[mag]'
ylabel = r'$\sigma V_r$ [km s$^{-1}$]'
zlabel = f'ruwe'

fig, ax = plt.subplots(figsize=figsize)

#
ax.scatter(x_multiple, y_multiple, facecolors='none', edgecolors='blue', s=pointsize, linewidth=1.5)
ax.scatter(x_single, y_single, facecolors='none', edgecolors='hotpink', s=pointsize*0.7, linewidth=1.5)
# ax.scatter(x_crit_23, y_crit_23, facecolors='none', edgecolors='orange', s=pointsize, linewidth=1.5)
plt.axhline(y=10, color='red', linestyle='--', lw=2.5, zorder=2)

# =============================================================================
# CUSTOM
# =============================================================================

ax.tick_params(axis='x', labelsize=tickssize, direction='in',
                  top=True, labeltop=False, which='both', labelbottom=True)
ax.tick_params(axis='y', labelsize=tickssize, direction='in',
                  right=True, labelright=False, which='both')
ax.tick_params('both', direction = 'in', length = 10, width = 1.5, which = 'major')
ax.tick_params('both', direction = 'in', length = 5, width = 0.5, which = 'minor')
ax.xaxis.set_tick_params(which='minor', bottom=True, top=True)
ax.minorticks_on()
ax.set_xlabel(xlabel, size=labelsize)
ax.set_ylabel(ylabel, size=labelsize)
# plt.gca().invert_yaxis()
# ax.set_xlim([0, 360])
ax.set_ylim([0.1, 100])
ax.set_yscale('log')

# =============================================================================
# OUTPUT
# =============================================================================

plt.savefig('Output/plot_G_rv.pdf', dpi=900, bbox_inches='tight')
plt.show()
