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
from matplotlib.ticker import ScalarFormatter

fpath = Path(matplotlib.get_data_path(), "/Users/ccifuentesr/Library/Fonts/MinionPro-MediumCapt.otf")
plt.rcParams["font.family"] = "Georgia"
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['axes.linewidth'] = 1.5

# =============================================================================
# DATA
# =============================================================================

input_file = 'cif03.v30'
df = pd.read_csv('Data/'+input_file+'.csv', sep=",", header=0)

filter_single = df['Component'] == '-'
filter_multiple = df['Component'] != '-'
filter_candidates = (df['ipd_gof_harmonic_amplitude'] > 0.1) & (df['ruwe'] > 1.4)

ipd_single = df[filter_single]['ipd_gof_harmonic_amplitude']
ipd_multiple = df[filter_multiple]['ipd_gof_harmonic_amplitude']
ruwe_single = df[filter_single]['ruwe']
ruwe_multiple = df[filter_multiple]['ruwe']

y_single = ipd_single
x_single = ruwe_single
y_multiple = ipd_multiple
x_multiple = ruwe_multiple

ipd_single_candidates = df[filter_single & filter_candidates]['ipd_gof_harmonic_amplitude']
ruwe_candidates = df[filter_single & filter_candidates]['ruwe']

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

ylabel = r'ipd_gof_harmonic_amplitude'
xlabel = r'ruwe'

fig, ax = plt.subplots(figsize=figsize)

ax.scatter(x_single, y_single, s=pointsize*0.7, linewidth=1.5, facecolors='none', edgecolors='hotpink', alpha=1)
# ax.scatter(ruwe_candidates, ipd_single_candidates, s=pointsize, linewidth=1.5,  facecolors='none', edgecolors='orange', zorder=1)
ax.scatter(x_multiple, y_multiple, s=pointsize, linewidth=1.5,  facecolors='none', edgecolors='blue', zorder=0)
plt.axhline(y=0.1, color='red', linestyle='--', lw=2.5, zorder=2)
plt.axvline(x=1.4, color='red', linestyle='--', lw=2.5, zorder=2)
#
ax.tick_params(axis='x', labelsize=tickssize, direction='in',
                  top=True, labeltop=False, which='both', labelbottom=True)
ax.tick_params(axis='y', labelsize=tickssize, direction='in',
                  right=True, labelright=False, which='both')
ax.tick_params('both', direction = 'in', length = 10, width = 1.5, which = 'major')
ax.tick_params('both', direction = 'in', length = 5, width = 0.5, which = 'minor')
ax.xaxis.set_tick_params(which='minor', bottom=True, top=True)
ax.minorticks_on()
ax.xaxis.set_major_formatter(ScalarFormatter())
ax.set_xticks([1, 2, 10, 20])
ax.set_xlabel(xlabel, size=labelsize)
ax.set_ylabel(ylabel, size=labelsize)
# plt.gca().invert_yaxis()
# ax.set_ylim([0, 51])
ax.set_xscale('log')
ax.xaxis.set_major_formatter(ScalarFormatter())
ax.set_xticks([1, 2, 10])
#
plt.savefig('Output/'+'plot_ruwe_ipd.pdf', dpi=900, bbox_inches='tight')
plt.show()
