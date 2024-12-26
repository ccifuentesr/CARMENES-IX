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

input_file = 'cif03.v06'
df = pd.read_csv('Data/'+input_file+'.csv', sep=",", header=0)

filter_rho = df['Bool_rho'] == True

Fe_H_phot = df['Fe_H_phot']
Fe_H_lit = df['Fe_H']

y = Fe_H_lit
x = Fe_H_phot

OC = y-x
xp = np.linspace(-2, 2, 1000)
yp_up = [+0.25 for i in xp]
yp_lo = [-0.25 for i in xp]

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

ylabel = r'[Fe/H]$_{\rm lit}$ [dex]'
xlabel = r'[Fe/H]$_{\rm phot}$ [dex]'

fig, ax = plt.subplots(2, 1, sharex='col',\
    gridspec_kw={'hspace': 0.1, 'wspace': 0.4, 'height_ratios': [10, 2], 'hspace': 0.03}, figsize=figsize)

ax[0].scatter(x, y, facecolors='none', edgecolors='blue', s=pointsize, zorder=1)
#
ax[0].plot(xp, xp, '-', c='grey', lw=2, zorder=0)
ax[0].plot(xp, xp+0.25, '--', c='grey', lw=1.5, zorder=0)
ax[0].plot(xp, xp-0.25, '--', c='grey', lw=1.5, zorder=0)
#
ax[1].scatter(x, OC, facecolors='none', edgecolors='blue', s=pointsize, zorder=1)
#
ax[1].axhline(y=0, c='grey', lw=2, ls='-', zorder=0)
ax[1].plot(xp, yp_up, c='grey', lw=2, ls='--', zorder=0)
ax[1].plot(xp, yp_lo, c='grey', lw=2, ls='--', zorder=0)

# =============================================================================
# CUSTOM
# =============================================================================

ax[0].tick_params(axis='x', labelsize=tickssize, direction='in',
                  top=True, labeltop=False, which='both', labelbottom=False)
ax[0].tick_params(axis='y', labelsize=tickssize, direction='in',
                  right=True, labelright=False, which='both')
ax[0].tick_params('both', direction = 'in', length = 10, width = 1.5, which = 'major')
ax[0].tick_params('both', direction = 'in', length = 5, width = 0.5, which = 'minor')
ax[0].xaxis.set_tick_params(which='minor', bottom=True, top=True)
ax[0].minorticks_on()
ax[0].set_xlabel('', size=labelsize)
ax[0].set_ylabel(ylabel, size=labelsize)
# ax[0].axes.xaxis.set_ticklabels([])
ax[0].set_xlim([-1.2, 1.2])
ax[0].set_ylim([-1.2, 1.2])
#
ax[1].set_ylabel('O-C\n [dex]', size=labelsize*0.8)
ax[1].tick_params(axis='x', labelsize=tickssize, direction='in',
                  top=True, labeltop=False, which='both', labelbottom=True)
ax[1].tick_params(axis='y', labelsize=tickssize*0.8, direction='in',
                  right=True, labelright=False, which='both')
ax[1].tick_params('both', direction = 'in', length = 10, width = 1.5, which = 'major')
ax[1].tick_params('both', direction = 'in', length = 5, width = 0.5, which = 'minor')
ax[1].set_xlabel(xlabel, size=labelsize)
ax[1].set_ylim(-1.2, 1.2)

#
plt.savefig('Output/'+'plot_metallicities.pdf', dpi=900, bbox_inches='tight')
plt.show()
