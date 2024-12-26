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

input_file = 'cif03.v14'
df = pd.read_csv('Data/'+input_file+'.csv', sep=",", header=1)

filter_rho = (df['Bool_rho'] == True)
filter_in = (np.abs(1 - df['WDS_sep2']/df['rho01']) < 0.105)
filter_out = (np.abs(1 - df['WDS_sep2']/df['rho01']) >= 0.105)
filter_true_gaia = (df['rho01'] > 0.196) # Highest resolution from Gaia in the sample

rho_WDS = df[filter_rho & filter_in & filter_true_gaia]['WDS_sep2']
rho_Gaia = df[filter_rho & filter_in & filter_true_gaia]['rho01']
rho_WDS_out = df[filter_rho & filter_out]['WDS_sep2']
rho_Gaia_out = df[filter_rho & filter_out]['rho01']

OC = rho_WDS - rho_Gaia
OC_out = rho_WDS_out - rho_Gaia_out
xp = np.linspace(0, 10000, 20000)
yp_up = [i*0.10 for i in xp]
yp_lo = [-i*0.10 for i in xp]
 
# breakpoint()

print(rho_WDS_out/rho_Gaia_out)

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

ylabel = r'$\rho_{\rm WDS}$'f'[arcsec]'
xlabel = r'$\rho_{Gaia}$'f'[arcsec]'

fig, ax = plt.subplots(2, 1, sharex='col',\
    gridspec_kw={'hspace': 0.1, 'wspace': 0.4, 'height_ratios': [10, 2], 'hspace': 0.03}, figsize=figsize)

#
ax[0].scatter(rho_Gaia, rho_WDS, facecolors='none', edgecolors='blue', s=pointsize, zorder=1)
ax[0].scatter(rho_Gaia_out, rho_WDS_out, facecolors='none', edgecolors='magenta', s=pointsize, zorder=2)
#
ax[0].plot(xp, xp, '-', c='grey', lw=2, zorder=0)
ax[0].plot(xp, xp*1.10, '--', c='grey', lw=1.5, zorder=0)
ax[0].plot(xp, xp*0.90, '--', c='grey', lw=1.5, zorder=0)
#
ax[1].scatter(rho_Gaia, OC, facecolors='none', edgecolors='blue', s=pointsize, zorder=1)
ax[1].scatter(rho_Gaia_out, OC_out, facecolors='none', edgecolors='magenta', s=pointsize, zorder=1)
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
ax[0].axes.xaxis.set_ticklabels([])
ax[0].set_xlim([0.1, 6000])
ax[0].set_ylim([0.1, 6000])
ax[0].set_xscale('log')
ax[0].set_yscale('log')
#
ax[1].set_ylabel('O-C\n [arcsec]', size=labelsize*0.8)
ax[1].tick_params(axis='x', labelsize=tickssize, direction='in',
                  top=True, labeltop=False, which='both', labelbottom=True)
ax[1].tick_params(axis='y', labelsize=tickssize*0.8, direction='in',
                  right=True, labelright=False, which='both')
ax[1].tick_params('both', direction = 'in', length = 10, width = 1.5, which = 'major')
ax[1].tick_params('both', direction = 'in', length = 5, width = 0.5, which = 'minor')
ax[1].set_xlabel(xlabel, size=labelsize)
ax[1].set_xlim([0.1, 6000])
ax[1].set_ylim(-1, 1)

#
plt.savefig('Output/'+'plot_rho_WDS.pdf', dpi=900, bbox_inches='tight')
plt.show()
