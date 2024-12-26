import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import rc
from pathlib import Path
from scipy.optimize import curve_fit

fpath = Path(matplotlib.get_data_path(), "/Users/ccifuentesr/Library/Fonts/MinionPro-MediumCapt.otf")
plt.rcParams["font.family"] = "Georgia"
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['axes.linewidth'] = 1.5

# =============================================================================
# DATA
# =============================================================================

input_file = 'cif03.v12'
output_file = 'hist_s'

df = pd.read_csv('Data/'+input_file+'.csv', sep=",", header=1)
df_nea = pd.read_csv('Data/nea_exoplanets.csv', sep=",", header=0)
df_dhi10 = pd.read_csv('Data/Dhital10.csv', sep=",", header=0) # Dhital et al. (2010)

s01_au = df['s01']
s02_au = df['s02']
concatenated_s = pd.concat([s01_au, s02_au])
pl_orbsmax = df_nea['pl_orbsmax']
s_dhi10_au = df_dhi10['s_au']

x1 = concatenated_s
x2 = pl_orbsmax

# =============================================================================
# PLOT
# =============================================================================

figsize = (15, 10)
linewidth = 3
elinewidth = 2
tickssize = 22
labelsize = 22
legendsize = 18

xlabel = r'$s$'f' [au]'
ylabel = f'Number of pairs (normalized)'

fig, ax = plt.subplots(figsize=figsize)


hist, bins = np.histogram(x1, bins=np.logspace(np.log10(0.005), np.log10(100000), 28))
hist = hist / np.max(hist)
plt.hist(bins[:-1], bins, weights=hist, histtype='step', color='blue', linewidth=linewidth, zorder=3)
hist, bins = np.histogram(x2, bins=np.logspace(np.log10(x2.min()), np.log10(x2.max()), 30))
hist = hist / np.max(hist)
plt.hist(bins[:-1], bins, weights=hist, histtype='step', color='DodgerBlue', linewidth=linewidth, alpha=0.5)

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
ax.set_xscale('log')

# =============================================================================
# OUT
# =============================================================================

plt.savefig('Output/'+output_file+'.pdf', dpi=900, bbox_inches='tight')
plt.show()
