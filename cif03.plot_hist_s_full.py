import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import rc
from pathlib import Path
from scipy.optimize import curve_fit
from scipy.stats import ks_2samp
from sklearn.metrics import r2_score, mean_squared_error
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.ticker import ScalarFormatter

fpath = Path(matplotlib.get_data_path(), "/Users/ccifuentesr/Library/Fonts/MinionPro-MediumCapt.otf")
plt.rcParams["font.family"] = "Georgia"
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['axes.linewidth'] = 1.5

# =============================================================================
# DATA
# =============================================================================

input_file = 'cif03.v19'
df = pd.read_csv('Data/' + input_file + '.csv', sep=",", header=1)
df_nea = pd.read_csv('Data/nea_exoplanets.csv', sep=",", header=0)

s01_au = df['s01']
s02_au = df['s02']
concatenated_s = np.sort(pd.concat([s01_au, s02_au]))

x = concatenated_s
x_cum = np.arange(1, len(x) + 1)

# =============================================================================
# FITTING
# =============================================================================

def opik(s, n0):
    s0 = np.min(s)
    return (n0 * (np.log10(s) - np.log10(s0)))

range1 = (x >= 7) & (x <= 80)
range2 = (x >= 80) & (x <= 3000)
x1, y1 = x[range1], x_cum[range1]
x2, y2 = x[range2], x_cum[range2]

# Fit both ranges
popt1, pcov1 = curve_fit(opik, x1, y1, maxfev=200000)
popt2, pcov2 = curve_fit(opik, x2, y2, maxfev=200000)

x_fit1 = np.linspace(7, 80, len(y1))
x_fit2 = np.linspace(80, 3000, len(y2))
y_fit1 = opik(x_fit1, *popt1)
y_fit2 = opik(x_fit2, *popt2)

# Coefficient R2
residuals1 = y1 - y_fit1
ss_res1 = np.sum(residuals1**2)
ss_tot1 = np.sum((y1 - np.mean(y1))**2)
r_squared1 = 1 - (ss_res1 / ss_tot1)
#
residuals2 = y2 - y_fit2
ss_res2 = np.sum(residuals2**2)
ss_tot2 = np.sum((y2 - np.mean(y2))**2)
r_squared2 = 1 - (ss_res2 / ss_tot2)

print(f'R2: {r_squared1}')
print(f'R2: {r_squared2}')
print(f'Range 1: a = {popt1[0]}')
print(f'Range 2: a = {popt2[0]}')

# =============================================================================
# PLOT
# =============================================================================

figsize = (12, 20)
pointsize = 80
linewidth = 2
tickssize = 22
labelsize = 22
legendsize = 18
cblabsize = 18
empty = '$\u25EF$'

xlabel = r'$s$'f' [au]'
ylabel1 = r'N($s$)' #f'Cumulative number of pairs'
ylabel2 = f'Number of pairs'

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=figsize)

# First Plot
ax1.scatter(x, x_cum, facecolors='none', edgecolors='blue', s=pointsize)
ax1.plot(x_fit1, opik(x_fit1, popt1[0]), 'r--', lw=linewidth * 1.4)
ax1.plot(x_fit2, opik(x_fit2, popt2[0]), 'r--', lw=linewidth * 1.4)
ax1.axvline(7, color='grey', linestyle='dashed', linewidth=linewidth, zorder=0)
ax1.axvline(80, color='grey', linestyle='dashed', linewidth=linewidth, zorder=0)
ax1.axvline(3000, color='grey', linestyle='dashed', linewidth=linewidth, zorder=0)
ax1.text(9, 0, "Range 1", fontsize=labelsize, color='black')
ax1.text(200, 0, "Range 2", fontsize=labelsize, color='black')

# Inline plot
ax_inset = inset_axes(ax1, width="30%", height="30%", loc='upper left', borderpad=2)
zoom_x = (x >= 50) & (x <= 140)
zoom_y = (x_cum[zoom_x] >= 3800) & (x_cum[zoom_x] <= 490) 
ax_inset.scatter(x[zoom_x], x_cum[zoom_x], facecolors='blue', s=pointsize*3)
ax_inset.plot(x_fit1, opik(x_fit1, popt1[0]), 'r--', lw=linewidth * 1.6)
ax_inset.plot(x_fit2, opik(x_fit2, popt2[0]), 'r--', lw=linewidth * 1.6)
ax_inset.axvline(80, color='grey', linestyle='dashed', linewidth=linewidth, zorder=0)
#
ax_inset.set_xlim(50, 140)
ax_inset.set_ylim(420, 490)
ax_inset.tick_params(axis='x', labelsize=tickssize*0.8, direction='in', top=True, labeltop=False, which='both', labelbottom=True)
ax_inset.tick_params(axis='y', labelsize=tickssize*0.8, direction='in', right=True, labelright=False, which='both')
ax_inset.tick_params('both', direction='in', length=10, width=1.5, which='major')
ax_inset.tick_params('both', direction='in', length=5, width=0.5, which='minor')
ax_inset.minorticks_on()
ax_inset.set_yticks([])  
ax_inset.set_xscale('log')
ax_inset.xaxis.set_major_formatter(ScalarFormatter())
ax_inset.set_xticks([60, 80, 100])

from matplotlib.patches import Rectangle
ax_inset.add_patch(Rectangle((0, 0), 1, 1, transform=ax_inset.transAxes, color='white', alpha=0.5))

# Customization for first plot
ax1.tick_params(axis='x', labelsize=tickssize, direction='in', top=True, labeltop=False, which='both', labelbottom=True)
ax1.tick_params(axis='y', labelsize=tickssize, direction='in', right=True, labelright=False, which='both')
ax1.tick_params('both', direction='in', length=10, width=1.5, which='major')
ax1.tick_params('both', direction='in', length=5, width=0.5, which='minor')
ax1.minorticks_on()
ax1.set_xlabel(xlabel, size=labelsize)
ax1.set_ylabel(ylabel1, size=labelsize)
ax1.set_xlim([0.005, 210000])
ax1.set_xscale('log')

# Second Plot
pl_orbsmax = df_nea['pl_orbsmax']

hist, bins = np.histogram(x, bins=np.logspace(np.log10(0.005), np.log10(100000), 20))
# hist = hist / np.max(hist)
ax2.hist(bins[:-1], bins, weights=hist, histtype='step', color='blue', linewidth=linewidth, zorder=3)
# hist_nea, bins_nea = np.histogram(pl_orbsmax, bins=np.logspace(np.log10(pl_orbsmax.min()), np.log10(pl_orbsmax.max()), 40))
# hist_nea = hist_nea / np.max(hist_nea) * np.max(hist)
# ax2.hist(bins_nea[:-1], bins_nea, weights=hist_nea, histtype='step', color='DodgerBlue', linewidth=linewidth, alpha=0.5)
ax2.axvline(7, color='grey', linestyle='dashed', linewidth=linewidth, zorder=0)
ax2.axvline(80, color='grey', linestyle='dashed', linewidth=linewidth, zorder=0)
ax2.axvline(3000, color='grey', linestyle='dashed', linewidth=linewidth, zorder=0)
# ax2.axvline(4.3, color='salmon', linestyle='dashed', linewidth=linewidth*1.2, zorder=0)

# Customization for second plot
ax2.tick_params(axis='x', labelsize=tickssize, direction='in', top=True, labeltop=False, which='both', labelbottom=True)
ax2.tick_params(axis='y', labelsize=tickssize, direction='in', right=True, labelright=False, which='both')
ax2.tick_params('both', direction='in', length=10, width=1.5, which='major')
ax2.tick_params('both', direction='in', length=5, width=0.5, which='minor')
ax2.minorticks_on()
ax2.set_xlabel(xlabel, size=labelsize)
ax2.set_ylabel(ylabel2, size=labelsize)
ax2.set_xscale('log')
ax2.set_xlim([0.005, 210000])

# =============================================================================
# OUT
# =============================================================================

plt.savefig('Output/hist_s_combined.pdf', dpi=900, bbox_inches='tight')
plt.show()
