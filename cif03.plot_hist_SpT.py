import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib import rc
from matplotlib.ticker import FormatStrFormatter, ScalarFormatter
from pathlib import Path

fpath = Path(matplotlib.get_data_path(), "/Users/ccifuentesr/Library/Fonts/MinionPro-MediumCapt.otf")
plt.rcParams["font.family"] = "Georgia"
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['axes.linewidth'] = 1.5

# =============================================================================
# DATA
# =============================================================================

filename = 'Data/cif03.v30.csv'

df = pd.read_csv(filename, sep=",", header=0)

filter_Karmn = df['Karmn'].notnull()

SpTnum = df[filter_Karmn]['SpTnum']
SpTypes = ['M0', 'M1', 'M2', 'M3', 'M4', 'M5', 'M6', 'M7', 'M8', 'M9']

star_ranges = [SpTnum[(SpTnum >= n) & (SpTnum <= n+0.5)] for n in range(70, 80, 1)]
sample_size = [len(star) for star in star_ranges]

bins = np.arange(len(sample_size))  # Create bins based on the length of the data
bin_centers = bins + 0.5

# =============================================================================
# PLOT
# =============================================================================

figsize = (12, 10)
pointsize = 60
linewidth = 2
elinewidth = 2
tickssize = 22
labelsize = 26
legendsize = 18
cblabsize = 18
bins_size = np.arange(70.5, 78.5, 1)
fig, ax = plt.subplots(figsize=figsize)
cm = plt.cm.get_cmap('magma_r')

xlabel = f'Spectral type'
ylabel = f'Number of stars'

colors = cm(np.linspace(0.1, 1, len(sample_size)))

bars = ax.bar(bins, sample_size, color=colors)

for bar, color in zip(bars, colors):
    bar.set_facecolor(color)

# =============================================================================
# CUSTOM
# =============================================================================

ax.set_xlabel(xlabel, size=labelsize, font=fpath)
ax.set_ylabel(ylabel, size=labelsize, font=fpath)
ax.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
ax.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
ax.set_xticklabels(SpTypes)

ax.tick_params(axis='x', labelsize=tickssize+2, direction='in',
               top=False, labeltop=False, which='both')
ax.tick_params(axis='y', labelsize=tickssize, direction='in',
               right=False, labelright=False, which='both')
ax.tick_params('both', direction = 'in', length = 10, width = 1.5, which = 'major')
ax.tick_params('both', direction = 'in', length = 5, width = 0.5, which = 'minor')
ax.minorticks_on()
ax.xaxis.set_tick_params(which='minor', bottom=False, top=False)
ax.set_xticks(bins)

plt.setp(ax.xaxis.get_majorticklabels(), rotation=0, ha="center")
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
#
plt.ylim(0.5, 1e3)
ax.set_yscale('log')
#

# =============================================================================
# OUT
# =============================================================================

output_file = 'hist_SpT'
plt.savefig('Output/'+output_file+'.pdf', bbox_inches='tight')
plt.show()