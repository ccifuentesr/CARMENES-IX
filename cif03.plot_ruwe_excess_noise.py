import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
from pathlib import Path

# Configure the font and plot settings
fpath = Path(matplotlib.get_data_path(), "/Users/ccifuentesr/Library/Fonts/MinionPro-MediumCapt.otf")
plt.rcParams["font.family"] = "Georgia"
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['axes.linewidth'] = 1.5

# =============================================================================
# DATA
# =============================================================================

input_file = 'cif03.v30'
df = pd.read_csv('Data/' + input_file + '.csv', sep=",", header=0)

# Filter data
filter_single = df['Type'] == 'Single'
filter_single_candidate = df['Type'] == 'Candidate'
filter_multiple = df['Type'] == 'Multiple'
filter_multiple_candidate = df['Type'] == 'Multiple+'

aen_single = df[filter_single]['astrometric_excess_noise']
aen_single_candidate = df[filter_single_candidate]['astrometric_excess_noise']
aen_multiple = df[filter_multiple]['astrometric_excess_noise']
aen_multiple_candidate = df[filter_multiple_candidate]['astrometric_excess_noise']
ruwe_single = df[filter_single]['ruwe']
ruwe_single_candidate = df[filter_single_candidate]['ruwe']
ruwe_multiple = df[filter_multiple]['ruwe']
ruwe_multiple_candidate = df[filter_multiple_candidate]['ruwe']

# Define variables
y_single = aen_single
x_single = ruwe_single
y_single_candidate = aen_single_candidate
x_single_candidate = ruwe_single_candidate
y_multiple = aen_multiple
x_multiple = ruwe_multiple
y_multiple_candidate = aen_multiple_candidate
x_multiple_candidate = ruwe_multiple_candidate

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

xlabel = r'ruwe'
ylabel = r'astrometric_excess_noise [mas]'

fig, ax = plt.subplots(figsize=figsize)

# Scatter plot
ax.scatter(x_single, y_single, facecolors='none', s=pointsize*0.7, linewidth=1.5, edgecolors='hotpink', zorder=10)
ax.scatter(x_single_candidate, y_single_candidate, s=pointsize*0.7, linewidth=1.5,  facecolors='none', edgecolors='hotpink', zorder=10)
ax.scatter(x_multiple, y_multiple, s=pointsize, linewidth=1.5, facecolors='none', edgecolors='blue')
ax.scatter(x_multiple_candidate, y_multiple_candidate, s=pointsize,linewidth=1.5, facecolors='none', edgecolors='blue')
ax.axvline(2, color='red', linestyle='dashed', linewidth=linewidth, zorder=20)
ax.axhline(0.2, color='red', linestyle='dashed', linewidth=linewidth, zorder=20)

# Adjust the axes
ax.set_xscale('log')
ax.set_yscale('log')
ax.xaxis.set_major_formatter(ScalarFormatter())
ax.set_xticks([1, 2, 10, 20])
ax.set_ylim(0.005, 12)
ax.set_xlabel(xlabel, size=labelsize)
ax.set_ylabel(ylabel, size=labelsize)
ax.tick_params(axis='x', labelsize=tickssize, direction='in',
               top=True, which='both', labelbottom=True)
ax.tick_params(axis='y', labelsize=tickssize, direction='in',
               right=True, which='both')
ax.tick_params('both', direction = 'in', length = 10, width = 1.5, which = 'major')
ax.tick_params('both', direction = 'in', length = 5, width = 0.5, which = 'minor')

# =============================================================================
# SAVE
# =============================================================================

plt.savefig('Output/' + 'plot_ruwe_excess.pdf', dpi=900, bbox_inches='tight')
plt.show()
