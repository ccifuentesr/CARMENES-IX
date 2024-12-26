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

input_file = 'cif03.v20'
df = pd.read_csv('Data/'+input_file+'.csv', sep=",", header=1)

d_pc = df['d_pc']
s_au = df['s01']
MG_mag = df['MG_mag']
dp = np.linspace(0, 1e7, 10000)
rho_values = [0.001, 0.01, 0.1, 1, 10, 100, 1000]

# =============================================================================
# PLOT
# =============================================================================

figsize = (12, 10)
pointsize = 60
linewidth = 2
elinewidth = 2
tickssize = 22
labelsize = 24
legendsize = 18
cblabsize = 18

xlabel = r'$d$'f' [pc]'
ylabel = r'$s$'f' [au]'
lines = []

fig, ax = plt.subplots(figsize=figsize)
ax.scatter(d_pc, s_au, s=pointsize, linewidth=1.5, facecolors='none', edgecolors='blue',zorder=10)

for rho in rho_values:
    s = rho * dp
    ax.plot(dp, s, linestyle='--', color='grey', label=str(rho)+'"', zorder=0)
    # ax.plot(dp, 0.196*dp, linestyle='--', color='salmon', label=str(rho)+'"', zorder=0)
    ax.text(2, rho*1.8, str(rho)+'"', fontsize=labelsize*0.8, color='grey', rotation=12, backgroundcolor='white')
    # ax.text(1.2, 0.22, '0.196"', fontsize=labelsize*0.8, color='salmon', rotation=10, backgroundcolor='white')
#
ax.tick_params(axis='x', labelsize=tickssize, direction='in',
                  top=True, labeltop=False, which='both', labelbottom=True)
ax.tick_params(axis='y', labelsize=tickssize, direction='in',
                  right=True, labelright=False, which='both')
ax.tick_params('both', direction = 'in', length = 10, width = 1.5, which = 'major')
ax.tick_params('both', direction = 'in', length = 5, width = 0.5, which = 'minor')
ax.xaxis.set_tick_params(which='minor', bottom=True, top=True)
ax.yaxis.set_tick_params(which='minor', left=True, right=True)  
ax.minorticks_on()
ax.set_xlabel(xlabel, size=labelsize)
ax.set_ylabel(ylabel, size=labelsize)
ax.set_ylim([2E-3, np.max(s_au)*1.2])
ax.set_xlim([1, 100])
ax.set_xscale('log')
ax.set_yscale('log')
ax.xaxis.set_major_formatter(ScalarFormatter())
ax.set_xticks([1, 10, 100])
#
plt.savefig('Output/'+'plot_s_d.pdf', dpi=900, bbox_inches='tight')
plt.show()