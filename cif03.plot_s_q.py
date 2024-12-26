import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc
from pathlib import Path
from matplotlib.colors import Normalize
from matplotlib.cm import get_cmap
import matplotlib.cm as cm

fpath = Path(matplotlib.get_data_path(), "/Users/ccifuentesr/Library/Fonts/MinionPro-MediumCapt.otf")
plt.rcParams["font.family"] = "Georgia"
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['axes.linewidth'] = 1.5

# =============================================================================
# DATA
# =============================================================================

input_file = 'cif03.v11'
df = pd.read_csv('Data/'+input_file+'.csv', sep=",", header=0)
filter = (df['SpTnum'] >= 70) & (df['SpTnum'] < 80)

x = df[filter]['s01']
y = df[filter]['q']

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
empty = '$\u25EF$'

xlabel = r'$s$ {au}'
ylabel = r'$q$'

fig, ax = plt.subplots(figsize=figsize)
#
ax.scatter(x, y, c='blue', marker=empty, s=pointsize, zorder=1)

ax.tick_params(axis='x', labelsize=tickssize, direction='in',
                  top=True, labeltop=False, which='both', labelbottom=True)
ax.tick_params(axis='y', labelsize=tickssize, direction='in',
                  right=True, labelright=False, which='both')
ax.tick_params('both', direction = 'in', length = 10, width = 1.5, which = 'major')
ax.tick_params('both', direction = 'in', length = 5, width = 0.5, which = 'minor')
ax.xaxis.set_tick_params(which='minor', bottom=True, top=True)
ax.minorticks_on()
# ax.set_xlim([0.04, 0.95])
ax.set_xlabel(xlabel, size=labelsize)
ax.set_ylabel(ylabel, size=labelsize)
# ax.invert_xaxis()

# =============================================================================
# SAVE
# =============================================================================

plt.savefig('Output/'+'plot_M_s.pdf', dpi=900, bbox_inches='tight')
plt.show()
