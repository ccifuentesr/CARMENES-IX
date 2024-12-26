import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc
from pathlib import Path
from matplotlib.colors import Normalize
import matplotlib.cm as cm
from matplotlib.cm import get_cmap

fpath = Path(matplotlib.get_data_path(), "/Users/ccifuentesr/Library/Fonts/MinionPro-MediumCapt.otf")
plt.rcParams["font.family"] = "Georgia"
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['axes.linewidth'] = 1.5

# =============================================================================
# DATA
# =============================================================================

input_file = 'cif03.v20'
df = pd.read_csv('Data/'+input_file+'.csv', sep=",", header=1)
filter = (df['SpTnum'] >= 70) & (df['SpTnum'] < 80)

x = df[filter]['Mass_A']
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

xlabel = r'$\mathcal{M}_A$ [$\mathcal{M}_\odot$]'

ylabel = r'$q$' #r'$q$

fig, ax = plt.subplots(figsize=figsize)
#
ax.scatter(x, y, c='blue', marker=empty, s=pointsize, zorder=1)
ax.axhline(y=1, color='salmon', linestyle='--', linewidth=linewidth, zorder=1)

num_shades = 10

cmap = get_cmap('magma')
norm = Normalize(vmin=0, vmax=num_shades)

# Method 1
# gradient = np.linspace(0, 1, 256)
# gradient = np.vstack((gradient, gradient))
# ymin, ymax = 0, 1
# ax.imshow(gradient, aspect='auto', extent=[0.077, 0.685, ymin, ymax], cmap=cmap, alpha = 0.5)
# ax.axvspan(0, 10, ymin=ymin, ymax=ymax, facecolor='none')

# Method 2
for i in range(num_shades):
    x0 = 0.077 + i * (0.685 - 0.077) / num_shades
    x1 = 0.077 + (i + 1) * (0.685 - 0.077) / num_shades
    color = cmap(norm(i))
    alpha = 0.3 #* (1 - abs(2 * (i / num_shades - 0.5)))   # Adjust this value if needed
    ax.axvspan(x0, x1, color=color, alpha=alpha, zorder=0)

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
ax.set_xlim([0.04, 0.95])
ax.set_xlabel(xlabel, size=labelsize)
ax.set_ylabel(ylabel, size=labelsize)
ax.invert_xaxis()

# =============================================================================
# SAVE
# =============================================================================

plt.savefig('Output/'+'plot_M_q.pdf', dpi=900, bbox_inches='tight')
plt.show()
