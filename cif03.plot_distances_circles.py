import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc
from pathlib import Path
from matplotlib.colors import Normalize
from matplotlib.patches import Circle
from matplotlib.colors import BoundaryNorm

fpath = Path(matplotlib.get_data_path(), "/Users/ccifuentesr/Library/Fonts/MinionPro-MediumCapt.otf")
plt.rcParams["font.family"] = "Georgia"
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['axes.linewidth'] = 1.5

# =============================================================================
# DATA
# =============================================================================

input_file = 'cif03.v13'
df = pd.read_csv('Data/'+input_file+'.csv', sep=",", header=1)

filter_karmn = (df['Karmn'].notnull()) & (df['d_pc'] < 33.4)

ra = df[filter_karmn]['ra']  
dec = df[filter_karmn]['dec']  
d_pc = df[filter_karmn]['d_pc']  
d_pc_norm = d_pc/d_pc.max()*33.4
MG_mag = df[filter_karmn]['MG_mag']
SpTnum = df[filter_karmn]['SpTnum']

theta = np.linspace(0, 2*np.pi, 100)
r = np.linspace(0, 35, 100)
R, Theta = np.meshgrid(r, theta)
color_data = R

# =============================================================================
# PLOT
# =============================================================================

plot_cmap = True
plot_stars = True
plot_circles = True

figsize = (10, 10)
pointsize = 20
tickssize = 16
labelsize = 20

fig, ax = plt.subplots(figsize=figsize, subplot_kw={'projection': 'polar'})

if plot_cmap == True:
    c = ax.pcolormesh(Theta, R, color_data, shading='auto', cmap='magma', vmin=0, vmax=35)
    ax.tick_params(axis='x', labelsize=labelsize*0.8, pad=10)
    ax.tick_params(axis='y', labelsize=labelsize*0.8, pad=10)
    ax.grid(linestyle='')
    ax.set_rmax(35)
    
if plot_circles == True:
    radii = [10, 20, 30]
    for radius in radii:
        ax.plot(theta, np.full_like(theta, radius), linestyle='dashed', color='white', lw=2)
        ax.text(0, radius, str(radius) + " pc", ha='center', va='center', fontsize=14, color='black',
                bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.5'))

if plot_stars == True:
    x = ra
    y = d_pc
    ax.scatter(0, 0, c='orange', marker='o', s=280, zorder=99)
    ax.scatter(x, y, c='none', marker='o', zorder=1, s=pointsize, edgecolors='white')
    ax.set_rlabel_position(20)
    
ax.set_yticklabels([])

plt.savefig('Output/'+'plot_distances.pdf', dpi=900, bbox_inches='tight')
plt.show()