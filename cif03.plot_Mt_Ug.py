import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.ticker import FormatStrFormatter
from matplotlib import rc
from pathlib import Path

fpath = Path(matplotlib.get_data_path(), "/Users/ccifuentesr/Library/Fonts/MinionPro-MediumCapt.otf")
plt.rcParams["font.family"] = "Georgia"
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['axes.linewidth'] = 1.5

# =============================================================================
# DATA
# =============================================================================

input_file = 'cif03.v20'
df = pd.read_csv('Data/'+input_file+'.csv', sep=",", header=1)
df = df[df['ID_system'] != 0]

Ug = df['Ug_J']
eUg = df['eUg_J']
Mt_Msol = df['Mtot']
eMt_Msol = df['eMtot']
s01 = df['s01']
id_system = df['ID_system']

plot = 'Ug'

# Group by 'ID_system' and take the first non-null values of 'Ug_J' and 'Mtot' for each system
grouped = df.groupby('ID_system').agg({
    'Ug_J': 'first',   # Get the first non-null value of Ug_J
    'Mtot': 'first',   # Get the first non-null value of Mtot
    'eUg_J': 'first',  # Get the first non-null value of eUg_J (if needed)
    'eMtot': 'first',  # Get the first non-null value of eMtot (if needed)
    's01': 'first',    # Get the first non-null value of s01 (if plotting s)
    'M_Msol': 'first',    # Get the first non-null value of s01 (if plotting s)
}).reset_index()

# Print the results as requested
if plot == 'Ug':
    print(grouped[['ID_system', 'Ug_J']])
elif plot == 's':
    print(grouped[['ID_system', 's01']])

# Prepare the data for plotting
x, ex = grouped['Mtot'], grouped['eMtot']

# Manually aded <1E33 Ug
Ug_added = np.array([1.742E+33, 2.332E+33, 9.221E+33, 6.964E+33, 5.198E+33, 3.550E+33, 2.615E+33, 6.480E+33])
Mtot_added = np.array([0.8196, 0.8027, 1.1843, 1.8191, 1.8191, 2.1694, 0.8646, 0.7386])
s01_added = np.array([60679.5365, 74092.5609, 62796.4854, 75945.4450, 75899.0502, 815.62684, 59342.6892, 36823.0677])

if plot == 'Ug':
    filename = 'plot_Mt_Ug.pdf'
    ylabel = r'$|U^{\ast}_g|$'f' [J]'
    y = grouped['Ug_J']
    z = np.log10(grouped['s01'])*35

if plot == 's':
    filename = 'plot_Mt_s.pdf'
    ylabel = r'$s$ [au]'
    y = grouped['s01']
    z = grouped['Ug_J']

def s_Ga(Mtot, t):
    s_pc = 1.212*Mtot/t
    s_au = s_pc * 206265
    return s_au

xp = np.linspace(0, 3.5, 100)

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

xlabel = r'$\mathcal{M}_T$ [$\mathcal{M}_\odot$]'

fig, ax = plt.subplots(figsize=figsize)

limit = 1E34 # J
#
if plot == 'Ug':
    ax.scatter(Mtot_added, Ug_added, facecolors='none', edgecolors='red', s=np.log10(s01_added)*35, zorder=0)
    ax.scatter(x[y > limit], y[y > limit], facecolors='none', edgecolors='blue', s=z[y > limit], zorder=0)
    ax.scatter(x[y < limit], y[y < limit], s=z[y < limit], facecolors='none', edgecolors='red', zorder=0)
    ax.axhline(limit, color='red', linestyle='dashed', linewidth=linewidth, zorder=9)

if plot == 's':
    ax.scatter(x[z > limit], y[z > limit], facecolors='none', edgecolors='blue', s=pointsize, zorder=1)
    # ax.scatter(x[z < limit], y[z < limit], facecolors='none', edgecolors='red', s=pointsize, zorder=0)
    ax.scatter(Mtot_added, s01_added, facecolors='none', edgecolors='red', s=pointsize, zorder=0)
    ax.axhline(1E4, color='red', linestyle='dashed', linewidth=linewidth, zorder=99) #1E4 au
    ax.plot(xp, s_Ga(xp, 1), ls='--', lw=linewidth, zorder=0, color='gray')
    ax.plot(xp, s_Ga(xp, 10), ls='--', lw=linewidth*1.5, zorder=0, color='gray')
    ax.plot(xp, s_Ga(xp, 100), ls='--', lw=linewidth*2, zorder=0, color='gray')
    ax.text(2.5, 6e5, '1 Ga', fontsize=20, color='k', rotation=-4)
    ax.text(2.5, 6e4, '10 Ga', fontsize=20, color='k', rotation=-4)
    ax.text(2.5, 6e3, '100 Ga', fontsize=20, color='k', rotation=-4)

# =============================================================================
# CUSTOM
# =============================================================================

ax.tick_params(axis='x', labelsize=tickssize, direction='in',
                  top=True, labeltop=False, which='both', labelbottom=True)
ax.tick_params(axis='y', labelsize=tickssize, direction='in',
                  right=False, labelright=False, which='both')
ax.tick_params('both', direction = 'in', length = 10, width = 1.5, which = 'major')
ax.tick_params('both', direction = 'in', length = 5, width = 0.5, which = 'minor')
# plt.xticks(xticks, labels=[str(tick) for tick in xticks])
# plt.yticks(yticks, labels=[str(tick) for tick in yticks])  
ax.xaxis.set_tick_params(which='minor', bottom=True, top=True)
ax.minorticks_on()
ax.set_xlabel(xlabel, size=labelsize)
ax.set_ylabel(ylabel, size=labelsize)
if plot == 'Ug':
    ax.set_xlim([0.1, 3.0])
    # ax.set_ylim([1e33, 10e37])
    ax.set_yscale('log')
    plt.gca().invert_yaxis()
if plot == 's':
    ax.set_xlim([0.1, 3.0])
    ax.set_ylim([0.5, 1E6])
    ax.set_yscale('log')
    plt.gca().invert_yaxis()
# ax.xaxis.set_major_formatter(FormatStrFormatter('%1.0f'))
# ax.xaxis.set_major_formatter(mticker.ScalarFormatter())
# ax.xaxis.get_major_formatter().set_scientific(False)
ax.set_xticks([1, 2, 3])
ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())

# =============================================================================
# SAVE
# =============================================================================

# plt.show()
plt.savefig('Output/'+filename, dpi=600, bbox_inches='tight')