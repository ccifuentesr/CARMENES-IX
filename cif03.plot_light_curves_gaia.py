import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
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

input_file = 'cif03.v06_variability'
df = pd.read_csv('Data/'+input_file+'.csv', sep=",", header=0)
ID_star = df['ID_star']

# =============================================================================
# PLOT
# =============================================================================

figsize = (12, 10)
pointsize = 100
linewidth = 2
elinewidth = 2
tickssize = 22
labelsize = 22
legendsize = 18
cblabsize = 18

xlabel = f'JD - 2455197.5'
ylabel = r'$G_\lambda$' f'[mag]'

for n in range(2650):
    fig, ax = plt.subplots(figsize=figsize)
    mask_in = (df['ID_star'] == n) & (df['BPrVFlag'] == 0) & (df['RPrVFlag'] == 0)  & (df['GrVFlag'] == 0)
    mask_out_BP = (df['ID_star'] == n) & (df['BPrVFlag'] == 1)
    mask_out_RP = (df['ID_star'] == n) & (df['RPrVFlag'] == 1)
    mask_out_G = (df['ID_star'] == n) & (df['GrVFlag'] == 1)
    #
    if df[mask_in]['RPmag'].empty:
        plt.close(fig)
        continue
    #
    ax.scatter(df[mask_in]['TimeG'], df[mask_in]['Gmag'], facecolors='none', edgecolors='green', s=pointsize)
    ax.scatter(df[mask_in]['TimeBP'], df[mask_in]['BPmag'], facecolors='none', edgecolors='blue', s=pointsize)
    ax.scatter(df[mask_in]['TimeRP'], df[mask_in]['RPmag'], facecolors='none', edgecolors='red', s=pointsize)
    ax.scatter(df[mask_out_RP]['TimeRP'], df[mask_out_RP]['RPmag'], s=pointsize, marker='x', color='red')
    ax.scatter(df[mask_out_BP]['TimeBP'], df[mask_out_BP]['BPmag'], s=pointsize, marker='x', color='blue')
    ax.scatter(df[mask_out_G]['TimeG'], df[mask_out_G]['Gmag'], s=pointsize, marker='x', color='green')
    #
    ax.tick_params(axis='x', labelsize=tickssize, direction='in', top=True, labeltop=False, which='both', labelbottom=True)
    ax.tick_params(axis='y', labelsize=tickssize, direction='in', right=True, labelright=False, which='both')
    ax.tick_params('both', direction = 'in', length = 10, width = 1.5, which = 'major')
    ax.tick_params('both', direction = 'in', length = 5, width = 0.5, which = 'minor')
    ax.xaxis.set_tick_params(which='minor', bottom=True, top=True)
    ax.minorticks_on()
    formatter = FuncFormatter(lambda y, pos: f'{y:.1f}')
    plt.gca().yaxis.set_major_formatter(formatter)
    ax.set_xlabel(xlabel, size=labelsize)
    ax.set_ylabel(ylabel, size=labelsize)
    ax.set_ylim([min(df[mask_in]['RPmag']) - .5, .5 + max(df[mask_in]['BPmag'])])
    ax.text(0.5, 0.95, df[mask_in]['Name'].iloc[0], transform=ax.transAxes, fontsize=labelsize, va='top', ha='center', bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.5'))
    ax.invert_yaxis()
    #
    plt.savefig(f'Output/plot_light_curve_{n}.pdf', dpi=900, bbox_inches='tight')
    plt.close()
