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

input_file = 'cif03.v06_variability_v01'
df = pd.read_csv('Data/'+input_file+'.csv', sep=",", header=0)
ID_star = df['ID_star']
Name = df['Name']
filter_clean_G = df['GrVFlag'] == 0
filter_clean_B = df['BPrVFlag'] == 0
filter_clean_R = df['RPrVFlag'] == 0
filter_single = df['Type'] == 'Single'
filter_single_ast = df['Type'] == 'Single*'
filter_multiple = (df['Type'] == 'Multiple') | (df['Type'] == 'Multiple*')

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

xlabel = f'$\sigma$$_R$$_P$' f'[mag]'
ylabel = f'$\sigma_G$' f'[mag]'

ID_star_unique = sorted(set(ID_star))

def compute_std(df, ID_star):
    std_results = {}
    for star_id in ID_star_unique:
        if star_id not in std_results:
            Gmag_values_single = df[(df['ID_star'] == star_id) & filter_clean_G & filter_single]['Gmag']
            BPmag_values_single = df[(df['ID_star'] == star_id) & filter_clean_B & filter_single]['BPmag']
            RPmag_values_single = df[(df['ID_star'] == star_id) & filter_clean_R & filter_single]['RPmag']
            Gmag_values_single_ast = df[(df['ID_star'] == star_id) & filter_clean_G & filter_single_ast]['Gmag']
            BPmag_values_single_ast = df[(df['ID_star'] == star_id) & filter_clean_B & filter_single_ast]['BPmag']
            RPmag_values_single_ast = df[(df['ID_star'] == star_id) & filter_clean_R & filter_single_ast]['RPmag']
            Gmag_values_multiple = df[(df['ID_star'] == star_id) & filter_clean_G & filter_multiple]['Gmag']
            BPmag_values_multiple = df[(df['ID_star'] == star_id) & filter_clean_B & filter_multiple]['BPmag']
            RPmag_values_multiple = df[(df['ID_star'] == star_id) & filter_clean_R & filter_multiple]['RPmag']
            if True:# len(Gmag_values_single) > 1 and len(BPmag_values_single) > 1 and len(RPmag_values_single) > 1 and \
            # len(Gmag_values_single_ast) > 1 and len(BPmag_values_single_ast) > 1 and len(RPmag_values_single_ast) > 1 and \
            # len(Gmag_values_multiple) > 1 and len(BPmag_values_multiple) > 1 and len(RPmag_values_multiple) > 1:
                std_results[star_id] = {
                    'Gmag_single_std': np.std(Gmag_values_single),
                    'BPmag_single_std': np.std(BPmag_values_single),
                    'RPmag_single_std': np.std(RPmag_values_single),
                    'Gmag_single_ast_std': np.std(Gmag_values_single_ast),
                    'BPmag_single_ast_std': np.std(BPmag_values_single_ast),
                    'RPmag_single_ast_std': np.std(RPmag_values_single_ast),
                    'Gmag_multiple_std': np.std(Gmag_values_multiple),
                    'BPmag_multiple_std': np.std(BPmag_values_multiple),
                    'RPmag_multiple_std': np.std(RPmag_values_multiple)
                }
    return std_results

results = compute_std(df, ID_star)

Gmag_single_std_values = [results[star_id]['Gmag_single_std'] for star_id in results]
BPmag_single_std_values = [results[star_id]['BPmag_single_std'] for star_id in results]
RPmag_single_std_values = [results[star_id]['RPmag_single_std'] for star_id in results]
Gmag_single_ast_std_values = [results[star_id]['Gmag_single_ast_std'] for star_id in results]
BPmag_single_ast_std_values = [results[star_id]['BPmag_single_ast_std'] for star_id in results]
RPmag_single_ast_std_values = [results[star_id]['RPmag_single_ast_std'] for star_id in results]
Gmag_multiple_std_values = [results[star_id]['Gmag_multiple_std'] for star_id in results]
BPmag_multiple_std_values = [results[star_id]['BPmag_multiple_std'] for star_id in results]
RPmag_multiple_std_values = [results[star_id]['RPmag_multiple_std'] for star_id in results]

fig, ax = plt.subplots(figsize=figsize)
#
ax.scatter(RPmag_single_std_values, Gmag_single_std_values, facecolors='none', edgecolors='lightblue', s=pointsize, linewidths=2)
ax.scatter(RPmag_single_ast_std_values, Gmag_single_ast_std_values, facecolors='none', edgecolors='tomato', s=pointsize, linewidths=2)
ax.scatter(RPmag_multiple_std_values, Gmag_multiple_std_values, facecolors='none', edgecolors='royalblue', s=pointsize, linewidths=1.5)
ax.plot([0, 1], [0,1], linestyle='--', color='grey')
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
ax.set_xlim(0.002, .5)
ax.set_ylim(0.002, .5)
ax.set_xscale('log')
ax.set_yscale('log')
plt.savefig(f'Output/plot_light_curve_sigma_RP.pdf', dpi=900, bbox_inches='tight')
plt.show()