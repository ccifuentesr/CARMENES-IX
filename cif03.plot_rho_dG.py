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

input_file = 'cif03.v16'
df = pd.read_csv('Data/'+input_file+'.csv', sep=",", header=1)
df_filtered = df[df['System'] != ''].copy()

df_WDS = df[(df['WDS_sep2'] < 999) & (df['System'] == '(AB)')]
rho_WDS = df_WDS['WDS_sep2']
dG_WDS = abs(df_WDS['WDS_mag1'] - df_WDS['WDS_mag2'])

results = []

# Group by 'ID_system' and calculate deltaG by comparing to the first item in each system
for system_id, group in df_filtered.groupby('ID_system'):
    # Skip 'ID_system' with value 0
    if system_id == 0:
        continue  # Skip this iteration

    if len(group) > 1:
        first_star_mag = group.iloc[0]['G_mag']
        first_star_system = group.iloc[0]['System']
        print(f"Processing system_id: {np.int64(system_id)} - {first_star_system}")

        # Convert first_star_system to string safely, handle cases where it's NaN or None
        if pd.notna(first_star_system):  # Only proceed if it's not NaN or None
            first_star_system = str(first_star_system)  # Ensure it's a string
        else:
            first_star_system = ''  # If it's NaN, treat it as an empty string

        # Initialize variables to store deltaG values
        deltaG_gaia2M = None
        deltaG_gaia = None

        # Compute deltaG for all subsequent stars in the same system
        for idx, row in group.iterrows():
            if idx != group.index[0]:  # Skip the first star
                deltaG = abs(row['G_mag'] - first_star_mag)

                # Assign deltaG to different variables based on the first star system type
                star_systems = ('A+B', '(AB)+C')#, 'Aab+B+C', 'Aab+B', '(AB)+C', 'A+(BC)', '(AB)+C+D', '(AB)+(CD)', 'A+Bab', 'Aab+Bab', 'A+Bab', 'AB+C+(DE)+F', 'A+B+C', 'Aab+(BC)', 'Aab+BC', 'Aab+B+C', '(AabBab)+Cab')
                if first_star_system in star_systems:
                    deltaG_gaia2M = deltaG
                elif first_star_system == 'AB':  # Gaia
                    deltaG_gaia = deltaG
                
                # Add the result to the list with appropriate deltaG
                results.append({
                    'ID_system': system_id,
                    'rho01': row['rho01'],  # Use the rho01 of the current star
                    'deltaG_gaia2M': deltaG_gaia2M,
                    'deltaG_gaia': deltaG_gaia
                })


df_results = pd.DataFrame(results)

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

xlabel = r'$\rho$'f' [arcsec]'
ylabel = r'$\Delta m$'f' [mag]'

fig, ax = plt.subplots(figsize=figsize)

ax.scatter(df_results['rho01'], df_results['deltaG_gaia2M'], s=pointsize*2, linewidth=1.5, facecolors='none', edgecolors='orange',zorder=10)
ax.scatter(df_results['rho01'], df_results['deltaG_gaia'], s=pointsize*2, linewidth=1.5, facecolors='none', edgecolors='dodgerblue',zorder=10)
ax.scatter(rho_WDS, dG_WDS, s=pointsize*2, linewidth=1.5, facecolors='none', edgecolors='blueviolet',zorder=0)

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
# ax.set_ylim([2E-3, np.max(s_au)*1.2])
# ax.set_xlim([1, 100])
ax.set_xscale('log')
plt.gca().invert_yaxis()

plt.savefig('Output/'+'plot_deltaG_rho.pdf', dpi=600, bbox_inches='tight')
plt.show()
