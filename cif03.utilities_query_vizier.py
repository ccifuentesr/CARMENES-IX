import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from astroquery.vizier import Vizier
import astropy.units as u
import matplotlib.pyplot as plt

# =============================================================================
# DATA
# =============================================================================

file = 'cif03.v06'
df = pd.read_csv('Data/'+file+'.csv', sep=",", header=0, nrows=2)

vizier = Vizier(columns=['g_transit_time', 'g_transit_mag'])
catalog_name = 'I/355/epphot'

def query_vizier(ra, dec, catalog_name):
    """Queries Vizier catalog <catalog_name> using coordinates
    
    Args:
        ra (float): Right Ascension
        dec (float): Declination
    
    Returns:
        result (astropy.table.Table): Query results
    """
    coord = SkyCoord(ra, dec, unit=(u.deg, u.deg), frame='icrs')
    result = vizier.query_region(coord, radius=1*u.arcsec, catalog=catalog_name)
    if result:
        return result[0]
    else:
        return None

# =============================================================================
# ACTION
# =============================================================================

print(query_vizier(49.442359174881400,+25.25015392509530, catalog_name).colnames)

results = []

for idx, row in df.iterrows():
    id_star = row['ID_star']
    ra = row['ra']
    dec = row['dec']
    result = query_vizier(1.3009212, 45.7858943, catalog_name)
    if result is not None:
        result['ID_star'] = id_star  # Add ID_star to the result table
        results.append(result)

# Concatenate all results into a single DataFrame
if results:
    all_results = pd.concat([res.to_pandas() for res in results], ignore_index=True)
else:
    all_results = pd.DataFrame()

# =============================================================================
# PLOTTING
# =============================================================================

# Ensure the catalog contains '_tab8_5' and '_tab8_9' columns
if '_tab8_5' in all_results.columns and '_tab8_9' in all_results.columns:
    plt.figure(figsize=(10, 6))
    plt.scatter(all_results['_tab8_5'], all_results['_tab8_9'], s=10, alpha=0.7)
    plt.xlabel('_tab8_5')
    plt.ylabel('_tab8_9')
    plt.title('_tab8_5 vs _tab8_9')
    plt.gca().invert_yaxis()  # Typically magnitudes are plotted with a reversed y-axis
    plt.grid(True)
    print(all_results['_tab8_9'])
    plt.show()
else:
    print("The required columns '_tab8_5' and '_tab8_9' are not present in the results.")