import numpy as np
import pandas as pd
from astroquery.gaia import Gaia

# =============================================================================
# DATA
# =============================================================================

input_file = 'cif03.v06_subset'
output_file = input_file + '_out.csv'
df = pd.read_csv('Data/' + input_file + '.csv', sep=",", header=0, dtype={'GaiaDR3_id': str})

# Handle non-finite values
df['GaiaDR3_id'] = pd.to_numeric(df['GaiaDR3_id'], errors='coerce')
df = df.dropna(subset=['GaiaDR3_id'])
df['GaiaDR3_id'] = df['GaiaDR3_id'].astype(int)

print(df['GaiaDR3_id'])

def query_gaia(source_id):
    """Queries in Gaia Archive using ADQL syntax
    
    Args:
        source_id (int): unique Gaia EDR3 identifier
    
    Returns:
        r (astropy.table.Table): Query results
    """
    Gaia.login(user='ccifuent', password='botzat-zYpzog51')
    query = f"""
    SELECT source_id, type_best_classification 
    FROM gaiadr3.vari_cepheid 
    WHERE source_id = {source_id}
    """
    job = Gaia.launch_job(query)
    results = job.get_results()
    return results

# =============================================================================
# ACTION
# =============================================================================

gaia_ids, var_types = [], []
for idx, row in df.iterrows():
    source_id = row['GaiaDR3_id']
    try:
        result = query_gaia(source_id)
        gaia_ids.append(source_id)
        if len(result) > 0:
            var_types.append(result['type_best_classification'][0])
        else:
            var_types.append(np.nan)
    except Exception as e:
        print(f"Error querying Gaia for source ID {source_id}: {e}")
        gaia_ids.append(np.nan)
        var_types.append(np.nan)

# =============================================================================
# WRITE OUT
# =============================================================================

output = 'y'

if output == 'y':
    df_output = pd.DataFrame({'gaia_id': gaia_ids, 'var_type': var_types})
    df_output.to_csv('Output/' + output_file, sep=',', encoding='utf-8', index=False)
else:
    pass
