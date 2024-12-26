import numpy as np
import pandas as pd
from astroquery.simbad import Simbad
import re
from astropy.coordinates import SkyCoord
import astropy.units as u

# =============================================================================
# DATA
# =============================================================================

input_file = 'cif03.v12'
output_file = f'{input_file}_out.csv'
df = pd.read_csv(f'Data/{input_file}.csv', sep=",", header=1)

ID_star = df['ID_star']
Name = df['Name']
ra = df['ra']
dec = df['dec']

# =============================================================================
# FUNCTION
# =============================================================================

def query_simbad(object, field):
    customSimbad = Simbad()
    customSimbad.add_votable_fields(field)
    return customSimbad.query_object(object)

def extract_identifier(name, string):
    search = Simbad.query_objectids(name)
    for id in search:
        if string in str(id):
            return re.sub(f'{string} ', '', re.sub(' +', ' ', id[0]))
    return ''

def get_default_name(ra, dec):
    coords = SkyCoord(ra=ra, dec=dec, unit=(u.deg, u.deg), frame='icrs')
    customSimbad = Simbad()
    result = customSimbad.query_region(coords, radius=3 * u.arcsec)
    return result['MAIN_ID'][0] if result is not None and len(result) > 0 else None

# =============================================================================
# ACTION
# =============================================================================

obtain_names = True
obtain_default_name = False

if obtain_names:
    identifiers = ['V*']#'GJ' , 'G ', 'HD', 'NAME', 'LP ', 'Wolf', 'Ross', 'BD', 'LTT', 'RX', 'PM ', 'PM J', 'Gaia DR2', 'Gaia EDR3', 'WDS', '2MASS J', 'Karmn', '2MASS', 'WISEA']
    name_lists = {id: [] for id in identifiers}

    for name in Name:
        for id in identifiers:
            try:
                name_lists[id].append(extract_identifier(name, id))
            except (EOFError, TypeError, ConnectionError):
                name_lists[id].append('-')

if obtain_default_name:
    Name_default = []
    for i in range(len(df)):
        try:
            Name_default.append(get_default_name(ra[i], dec[i]))
        except (EOFError, TypeError, ConnectionError):
            Name_default.append('')

# =============================================================================
# WRITE OUT
# =============================================================================

if obtain_names == True:
    df_params = {**{'ID_star': ID_star, 'Name': Name}, **name_lists}
    df_output = pd.DataFrame(data=df_params)
    df_output.to_csv(f'Output/{output_file}', sep=',', encoding='utf-8')

if obtain_default_name == True:
    df_params = {'ID_star': ID_star, 'Name_default': Name_default}
    df_output = pd.DataFrame(data=df_params)
    df_output.to_csv(f'Output/{output_file}', sep=',', encoding='utf-8')
