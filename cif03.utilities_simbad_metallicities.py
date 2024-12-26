import numpy as np
import pandas as pd
import astroquery
from astroquery.gaia import Gaia
from astroquery.simbad import Simbad
import re
import csv
astroquery.simbad.conf.server

# =============================================================================
# DATA
# =============================================================================

input_file = 'cif03.v05'
output_file = input_file + '_out.csv'
df = pd.read_csv('Data/'+input_file+'.csv', sep=",", header=0)

Name = df['Name']
ID_star = df['ID_star']

# =============================================================================
# FUNCTION
# =============================================================================

def Simbad_name(object):
    """
    Gives the default identifier in Simbad for an object.
    """
    customSimbad = Simbad()
    customSimbad.add_votable_fields('fe_h')
    result = customSimbad.query_object(object)
    return result['Fe_H_Fe_H'][0], result['Fe_H_bibcode'][0]

# =============================================================================
# ACTION
# =============================================================================

Fe_H, Fe_H_bib = [], []

for i in range(len(df)):
    try:
        Fe_H.append(Simbad_name(Name[i])[0])
        Fe_H_bib.append(Simbad_name(Name[i])[1])
    except EOFError:
        Fe_H.append('')
        Fe_H_bib.append('')
    except TypeError:
        Fe_H.append('')
        Fe_H_bib.append('')
    except ConnectionError:
        Fe_H.append('')
        Fe_H_bib.append('')

# =============================================================================
# WRITE OUT
# =============================================================================

df_out = pd.DataFrame(data={'ID_star': ID_star, 'Fe_H': Fe_H, 'Fe_H_bib': Fe_H_bib})
df_out.to_csv('Output/'+output_file, sep=',', encoding='utf-8')
