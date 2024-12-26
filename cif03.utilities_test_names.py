import numpy as np
import pandas as pd
import astropy.coordinates as coord
import astropy.units as u
import re
import sys
#!{sys.executable} -m pip install termcolor
from termcolor import colored
from astroquery.simbad import Simbad
from astropy.coordinates import SkyCoord 
from astropy import units as u

input_file = 'cif03.csv'
df = pd.read_csv(input_file, sep=',', header=1)

names = df['Name']
names = names
count_ok = 0
count_ko = 0
failed_names = []

def simbad_name(object_name):
    global count_ok, count_ko
    # Attempt to query object identifiers
    result_table_id = Simbad.query_objectids(object_name)
    if result_table_id is not None:
        print(colored(f"\n{object_name}\n", 'green'), result_table_id)
        count_ok += 1
    else: 
        print(colored(f"\nFailed for {object_name}", 'red'))
        count_ko += 1
        failed_names.append(object_name)
    return print(colored(f"\nValid names: {count_ok}\n", 'green'), colored(f"Invalid names: {count_ko}\n", 'red'), colored(f"Failed names: {failed_names}", 'yellow'))

def search_identifier(identifier, object_name):
    result_table = Simbad.query_objectids(object_name)
    matches = [row for row in result_table['ID'] if re.search(identifier, row)]
    if matches:
        for match in matches:
            print(colored(f"Found {identifier} identifier for {object_name}:", 'green'), f'{match}')
    else:
        print(colored(f"No {identifier} identifier for {object_name}", 'red'))

for name in names:
    simbad_name(name)
    # print(colored(f"Star name: {name}", 'green'))
    # identifier = 'GJ'
    # search_identifier(identifier, name)
    