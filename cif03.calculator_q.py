import numpy as np
import pandas as pd

# =============================================================================
# DATA
# =============================================================================

input_file = 'cif03.v10'
output_file = input_file + '_out.csv'
df = pd.read_csv('Data/' + input_file + '.csv', sep=",", header=0)

# Filter rows where 'ID_system' is not '0'
df_filtered = df[df['ID_system'] != '0']

# # =============================================================================
# # FUNCTIONS
# # =============================================================================

# def append_masses(i, j):
#     Mass_A.append(df['M_Msol'][j])
#     eMass_A.append(df['eM_Msol'][j])
#     Mass_B.append(df['M_Msol'][i])
#     eMass_B.append(df['eM_Msol'][i])

# def calculate_total_mass(indexes):
#     masses = df.loc[indexes, 'M_Msol']
#     errors = df.loc[indexes, 'eM_Msol']
#     Mass_total.append(masses.sum())
#     eMass_total.append(np.sqrt((errors**2).sum()))

# # =============================================================================
# # COMPUTE
# # =============================================================================

# Mass_A, eMass_A = [], []
# Mass_B, eMass_B = [], []
# Mass_total, eMass_total = [], []

# for i in range(3, len(df) - 3):
#     if not pd.isnull(df['System'][i]):
#         Mass_A.append(df['M_Msol'][i])
#         eMass_A.append(df['eM_Msol'][i])
#         Mass_B.append(empty)
#         eMass_B.append(empty)
#     elif ID_system[i] == ID_system[i - 1] and not pd.isnull(df['System'][i - 1]):
#         append_masses(i, i - 1)
#     elif ID_system[i] == ID_system[i - 2] and not pd.isnull(df['System'][i - 2]):
#         append_masses(i, i - 2)
#     else:
#         Mass_A.append(empty)
#         eMass_A.append(empty)
#         Mass_B.append(empty)
#         eMass_B.append(empty)

# for i in range(3, len(df) - 3):
#     if ID_system[i] == ID_system[i + 1] == ID_system[i + 2]:
#         calculate_total_mass([i, i + 1, i + 2])
#     elif ID_system[i - 1] == ID_system[i] == ID_system[i + 1]:
#         calculate_total_mass([i - 1, i, i + 1])
#     elif ID_system[i - 2] == ID_system[i - 1] == ID_system[i]:
#         calculate_total_mass([i - 2, i - 1, i])
#     elif ID_system[i] == ID_system[i + 1] != ID_system[i + 2] != ID_system[i - 1]:
#         calculate_total_mass([i, i + 1])
#     elif ID_system[i - 1] == ID_system[i] != ID_system[i - 2] != ID_system[i + 1] != ID_system[i - 2]:
#         calculate_total_mass([i - 1, i])
#     else:
#         Mass_total.append(empty)
#         eMass_total.append(empty)

# ID_star = df['ID_star'][3:len(df) - 3].tolist()

# df_results = {
#     'ID_star': ID_star,
#     'Mass_A': Mass_A,
#     'eMass_A': eMass_A,
#     'Mass_B': Mass_B,
#     'eMass_B': eMass_B,
#     'Mass_tot': np.round(Mass_total, 5),
#     'eMass_tot': np.round(eMass_total, 5)
# }

# # =============================================================================
# # OUT
# # =============================================================================

# save_csv = 'yes'
# output_full = 'no'

# if save_csv == 'yes':
#     df_append = pd.DataFrame(data=df_results)
#     if output_full == 'yes':
#         output = pd.concat([df, df_append], axis=1)
#     else:
#         output = df_append
#     output.to_csv('Output/' + output_file, sep=',', encoding='utf-8', index=False)

# =============================================================================
# COMPUTE
# =============================================================================

# Group by 'ID_system' and calculate the sum of masses and the combined uncertainty
grouped = df_filtered.groupby('ID_system').agg(
    Mtot=('M_Msol', 'sum'),
    eMtot=('eM_Msol', lambda x: np.sqrt((x**2).sum()))
).reset_index()

# Round the computed values to 4 decimals
grouped['Mtot'] = grouped['Mtot'].round(4)
grouped['eMtot'] = grouped['eMtot'].round(4)

# Merge the results with the original DataFrame to get all 'ID_star' rows
df_result = pd.merge(df, grouped, on='ID_system', how='left')

# Only keep Mtot and eMtot values where 'Component' starts with 'A', else set them to NaN
df_result.loc[~df_result['Component'].str.startswith('A'), ['Mtot', 'eMtot']] = np.nan

# Also set Mtot and eMtot to NaN where ID_system is '0'
df_result.loc[df_result['ID_system'] == '0', ['Mtot', 'eMtot']] = np.nan

# Select only the required columns
df_output = df_result[['ID_star', 'Mtot', 'eMtot']]

# =============================================================================
# OUT
# =============================================================================

save_csv = 'yes'

if save_csv == 'yes':
    df_output.to_csv('Output/' + output_file, sep=',', encoding='utf-8', index=False)