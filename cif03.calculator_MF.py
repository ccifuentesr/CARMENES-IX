import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc
from pathlib import Path
from statsmodels.stats.proportion import proportion_confint

# =============================================================================
# DATA AND FILTERS
# =============================================================================

input_file = 'cif03.v30'
df = pd.read_csv('Data/'+input_file+'.csv', sep=",", header=0)


# filter_karmn = (pd.isnull(df['Karmn']) == False)       # Only M dwarfs
filter_karmn = (pd.isnull(df['Karmn']) == False) & ((df['SpTnum'] >= 75))       # Only M dwarfs < M5.5
# filter_karmn = (pd.isnull(df['Karmn']) == False) & (df['Category'] == 3)       # Only M dwarfs with categories
# filter_karmn = (df['ID_star'] == df['ID_star'])         # All primaries
filter_primaries = (pd.isnull(df['System']) == False)   # Only components 'A' (primaries)
filter_candidate = (df['Type'] == 'Candidate')          # Only candidates to binaries from single 
filter_candidate_multiple = (df['Type'] == 'Multiple+') # Only candidates to binaries from multiple 
filter_double = (df['Class'] == 'Double')
filter_triple = (df['Class'] == 'Triple')
filter_quadruple = (df['Class'] == 'Quadruple')
filter_quintuple = (df['Class'] == 'Quintuple')
filter_sextuple = (df['Class'] == 'Sextuple')

# =============================================================================
# MF & UNCERTAINTIES
# =============================================================================

n_all =  len(df[filter_karmn])
n_multiples = len(df[filter_karmn & filter_primaries])
n_doubles = len(df[filter_karmn & filter_double & filter_primaries])
n_doubles_cand = len(df[filter_karmn & filter_candidate])
n_triples = len(df[filter_karmn & filter_triple & filter_primaries])
n_triples_cand = len(df[filter_karmn & filter_triple & filter_candidate_multiple])
n_quadruples = len(df[filter_karmn & filter_quadruple & filter_primaries])
n_quintuples = len(df[filter_karmn & filter_quintuple & filter_primaries])
n_sextuples = len(df[filter_karmn & filter_sextuple & filter_primaries])
n_with_candidates = n_multiples + n_doubles_cand
MF = n_multiples/n_all
MF_cand = n_with_candidates/n_all
CSF = (n_doubles + n_triples*2 + n_quadruples*3 + n_quintuples*4)/n_all
CSF_cand = (n_doubles + n_doubles_cand + n_triples*2 + n_triples_cand*2 + n_quadruples*3 + n_quintuples*4)/n_all

confint_low_MF, confint_upp_MF = proportion_confint(n_multiples, n_all, alpha=0.05, method='wilson')
confint_low_cand_MF, confint_upp_cand_MF = proportion_confint(n_with_candidates, n_all, alpha=0.05, method='wilson')
confint_low_CSF, confint_upp_CSF = proportion_confint((n_doubles + n_triples*2 + n_quadruples*3 + n_quintuples*4), n_all, alpha=0.05, method='wilson')
confint_low_cand_CSF, confint_upp_cand_CSF = proportion_confint(n_with_candidates, n_all, alpha=0.05, method='wilson')
#
latex_MF = f"MF $= {100*MF:.1f}^{{+{100*(confint_upp_MF - MF):.1f}}}_{{-{100*(MF - confint_low_MF):.1f}}}$"
latex_MF_cand = f"MF+ $= {100*MF_cand:.1f}^{{+{100*(confint_upp_cand_MF - MF_cand):.1f}}}_{{-{100*(MF_cand - confint_low_cand_MF):.1f}}}$"
latex_CSF = f"CSF $= {CSF:.3f}^{{+{(confint_upp_CSF - CSF):.3f}}}_{{-{(CSF - confint_low_CSF):.3f}}}$"
latex_CSF_cand = f"CSF+ $= {CSF_cand:.3f}^{{+{(confint_upp_cand_CSF - CSF_cand):.3f}}}_{{-{(CSF_cand - confint_low_cand_CSF):.3f}}}$"

print(f'MF = {MF*100:.2f}%')
print(f'95% CI = [{confint_low_MF*100:.2f}%, {confint_upp_MF*100:.2f}%]')
print(f'MF+ = {MF_cand*100:.2f}%')
print(f'95% CI = [{confint_low_cand_MF*100:.2f}%, {confint_upp_cand_MF*100:.2f}%]')
print(f'CSF = {CSF:.4f}')
print(f'95% CI = [{confint_low_CSF:.4f}, {confint_upp_CSF:.4f}]')
print(f'CSF+ = {CSF_cand:.4f}')
print(f'95% CI = [{confint_low_cand_CSF:.4f}, {confint_upp_cand_CSF:.4f}]')
print('---\n')
print(latex_MF)
print(latex_MF_cand)
print(latex_CSF)
print(latex_CSF_cand)
print('---\n')
print(f'Total = {len(df):.0f}')
print(f'Primaries = {n_multiples:.0f}')
print(f'Single = {len(df[filter_karmn & (df['Type'] == 'Single')]):.0f}')
print(f'Candidate = {len(df[filter_karmn & (df['Type'] == 'Candidate')]):.0f}')
print(f'Multiple+ = {len(df[filter_karmn & (df['Type'] == 'Multiple+')]):.0f}')
print(f'Multiple = {len(df[filter_primaries & filter_karmn & (df['Type'] == 'Multiple')]):.0f}')
print('---\n')
print(f'Double = {n_doubles:.0f}')
print(f'Triple = {n_triples:.0f}')
print(f'Quadruple = {n_quadruples:.0f}')
print(f'Quintuple = {n_quintuples:.0f}')
print(f'Sextuple = {n_sextuples:.0f}')
