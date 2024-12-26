import numpy as np
import pandas as pd
import uncertainties
from uncertainties.umath import *

# =============================================================================
# DATA
# =============================================================================

input_file = 'cif03.v05'
output_file = input_file + '_out.csv'
df = pd.read_csv('Data/'+input_file+'.csv', sep=",", header=0)

ID_star = df['ID_star']
parallax = df['parallax']
parallax_error = df['parallax_error']
G_mag = df['G_mag']
eG_mag = df['eG_mag']
RP_mag = df['RP_mag']
eRP_mag = df['eRP_mag']
BP_mag = df['BP_mag']
eBP_mag = df['eBP_mag']
W1_mag = df['W1mag']
eW1_mag = df['e_W1mag']
W2_mag = df['W2mag']
eW2_mag = df['e_W2mag']

# =============================================================================
# FUNCTION
# =============================================================================

def metallicity_duque(G, eG, BP, eBP, RP, eRP, W1, eW1, W2, eW2, pi_mas, epi_mas):
    """
    Estimation of metallicity as described in Duque-Arribas et al. (2023) 
    https://ui.adsabs.harvard.edu/abs/2023ApJ...944..106D/abstract
    """
    parallax = uncertainties.ufloat(pi_mas, epi_mas)
    G_mag = uncertainties.ufloat(G, eG)
    BP_mag = uncertainties.ufloat(BP, eBP)
    RP_mag = uncertainties.ufloat(RP, eRP)
    W1_mag = uncertainties.ufloat(W1, eW1)
    W2_mag = uncertainties.ufloat(W2, eW2)
    MG_mag = G_mag - 5*log10(1000/parallax) + 5
    #
    a = uncertainties.ufloat(5.87, 0.44)
    b = uncertainties.ufloat(-1.40, 0.12)
    c = uncertainties.ufloat(5.08, 0.24)
    d = uncertainties.ufloat(-2.45, 0.14)
    e = uncertainties.ufloat(-0.788, 0.052)
    f = uncertainties.ufloat(0.1094, 0.0071)
    X = W1_mag - W2_mag
    Y = BP_mag - RP_mag
    Z = MG_mag
    Fe_H = a + b*X + c*Y + d*Z + e*Y**2 + f*Z**2
    return Fe_H

# =============================================================================
# ACTION
# =============================================================================

Fe_H = [np.round(metallicity_duque(G_mag[i], eG_mag[i], BP_mag[i], eBP_mag[i], RP_mag[i], eRP_mag[i], W1_mag[i], eW1_mag[i], W2_mag[i], eW2_mag[i], parallax[i], parallax_error[i]).n, 4) for i in range(len(df))]
eFe_H = [np.round(metallicity_duque(G_mag[i], eG_mag[i], BP_mag[i], eBP_mag[i], RP_mag[i], eRP_mag[i], W1_mag[i], eW1_mag[i], W2_mag[i], eW2_mag[i], parallax[i], parallax_error[i]).s, 4) for i in range(len(df))]

# =============================================================================
# WRITE OUT
# =============================================================================

df_append = pd.DataFrame({'ID_star': df['ID_star'], 'Fe_H': Fe_H, 'eFe_H': eFe_H})

save_csv = True

if save_csv == True:
    df_append = pd.DataFrame(data=df_append)
    df_append.to_csv('Output/' + output_file, sep=',', encoding='utf-8')
