import numpy as np
import uncertainties
from uncertainties.umath import *
import astropy.constants as const
import pandas as pd

# =============================================================================
# DATA
# =============================================================================

input_file = 'cif03.v00'
output_file = input_file + '_out.csv'
df = pd.read_csv('Data/'+input_file+'.csv', sep=",", header=0)

L_Lsol = [uncertainties.ufloat(df['L_Lsol'][i], df['eL_Lsol'][i]) for i in range(len(df))]
Teff_K = [uncertainties.ufloat(df['Teff_K'][i], df['eTeff_K'][i]) for i in range(len(df))]
# parallax = [uncertainties.ufloat(df['parallax'][i], df['parallax_error'][i]) for i in range(len(df))]
# rho_arcsec = df['rho_arcsec']

write = True

# =============================================================================
# CONSTANTS
# =============================================================================

Mterra = 5.97237*1e24 # kg (exact)
Rterra = 6.3781*1e6 # m (exact)
Lsol = 3.828*1e26 # W (exact)
Rsol = 6.957*1e8 # m (exact)
Msol = 1.98847*1E30  # kg (exact)
au_m = 149597870700 # m (exact)
G = uncertainties.ufloat(6.67430*1e-11, 0.00015*1e-11) # m3 kg−1 s−2
GM = 1.3271244*1e20 # m3 s−2 (exact)
sigma = 5.670374419*1e-8 # W m-2 K-4 (exact)
Ssol = 1361 # W m-2

# =============================================================================
# FUNCTIONS
# =============================================================================

def Radius_SB(L_Lsol, Teff_K):
    """Stellar radius and its error from the Stefan–Boltzmann law under the black body approximation.

    Args:
        Lbol (float): Bolometric luminosity in solar units.
        Teff (float): Effective temperature in Kelvin.

    Returns:
        float: Stellar radius in solar units.
        float: Stellar radius error in solar units.

    Nominal solar values from the IAU B.3 resolution
    on recommended nominal conversion constants for selected solar and planetary properties:
    https://www.iau.org/static/resolutions/IAU2015_English.pdf

    Nominal solar luminosity: 3.828 x 10+26 W (exact)
    Nominal solar radius: 6.957 x 10+8 m (exact)

    Stefan-Boltzmann constant value from 2018 CODATA recommended values:
    https://physics.nist.gov/cuu/pdf/wall_2018.pdf

    Stefan-Boltzman constant, sigma: 5.670 374 419 x 10-8 W m-2 K-4 (exact)

    """
    R_Rsol = 1/Rsol * sqrt(L_Lsol*Lsol/(4*np.pi*sigma*Teff_K**4))
    return(R_Rsol)

def Teff_SB(L_Lsol, R_Rsol):
    """Effective temperature and its error from the Stefan–Boltzmann law under the black body approximation.

    Args:
        Lbol (float): Bolometric luminosity in solar units.
        Radius (float): Stellar radius in solar units.

    Returns:
        float: Effective temperature in Kelvin.
        float: Effective temperature error in Kelvin.
    """
    Teff_K = (L_Lsol*Lsol / (R_Rsol*Rsol)**2 / (4*np.pi*sigma))**(1/4)
    return(Teff_K)

def Mass_sch19(R_Rsol):
    """Stellar mass and its error from the empirical relation by Schweitzer et al. 2019
    (2019A&A...625A..68S), based on masses and radii of eclipsing binaries.

    Args:
        R_Rsol (float): Stellar radius in solar units.
        eR_Rsol (float): Stellar radius uncertainty in solar units.

    Returns:
        float: Stellar mass in solar units.
        float: Stellar mass error in solar units.

    (See Equation 6 in Schweitzer et al. 2019 and references therein).
    """
    a = uncertainties.ufloat(-0.024048024, 0.007592668)
    b = uncertainties.ufloat(1.0552427, 0.017044148)
    M_Msol = a + b * R_Rsol
    return(M_Msol)

def Porb(rho_arcsec, parallax, M1_Msol, M2_Msol):
    """Calculates the orbital period in years.

    Args:
        rho (float): angular separation in arcsec.
        parallax (float): parallax in milliarcseconds.
        M1_Msol (float): Mass of the primary in solar units.
        M2_Msol (float): Mass of the secondary in solar units.

    Returns:
        float: orbital period in days.

    Nominal solar values from the IAU B2 (AU) and B3 (GM) resolutions.
    Gravitational constant value from 2018 CODATA recommended values.
    """
    d_pc = 1000/parallax
    s_au = rho_arcsec * d_pc # au
    s_m = s_au * au_m # m
    mu = GM * (M1_Msol + M2_Msol) # m3 s-2
    Porb_s = 2*np.pi * sqrt(s_m**3/mu) # s
    Porb_d = Porb_s/86400 # d
    Porb_a = Porb_d/365.25 # yr
    return(Porb_a)

# =============================================================================
# WRITE OUT
# =============================================================================

if write:
    Radius = [Radius_SB(L_Lsol[i], Teff_K[i]) for i in range(len(df))]
    Mass = [Mass_sch19(Radius[i]) for i in range(len(df))]
    # Porb = [Porb(rho_arcsec[i], parallax[i], Mass[i]) for i in range(len(df))]

    R_Rsol = [Radius[i].n for i in range(len(df))]
    eR_Rsol = [Radius[i].s for i in range(len(df))]
    M_Msol = [Mass[i].n for i in range(len(df))]
    eM_Msol = [Mass[i].s for i in range(len(df))]
    # Porb_d = [Porb[i].n for i in range(len(df))]
    # ePorb_d = [Porb[i].s for i in range(len(df))]

    df_append = pd.DataFrame({'ID_star': df['ID_star'], 'R_Rsol': R_Rsol, 'eR_Rsol': eR_Rsol,\
    'M_Msol': M_Msol, 'eM_Msol': eM_Msol})#, 'Porb_d': Porb_d, 'ePorb_d': ePorb_d})

    save_csv = True
    output_full = False

    if save_csv:
        df_append = pd.DataFrame(data=df_append)
        if output_full:
            output = pd.concat([df, df_append], axis=1)
        else:
            output = df_append
        output.to_csv('Output/' + output_file, sep=',', encoding='utf-8')

L_Lsol_test = uncertainties.ufloat(0.3249495600,0.0020070933)
Teff_K_test = uncertainties.ufloat(4500, 50)

print(Radius_SB(L_Lsol_test, Teff_K_test))
print(Mass_sch19(Radius_SB(L_Lsol_test, Teff_K_test)))