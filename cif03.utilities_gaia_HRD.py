import numpy as np
import pandas as pd
import uncertainties
from uncertainties.umath import *
import matplotlib.pyplot as plt
import astroquery
from astroquery.gaia import Gaia
from astroquery.simbad import Simbad
from astropy.coordinates import SkyCoord
import astropy.units as u
from PyAstronomy import pyasl
import re
import matplotlib
astroquery.simbad.conf.server

# =============================================================================
# DATA
# =============================================================================

Gaia.login(user='ccifuent', password='nOtmY*pAsSw0rd')
job = Gaia.launch_job(
    "SELECT parallax, phot_g_mean_mag, phot_bp_mean_mag, phot_rp_mean_mag "
    "FROM gaiadr3.gaia_source WHERE parallax_over_error >= 5 AND random_index "
    "BETWEEN 0 AND 9999999999999999999")
r = job.get_results()

x = [r[i]['phot_bp_mean_mag']-r[i]['phot_rp_mean_mag'] for i in range(len(r))]
y = [r[i]['phot_g_mean_mag']-5*np.log10(1000/r[i]['parallax']) for i in range(len(r))]

# =============================================================================
# PLOT
# =============================================================================

figsize = (12, 10)
pointsize = 8
linewidth = 2
elinewidth = 2
tickssize = 22
labelsize = 22
legendsize = 18
cblabsize = 18

ylabel = r'$M_G$ [mag]'
xlabel = r'$G_{BP}-G_{RP}$ [mag]'

fig, ax = plt.subplots(figsize=figsize)
cmap = plt.get_cmap('magma_r')

sc = ax.scatter(x, y, c=x, cmap=cmap, s=pointsize)

# =============================================================================
# CUSTOM
# =============================================================================

ax.tick_params(axis='x', labelsize=tickssize, direction='in',
                  top=True, labeltop=False, which='both', labelbottom=True)
ax.tick_params(axis='y', labelsize=tickssize, direction='in',
                  right=True, labelright=False, which='both')
ax.tick_params('both', direction = 'in', length = 10, width = 1.5, which = 'major')
ax.tick_params('both', direction = 'in', length = 5, width = 0.5, which = 'minor')
ax.xaxis.set_tick_params(which='minor', bottom=True, top=True)
ax.minorticks_on()
ax.set_xlabel(xlabel, size=labelsize)
ax.set_ylabel(ylabel, size=labelsize)
plt.gca().invert_yaxis()
sc.set_clim(vmin=min(x), vmax=max(x))

# =============================================================================
# OUTPUT
# =============================================================================

# plt.savefig('Output/plot_Gaia_HDR.pdf', dpi=900, bbox_inches='tight')
plt.show()
