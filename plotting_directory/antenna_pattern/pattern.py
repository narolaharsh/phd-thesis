import numpy as np
import matplotlib.pyplot as plt
import bilby
import scienceplots

from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy import units as u
import ligo.skymap.plot
from ligo.skymap.io.fits import read_sky_map



from astropy.coordinates import SkyCoord
import astropy.units as u

plt.style.use(['science'])


plt.rcParams['axes.labelsize'] = 11.
plt.rc('xtick',labelsize=11)
plt.rc('ytick',labelsize=11)


red = (0.80, 0.36, 0.36)
green = (0.62, 0.79, 0.62)
blue = (0.27, 0.51, 0.71)


ifos = bilby.gw.detector.InterferometerList(['V1'])


N = int(1e5)
ra = np.random.uniform(-np.pi, np.pi, N)
dec = np.random.uniform(-np.pi/2, np.pi/2, N)
psi = 0
tc = 1187008882.4


ra_deg = np.rad2deg(ra) *u.deg
dec_deg = np.rad2deg(dec) *u.deg

plus = np.zeros(N)
cross = np.zeros(N)
gridsize = 40
cmap = 'copper'
for ii in range(N):
    plus[ii] = ifos[0].antenna_response(ra[ii], dec[ii], tc, psi, 'plus')
    cross[ii] =  ifos[0].antenna_response(ra[ii], dec[ii], tc, psi, 'cross')

fig, axes = plt.subplots(1, 3, sharex=True, sharey=True, gridspec_kw={'wspace':.05}, figsize = (12, 3))
ax = axes[0]
ax.hexbin(ra_deg, dec_deg, C = plus, cmap = cmap, gridsize = gridsize)
ax.set_xticks([-180, -90, 0, 90, 180])
ax.set_yticks([-90, -45, 0, 45, 90])
ax.set_xlabel("Right Ascension [deg]")
ax.set_ylabel("Declination [deg]")
ax.set_title('Sensitivity to plus polarization $(F_+)$')

ax = axes[1]
ax.hexbin(ra_deg, dec_deg, C = cross, cmap = cmap, gridsize = gridsize)
ax.set_xticks([-180, -90, 0, 90, 180])
ax.set_yticks([-90, -45, 0, 45, 90])
ax.set_xlabel("Right Ascension [deg]")
#ax.set_ylabel("Declination [deg]")
ax.set_title('Sensitivity to cross polarization $(F_{\\times})$')


ax = axes[2]
total = (plus**2 + cross**2)**0.5
data = ax.hexbin(ra_deg, dec_deg, C = total, cmap = cmap, gridsize = gridsize)
ax.set_xticks([-180, -90, 0, 90, 180])
ax.set_yticks([-90, -45, 0, 45, 90])
ax.set_xlabel("Right Ascension [deg]")
#ax.set_ylabel("Declination [deg]")
ax.set_title('$\\sqrt{F^2_+ + F^2_{\\times}}$')


fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.81, 0.15, 0.015, 0.7])
fig.colorbar(data, cax=cbar_ax)
fig.savefig('../../figures/joint_pattern.pdf')
fig.savefig('joint_pattern.pdf')

