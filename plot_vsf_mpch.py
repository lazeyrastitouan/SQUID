import matplotlib.pyplot as plt
from numpy import*
from astropy import cosmology
import astropy.units as u
import matplotlib
import argparse
import sys
import astropy 

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

matplotlib.rc('xtick', labelsize=16)
matplotlib.rc('ytick', labelsize=16)

def deg2mpch(thetasq, z , dA):
	'''
	Convert circular area in degree square to radius in mpc/h 
	Parameters: 
	thetasq : area in degree square
	z       : redshift
	da      : angular diameter distance
	'''
	theta = sqrt(thetasq/pi) #in degrees
	theta = theta * pi/180.  #in radians
	return  theta * dA * (1.+z) * h

zlo = 0.3
h   = 0.6831
#BOSS DR12 cosmology w/ voids parameters (Contarini et al. 2022)
cosmo  = cosmology.Planck18.clone(H0=100.*h, Om0=0.29)
cosmo.meta['sigma8']=0.79
print(cosmo)
print(cosmo.meta)

d_A03 = cosmo.angular_diameter_distance(zlo) #Mpc
cat03 = loadtxt("R_Nvoids_sNvoids_BOSSDR12_zlo_cumul.dat").T
dr11  = loadtxt("R_Nvoids_sNvoids_BOSSDR11_LOWZ_cumul.dat").T
R03   = cat03[0] #cMpc/h
R11   = dr11[0] #cmpc/h

#Merger localisation from ET projection
#z=0.3
zfix03            = '03'
tab9003           = deg2mpch(array([12.6]),zlo,d_A03.value) #mpc/h
localisation9003  = array([mean(tab9003)]*2)
tab0505           = deg2mpch(array([36.5]),zlo,d_A03.value) #mpc/h
localisation0505  = array([tab0505[0]]*2)

print(localisation9003)
print(localisation0505)

#prepare figure
plt.figure()
plt.yscale('log')
plt.xlabel(r'$R_v$ [Mpc/$h$]', fontsize=17)
plt.ylabel(r'$N_{\textrm{voids}} \,(> R_v)$',fontsize = 17)

plt.text(7,1,r'$z=0.3$', fontsize=14)

plt.ylim(0.7,5000.)
plt.xlim(5,100)
plt.plot(R11,dr11[1],marker='o',ls='--',color='blue',label=r'BOSS DR11')
plt.errorbar(R11,dr11[1],xerr=0.0,yerr=dr11[2],ls='--',fmt='bo')
plt.plot(R03,cat03[1],marker='o',ls='--',color='teal',label=r'BOSS DR12')
plt.errorbar(R03,cat03[1],xerr=0.0,yerr=cat03[2],ls='--',marker='o',fmt='teal')

#finalise and save figure
#plt.plot(localisation5003,[1e-1,1e5],color='k',ls='--',label=r'50\%  area')
plt.plot(localisation9003,[1e-1,1e5],color='k',ls='--',label=r'1-1 $M_\odot$ 90\% loc.')
plt.plot(localisation0505,[1e-1,1e5],color='grey',ls=':',label=r'0.5-0.5 $M_\odot$ 90\% loc.')
#plt.fill_betweenx([1e-1,1e5],localisation5003-slocalisation5003,localisation5003+slocalisation5003,color='k',alpha=0.2)
#plt.fill_betweenx([1e-1,1e5],localisation9003-slocalisation9003,localisation9003+slocalisation9003,color='k',alpha=0.2)
plt.legend(loc='upper right',ncol=1, prop={'size':11}, frameon=False)
plt.savefig("Nvoids_of_R_detection_area_mpc_z0p3.pdf",bbox_inches='tight')
plt.show()
