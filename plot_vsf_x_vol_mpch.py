import matplotlib.pyplot as plt
from numpy import*
from astropy import cosmology
import astropy.units as u
import matplotlib
import argparse

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
R03   = cat03[0]/(1.+zlo) #Mpc/h
R11   = dr11[0]/(1.+zlo)
V03   = 4./3.*pi*R03**3 
V11   = 4./3.*pi*R11**3

#Merger localisation from ET projection
#z=0.3
thetas = array([12.557,36.49,104.0]) #deg^2
locs   = deg2mpch(thetas, zlo, d_A03.value)  
print(locs)

#prepare figure
plt.figure()
plt.yscale('log')
plt.xlabel(r'$R_v$ [Mpc$/h$]', fontsize=17)
plt.ylabel(r'$N_{\textrm{voids}} \,(> R_v) \times V(R_v)$',fontsize = 17)
plt.xlim(0,100.)
plt.ylim(1e6,1e9)
plt.plot(R11,dr11[1]*V11,marker='o',ls='--',color='blue',label=r'BOSS DR11')
plt.errorbar(R11,dr11[1]*V11,xerr=0.0,yerr=dr11[2],ls='--',marker='o',fmt='blue')
plt.plot(R03,cat03[1]*V03,marker='o',ls='--',color='teal',label=r'BOSS DR12')
plt.errorbar(R03,cat03[1]*V03,xerr=0.0,yerr=cat03[2],ls='--',marker='o',fmt='teal')
plt.plot([locs[0],locs[0]],[1e6,1e10],color='orange',ls=':',label=r'1 - 1 $M_\odot$')
plt.plot([locs[1],locs[1]],[1e6,1e10],color='red',ls=':',label=r'0.5 - 0.5 $M_\odot$')
plt.plot([locs[2],locs[2]],[1e6,1e10],color='yellow',ls=':',label=r'0.25 - 0.25 $M_\odot$')
plt.legend(loc='upper right',prop={'size':13},ncols=1, frameon=False)
plt.savefig("NvoidsVolume_of_R_mpc.pdf",bbox_inches='tight')
plt.show()
