import matplotlib.pyplot as plt
from numpy import*
from astropy import cosmology
import astropy.units as u
import matplotlib
import argparse

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

matplotlib.rc('xtick', labelsize=14)
matplotlib.rc('ytick', labelsize=14)

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
zhi = 0.55
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

d_A05 = cosmo.angular_diameter_distance(zhi) #Mpc
cat05 = loadtxt("R_Nvoids_sNvoids_BOSS_zhi.dat").T
R05   = cat05[0] #cMpc/h

fe    = [1.,0.75,0.5] #fraction of radius within which void is empty

#Merger localisation from ET projection
#z=0.3
zfix03            = '03'
tab5003           = array([3.0,3.0,3.0,3.0]) #1Msun
tab9003           = array([8.0,11.0,10.0,10.0])
tab5003           = deg2mpch(tab5003,zlo,d_A03.value) #mpc/h
tab9003           = deg2mpch(tab9003,zlo,d_A03.value)
localisation5003  = array([mean(tab5003)]*2)
localisation9003  = array([mean(tab9003)]*2)
slocalisation5003 = array([std(tab5003)]*2) 
slocalisation9003 = array([std(tab9003)]*2)
#z=0.55
zfix05            = '055'
tab5005           = array([0.482,0.482,0.482,0.482]) #10Msun tab50 = array([8.0,17.0,52.0,404.0])
tab9005           = array([1.61,1.61,1.61,1.61])
tab5005           = deg2mpch(tab5005,zhi,d_A05.value) #mpc/h
tab9005           = deg2mpch(tab9005,zhi,d_A05.value)
localisation5005  = array([mean(tab5005)]*2) 
localisation9005  = array([mean(tab9005)]*2) 
slocalisation5005 = array([std(tab5005)]*2) 
slocalisation9005 = array([std(tab9005)]*2) 

#prepare figure
colors = ['blue','teal','cyan']
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, sharey=True, figsize=(13,4))
plt.subplots_adjust(hspace=0.1)
ax1.set_yscale('log')
ax2.set_yscale('log')
ax3.set_yscale('log')
ax1.set_xlabel(r'$R_v$ [Mpc/$h$]', fontsize=14)
ax1.set_ylabel(r'$N_{\textrm{voids}} \,(> R_v)$',fontsize = 14)
ax2.set_xlabel(r'$R_v$ [Mpc/$h$]', fontsize=14)
ax3.set_xlabel(r'$R_v$ [Mpc/$h$]', fontsize=14)
#remove y label for right plot
for ax in fig.get_axes():
    ax.label_outer()
ax1.title.set_text(r'BOSS DR11, z='+str(zlo))
ax2.title.set_text(r'BOSS DR12, z='+str(zlo))
ax3.title.set_text(r'BOSS DR12, z='+str(zhi))

ax3.text(70,8e2,r'2 x 10 $M_\odot$', fontsize=12)
ax2.text(70,8e2,r'2 x 1 $M_\odot$', fontsize=12)
ax1.text(70,8e2,r'2 x 1 $M_\odot$', fontsize=12)

#ax1.set_xlim(-1.0,23.)
#ax2.set_xlim(-1.0,23.)
#ax3.set_xlim(-1.0,13.)
for ax in fig.get_axes():
	ax.set_ylim(0.7,20000.)
#compute theta and plot
for f in fe:
	i           = fe.index(f)	
	if f == 1.:
		ax1.plot(R11,dr11[1],marker='o',ls='--',color=colors[i])
		ax1.errorbar(R11,dr11[1],xerr=0.0,yerr=dr11[2],ls='--',fmt=colors[i])
		ax2.plot(R03,cat03[1],marker='o',ls='--',color=colors[i],label=r'R')
		ax2.errorbar(R03,cat03[1],xerr=0.0,yerr=cat03[2],ls='--',fmt=colors[i])
		ax3.plot(R05,cat05[1],marker='o',ls='--',color=colors[i])
		ax3.errorbar(R05,cat05[1],xerr=0.0,yerr=cat05[2],ls='--',fmt=colors[i])
	else:
		ax1.plot(R11*f,dr11[1],marker='o',ls='--',color=colors[i],markerfacecolor='none', markeredgecolor=colors[i])
		ax1.errorbar(R11*f,dr11[1],xerr=0.0,yerr=dr11[2],ls='--',fmt=colors[i])
		ax3.plot(R05*f,cat05[1],marker='o',ls='--',color=colors[i],markerfacecolor='none', markeredgecolor=colors[i])
		ax3.errorbar(R05*f,cat05[1],xerr=0.0,yerr=cat05[2],ls='--',fmt=colors[i])
		ax2.plot(R03*f,cat03[1],marker='o',ls='--',color=colors[i],markerfacecolor='none', markeredgecolor=colors[i],label=r''+str(f)+' R')
		ax2.errorbar(R03*f,cat03[1],xerr=0.0,yerr=cat03[2],ls='--',fmt=colors[i])

#finalise and save figure
ax1.plot(localisation5003,[1e-1,1e5],color='k',ls='--',label=r'50\%  area')
ax1.plot(localisation9003,[1e-1,1e5],color='k',ls=':',label=r'90\% area')
ax2.plot(localisation5003,[1e-1,1e5],color='k',ls='--')
ax2.plot(localisation9003,[1e-1,1e5],color='k',ls=':')
ax3.plot(localisation5005,[1e-1,1e5],color='k',ls='--')
ax3.plot(localisation9005,[1e-1,1e5],color='k',ls=':')
ax1.fill_betweenx([1e-1,1e5],localisation5003-slocalisation5003,localisation5003+slocalisation5003,color='k',alpha=0.2)
ax1.fill_betweenx([1e-1,1e5],localisation9003-slocalisation9003,localisation9003+slocalisation9003,color='k',alpha=0.2)
ax2.fill_betweenx([1e-1,1e5],localisation5003-slocalisation5003,localisation5003+slocalisation5003,color='k',alpha=0.2)
ax2.fill_betweenx([1e-1,1e5],localisation5003-slocalisation5003,localisation5003+slocalisation5003,color='k',alpha=0.2)
ax1.legend(loc='upper right',prop={'size':12}, frameon=False)
ax2.legend(loc='upper right',prop={'size':12}, frameon=False)
ax3.legend(loc='upper right',prop={'size':12}, frameon=False)
plt.savefig("Nvoids_of_R_detection_area_mpc.pdf",bbox_inches='tight')
plt.show()
