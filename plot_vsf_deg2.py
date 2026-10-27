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

zlo = 0.3
zhi = 0.55
h   = 0.6831
#MultiDark cosmology
#cosmo  = cosmology.FlatLambdaCDM(H0=100.*h, Om0=0.307115, Ob0=0.048206)#, Ode0=0.692885)
#Planck + BOSS DR12 cosmology (Philcox and Ivanov 2022)
#cosmo  = cosmology.Planck18.clone(H0=100.*h, Om0=0.32)
#BOSS DR12 cosmology w/ voids parameters (Contarini et al. 2022)
cosmo  = cosmology.Planck18.clone(H0=100.*h, Om0=0.29)
cosmo.meta['sigma8']=0.79
print(cosmo)
print(cosmo.meta)

d_A03 = cosmo.angular_diameter_distance(zlo) #Mpc
cat03 = loadtxt("R_Nvoids_sNvoids_BOSSDR12_zlo_cumul.dat").T
dr11  = loadtxt("R_Nvoids_sNvoids_BOSSDR11_LOWZ_cumul.dat").T
R03   = cat03[0]/h/(1.+zlo) #Mpc
R11   = dr11[0]/h/(1.+zlo)
V03   = 4./3.*pi*R03**3
V11   = 4./3.*pi*R11**3

d_A05 = cosmo.angular_diameter_distance(zhi) #Mpc
cat05 = loadtxt("R_Nvoids_sNvoids_BOSS_zhi.dat").T
R05   = cat05[0]/h/(1.+zhi) #Mpc
V05   = 4./3.*pi*R05**3

fe    = [1.]#,4./3.,2.] #fraction of radius within which void is empty

#Merger localisation from ET projection
#z=0.3
zfix03  = '03'
tab5003 = array([3.0,3.0,3.0,3.0]) #1Msun
tab9003 = array([8.0,11.0,10.0,10.0])
localisation5003  = array([mean(tab5003)]*2)
localisation9003  = array([mean(tab9003)]*2) #average degree^2
slocalisation5003 = array([std(tab5003)]*2) #average degree^2
slocalisation9003 = array([std(tab9003)]*2) #std dev degree^2
#z=0.55
zfix05  = '055'
tab5005 = array([0.482,0.482,0.482,0.482]) #10Msun tab50 = array([8.0,17.0,52.0,404.0])
tab9005 = array([1.61,1.61,1.61,1.61])
localisation5005  = array([mean(tab5005)]*2) 
localisation9005  = array([mean(tab9005)]*2) #average degree^2
slocalisation5005 = array([std(tab5005)]*2) #average degree^2
slocalisation9005 = array([std(tab9005)]*2) #std dev degree^2

#prepare figure
colors = ['blue','teal','cyan']
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, sharey=True, figsize=(13,4))
plt.subplots_adjust(hspace=0.1)
ax1.set_yscale('log')
ax2.set_yscale('log')
ax3.set_yscale('log')
ax1.set_xlabel(r'$R_v$ [deg$^2$]', fontsize=14)
ax1.set_ylabel(r'$N_{\textrm{voids}} \,(> R_v) \times V(R_v)$',fontsize = 14)
ax2.set_xlabel(r'$R_v$ [deg$^2$]', fontsize=14)
ax3.set_xlabel(r'$R_v$ [deg$^2$]', fontsize=14)
#remove y label for right plot
for ax in fig.get_axes():
    ax.label_outer()
ax1.title.set_text(r'BOSS DR11, z='+str(zlo))
ax2.title.set_text(r'BOSS DR12, z='+str(zlo))
ax3.title.set_text(r'BOSS DR11, z='+str(zhi))

ax3.text(7.5,2000,r'2 x 10 $M_\odot$', fontsize=12)
ax2.text(16,200,r'2 x 1 $M_\odot$', fontsize=12)
ax1.text(16,350,r'2 x 1 $M_\odot$', fontsize=12)

ax1.set_xlim(-1.0,130.)
ax2.set_xlim(-1.0,130.)
ax3.set_xlim(-1.0,55.)
for ax in fig.get_axes():
	ax.set_ylim(2e6,3e9)
#compute theta and plot
for f in fe:
	i           = fe.index(f)
	theta03     = arctan(R03/f/d_A03.value) * 180./pi #degrees
	void_area03 = pi*(theta03)**2. #degree^2 assuming spherical void
	theta11     = arctan(R11/f/d_A03.value) * 180./pi
	void_area11 = pi*(theta11)**2.
	theta05     = arctan(R05/f/d_A05.value) * 180./pi
	void_area05 = pi*(theta05)**2.
	
	if f == 1.:
		ax1.plot(void_area11,dr11[1]*V11,marker='o',ls='--',color=colors[i])
		ax1.errorbar(void_area11,dr11[1]*V11,xerr=0.0,yerr=dr11[2],ls='--',fmt=colors[i])
		ax2.plot(void_area03,cat03[1]*V03,marker='o',ls='--',color=colors[i],label=r'R')
		ax2.errorbar(void_area03,cat03[1]*V03,xerr=0.0,yerr=cat03[2],ls='--',fmt=colors[i])
		ax3.plot(void_area05,cat05[1]*V05,marker='o',ls='--',color=colors[i])
		ax3.errorbar(void_area05,cat05[1]*V05,xerr=0.0,yerr=cat05[2],ls='--',fmt=colors[i])
	else:
		ax1.plot(void_area11,dr11[1],marker='o',ls='--',color=colors[i],markerfacecolor='none', markeredgecolor=colors[i])
		ax1.errorbar(void_area11,dr11[1],xerr=0.0,yerr=dr11[2],ls='--',fmt=colors[i])
		ax3.plot(void_area05,cat05[1],marker='o',ls='--',color=colors[i],markerfacecolor='none', markeredgecolor=colors[i])
		ax3.errorbar(void_area05,cat05[1],xerr=0.0,yerr=cat05[2],ls='--',fmt=colors[i])
		ax2.plot(void_area03,cat03[1],marker='o',ls='--',color=colors[i],markerfacecolor='none', markeredgecolor=colors[i],label=r'R/'+str(int(f)))
		ax2.errorbar(void_area03,cat03[1],xerr=0.0,yerr=cat03[2],ls='--',fmt=colors[i])

#finalise and save figure
#ax1.plot(localisation5003,[1e-1,1e5],color='k',ls='--',label=r'50\%  area')
ax1.plot(localisation9003,[1e6,1e10],color='k',ls=':',label=r'90\% area')
#ax2.plot(localisation5003,[1e-1,1e5],color='k',ls='--')
ax2.plot(localisation9003,[1e6,1e10],color='k',ls=':')
#ax3.plot(localisation5005,[1e-1,1e5],color='k',ls='--')
ax3.plot(localisation9005,[1e6,1e10],color='k',ls=':')
ax1.fill_betweenx([1e-1,1e5],localisation5003-slocalisation5003,localisation5003+slocalisation5003,color='k',alpha=0.2)
ax1.fill_betweenx([1e-1,1e5],localisation9003-slocalisation9003,localisation9003+slocalisation9003,color='k',alpha=0.2)
ax2.fill_betweenx([1e-1,1e5],localisation5003-slocalisation5003,localisation5003+slocalisation5003,color='k',alpha=0.2)
ax2.fill_betweenx([1e-1,1e5],localisation5003-slocalisation5003,localisation5003+slocalisation5003,color='k',alpha=0.2)
ax1.legend(loc=0,prop={'size':12}, frameon=False)
#ax2.legend(loc=0,prop={'size':12}, frameon=False)
#ax3.legend(loc=0,prop={'size':12}, frameon=False)
plt.savefig("NvoidsVolume_of_R.pdf",bbox_inches='tight')
plt.show()
