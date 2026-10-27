import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from plot_shaded import plot_upper_limit

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

matplotlib.rc('xtick', labelsize=14)
matplotlib.rc('ytick', labelsize=14)


#[0]sigma, [1]Nhod full, [2]Ngal full, [3]err(Nfull), [4]Nhod frac, [5]err(Nhodfrac), [6]Ngalfrac, [7]err(Ngalfrac), [8]Npbhfull, [9]err(Npbhfull), [10]Npbhfrac, [11]err(Npbhfrac)
rates_m05 = np.loadtxt("sigma_Nabh_Nbns_Npbh_Mmin0.5_Rmin49.8_ET_2L.txt").T
rates_m1  = np.loadtxt("sigma_Nabh_Nbns_Npbh_Mmin1.0_Rmin29.3_ET_2L.txt").T

years   = 10
yfactor = 19.53/3.61*8.02/23.22

fig, axs = plt.subplots(1,2,figsize=(12,4),sharey=True)
plt.subplots_adjust(wspace=0.07, right=0.72)
axs[0].set_ylabel(r'BOSS-like \\ \# event / 10yr',fontsize = 14)

sigma         = np.array([0,0.3,0.6,0.8,1,1.2,1.5])
hod_frac_05   = [(rates_m05[4,0]+rates_m05[5,0])*years]*len(sigma)
astro_full_05 = [(rates_m05[1,0]+rates_m05[3,0])*years]*len(sigma)
gal_frac_05   = [(rates_m05[6,0]+rates_m05[7,0])*years]*len(sigma)
hod_frac_1    = [(rates_m1[4,0]+rates_m1[5,0])*years]*len(sigma)
astro_full_1  = [(rates_m1[1,0]+rates_m1[3,0])*years]*len(sigma)
gal_frac_1    = [(rates_m1[6,0]+rates_m1[7,0])*years]*len(sigma)

#axs[0]: Mmin = 0.5 
axs[0].plot(sigma,hod_frac_05,'b--',label=r'ABH HOD void frac')
axs[0].errorbar(sigma,hod_frac_05,xerr=0.0,yerr=[rates_m05[5,0]*years]*len(sigma),uplims=True,capsize=2,fmt='b.')
axs[0].plot(sigma,astro_full_05,'b-',label=r'ABH full voids')
axs[0].errorbar(sigma,astro_full_05,xerr=0.0,yerr=[rates_m05[3,0]*years]*len(sigma),uplims=True,capsize=2,fmt='b.')
axs[0].plot(sigma+0.02,gal_frac_05,color='cyan',ls='--',label=r'ABH Obs void frac')
axs[0].errorbar(sigma+0.02,gal_frac_05,xerr=0.0,yerr=[rates_m05[7,0]*years]*len(sigma),uplims=True,capsize=2,marker='.',ls='--',fmt='cyan')
axs[0].plot(rates_m05[0],rates_m05[10]*years,'r--',label=r'PBH void frac')
axs[0].errorbar(rates_m05[0],rates_m05[10]*years,xerr=0.0,yerr=rates_m05[11]*years,fmt='r.')
axs[0].plot(rates_m05[0],rates_m05[8]*years,'r-',label=r'PBH full voids')
axs[0].errorbar(rates_m05[0],rates_m05[8]*years,xerr=0.0,yerr=rates_m05[9]*years,fmt='r.')
#axs[1]: Mmin = 1
axs[1].plot(sigma,hod_frac_1,'b--',label=r'ABH HOD void frac')
axs[1].errorbar(sigma,hod_frac_1,xerr=0.0,yerr=[rates_m1[5,0]*years]*len(sigma),uplims=True,capsize=2,fmt='b.')
axs[1].plot(sigma,astro_full_1,'b-',label=r'ABH full voids')
axs[1].errorbar(sigma,astro_full_1,xerr=0.0,yerr=[rates_m1[3,0]*years]*len(sigma),uplims=True,capsize=2,fmt='b.')
axs[1].plot(sigma+0.02,gal_frac_1,color='cyan',ls='--',label=r'ABH Obs void frac')
axs[1].errorbar(sigma+0.02,gal_frac_1,xerr=0.0,yerr=[rates_m1[7,0]*years]*len(sigma),uplims=True,capsize=2,marker='.',ls='--',fmt='cyan')
axs[1].plot(rates_m1[0],rates_m1[10]*years,'r--',label=r'PBH void frac')
axs[1].errorbar(rates_m1[0],rates_m1[10]*years,xerr=0.0,yerr=rates_m1[11]*years,fmt='r.')
axs[1].plot(rates_m1[0],rates_m1[8]*years,'r-',label=r'PBH full voids')
axs[1].errorbar(rates_m1[0],rates_m1[8]*years,xerr=0.0,yerr=rates_m1[9]*years,fmt='r.')

#axs[0].set_title(r"$M_{min} = 0.5 M_\odot$, $R_{min} = 49.8$ Mpc",fontsize=14)
#axs[1].set_title(r"$M_{min} = 1 M_\odot$, $R_{min} = 25.7$ Mpc",fontsize=14)
axs[0].text(0.3,0.022,r"$M_{min} = 0.5 M_\odot$",fontsize=14)
axs[0].text(0.3,0.013,r"$R_{min} = 49.8$ Mpc",fontsize=14)
axs[1].text(0.3,0.022,r"$M_{min} = 1 M_\odot$",fontsize=14)
axs[1].text(0.3,0.013,r"$R_{min} = 25.7$ Mpc",fontsize=14)


plot_upper_limit(ax, cdf_gal[0], merger_rate_bns_hod_voids+err_merger_rate_bns_hod_voids[1], 'green', r'BNS sim')

for ax in axs:
	ax.set_yscale('log')
	ax.set_xlabel(r'$\sigma$',fontsize = 16)
	ax.set_xlim(0.27,1.23)
	ax.set_ylim(0.01,9)

# Secondary DESI-like y axes
secax0 = axs[0].secondary_yaxis('right',functions=(lambda y: y * yfactor,lambda y: y / yfactor))
secax1 = axs[1].secondary_yaxis('right',functions=(lambda y: y * yfactor,lambda y: y / yfactor))

# Ticks on ax[0], but no tick labels
secax0.tick_params(axis='y',which='both',labelright=False)
# Ticks + labels on ax[1]
secax1.set_ylabel(r'DESI-like \\ \# event / 10yr',fontsize=14,rotation=270,labelpad=20)

# Use a single shared legend, since both panels show the same line types.
# The handles/labels are taken from axs[0] and the legend is placed
# to the right of axs[1].
handles, labels = axs[0].get_legend_handles_labels()
axs[1].legend(
	handles,
	labels,
	loc='upper center',
	bbox_to_anchor=(-0.2, 1.2),
	ncol=3,
	frameon=False,
	prop={'size':12}
)

plt.savefig("Number_of_events.pdf",bbox_inches='tight')
plt.show()
