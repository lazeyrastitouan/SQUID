import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import argparse
from void_factors import*
from rate_pfav import rate_pfav
from plot_fancy import plot_upper_limit
from plot_fancy import add_y_jump
from matplotlib.patches import Patch

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

matplotlib.rc('xtick', labelsize=16)
matplotlib.rc('ytick', labelsize=16)

config            = 'ET_2L'
mmin              = 1
Rvs               = [85.01,65.01,65.01]
loc               = 29.3
nr                = len(Rvs)
sigmas            = [0.3,0.6,0.8]
ns                = len(sigmas)
gal_void_frac     = 0.09
err_gal_void_frac = np.array([0.02,0.035]) #8046./120606. #~7% corresponding to Pan et al. SDSSDR7 (11% for other sample)
cdm_void_frac     = 0.15
err_cdm_void_frac = np.array([0.05,0.05])

#prepare figure
fig, axs  = plt.subplots(1,ns,figsize=(12,5), squeeze=False, layout='constrained') #squeeze=False so indexing works even for 1 panel, constrained so panels and suptitle don't overlap
title = r"ET 2L, $z=0.3$, $M_{\rm{min}}$ = "+str(mmin)+'$M_\odot$' #shown as the legend title (a suptitle would overlap the 'outside upper' legend)
#broken y axis for selected panels (column index -> limits). Visible: [ymin,ylow] and [yhigh,ymax]; ylow..yhigh is cut out.
y_jumps = {1: dict(ymin=1e-13, ylow=1e-10, yhigh=1e-6, ymax=2e-2, gap=0.5),
           2: dict(ymin=5e-30, ylow=1e-26, yhigh=1e-6, ymax=2e-2, gap=0.5)}


for ir, Rv in enumerate(Rvs):
	sigma = sigmas[ir]
	#void filling factor, same for all components
	fill     = np.load('volume_filling_observed_z03.npz')
	mask     = any(isclose(fill['R'][:,None],Rv,atol=0.003),axis=1)
	if len(fill['R'][mask]) > 1:
		exit("Your tolerance for 'R' in the filling factor is too large and you have more 'R' values than void radii)")
	vfill    = fill['fv_of_R'][mask]
	errvfill = fill['err_fv_of_R'][mask]

	#CDFS
	cdf_gal      = np.loadtxt('CDFs/Rv_CDF_'+str(float(Rv))+'_HamDM_SutDM_SutH_RicDM_RicGal.txt',usecols=(0,3)).T
	cdf_gal_obs  = np.loadtxt("CDFs/Rv_CDF_"+str(float(Rv))+"_HamDM_SutDM_SutH_RicDM_RicGal.txt",usecols=(0,5)).T
	cdf_dm       = np.loadtxt("CDFs/Rv_CDF_"+str(float(Rv))+"_HamDM_SutDM_SutH_RicDM_RicGal.txt",usecols=(0,1)).T

	#--- Merger rates*Pfav (i.e only localisation criterion)
	#ABH #even for astro merger rates depend on mmin and sigma as changes Mmax
	merger_rate_abh_all     = rate_pfav(sigma=sigma, m_min=mmin, network=config, channel='BBH',rate_case='median')(Rv)
	up                      = rate_pfav(sigma=sigma, m_min=mmin, network=config, channel='BBH',rate_case='hi')(Rv)
	down                    = rate_pfav(sigma=sigma, m_min=mmin, network=config, channel='BBH',rate_case='lo')(Rv)
	err_merger_rate_abh_all = np.array([merger_rate_abh_all-down,up-merger_rate_abh_all]) 

	#BNS
	up                      = rate_pfav(sigma=sigma, m_min=mmin, network=config, channel='BNS',rate_case='hi')(Rv)
	down                    = rate_pfav(sigma=sigma, m_min=mmin, network=config, channel='BNS',rate_case='lo')(Rv)
	merger_rate_bns_all     = (up+down)/2 #no median from catalogs so estimate it like this
	err_merger_rate_bns_all = np.array([merger_rate_bns_all-down,up-merger_rate_bns_all])

	#EARLY PBH
	merger_rate_pbh_all     = rate_pfav(sigma=sigma, m_min=mmin, network=config, channel='PBH')(Rv)
		
	
	#--- multiply by voids factors
	#ABH sims
	merger_rate_abh_hod_voids     = RinV(merger_rate_abh_all, gal_void_frac, cdf_gal[1], vfill)
	err_merger_rate_abh_hod_voids = err_RinV(merger_rate_abh_all, gal_void_frac,cdf_gal[1], vfill, err_merger_rate_abh_all, err_gal_void_frac, errvfill)
	#ABH obs
	merger_rate_abh_obs_voids     = RinV(merger_rate_abh_all, gal_void_frac, cdf_gal_obs[1], vfill)
	err_merger_rate_abh_obs_voids = err_RinV(merger_rate_abh_all, gal_void_frac, cdf_gal_obs[1], vfill, err_merger_rate_abh_all, err_gal_void_frac, errvfill)

	#BNS sims
	merger_rate_bns_hod_voids     = RinV(merger_rate_bns_all, gal_void_frac, cdf_gal[1], vfill)
	err_merger_rate_bns_hod_voids = err_RinV(merger_rate_bns_all, gal_void_frac, cdf_gal[1], vfill,err_merger_rate_bns_all, err_gal_void_frac, errvfill)
	#BNS obs
	merger_rate_bns_obs_voids     = RinV(merger_rate_bns_all, gal_void_frac, cdf_gal_obs[1], vfill)
	err_merger_rate_bns_obs_voids = err_RinV(merger_rate_bns_all, gal_void_frac, cdf_gal_obs[1], vfill, err_merger_rate_bns_all, err_gal_void_frac, errvfill)

	#EARLY PBH
	merger_rate_pbh_voids         = RinV(merger_rate_pbh_all, cdm_void_frac, cdf_dm[1], vfill)
	err_merger_rate_pbh_voids     = err_RinV(merger_rate_pbh_all, cdm_void_frac, cdf_dm[1], vfill, np.array([0,0]), err_cdm_void_frac, errvfill)


	#--- PLOT
	ax = axs[0,ir]
	#Early time pbh
	ax.fill_between(cdf_dm[0], merger_rate_pbh_voids+err_merger_rate_pbh_voids[1], merger_rate_pbh_voids-err_merger_rate_pbh_voids[0], color='red', alpha=0.5, label=r'PBH')
	iloc = find_nearest(loc/Rv,cdf_dm[0])
	ax.plot([loc/Rv,loc/Rv], [merger_rate_pbh_voids[iloc]+err_merger_rate_pbh_voids[1,iloc], merger_rate_pbh_voids[iloc]-err_merger_rate_pbh_voids[0,iloc]],'k-')

	#ABH sims
	ax.fill_between(cdf_gal[0], merger_rate_abh_hod_voids+err_merger_rate_abh_hod_voids[1], merger_rate_abh_hod_voids-err_merger_rate_abh_hod_voids[0], color='blue', alpha=0.5, label=r'ABH sim')#, Sutter+ 2014')
	#ABH obs
	ax.fill_between(cdf_gal_obs[0], merger_rate_abh_obs_voids+err_merger_rate_abh_obs_voids[1], merger_rate_abh_obs_voids-err_merger_rate_abh_obs_voids[0], color='cornflowerblue', alpha=0.5, label=r'ABH obs')#, Ricciardelli+ 2014')
	
	#BNS sims
	plot_upper_limit(ax, cdf_gal[0], merger_rate_bns_hod_voids+err_merger_rate_bns_hod_voids[1], 'green', r'BNS sim')
	#ax.fill_between(cdf_gal[0], merger_rate_bns_hod_voids+err_merger_rate_bns_hod_voids[1], merger_rate_bns_hod_voids-err_merger_rate_bns_hod_voids[0], color='green', alpha=0.5, label=r'BNS sim')
	#BNS obs
	plot_upper_limit(ax, cdf_gal_obs[0], merger_rate_bns_obs_voids+err_merger_rate_bns_obs_voids[1], 'mediumseagreen', r'BNS obs')
	#ax.fill_between(cdf_gal_obs[0], merger_rate_bns_obs_voids+err_merger_rate_bns_obs_voids[1], merger_rate_bns_obs_voids-err_merger_rate_bns_obs_voids[0], color='mediumseagreen', alpha=0.5, label=r'BNS obs')
		
	#ymin, ymax = axs[ir,ids].get_ylim()
	ax.text(0.55,0.075,r'$\sigma=$ '+str(sigma),transform=ax.transAxes,fontsize=13,va='top', ha='left') #puts the text at top left regardless of the axis
	ax.text(0.55,0.15,r'$R_v=$ '+str(round(Rv))+' Mpc/$h$',transform=ax.transAxes,fontsize=13,va='top', ha='left')
	#ax.text(0.05,0.77,r'$z = 0.3$',transform=ax.transAxes,fontsize=13,va='top', ha='left')
	ax.set_xlim(0,1)
	if ir in y_jumps:
		add_y_jump(ax, **y_jumps[ir])
	else:
		ax.set_yscale('log')
		ax.set_ylim(1e-6,2e-2)
	#ax.legend(loc='lower right',prop={'size':12},frameon=False,ncols=2)
	#axis labels: x only on the last row, y only on the first column
	ax.set_xlabel(r'$r/R_v$', fontsize=18)
	if ir == 0:
		ax.set_ylabel(r'Merger rate [yr$^{-1}$Gpc$^{-3}]$',fontsize=18,labelpad=15)

#--- Single legend for all panels, 3 columns: PBH | ABH | BNS
#matplotlib fills legend columns top to bottom, so PBH gets an invisible dummy entry
#below it to keep ABH (sim, obs) and BNS (sim, obs) in their own columns.
h, l  = axs[0,0].get_legend_handles_labels()
hd    = dict(zip(l, h))
order = [('PBH', hd['PBH']), ('ABH sim', hd['ABH sim']), ('ABH obs', hd['ABH obs']), ('BNS sim', hd['BNS sim']), ('BNS obs', hd['BNS obs'])]
fig.legend([o[1] for o in order], [o[0] for o in order], loc='outside upper center',ncols=5, frameon=False, fontsize=14,title=title, title_fontsize=18)


fout='Merger_rates_observedVSF_z0p3_ABH_PBH_BNS_'+config+'_Mmin'+str(mmin)+'_selectedpanels.pdf'
print("Saving figure to " + fout)
plt.savefig(fout,bbox_inches='tight')
plt.show()
