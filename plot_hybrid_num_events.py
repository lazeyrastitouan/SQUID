from numpy import*
import matplotlib.pyplot as plt
import matplotlib

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

matplotlib.rc('xtick', labelsize=16)
matplotlib.rc('ytick', labelsize=16)

A = True
B = True

if A:
	events05 = loadtxt("hybrid_number_events_of_R_0.5Msol.txt").T
	events1  = loadtxt("hybrid_number_events_of_R_1.0Msol.txt").T
	
	plt.figure(1)
	plt.yscale('log')
	plt.plot(events05[0],events05[1],color='blue',label=r'Astro sims')
	plt.errorbar(events05[0],events05[1],xerr=0.0,yerr=[events05[3],events05[2]],marker='.',fmt='blue')
	plt.plot(events05[0],events05[4],color='cyan',label=r'Astro obs')
	plt.errorbar(events05[0],events05[4],xerr=0.0,yerr=[events05[5],events05[6]],marker='.',fmt='cyan')
	plt.plot(events05[0],events05[7],color='red',label=r'PBH')
	plt.errorbar(events05[0],events05[7],xerr=0.0,yerr=events05[8],marker='.',fmt='red')
	plt.plot(events1[0],events1[1],color='blue',ls='--')
	plt.errorbar(events1[0],events1[1],xerr=0.0,yerr=[events1[3],events1[2]],ls='--',marker='.',fmt='blue')
	plt.plot(events1[0],events1[4],color='cyan',ls='--')
	plt.errorbar(events1[0],events1[4],xerr=0.0,yerr=[events1[5],events1[6]],ls='--',marker='.',fmt='cyan')
	plt.plot(events1[0],events1[7],color='red',ls='--')
	plt.errorbar(events1[0],events1[7],xerr=0.0,yerr=events1[8],ls='--',marker='.',fmt='red')
	plt.plot([1,2],[-1,-1],'k--',label=r'$R_v >$ loc(1$M_\odot$)') 
	plt.plot([1,2],[-1,-1],'k-',label=r'$R_v >$ loc(0.5$M_\odot$)') 
	plt.xlim(25,95)
	plt.ylim(1e-2,10)
	plt.legend(prop={'size':14},frameon=False)
	plt.xlabel(r"$\geq$ Scale [Mpc/h]",fontsize=16)
	plt.ylabel(r"\# events/10 yr",fontsize=16)
	plt.savefig("hybrid_number_of_events_scale_dep.pdf",bbox_inches='tight')
	plt.show()

if B:
	
	events1_30 = load("hybrid_number_events_radial_profile_Rv30_1.0Msol.npz")
	events1    = load("hybrid_number_events_radial_profile_Rv60_1.0Msol.npz")
	events05   = load("hybrid_number_events_radial_profile_Rv60_0.5Msol.npz")
	
	fig, axs  = plt.subplots(1,3,figsize=(16,4),sharey=True)
	plt.subplots_adjust(wspace=0.09)
	axs[0].set_ylabel(r'events / 10yr',fontsize=14)
	
	axs[0].fill_between(events1_30['r'],events1_30['Nhod']-events1_30['err_Nhod'][0],events1_30['Nhod']+events1_30['err_Nhod'][1],color='blue',alpha=0.5)
	axs[0].fill_between(events1_30['r'],events1_30['Ngal']-events1_30['err_Ngal'][0],events1_30['Ngal']+events1_30['err_Ngal'][1],color='cyan',alpha=0.5)
	axs[0].fill_between(events1_30['r'],events1_30['Npbh']-events1_30['err_Npbh'][0],events1_30['Npbh']+events1_30['err_Npbh'][1],color='red',alpha=0.5)
	axs[0].plot([25.7/30,25.7/30],[1e-4,2],'k:')
	axs[0].text(0.03,0.6,r'$R_v$ = 30 Mpc/h',fontsize=13)
	axs[0].text(0.03,0.35,r'$M_{min}$ = 1$M_{\odot}$',fontsize=13)
	
	axs[1].fill_between(events1['r'],events1['Nhod']-events1['err_Nhod'][0],events1['Nhod']+events1['err_Nhod'][1],color='blue',alpha=0.5)
	axs[1].fill_between(events1['r'],events1['Ngal']-events1['err_Ngal'][0],events1['Ngal']+events1['err_Ngal'][1],color='cyan',alpha=0.5)
	axs[1].fill_between(events1['r'],events1['Npbh']-events1['err_Npbh'][0],events1['Npbh']+events1['err_Npbh'][1],color='red',alpha=0.5)
	axs[1].plot([25.7/60,25.7/60],[1e-4,2],'k:')
	axs[1].text(0.03,0.6,r'$R_v$ = 60 Mpc/h',fontsize=13)
	axs[1].text(0.03,0.35,r'$M_{min}$ = 1$M_{\odot}$',fontsize=13)
	
	axs[2].fill_between(events05['r'],events05['Nhod']-events05['err_Nhod'][0],events05['Nhod']+events05['err_Nhod'][1],color='blue',alpha=0.5,label=r'Astro sims')
	axs[2].fill_between(events05['r'],events05['Ngal']-events05['err_Ngal'][0],events05['Ngal']+events05['err_Ngal'][1],color='cyan',alpha=0.5,label=r'Astro obs')
	axs[2].fill_between(events05['r'],events05['Npbh']-events05['err_Npbh'][0],events05['Npbh']+events05['err_Npbh'][1],color='red',alpha=0.5,label=r'PBH')
	axs[2].plot([49.8/60,49.8/60],[1e-4,2],'k:')
	axs[2].text(0.03,0.6,r'$R_v$ = 60 Mpc/h',fontsize=13)
	axs[2].text(0.03,0.35,r'$M_{min}$ = 0.5$M_{\odot}$',fontsize=13)
	
	for ax in axs:
		ax.set_yscale('log')
		ax.set_xlabel(r'$r/R_v$',fontsize=14)
		ax.set_xlim(0,1)
		ax.set_ylim(2e-4,1.0)
	axs[2].legend(loc='lower left',frameon=False,prop={'size':13})
	plt.savefig("rprofile_hybrid_number_of_events_fixedRv.pdf",bbox_inches='tight')
	plt.show()

	
