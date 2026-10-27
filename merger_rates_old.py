from numpy import*
import matplotlib.pyplot as plt
import matplotlib
import argparse
from void_factors import* 

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

matplotlib.rc('xtick', labelsize=14)
matplotlib.rc('ytick', labelsize=14)

parser = argparse.ArgumentParser()
parser.add_argument('mmin',nargs=1,type=float, help='mmin')
parser.add_argument('sid',nargs=1,type=int, help='sigma index')
args     = parser.parse_args()

#void filling factor, same for all components
Rv       = array([60])
rs       = len(Rv)
fill     = load('volume_filling_observed_z03.npz')
mask     = any(isclose(fill['R'][:,None],Rv[None,:],atol=0.003),axis=1)
if len(fill['R'][mask]) > len(Rv):
	exit("Your tolerance for 'R' in the filling factor is too large and you have more 'R' values than void radii)")
else:
	print(fill['R'][mask])
vfill    = fill['fv_of_R'][mask]
errvfill = fill['err_fv_of_R'][mask]

#ASTRO
merger_rate_gal_all       = 35.1*0.003 #46*0.003 #100. #mergers/year/Gpc^-3 z=0.3
err_merger_rate_gal_all   = array([11.8,19.3])*0.003 
gal_void_frac             = 0.09
err_gal_void_frac         = array([0.02,0.035]) #8046./120606. #~7% corresponding to Pan et al. SDSSDR7 (11% for other sample)

#PBH
#pbh fakhry
merger_rate_pbh = 5e-5*1.5 #Optimistic case, assume 2 1Msun pbh (factor 1.5 z dependence from z=0 to 0.3) Einasto profile, f=0.01

#Pippa, late time mergers
Rlate70 = [8.4e-4,8.4e-4,9.1e-4]
Rlate55 = [8e-4,8e-4,8.7e-4]

#Early time pbh
mmin              = args.mmin[0] 
ind_s             = args.sid[0]

im = int(mmin)
if mmin == 0.5:
	merger_early_pbh  = genfromtxt("REBM0pt5_z015-043.txt")[:,1:4]*0.01 #rows sigma, cols mass 
elif mmin == 1:
	
	merger_early_pbh  = genfromtxt("REBM1_z015-043.txt")[:,1:4]*0.01 #rows sigma, cols mass
sigma             = [1.2,1,0.8,0.6,0.3]
sigma_str         = ['1p2','1','0p8','0p6','0p3']
masses            = ['1-1','5-5','10-10']
cdm_void_frac     = 0.15
err_cdm_void_frac = array([0.05,0.05])

#FIGURE
fig, axs = plt.subplots(1,rs)#,figsize=(12,4))
if rs == 1:
	axs = [axs]
axs[0].set_ylabel(r'Merger rate [yr$^{-1}$Gpc$^{-3}]$',fontsize=14)
#colors = [['darkblue','cornflowerblue','teal'],['orchid','magenta','purple'],['darkgreen','green','lime'],['green','lime','orange'],['red','darkred','khaki']]
colors = [['orange','salmon','tomato','red','darkred'],['orange','salmon','tomato','red','darkred'],['darkgreen','lime','khaki','teal','green']]


### ---------- Rv LOOP ---------- ###
for i in range(rs):
	#figure
	axs[i].set_yscale('log')
	if mmin == 0.5:
		axs[i].set_title(r'$0.15 \leq z \leq  0.43$, $M_{peak} = 5 M_\odot$, $M_{min}= 0.5  M_\odot$',fontsize=14)
		axs[i].text(0.05,8e-3,r'$R_v=$ '+str(Rv[i])+' Mpc/$h$',fontsize=13)
		axs[i].set_ylim(4e-7,2e-2) #0.5 Msun
	elif mmin == 1:
		axs[i].set_title(r'$0.15 \leq z \leq  0.43$, $M_{peak} = 5 M_\odot$, $M_{min}= 1  M_\odot$',fontsize=14)
		axs[i].text(0.05,4e-3,r'$R_v=$ '+str(Rv[i])+' Mpc/$h$',fontsize=13)
		axs[i].set_ylim(1e-6,1e-2) #1 Msun
	axs[i].set_xlabel(r'$r/R_v$', fontsize=14)
	axs[i].set_xlim(0,1)
	if mmin==0.5:
		axs[i].plot([49.8/Rv[i],49.8/Rv[i]],[1e-8,2],'k:')
	elif mmin == 1:
		axs[i].plot([25.7/Rv[i],25.7/Rv[i]],[1e-8,2],'k:')                
	#plot Fakhry
	#axs[i].plot([0,1],[merger_rate_pbh,merger_rate_pbh],color='yellow', ls=':',label=r'Fakhry+ 2024, $f_{DM}=0.01$')

	#pbh late
	#if Rv[i]==60:
#		axs[i].plot([0,0.3333333,1],Rlate55,color='gold',ls='--',marker='',label=r'PBH late channel 1-1M$_\odot$')
#	elif Rv[i]==70:
#		axs[i].plot([0,0.3333333,1],Rlate70,color='gold',ls='--',marker='',label=r'PBH late channel 1-1M$_\odot$')

	#Early time pbh
	cdf_dm                    = loadtxt("Rv_CDF_"+str(Rv[i])+"_HamDM_SutDM_SutH_RicDM_RicGal.txt",usecols=(0,1)).T
	for s in [ind_s]:
		for m in [1]:
			merger_rate_pbh_voids     = RinV(merger_early_pbh[s,m],cdm_void_frac,cdf_dm[1],vfill[i])
			err_merger_rate_pbh_voids = err_RinV(merger_early_pbh[s,m],cdm_void_frac,cdf_dm[1],vfill[i],array([0,0]),err_cdm_void_frac,errvfill[i])

			axs[i].fill_between(cdf_dm[0],merger_rate_pbh_voids+err_merger_rate_pbh_voids[1],merger_rate_pbh_voids-err_merger_rate_pbh_voids[0],color=colors[im][s],alpha=0.5,label=r'PBH, $\sigma =$ '+str(sigma[s]))


	#astro sims
	cdf_gal                   = loadtxt('Rv_CDF_'+str(Rv[i])+'_HamDM_SutDM_SutH_RicDM_RicGal.txt',usecols=(0,3)).T
	merger_rate_hod_voids     = RinV(merger_rate_gal_all,gal_void_frac,cdf_gal[1],vfill[i])
	err_merger_rate_hod_voids = err_RinV(merger_rate_gal_all,gal_void_frac,cdf_gal[1],vfill[i],err_merger_rate_gal_all,err_gal_void_frac,errvfill[i]) 
	#fout=open('Merger_rate_HOD_Sutter+2014_R'+str(Rv[i])+'.txt','w')
	#for i in range(len(cdf_gal[0])):
	#	fout.write(str(cdf_gal[0,i])+"\t"+str(merger_rate_hod_voids[i])+"\t"+str(err_merger_rate_hod_voids[0,i])+"\t"+str(err_merger_rate_hod_voids[1,i])+"\n")
	#fout.close()
	#axs[i].plot(cdf_gal[0],merger_rate_hod_voids,'b--',label=r'Astro, HOD, Sutter+ 2014')
	axs[i].fill_between(cdf_gal[0],merger_rate_hod_voids+err_merger_rate_hod_voids[1],merger_rate_hod_voids-merger_rate_hod_voids[0],color='blue',alpha=0.5,label=r'Astro HOD')#, Sutter+ 2014')

	#astro obs
	cdf_gal_obs               = loadtxt("Rv_CDF_"+str(Rv[i])+"_HamDM_SutDM_SutH_RicDM_RicGal.txt",usecols=(0,5)).T 
	merger_rate_gal_voids     = RinV(merger_rate_gal_all, gal_void_frac, cdf_gal_obs[1], vfill[i])
	err_merger_rate_gal_voids = err_RinV(merger_rate_gal_all, gal_void_frac, cdf_gal_obs[1], vfill[i], err_merger_rate_gal_all, err_gal_void_frac,errvfill[i]) 
	#axs[i].plot(cdf_gal_obs[0],merger_rate_gal_voids,'b--',label=r'Mergers from astrophysical sources')
	axs[i].fill_between(cdf_gal_obs[0],merger_rate_gal_voids+err_merger_rate_gal_voids[1],merger_rate_gal_voids-merger_rate_gal_voids[0],color='cornflowerblue',alpha=0.5,label=r'Astro obs')#, Ricciardelli+ 2014')


	axs[i].legend(loc='lower right',prop={'size':12},frameon=False,ncols=2)
if mmin == 0.5:
	fout='Merger_rate_observedVSF_z0p15-0p43_Mpeak5_0p5Msol_s'+sigma_str[ind_s]+'.pdf'
elif mmin == 1:
	fout='Merger_rate_observedVSF_z0p15-0p43_Mpeak5_1Msol_s'+sigma_str[ind_s]+'.pdf'
print("Saving figure to " + fout)
plt.savefig(fout,bbox_inches='tight')
#plt.show()


#------------------------------------------------
#------------------------------------------------
#Hybrid merger rate: fix \mathcal(R) to anything, and then multiply by i) fcdm CDF(rho_m) and ii) fgal CDF(rho_g)
#fig2, axs2 = plt.subplots(1,rs,figsize=(12,4))
#axs2[0].set_ylabel(r'Merger rate [yr$^{-1}$Gpc$^{-3}]$',fontsize=14)
#
#for i in range(rs):
#	cdf_dm                     = loadtxt("Rv_CDF_"+str(Rv[i])+"_HamDM_SutDM_SutH_RicDM_RicGal.txt",usecols=(0,1)).T
#	merger_rate_hybrid_dm      = RinV(merger_rate_gal_all,cdm_void_frac,cdf_dm[1],vfill[i])
#	err_merger_rate_hybrid_dm  = err_RinV(merger_rate_gal_all,cdm_void_frac,cdf_dm[1],vfill[i],err_merger_rate_gal_all,err_cdm_void_frac,errvfill[i])
#
#	cdf_gal                    = loadtxt('Rv_CDF_'+str(Rv[i])+'_HamDM_SutDM_SutH_RicDM_RicGal.txt',usecols=(0,3)).T
#	merger_rate_hybrid_gal     = RinV(merger_rate_gal_all,gal_void_frac,cdf_gal[1],vfill[i])
#	err_merger_rate_hybrid_gal = err_RinV(merger_rate_gal_all,gal_void_frac,cdf_gal[1],vfill[i],err_merger_rate_gal_all,err_gal_void_frac,errvfill[i])
#	
#	axs2[i].set_yscale('log')
#	
#	if i == 0:
#		axs2[i].plot([49.8/Rv[i],49.8/Rv[i]],[7e-7,2e-3],'k:', label=r'$0.5-0.5M_\odot$ loc.')
#		axs2[i].plot([25.7/Rv[i],25.7/Rv[i]],[7e-7,2e-3],'k--', label=r'$1-1M_\odot$ loc.')
#		axs2[i].fill_between(cdf_dm[0],merger_rate_hybrid_dm+err_merger_rate_hybrid_dm[1],merger_rate_hybrid_dm-merger_rate_hybrid_dm[0],color='blue',alpha=0.5)
#		axs2[i].fill_between(cdf_gal[0],merger_rate_hybrid_gal+err_merger_rate_hybrid_gal[1],merger_rate_hybrid_gal-merger_rate_hybrid_gal[0],color='red',alpha=0.5)
#	else:
#		axs2[i].plot([49.8/Rv[i],49.8/Rv[i]],[7e-7,2e-3],'k:')
#		axs2[i].plot([25.7/Rv[i],25.7/Rv[i]],[7e-7,2e-3],'k--')
#		axs2[i].fill_between(cdf_dm[0],merger_rate_hybrid_dm+err_merger_rate_hybrid_dm[1],merger_rate_hybrid_dm-merger_rate_hybrid_dm[0],color='blue',alpha=0.5,label=r'$f_{{dm}}\times$CDF$(\rho_{dm})$')
#		axs2[i].fill_between(cdf_gal[0],merger_rate_hybrid_gal+err_merger_rate_hybrid_gal[1],merger_rate_hybrid_gal-merger_rate_hybrid_gal[0],color='red',alpha=0.5,label=r'$f_{g}\times$CDF$(\rho_g)$')
#
#	axs2[i].text(0.015,5e-4,r'$R_v=$ '+str(Rv[i])+' Mpc/$h$',fontsize=13)
#	axs2[i].legend(loc='lower right',prop={'size':12},frameon=False,facecolor='white',framealpha=1,ncols=2)
#	axs2[i].set_ylim(2e-7,1e-3) 
#	axs2[i].set_xlim(0,1) 
#
#	
#fout2='Merger_rate_hybrid_observedVSF_z0p3.pdf'
#print("Saving figure to " + fout2)
#plt.savefig(fout2,bbox_inches='tight')
#plt.show()

