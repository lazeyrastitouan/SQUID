from numpy import*
import matplotlib.pyplot as plt
import matplotlib
import argparse
import re, pathlib

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

matplotlib.rc('xtick', labelsize=14)
matplotlib.rc('ytick', labelsize=14)


def RinV(R,fv,cdf,Vfill):
	'''
	Multiply total merger rate R by void variables: fraction of component in voids fv, cumulative density function normalized to 1 at R_void, void filling factor
	'''
	return R * fv * cdf * Vfill

def err_RinV(R,fv,cdf,Vfill,errR,errfv,errVfill):
        '''
        Error propagation for RinV (2D as asymetrical error)
        errR, errfv have to be 2D arrays
        '''
        sRinV_lo = sqrt((errR[0]/R)**2.+(errfv[0]/fv)**2.+(errVfill/Vfill)**2.) * RinV(R,fv,cdf,Vfill)
        sRinV_hi = sqrt((errR[1]/R)**2.+(errfv[1]/fv)**2.+(errVfill/Vfill)**2.) * RinV(R,fv,cdf,Vfill)

        return vstack((sRinV_lo,sRinV_hi))

def find_nearest(array, value):
    idx = (abs(array - value)).argmin()
    return idx

### ---------- SETUP ---------- ### 
#Survey volume BOSS DR12 LOWZ 2.3 Gpc^3
vol_survey = 3.61 #Gpc^3 BOSSDR11 LOWZ (Nadathur 2016)
vol_tot    = 19.53 #Gpc^3 Full sky volume for the redshift slice corresponding to DR11 LOWZ
#void filling factor, same for all components
fill      = load('volume_filling_observed_z03.npz')
vfill     = fill['fv_of_R']
err_vfill = fill['err_fv_of_R'] 
Rv        = fill['R']
rs        = len(Rv)
Rmins     = [49.8,25.7] #minimum size for localisation 0.5-0.5, 1-1 respectively

#ASTRO
N_gal_all         = 552.04*0.003 #events/year 
err_N_gal_all     = array([13.4,22.4])*0.003 #-,+ 
gal_void_frac     = 0.09
err_gal_void_frac = array([0.02,0.035]) #8046./120606. #~7% corresponding to Pan et al. SDSSDR7 (11% for other sample)

#For arly time pbh
sigma             = [1.2,1.0,0.8,0.6,0.3]
swanted           = [0,1,2,3,4]
ls                = len(swanted)
#masses            = ['1-1','5-5','10-10']
cdm_void_frac     = 0.15
err_cmd_void_frac = array([0.05,0.05])

#Since the CDF might not have been computed for exactly the same Rs as the VSF, need to find the closest R value of CDF (assum CDF don't vary much with R)
cdfpath     = pathlib.Path("/Users/titouanlazeyras/Desktop/BICOCCA/SQUID/CDFs/")
R_available = []
pattern     = re.compile("Rv_CDF_([0-9.]+)_HamDM")
for file in cdfpath.glob('*.txt'):
	match = pattern.search(file.name)
	if match:
		R_available.append(float(match.group(1)))
R_available = array(R_available)

### ---------- INTEGRATION OVER Rv ---------- ###
#i)  Take full voids, CDF=1
#ii) Take fraction of voids such that f*Rv = Rmin cMpc
for Rmin in Rmins:
	if Rmin == 49.8:
		Mmin = 0.5
		N_early_pbh = genfromtxt("Nevent0pt5_z015-043.txt")[:,1]*0.01 #rows sigma, cols mass. We load all sigmas but only mergers for 1-1 Msun
	if Rmin == 25.7:
		Mmin = 1
		N_early_pbh = genfromtxt("Nevent1_z015-043.txt")[:,1]*0.01 #rows sigma, cols mass

	N_hod_full, err_N_hod_full, N_hod_frac, err_N_hod_frac = [],[],[],[] #astro sims
	N_gal_full, err_N_gal_full, N_gal_frac, err_N_gal_frac = [],[],[],[] #astro obs
	N_pbh_full, err_N_pbh_full, N_pbh_frac, err_N_pbh_frac = [[] for x in range(ls)],[[] for x in range(ls)],[[] for x in range(ls)],[[] for x in range(ls)] #pbh early down to 0.1
	for i in range(rs):
		if Rv[i]< Rmin:
			continue
	
		#i)  Take full voids, CDF=1
		N_hod_full.append(RinV(N_gal_all, gal_void_frac, 1, vfill[i]))
		err_N_hod_full.append(err_RinV(N_gal_all, gal_void_frac, 1, vfill[i], err_N_gal_all, err_gal_void_frac, err_vfill[i]))
		N_gal_full.append(RinV(N_gal_all, gal_void_frac, 1, vfill[i]))
		err_N_gal_full.append(err_RinV(N_gal_all, gal_void_frac, 1, vfill[i], err_N_gal_all, err_gal_void_frac, err_vfill[i]))	
		for s in swanted:
			si = swanted.index(s)
			N_pbh_full[si].append(RinV(N_early_pbh[s], cdm_void_frac, 1, vfill[i]))
			err_N_pbh_full[si].append(err_RinV(N_early_pbh[s], cdm_void_frac, 1, vfill[i], array([0,0]), err_cmd_void_frac, err_vfill[i])) #no error on N_pbh calculated analytically by Pippa

		#ii) Take fraction of voids such that f*Rv = Rmin cMpc 
		#CDFs computed more sparsely than Vfill, lets find the CDF that maches Rv[i]
		idcdf     = find_nearest(R_available,Rv[i])
		cdfs      = loadtxt('CDFs/Rv_CDF_'+str(R_available[idcdf])+'_HamDM_SutDM_SutH_RicDM_RicGal.txt',usecols=(0,1,3,5)).T #Rfrac, DM, gal_hod, gal_obs
		#find fraction at which Rv*f = Rmin
		idx       = find_nearest(cdfs[0]*Rv[i],Rmin)

		N_hod_frac.append(RinV(N_gal_all, gal_void_frac, cdfs[2,idx], vfill[i]))
		err_N_hod_frac.append(err_RinV(N_gal_all, gal_void_frac, cdfs[2,idx], vfill[i], err_N_gal_all, err_gal_void_frac, err_vfill[i]))
		N_gal_frac.append(RinV(N_gal_all, gal_void_frac, cdfs[3,idx], vfill[i]))
		err_N_gal_frac.append(err_RinV(N_gal_all, gal_void_frac, cdfs[3,idx], vfill[i], err_N_gal_all, err_gal_void_frac, err_vfill[i]))
		for s in swanted:
			si = swanted.index(s)
			N_pbh_frac[si].append(RinV(N_early_pbh[s], cdm_void_frac, cdfs[1,idx], vfill[i]))
			err_N_pbh_frac[si].append(err_RinV(N_early_pbh[s], cdm_void_frac, cdfs[1,idx], vfill[i], array([0,0]), err_cmd_void_frac, err_vfill[i]))

	#trick needed manipulation to get correct shape
	err_N_hod_full = array(err_N_hod_full).reshape(shape(err_N_hod_full)[0],shape(err_N_hod_full)[1])	
	err_N_hod_frac = array(err_N_hod_frac).reshape(shape(err_N_hod_frac)[0],shape(err_N_hod_frac)[1])
	err_N_gal_frac = array(err_N_gal_frac).reshape(shape(err_N_gal_frac)[0],shape(err_N_gal_frac)[1])
	err_N_pbh_full = array(err_N_pbh_full).reshape(shape(err_N_pbh_full)[0],shape(err_N_pbh_full)[1],shape(err_N_pbh_full)[2])    
	err_N_pbh_frac = array(err_N_pbh_frac).reshape(shape(err_N_pbh_frac)[0],shape(err_N_pbh_frac)[1],shape(err_N_pbh_frac)[2])

	#integral
	filter = Rv>=Rmin
	f = open("sigma_Nhod_Ngal_NpbhMmin"+str(Mmin)+"_Rmin"+str(Rmin)+".txt",'w')
	f.write("#sigma\tNhod full\tNgal full\terr(Nfull)\tNhod frac\terr(Nhodfrac)\tNgalfrac\terr(Ngalfrac)\tNpbhfull\terr(Npbhfull)\tNpbhfrac\terr(Npbhfrac)\n")
	n_events_hod_full = trapz(N_hod_full,Rv[filter])*vol_survey/vol_tot
	err_n_full        = trapz(err_N_hod_full[:,1],Rv[filter])*vol_survey/vol_tot #Same for hod and gal if full. We consider only the upper errorbar to plot upper limits
	n_events_gal_full = trapz(N_gal_full,Rv[filter])*vol_survey/vol_tot
	n_events_hod_frac = trapz(N_hod_frac,Rv[filter])*vol_survey/vol_tot
	err_n_hod_frac    = trapz(err_N_hod_frac[:,1],Rv[filter])*vol_survey/vol_tot
	n_events_gal_frac = trapz(N_gal_frac,Rv[filter])*vol_survey/vol_tot
	err_n_gal_frac    = trapz(err_N_gal_frac[:,1],Rv[filter])*vol_survey/vol_tot
	for i in range(ls):
		n_events_pbh_full = trapz(N_pbh_full[i],Rv[filter])*vol_survey/vol_tot
		n_events_pbh_frac = trapz(N_pbh_frac[i],Rv[filter])*vol_survey/vol_tot
		err_n_pbh_full    = trapz(err_N_pbh_full[i,:,1],Rv[filter])*vol_survey/vol_tot #These are symetric by construction
		err_n_pbh_frac    = trapz(err_N_pbh_frac[i,:,1],Rv[filter])*vol_survey/vol_tot
		f.write(str(sigma[i])+"\t"+str(n_events_hod_full)+"\t"+str(n_events_gal_full)+"\t"+str(err_n_full)+"\t"+str(n_events_hod_frac)+"\t"+str(err_n_hod_frac)+"\t"+str(n_events_gal_frac)+"\t"+str(err_n_gal_frac)+"\t"+str(n_events_pbh_full)+"\t"+str(err_n_pbh_full)+"\t"+str(n_events_pbh_frac)+"\t"+str(err_n_pbh_frac)+"\n")

	f.close()
