######
# GET NUMBER OF EVENT PER YEAR FROM VOIDS FOR A GIVEN Mmin. RATES COMPUTED EITHER AT FIXED z OR OVER A z SLICE
# ACCURATE INTEGRAL WITH ERRORS. CONSIDE ONLY UPPER LIMITS FOR ASTROPHYSICAL OBJECTS MERGER, ERRORBARS FOR PBH
#####


import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import argparse
import re, pathlib
from void_factors import*
import rate_pfav as rpf

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

matplotlib.rc('xtick', labelsize=14)
matplotlib.rc('ytick', labelsize=14)

def integrate_with_error(Rg, yg, Rv, t1, s1, t2, t3, s3):
	"""
	Exact result of the integral of R_mean * f * CDF * Vfill dRv	
	where R_mean is the mean rate over the redshift slice
	Arguments:
	Rv: voids radii values
	Rg, yg: the nodes and values of the rate used to interpolate between Rv values in rates_x_pfav.py 
	t1: cdm/gal frac, s1 = (s1_minus, s1_plus), both given as positive numbers its potentitally asymetric arror
	t2: CDF, no error
	t3: Vfill, s3 its error (symmetric)
    	Returns: 
	I, sI_minus, sI_plus.
	Usage:
	N, N_m, N_p = integrate_with_error(R.R_grid,R.values,Rv,void_frac,err_void_frac,CDF,Vfill,err_Vfill) 
	"""
	s1m, s1p = s1

	lo, hi = max(Rv[0], Rg[0]), Rv[-1]
	if hi <= lo:
		return 0.0, 0.0, 0.0

	i0, i1 = np.searchsorted(Rv, lo, 'right'), np.searchsorted(Rv, hi, 'left')
	j0, j1 = np.searchsorted(Rg, lo, 'right'), np.searchsorted(Rg, hi, 'left')
	knots = np.sort(np.concatenate(([lo], Rv[i0:i1], Rg[j0:j1], [hi])))

	g = np.interp(knots, Rg, yg, left=0.0, right=yg[-1])
	h = np.diff(knots)
	ga, gb = g[:-1], g[1:]

	c = np.zeros(len(knots))
	c[:-1] += h * (2*ga + gb) / 6.0
	c[1:]  += h * (ga + 2*gb) / 6.0

	n = len(Rv)
	idx = np.clip(np.searchsorted(Rv, knots, 'right') - 1, 0, n - 2)
	t = (knots - Rv[idx]) / (Rv[idx+1] - Rv[idx])
	w = (np.bincount(idx,     c * (1 - t), minlength=n)+ np.bincount(idx + 1, c * t,       minlength=n))

	J  = np.sum(w * t2 * t3)
	sJ = np.sqrt(np.sum((w * t2 * s3)**2))
	I  = t1 * J

	# I = t1 * J: asymmetric t1 error on each side, plus symmetric J error
	sIm = np.sqrt((J * s1m)**2 + (t1 * sJ)**2)
	sIp = np.sqrt((J * s1p)**2 + (t1 * sJ)**2)
	
	return I, sIm, sIp

### ---------- SETUP ---------- ###
parser = argparse.ArgumentParser()
parser.add_argument('-config', nargs='?', default='ET_2L', type=str, help='ET configuration')
parser.add_argument('-mmin', nargs='?', default=0.5, type=float, help='PBH mass function minimum mass')
args   = parser.parse_args()
config = args.config
mmin   = args.mmin

if mmin == 1:
	Rmin = 29.3
elif mmin == 0.5:
	Rmin = 49.8
else:
	exit("Invalid Mmin")

print("ET configuration: "+config)
print("Rmin "+str(Rmin)+", Mmin "+str(mmin))

#Survey volume BOSS DR12 LOWZ 2.3 Gpc^3
vol_survey = 3.61 #Gpc^3 BOSSDR11 LOWZ (Nadathur 2016)
vol_tot    = 19.53 #Gpc^3 Full sky volume for the redshift slice corresponding to DR11 LOWZ
#void filling factor, same for all components
fill      = np.load('volume_filling_observed_z03.npz')
#select relevant part for integral over Rv, i.e. where Rv>=Rmin
filter    = fill['R']>=Rmin
vfill     = fill['fv_of_R'][filter]
err_vfill = fill['err_fv_of_R'][filter] 
Rv        = fill['R'][filter]
rs        = len(Rv)

#ASTRO
gal_void_frac     = 0.09
err_gal_void_frac = np.array([0.02,0.035]) #8046./120606. #~7% corresponding to Pan et al. SDSSDR7 (11% for other sample)

#PBH
sigmas            = [1.2,0.8,0.6,0.3]
ls                = len(sigmas)
cdm_void_frac     = 0.15
err_cdm_void_frac = np.array([0.05,0.05])

#Since the CDF might not have been computed for exactly the same Rs as the VSF, need to find the closest R value of CDF (assum CDF don't vary much with R)
cdfpath     = pathlib.Path("/Users/titouanlazeyras/Desktop/BICOCCA/SQUID/CDFs/")
R_available = []
pattern     = re.compile("Rv_CDF_([0-9.]+)_HamDM")
for file in cdfpath.glob('*.txt'):
	match = pattern.search(file.name)
	if match:
		R_available.append(float(match.group(1)))
R_available = np.array(R_available)
#lets find for each Rv the CDF value at fraction of voif f such that f*Rv = Rmin and put it into a table
CDF_DM, CDF_gal_HOD, CDF_gal_obs = [], [], []
for i in range(rs):
	idcdf = find_nearest(R_available,Rv[i])
	cdfs  = np.loadtxt('CDFs/Rv_CDF_'+str(R_available[idcdf])+'_HamDM_SutDM_SutH_RicDM_RicGal.txt',usecols=(0,1,3,5)).T #Rfrac, DM, gal_hod, gal_obs
	#find fraction at which Rv*f = Rmin
	idx   = find_nearest(cdfs[0]*Rv[i],Rmin)
	CDF_DM.append(cdfs[1,idx])
	CDF_gal_HOD.append(cdfs[2,idx])
	CDF_gal_obs.append(cdfs[3,idx])
CDF_DM      = np.array(CDF_DM)
CDF_gal_HOD = np.array(CDF_gal_HOD)
CDF_gal_obs = np.array(CDF_gal_obs)

#prepare file to write
fout = "sigma_Nabh_Nbns_Npbh_Mmin"+str(mmin)+"_Rmin"+str(Rmin)+'_'+config+".txt"
f    = open(fout,'w')
f.write("#sigma\tNabh full\tNbns full\tNpbh full\terr(Npbhfull)\tNabh hod frac\tNabh obs frac\tNbns hod frac\tNbns obs frac\tNpbhfrac\terr(Npbhfrac)\n")
# loop over sigma
for sigma in sigmas:
	### ---------- AVERAGE RxPfav OVER Z SLICE ---------- ###
	#ABH - we only consider upper limits as it is what matter for discriminating with PBHS
	merger_rate_abh_all = rpf.rate_pfav_slice(sigma=sigma, m_min=mmin, network=config, channel='BBH',rate_case='hi',m_c=1.0, z_min=0.15, z_max=0.43, nz=57)
	#BNS -  same as ABH
	merger_rate_bns_all = rpf.rate_pfav_slice(sigma=sigma, m_min=mmin, network=config, channel='BNS',rate_case='hi',m_c=1.0, z_min=0.15, z_max=0.43, nz=57)
	#EARLY PBH - no error on rates from all sky
	merger_rate_pbh_all = rpf.rate_pfav_slice(sigma=sigma, m_min=mmin, network=config, channel='PBH',rate_case=None,m_c=1.0, z_min=0.15, z_max=0.43, nz=57)

	### ---------- INTEGRATION OVER Rv ---------- ###
	#i)  Take full voids, CDF=1
	#ABH
	N_abh_full, N_abh_full_m, N_abh_full_p = integrate_with_error(merger_rate_abh_all.R_grid,merger_rate_abh_all.values,Rv,gal_void_frac,err_gal_void_frac,1,vfill,err_vfill)	
	test = np.trapz(merger_rate_abh_all(Rv)*gal_void_frac*vfill, Rv)
	#upper limit, multiplied by volume factor
	N_abh_full = (N_abh_full+N_abh_full_p)*vol_survey/vol_tot 
	#BNS
	N_bns_full, N_bns_full_m, N_bns_full_p = integrate_with_error(merger_rate_bns_all.R_grid,merger_rate_bns_all.values,Rv,gal_void_frac,err_gal_void_frac,1,vfill,err_vfill)
	N_bns_full = (N_bns_full+N_bns_full_p)*vol_survey/vol_tot
	#PBH
	N_pbh_full, N_pbh_full_m, N_pbh_full_p = integrate_with_error(merger_rate_pbh_all.R_grid,merger_rate_pbh_all.values,Rv,cdm_void_frac,err_cdm_void_frac,1,vfill,err_vfill)
	N_pbh_full   *= vol_survey/vol_tot
	N_pbh_full_m *= vol_survey/vol_tot

	#ii) Take fraction of voids such that f*Rv = Rmin cMpc
	#ABH HOD and obs
	N_abh_frac_hod, N_abh_frac_hod_m, N_abh_frac_hod_p = integrate_with_error(merger_rate_abh_all.R_grid,merger_rate_abh_all.values,Rv,gal_void_frac,err_gal_void_frac,CDF_gal_HOD,vfill,err_vfill)
	N_abh_frac_hod = (N_abh_frac_hod+N_abh_frac_hod_p)*vol_survey/vol_tot
	N_abh_frac_obs, N_abh_frac_obs_m, N_abh_frac_obs_p = integrate_with_error(merger_rate_abh_all.R_grid,merger_rate_abh_all.values,Rv,gal_void_frac,err_gal_void_frac,CDF_gal_obs,vfill,err_vfill)
	N_abh_frac_obs = (N_abh_frac_obs+N_abh_frac_obs_p)*vol_survey/vol_tot
	#BNS HOD and obs
	N_bns_frac_hod, N_bns_frac_hod_m, N_bns_frac_hod_p = integrate_with_error(merger_rate_bns_all.R_grid,merger_rate_bns_all.values,Rv,gal_void_frac,err_gal_void_frac,CDF_gal_HOD,vfill,err_vfill)
	N_bns_frac_hod = (N_bns_frac_hod+N_bns_frac_hod_p)*vol_survey/vol_tot
	N_bns_frac_obs, N_bns_frac_obs_m, N_bns_frac_obs_p = integrate_with_error(merger_rate_bns_all.R_grid,merger_rate_bns_all.values,Rv,gal_void_frac,err_gal_void_frac,CDF_gal_obs,vfill,err_vfill)
	N_bns_frac_obs = (N_bns_frac_obs+N_bns_frac_obs_p)*vol_survey/vol_tot
	#PBH 
	N_pbh_frac, N_pbh_frac_m, N_pbh_frac_p = integrate_with_error(merger_rate_pbh_all.R_grid,merger_rate_pbh_all.values,Rv,cdm_void_frac,err_cdm_void_frac,CDF_DM,vfill,err_vfill)
	N_pbh_frac   *= vol_survey/vol_tot
	N_pbh_frac_m *= vol_survey/vol_tot #symmetric

	### ---------- WRITE  ---------- ###
	f.write(str(sigma)+"\t"+str(N_abh_full)+"\t"+str(N_bns_full)+"\t"+str(N_pbh_full)+"\t"+str(N_pbh_full_m)+"\t"+str(N_abh_frac_hod)+"\t"+str(N_abh_frac_obs)+"\t"+str(N_bns_frac_hod)+"\t"+str(N_bns_frac_obs)+"\t"+str()+"\t"+str(N_pbh_frac)+"\t"+str(N_pbh_frac_m)+"\n")

f.close()
print("Wrote results to "+fout)
