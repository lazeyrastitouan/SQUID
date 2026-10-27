#########################################################################################################################################################
#																			#
#	COMPUTE NUMBER OF EVENTS FROM 															#
#	i)   ALL VOIDS 																	#
#	AND AS A FUNCTION OF SCALE IN 2 WAYS 														#
#	ii)  INSIDE VOIDS OF GIVEN RADII 														#
#	iii) ABOVE A GIVEN SCALE USING ALL VOIDS OF SUFFICIENT SIZE) 											#
#	WHEN FIXING THE FULL SKY NUMBER OF EVENTS AND ADAPTING THE MERGER RATES IN VOIDS TO THE VALUES NEEDED TO MATCH THis HYPOTETICAL OBSERVATION 	#
#	OBJECTIVE: SEE IF NUMBER FROM VOIDS DIFFERS SUFFICIENTLY TO DISCRIMINTE PBH FROM ASTRO BH AGNOSTICALLY, AND IF SCALE DEPENDENCE HELPS		#
#																			#
#########################################################################################################################################################


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
vol_survey = 3.6 #8.0 #Gpc^3 DESI BGS
vol_tot    = 19.53 #23.22 #Gpc^3 z=0.05 - 0.45
Ntarget    = 10000 #events per 10 year
Rtarget    = [[50,60,70,80,90],[30,40,50,60,70,80,90]]
fill              = load('volume_filling_observed_z03.npz')
vfill             = fill['fv_of_R']
err_vfill         = fill['err_fv_of_R'] 
Rv                = fill['R']
rs                = len(Rv)
Rmins             = [49.8,25.7]
Mmins             = [0.5,1.0]
gal_void_frac     = 0.09
err_gal_void_frac = array([0.02,0.035]) 
Pfav_gal          = 0.003
cdm_void_frac     = 0.15
err_cdm_void_frac = array([0.05,0.05])
Pfav_pbh          = 0.01

### --------- WHAT TO DO --------- ###
A = False #Fixed scale, all events from all voids larger than this scale, use CDF at scale
B = False  #For a given void radius, radial profile

#Since the CDF might not have been computed for exactly the same Rs as the VSF, need to find the closest R value of CDF (assum CDF don't vary much with R)
cdfpath     = pathlib.Path("/Users/titouanlazeyras/Desktop/BICOCCA/SQUID/CDFs/")
R_available = []
pattern     = re.compile("Rv_CDF_([0-9.]+)_HamDM")
for file in cdfpath.glob('*.txt'):
	match = pattern.search(file.name)
	if match:
		R_available.append(float(match.group(1)))
R_available = array(R_available)

print("Your target number of event per 10 yr is "+str(Ntarget))
for Rmin in Rmins:
	im = Rmins.index(Rmin)
	N_gal_in_void, err_N_gal_in_void, N_pbh_in_void, err_N_pbh_in_void = [],[],[],[]
	for i in range(rs):
		if Rv[i]< Rmin:
			continue
		#i) Predict number of events from all voids larger than localisation
		N_gal_in_void.append(RinV(Ntarget, gal_void_frac, 1, vfill[i])*Pfav_gal)
		err_N_gal_in_void.append(err_RinV(Ntarget, gal_void_frac, 1, vfill[i], array([0,0]), err_gal_void_frac, err_vfill[i])*Pfav_gal)
		N_pbh_in_void.append(RinV(Ntarget, cdm_void_frac, 1, vfill[i])*Pfav_pbh)
		err_N_pbh_in_void.append(err_RinV(Ntarget, cdm_void_frac, 1, vfill[i], array([0,0]), err_cdm_void_frac, err_vfill[i])*Pfav_pbh)

	err_N_gal_in_void = array(err_N_gal_in_void).reshape(shape(err_N_gal_in_void)[0],shape(err_N_gal_in_void)[1])
	err_N_pbh_in_void = array(err_N_pbh_in_void).reshape(shape(err_N_pbh_in_void)[0],shape(err_N_pbh_in_void)[1])

	#integral in full voids
	filter = Rv>=Rmin
	n_events_gal_full = trapz(N_gal_in_void,Rv[filter])*vol_survey/vol_tot
	err_n_gal_full_p  = trapz(err_N_gal_in_void[:,1],Rv[filter])*vol_survey/vol_tot #upper error
	err_n_gal_full_m  = trapz(err_N_gal_in_void[:,0],Rv[filter])*vol_survey/vol_tot #lower error
	n_events_pbh_full = trapz(N_pbh_in_void,Rv[filter])*vol_survey/vol_tot
	err_n_pbh_full    = trapz(err_N_pbh_in_void[:,1],Rv[filter])*vol_survey/vol_tot 

	print("Rmin (Mmin) = "+str(Rmin)+" ("+str(Mmins[im])+")")
	print("N events all voids\nPBH: "+str(n_events_pbh_full)+" +/- "+str(err_n_pbh_full))
	print("Astro : "+str(n_events_gal_full)+" + "+str(err_n_gal_full_p)+" - "+str(err_n_gal_full_m))
