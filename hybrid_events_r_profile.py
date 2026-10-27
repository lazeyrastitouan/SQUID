#########################################################################################################################################################
#																			#
#	COMPUTE NUMBER OF EVENTS AS A FUNCTION OF SCALE IN 2 WAYS 											#
#	i)  INSIDE VOIDS OF GIVEN RADII 														#
#	ii) ABOVE A GIVEN SCALE USING ALL VOIDS OF SUFFICIENT SIZE) 											#
#	BY FIXING THE NUMBER OF EVENTS SEEN FROM VOIDS AND ADAPTING THE MERGER RATES TO THE VALUES NEEDED TO MATCH THESE HYPOTETICAL OBSERVATIONS 	#
#	OBJECTIVE: SEE IF SCALE DEPENDENCE ALONE IS SUFFICIENT TO DISCRIMINTE PBH FROM ASTRO BH								#
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
vol_survey = 8.0 #Gpc^3 DESI BGS
Ntarget    = 10 #events per 10 years
Rtarget    = [[50,60,70,80,90],[30,40,50,60,70,80,90]]
fill              = load('volume_filling_observed_z03.npz')
vfill             = fill['fv_of_R']
err_vfill         = fill['err_fv_of_R'] 
Rv                = fill['R']
rs                = len(Rv)
err_merger_rate   = array([0,0])
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
B = True  #For a given void radius, radial profile

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
	if A:
		f1  = open("hybrid_number_events_of_R_"+str(Mmins[im])+"Msol.txt",'w')
	
	#i)  Find the merger rate for astro and early pbh such that N=10/10yr at CDF=1
	filter = Rv>=Rmin
	int_full = trapz(vfill[filter],Rv[filter])
	merger_rate_gal_full = Ntarget/(int_full*vol_survey*gal_void_frac*1.0*Pfav_gal)
	merger_rate_pbh_full = Ntarget/(int_full*vol_survey*cdm_void_frac*1.0*Pfav_pbh) 

	#ii) Compute merger rate for various fraction of voids corresponding to different scales to have the r profile
	for rt in Rtarget[im]:
		if A:
			merger_rate_gal_frac, err_merger_rate_gal_frac = [],[] 
			merger_rate_hod_frac, err_merger_rate_hod_frac = [],[]
			merger_rate_pbh_frac, err_merger_rate_pbh_frac = [],[]
			for i in range(rs):
				if Rv[i]< rt:
					continue
				#CDFs computed more sparsely than Vfill, lets find the CDF value that matches Rv[i] the closest
				idcdf = find_nearest(R_available,Rv[i])
				cdfs  = loadtxt('CDFs/Rv_CDF_'+str(R_available[idcdf])+'_HamDM_SutDM_SutH_RicDM_RicGal.txt',usecols=(0,1,3,5)).T #Rfrac, DM, gal_hod, gal_obs
				#find fraction at which Rv*f = rt
				idx   = find_nearest(cdfs[0]*Rv[i],rt)
		
				merger_rate_hod_frac.append(RinV(merger_rate_gal_full,gal_void_frac,cdfs[2,idx],vfill[i])*Pfav_gal)
				err_merger_rate_hod_frac.append(err_RinV(merger_rate_gal_full,gal_void_frac,cdfs[2,idx],vfill[i],err_merger_rate,err_gal_void_frac,err_vfill[i])*Pfav_gal)
				merger_rate_gal_frac.append(RinV(merger_rate_gal_full, gal_void_frac,cdfs[3,idx],vfill[i])*Pfav_gal)
				err_merger_rate_gal_frac.append(err_RinV(merger_rate_gal_full, gal_void_frac,cdfs[3,idx],vfill[i], err_merger_rate, err_gal_void_frac,err_vfill[i])*Pfav_gal)
				merger_rate_pbh_frac.append(RinV(merger_rate_pbh_full,cdm_void_frac,cdfs[1,idx],vfill[i])*Pfav_pbh)
				err_merger_rate_pbh_frac.append(err_RinV(merger_rate_pbh_full,cdm_void_frac,cdfs[1,idx],vfill[i],err_merger_rate,err_cdm_void_frac,err_vfill[i])*Pfav_pbh)

			err_merger_rate_hod_frac = array(err_merger_rate_hod_frac).reshape(shape(err_merger_rate_hod_frac)[0],shape(err_merger_rate_hod_frac)[1])
			err_merger_rate_gal_frac = array(err_merger_rate_gal_frac).reshape(shape(err_merger_rate_gal_frac)[0],shape(err_merger_rate_gal_frac)[1])
			err_merger_rate_pbh_frac = array(err_merger_rate_pbh_frac).reshape(shape(err_merger_rate_pbh_frac)[0],shape(err_merger_rate_pbh_frac)[1])
		
			filter   = Rv>=rt
			Nhod        = trapz(merger_rate_hod_frac,Rv[filter])*vol_survey
			err_Nhod_up = trapz(err_merger_rate_hod_frac[:,1],Rv[filter])*vol_survey #errors are asymetric
			err_Nhod_dn = trapz(err_merger_rate_hod_frac[:,0],Rv[filter])*vol_survey
			Ngal        = trapz(merger_rate_gal_frac,Rv[filter])*vol_survey
			err_Ngal_up = trapz(err_merger_rate_gal_frac[:,1],Rv[filter])*vol_survey
			err_Ngal_dn = trapz(err_merger_rate_gal_frac[:,0],Rv[filter])*vol_survey
			Npbh        = trapz(merger_rate_pbh_frac,Rv[filter])*vol_survey		
			err_Npbh    = trapz(err_merger_rate_pbh_frac[:,1],Rv[filter])*vol_survey #errors are symetric for pbh 

			f1.write(str(rt)+"\t"+str(Nhod)+"\t"+str(err_Nhod_up)+"\t"+str(err_Nhod_dn)+"\t"+str(Ngal)+"\t"+str(err_Ngal_up)+"\t"+str(err_Ngal_dn)+"\t"+str(Npbh)+"\t"+str(err_Npbh)+"\n")
		if B:
			idcdf = find_nearest(R_available,rt)
			cdfs  = loadtxt('CDFs/Rv_CDF_'+str(R_available[idcdf])+'_HamDM_SutDM_SutH_RicDM_RicGal.txt',usecols=(0,1,3,5)).T #Rfrac, DM, gal_hod, gal_obs
			idvf  = find_nearest(Rv,rt)
			merger_rate_hod     = merger_rate_gal_full*gal_void_frac*cdfs[2]*vfill[idvf]*Pfav_gal*vol_survey #since multiplied by volume this is actually a number of event per 10 year
			err_merger_rate_hod = err_RinV(merger_rate_gal_full,gal_void_frac,cdfs[2],vfill[idvf],err_merger_rate,err_gal_void_frac,err_vfill[idvf])*Pfav_gal*vol_survey
			merger_rate_gal     = merger_rate_gal_full*gal_void_frac*cdfs[3]*vfill[idvf]*Pfav_gal*vol_survey
			err_merger_rate_gal = err_RinV(merger_rate_gal_full,gal_void_frac,cdfs[3],vfill[idvf],err_merger_rate,err_gal_void_frac,err_vfill[idvf])*Pfav_gal*vol_survey
			merger_rate_pbh     = merger_rate_pbh_full*cdm_void_frac*cdfs[1]*vfill[idvf]*Pfav_pbh*vol_survey
			err_merger_rate_pbh = err_RinV(merger_rate_pbh_full,cdm_void_frac,cdfs[1],vfill[idvf],err_merger_rate,err_cdm_void_frac,err_vfill[idvf])*Pfav_pbh*vol_survey
			
			f2 = "hybrid_number_events_radial_profile_Rv"+str(rt)+"_"+str(Mmins[im])+"Msol"
			savez(f2,r=cdfs[0],Nhod=merger_rate_hod,err_Nhod=err_merger_rate_hod,Ngal=merger_rate_gal,err_Ngal=err_merger_rate_gal,Npbh=merger_rate_pbh,err_Npbh=err_merger_rate_pbh)
	if A:
		f1.close()
