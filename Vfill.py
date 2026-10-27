# ==========================================================================
# Theoretical size function of cosmic voids
# ==========================================================================
# to ensure compatibility in Python versions 2.x and 3.x
#from __future__ import print_function
 
import numpy as np
import os
import matplotlib.pyplot as plt
#import CosmoBolognaLib as cbl
import matplotlib

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

matplotlib.rc('xtick', labelsize=18)
matplotlib.rc('ytick', labelsize=18)

def vol(R):
    return 4./3.*np.pi*R**3.

'''
# define a cosmological model, using Planck 2018 parameters #
cosm = cbl.Cosmology(cbl.CosmologicalModel__Planck18_)
cosm.set_parameter(cbl.CosmologicalParameter__sigma8_,0.79)
cosm.set_parameter(cbl.CosmologicalParameter__Omega_matter_,0.29)
cosm.set_parameter(cbl.CosmologicalParameter__Omega_DE_,0.71)
cosm.set_parameter(cbl.CosmologicalParameter__Omega_neutrinos_,0.0)
s8 = cosm.value(cbl.CosmologicalParameter__sigma8_)

print(cosm.print_parameters())

# Minimum and maximum of effective void radii
R_min = 1e-3
R_max = 100.

# number of radii at which the size function is computed
n_val = 10000
 
# list of effective void radii with logarithmic binning
RR = np.logspace(np.log10(R_min), np.log10(R_max), n_val , endpoint=True)
# redshift of the samples: 0.3 for LOWZ, 0.55 for CMASS
zz = 0.3
h  = cosm.value(cbl.CosmologicalParameter__hh_)

#Nonlinear void density threshold and effective bias 
#Matching Contarini BOSSDR12 clean, Nadathur BOSSDR11 
vdn_params = [{'b_eff':1.36/s8, 'slope':1., 'offset':1.5, 'deltav_NL':-0.795, 'delta_c':1.686},{'b_eff':2./s8, 'slope':2.5, 'offset':0., 'deltav_NL':-0.795, 'delta_c':1.686},{'b_eff':1.,'slope':1.,'offset':0., 'deltav_NL':-0.795, 'delta_c':1.686}]

#BOSS DR12 volume Veff = 5.1 Gpc3 for CMASS and 2.3 Gpc3 for LOWZ 1.15 Gpc3 for DR11
volume   = 2.3e9*h**3 #(Mpc/h)^3
ds       = ['DR12','DR11','CDM']
BOSSDR12 = np.array(([38,42,47],[15,11,4]))
BOSSDR11 = np.array(([8.90625,10.625,13.75,17.65625,20.9375,24.84375,32.03125,37.96875,44.84375,53.28125,57.96875],[3.471958582498651e-9,1.085200856825879e-8,2.1881942729662105e-8,3.702692861806392e-8,4.041938092683189e-8,4.816522514798774e-8,6.084975646464324e-8,5.7395458821001434e-8,4.543099853147129e-8,3.199369555101389e-8,1.9468085826831354e-8]))

plt.figure(1)
plt.yscale('log')
#plt.plot(BOSSDR12[0],BOSSDR12[1])
plt.plot(BOSSDR11[0],BOSSDR11[1],'k.')
plt.xlabel(r'$R_\textrm{eff} \, [\textrm{Mpc}/h]$',fontsize=24, labelpad=14)
plt.ylabel(r'$\frac{\textrm{d}n}{\textrm{d}R} \; [(h/\textrm{Mpc})^4]$',fontsize=24, labelpad=14)
plt.text(20,1e-11,r'$z=0.3$',fontsize=16)
plt.ylim(5e-13,5e-3)
plt.xlim(0,100)
filter = (20.0 < RR) & (90 > RR)
'''
VSFfile = np.loadtxt("vsf_Nadathur_BOSSDR11_LOWZ_quad_fit.txt").T
RR      = VSFfile[0]
VSF     = VSFfile[1]
err     = VSFfile[2]
fout    = 'volume_filling_observed_z03.npz'
fvup    = []
fvuperr = []
#denominator
Ntot    = np.trapz(VSF*vol(RR),RR)
print(Ntot)
Ntoterr = np.trapz((VSF+err)*vol(RR),RR)-Ntot #error on Ntot is symmetrical
fv      = VSF*vol(RR)/Ntot
fverr   = fv*np.sqrt((err/VSF)**2+(Ntoterr/Ntot)**2.)
for j in range(len(RR)):
    '''VSF = np.empty(n_val)
    for i in range(n_val):
        #Sheth & van de Weygaert size function
        #print('%.6e' %cosm.size_function(RR[i], zz, "SvdW", b_eff), '(h/Mpc)^3 at R = %.3g' %RR[i], ' Mpc/h ')
        # Jennings +13 VSF
        #arguments of the function: const double RV, const double redshift, const std::string model, const double b_eff, double slope, double offset, const double deltav_NL, const double del_c, const std::string method_Pk, const bool store_output, const std::string output_root, const std::string interpType, const double k_max, const std::string input_file, const bool is_parameter_file
        VSF[i] = cosm.size_function(RR[i], zz, "Vdn", vdn_params[j]['b_eff'],vdn_params[j]['slope'],vdn_params[j]['offset'],vdn_params[j]['deltav_NL'],vdn_params[j]['delta_c']) 
    
    #Nvoids (const double min_r, const double max_r, const int num_bins, const double mean_z, const double Volume, const std::string model, const double b_eff, double slope, double offset, const double deltav_NL, const double del_c, const std::string method_Pk, const bool store_output, const std::string output_root, const std::string interpType, const double k_max, const std::string input_file, const bool is_parameter_file)
    #NV  = cosm.Nvoids(R_min,R_max,n_val,zz,volume,"Vdn",b_eff,slope,offset,deltav_NL,delta_c)

    #VSF=VSF*volume
    '''
   
    #now lets integrate from various minimum R to get the volume filling factor above a given radius
    integral = np.trapz(VSF[j:]*vol(RR[j:]),RR[j:])
    fvup.append(integral/Ntot)
    interr = np.trapz((VSF[j:]+err[j:])*vol(RR[j:]),RR[j:])-integral
    fvuperr.append(integral/Ntot*np.sqrt((interr/integral)**2.+(Ntoterr/Ntot)**2))

fvup    = np.array(fvup)
fvuperr = np.array(fvuperr)

plt.yscale('log')
plt.plot(RR,fv,'k--',label=r'filling factor at Rv')
plt.fill_between(RR,fv-fverr,fv+fverr,color='k',alpha=0.2)
plt.plot(RR,fvup,'b--',label=r'filling factor above Rv')
plt.fill_between(RR,fvup-fvuperr,fvup+fvuperr,color='b',alpha=0.2)
plt.xlabel(r'Rv [mpc/h]',fontsize=14)
plt.xlim(9.5,92)
plt.ylim(1e-5,2)

#np.savez(fout,R=RR,fv_of_R=fv,err_fv_of_R=fverr,fv_above_R=fvup,err_fv_above_R=fvuperr)

plt.legend(loc=0,prop={'size':14},frameon=False)
#plt.savefig('Vfill_observed_VSF_DR11_lowz.pdf',bbox_inches='tight')
plt.show()
