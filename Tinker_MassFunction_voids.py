from numpy import*
from astropy import cosmology
import matplotlib.pyplot as plt
import argparse
import matplotlib

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

matplotlib.rc('xtick', labelsize=18)
matplotlib.rc('ytick', labelsize=18)

#### NOTES ####
# To obtain the linear power spectrum to compute sigma, we assume a modified planck 2018 cosmology with OmegaM = 0.04, an enhanced baryon fraction of 0.18 and an open universe
# With h=0.6831 this gives Och^2 = 0.04*0.82*0.6831**2 = 0.015305320008 and Obh^2 = 0.04*0.18*0.6831**2 = 0.0033597043920000004
# We also have OL=0.71 in Planck, so Ok = 1.-0.71-0.04 = 0.25

#parser   = argparse.ArgumentParser(description='python Tinker_MassFunction_voids.py < z (float) >')
#parser.add_argument('red',nargs=1,type=str, help='Redshift (float)')
#args     = parser.parse_args()
#z        = args.red[0]
z = 0.0

#Cosmology, same as VSF (~Planck 2018) modified for being in void
Omv = 0.04 #0.29 0.04 #0.307
Om  = 0.29 #0.307
rho_crit = 27.75e10 #h^2 Msun/Mpc^3
rho_m    = rho_crit*Om #cosmo.Om0 #h^2 Msun/Mpc^3 
rho_b    = rho_crit*Omv 

#charger le PS et creer la masse
#load the PS and Create the mass
PStxt = loadtxt("Planck_2018_voids_Om0p04_OL0p71_Ok0p25_s80p2233_matterpower_z0.dat").T
#PStxt  = loadtxt("SQUID_cosmo_Om0p29_OL0p71_Ob0p0293_s80p79_matterpower_z0.dat").T
k = PStxt[0] #h/Mpc
PS = PStxt[1] #(Mpc/h)^3

Mlog = linspace(2.0,16.0,14000)
M = 10**Mlog #Msun/h    
R_M = (3.0*M/(4.0*pi*rho_m))**(1.0/3.0) #Mpc/h

Sk = shape(k)[0]
SM = shape(M)[0]

#mass function calculation for all R
#tophat filter for all R
W_M = zeros((SM, Sk),dtype = 'd')
for i in range(SM) :
        for l in range(Sk) :
                W_M[i,l] = 3.0/(k[l]*R_M[i])**3.0*(sin(k[l]*R_M[i])-k[l]*R_M[i]*cos(k[l]*R_M[i]))

#sigma calculation for All R_M
s_m_sq = trapz(1.0/(2.0*(pi)**2.0)*W_M**2.*PS*k**2.,k)
s_m    = sqrt(s_m_sq)

#dsigma/dM (derivative computed by hand, here just rcompute integral over k)Cc = 4.0*pi*rho_m/3.0
dlnSdlnM = gradient(log(s_m),Mlog)
dSdM     = gradient(s_m,M)

#mass function
#from Tinker et al. 2008, equations 2-8
A0 = 0.186
a0 = 1.47
b0 = 2.57
c  = 1.19
A  = A0*(1.+z)**(-0.14)
a  = a0*(1.+z)**(-0.06)
af = 0.00771203275022638 #for Delta = 200 af = 10.**(-1.2*(3./(4.*log10(200./75.))))
b  = b0*(1.+z)**af
f_s_m = A*((s_m/b)**(-a)+1.0)*exp(-c/s_m_sq)
dndM  = -f_s_m/s_m*rho_b/M*dSdM

ncumul = empty(SM)
for i in range(SM):
	dndMinv = dndM[::-1][0:i]
	ncumul[i] = sum(dndMinv)
	del dndMinv
ncumul = ncumul[::-1]

#Mass function taken from paper
Gott = loadtxt("MassFunctionVoidsGottlober2003.txt").T
Tink = loadtxt("/Users/titouanlazeyras/Desktop/BICOCCA/codes/Python_codes/MassFunction/Sigma_TinkerMassFunction_OfM").T 
lT   = len(Tink[0])
Tcumul=empty(lT)
for i in range(lT):
	Tinv = Tink[2][::-1][0:i]
	Tcumul[i] = sum(Tinv)
	del Tinv
Tcumul = Tcumul[::-1]

plt.loglog()
#plt.plot(Tink[0],Tink[2],'k:')
#plt.plot(M,M**2/rho_m*dndM,'b:')
#plt.plot(M,dndM,'b-')
plt.plot(M,M**2./rho_m*dndM,'b--')
plt.plot(Gott[0],Gott[1],'k--')
#plt.xlim(1e6,4e15)
#plt.ylim(2e-7,20.)
plt.xlabel(r'M [$M_\odot$/h]',fontsize=16)
plt.ylabel(r'$n(>M)$ [Mpc$^{-3}$ $h^{-3}$]',fontsize=16)
#plt.savefig('M_dndM_Tinker_voids_Om0p04_OL0p71_Ok0p25_s80p2233_z0p3.pdf',bbox_inches='tight')
plt.show()

#WRITE M dndM
#fout = open('M_sigma_dndM_Tinker_Om0p29_OL0p71_s80p79_z0.txt','w')
#for i in range(SM):
#	fout.write(str(M[i])+"\t"+str(s_m[i])+"\t"+str(dndM[i])+"\n")
#fout.close()
