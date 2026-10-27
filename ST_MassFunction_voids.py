from numpy import*
from astropy import cosmology
import matplotlib.pyplot as plt
import argparse
import matplotlib
from scipy.special import gamma
from scipy.special import hyp2f1
from scipy.integrate import quad
from scipy.optimize import curve_fit

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

matplotlib.rc('xtick', labelsize=18)
matplotlib.rc('ytick', labelsize=18)


def extrapolate(Mlog,a,b):
        return a*Mlog+ b

parser = argparse.ArgumentParser()
parser.add_argument('-red', nargs='?', default=0.3, type=float, help='redshift')
parser.add_argument('-dmv',nargs='?', default=-0.24, type=float, help='delta voids')
args   = parser.parse_args()

#### NOTES ####
#Follow Verza, Carbone, Renzi 2022

#Cosmology
#cosmo = cosmology.Planck18.clone(H0=100.*h, Om0=0.29)
#cosmo.meta['sigma8']=0.79
z        = args.red #0.3 or 0.55
a        = 1./(1.+z)
#Omv      = 0.174 #dmv = -0.4
#Omv      = 0.219 #dmv = -0.24 (R=70, f=1)
#Omv      = 0.189 #dmv = -0.39 (R=70, f=1/3)
#Omv      = 0.194 #dmv = -0.33 (R=55, f=1)
#Omv      = 0.150 #dmv = -0.48 (R=55, f=1/3) 
Om0      = 0.29
Om       = Om0/a**3 / (Om0/a**3 + 1. - Om0)
rho_crit = 27.75e10 #h^2 Msun/Mpc^3
rho_m    = rho_crit*Om #cosmo.Om0 #h^2 Msun/Mpc^3 
h        = 0.6831
Df       = lambda x : x*hyp2f1(1./3.,1.,11./6.,x**3.*(Om0-1)/Om0) #Silveira & waga 1994 corrected
dc       = 1.686 #standard one alsi in Massara and sheth 
dcv      = 3. * (12.*pi)**(2./3.)*(1. + 0.0123*log10(Om))/20.
dsc      = dcv * Df(1)/Df(a)
dmv      = args.dmv #-0.48 -0.33 -0.39 -0.24#Omv/Om-1. #nonlinear
dmvL     = dc*(1.-(1+dmv)**(-1./dc))/(Df(a)/Df(1)) #rescaled to 0, relation between linear and nonlinear density contrast without taking volume effect from Massara and sheth 2018 1811.03132  
#if dmv == -0.39 or dmv == -0.48:
	#inner parts of voids
#	qU = 1.14
#	pU = 0.02
#else:
	#outer parts
qU = 0.783 #From figure 3 of Verza+ 2022 
pU = 0.215
Av = 1. / (1. + gamma(0.5-pU)/(2.**pU*sqrt(pi)))

#charger le PS et creer la masse
#load the PS and Create the mass
#PStxt = loadtxt("Gottlober_2003_voids_"+PSsuffix+"_matterpower_z0.dat").T
PStxt = loadtxt("SQUID_cosmo_Om0p29_OL0p71_Ob0p0293_s80p79_matterpower_z0.dat").T
k     = PStxt[0] #h/Mpc
PS    = PStxt[1] #(Mpc/h)^3

Mlog = linspace(2.0,16.0,14000)
M    = 10**Mlog #Msun/h    
R_M  = (3.0*M/(4.0*pi*rho_m))**(1.0/3.0) #Mpc/h rho_m or rho_m(1+dmv)??

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

#dsigma/dM 
#dlnSdlnM = gradient(log(s_m),Mlog)
dSdM     = gradient(s_m,M)

#mass function
dndlnM = zeros(SM, dtype = 'd')
hipass = Mlog >=7.9
nut    = (dsc-dmvL*Df(1)/Df(a)) / s_m[hipass]
fnut   = sqrt(2./pi) * Av * (1. + (qU*nut**2.)**(-pU)) *sqrt(qU) * exp(-qU*nut**2./2.)
dnutdM = - (dc-dmvL*Df(1)/Df(a)) * (1./s_m_sq[hipass]) * dSdM[hipass]
dndlnM[hipass] = (1. + dmv) * rho_m * fnut * dnutdM

#extrapolate
midpass = (Mlog >= 7.9) & (Mlog<=8.1)
popt, pcov = curve_fit(extrapolate, Mlog[midpass],log10(dndlnM[midpass]),p0=[-1.,1e-3])
lopass = Mlog<7.9
dndlnM[lopass] = 10.**(popt[1] + popt[0]*Mlog[lopass])

'''
#from Gottlober et al. 2003, equation 8
A   = 0.322
p   = 0.3
a   = 0.707

dndM = - sqrt(2./pi) * A * (1. + (a*dc**2./s_m_sq)**(-p)) * sqrt(a) * rho_b/M * dc/s_m_sq * dSdM * exp(-a*dc**2./(2.*s_m_sq))


ncumul = empty(SM)
for i in range(SM):
	dndMinv = dndM[::-1][0:i]
	ncumul[i] = sum(dndMinv)
	del dndMinv
ncumul = ncumul[::-1]
'''

#Mass function taken from other sources and nonvoids
#Gott = loadtxt("MassFunctionVoidsGottlober2003.txt").T
Tink = loadtxt("/Users/titouanlazeyras/Desktop/BICOCCA/SQUID/Tinker_M_dndlnM_Om0p29_h0p6831_Mpc-3").T 
'''
lT   = len(Tink[0])
Tcumul=empty(lT)
for i in range(lT):
	Tinv = Tink[2][::-1][0:i]
	Tcumul[i] = sum(Tinv)
	del Tinv
Tcumul = Tcumul[::-1]
'''
'''
plt.loglog()
plt.plot(Tink[0],Tink[1],'k:')
plt.plot(M,dndlnM,'b:')
#plt.plot(M,M*ncumul,'b--')
#plt.plot(Gott[0],Gott[1],'k--')
plt.xlim(1e13,7e15)
#plt.ylim(1e-7,1e2)
plt.xlabel(r'M [$M_\odot$/h]',fontsize=16)
plt.ylabel(r'$dn/dlnM$',fontsize=16)
#plt.ylabel(r'$n(>M)$ [Mpc$^{-3}$ $h^{-3}$]',fontsize=16)
#plt.savefig('M_dndM_ST_voids_'+PSsuffix+'z0.pdf',bbox_inches='tight')
plt.show()
'''

#M dndlnM
fout = open('M_dndlnM_ST_voids_Mpc-3_h0p6831_dmv'+str(dmv)+'_z'+str(z)+'.txt','w')
for i in range(SM):
	fout.write(str(M[i])+"\t"+str(dndlnM[i])+"\n") #*1e9*h**3 to get Gpc^-3
fout.close()
print("Wrote mass function to M_dndlnM_ST_voids_Mpc-3_h0p6831_dmv"+str(dmv)+'_z'+str(z)+'.txt')
