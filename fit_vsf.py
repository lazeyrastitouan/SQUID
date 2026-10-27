import matplotlib.pyplot as plt
from numpy import*
from astropy import cosmology
import astropy.units as u
import matplotlib
import argparse
from scipy.optimize import curve_fit

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

matplotlib.rc('xtick', labelsize=16)
matplotlib.rc('ytick', labelsize=16)

def quadfit(R,a,b,c):
	return a*R*R + b*R + c

def quadratic_fit_uncertainty(R, pcov):
    J = vstack([R**2, R, ones_like(R)]).T #Jacobian
    return sqrt(einsum("ij,jk,ik->i", J, pcov, J))

'''
dr11      = loadtxt("R_Nvoids_sNvoids_BOSSDR11_LOWZ_cumul_raw.dat").T
R11       = dr11[0] #cmpc/h
VSFc      = dr11[1][::-1]
err11     = dr11[2][::-1]-VSFc
nbin      = len(R11)
VSF       = zeros(nbin)
errVSF    = zeros(nbin)
VSF[0]    = VSFc[0]
errVSF[0] = err11[0]

for i in range(nbin-1):
	VSF[i+1]    = VSFc[i+1]-VSFc[i] 
	errVSF[i+1] = sqrt(err11[i+1]**2.+err11[i]**2.)

VSF    = VSF[::-1]
errVSF = errVSF[::-1]
'''

#from dn/dR
dndRraw   = loadtxt("dn_dR_BOSS_DR11_LOWZ_Nadathur_raw.dat").T
R         = dndRraw[0]
deltaR    = zeros(len(R))
deltaR[0] = R[1]-R[0] #equal bins in log10(R)
for i in range(1,len(R)):
        deltaR[i] = R[i]-R[i-1]
dndR    = dndRraw[1]*deltaR
errsup  = dndRraw[2]*deltaR-dndR
errinf  = dndR-dndRraw[3]*deltaR
#symetrize erro as curve_fit doesn take asymetric onesr
errdndR = (errsup+errinf)/2.
logdndR = log10(dndR)
logerr  = errdndR/(dndR * log(10)) #approximate error in log10

#from cumulative VSF
dr11   = loadtxt("R_nvoids_snvoids_BOSSDR11_LOWZ.dat").T
filter = dr11[1] > 0
R11    = dr11[0][filter]
VSF    = log10(dr11[1][filter])
errVSF = dr11[2][filter]/(dr11[1][filter] * log(10)) #approximate error in log10

#combine both to maximize R range
maxval        = max(R)
filtercombine = R11 > maxval
Rfinal        = concatenate((R,R11[filtercombine]))
VSFfinal      = concatenate((logdndR,VSF[filtercombine]))
errorfinal    = concatenate((logerr,errVSF[filtercombine])) 
errorfinal    = where(errorfinal>0,errorfinal,0.11/log(10)) 

#fit + error on fit
popt, pcov  = curve_fit(quadfit,Rfinal,VSFfinal,sigma=errorfinal,absolute_sigma=True)
rmock       = linspace(Rfinal.min(), Rfinal.max(), 10000)
vsf_fit     = quadfit(rmock, popt[0], popt[1], popt[2])
vsf_fit_err = quadratic_fit_uncertainty(rmock, pcov)

#save
fout = open("vsf_Nadathur_BOSSDR11_LOWZ_quad_fit.txt",'w')
for i in range(len(rmock)):
	fout.write(str(rmock[i])+"\t"+str(10**vsf_fit[i])+"\t"+str(10**(vsf_fit[i]+vsf_fit_err[i])-10**vsf_fit[i])+"\n")
fout.close()

# plot
plt.figure()
plt.yscale('log')
plt.xlabel(r'$R_v$ [Mpc/$h$]', fontsize=17)
plt.ylabel(r'$N_{\textrm{voids}} \,(R_v)$',fontsize = 17)
plt.plot(R11,10**VSF,'b.',label=r'DR11 LOWZ')
plt.errorbar(R11,10**VSF,xerr=0.0,yerr=dr11[2][filter],fmt='b.')
plt.errorbar(R,dndR,xerr=0.0,yerr=errdndR,fmt='gx')
plt.errorbar(Rfinal,10**VSFfinal,xerr=0.0,yerr=errorfinal*10**VSFfinal*log(10),fmt='r.')
plt.plot(rmock, 10**vsf_fit, 'k--',label="quadratic fit")
plt.fill_between(rmock,10**(vsf_fit-vsf_fit_err),10**(vsf_fit+vsf_fit_err),color='black',alpha=0.3)
plt.legend(loc='upper right',ncol=1, prop={'size':11}, frameon=False)
plt.savefig("Nvoids_of_R_With_fit_mpc_z0p3.pdf",bbox_inches='tight')
plt.show()
