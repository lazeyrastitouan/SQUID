### compute the mean density (and thus Omega_m) in voids of various radius following the profile of Hamaus et al PRL 2014

from numpy import*
from scipy.interpolate import CubicSpline
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import matplotlib

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

matplotlib.rc('xtick', labelsize=14)
matplotlib.rc('ytick', labelsize=14)

def CDFR(r,Re,alpha,beta):
	return (r/Re)**alpha*exp((r/Re)**beta -1.)

def dprofile_Hamaus(dc,rs,Rv):
	'''
	Compute the density profile following Hamaus+ 20214 fit
	Rv is the void radius
	dc, rs are the fit's free parameters at Rv
	'''
	a  = -2.*rs/Rv + 4.
	if rs/Rv <= 0.91:
		b = 17.5*rs/Rv - 6.5
	else:
		b = -9.8*rs/Rv + 18.4

	r = linspace(0,Rv,100)
	rho = rhom * (dc * ((1. - (r/rs)**a)/(1. + (r/Rv)**b)) + 1.) #density
	
	return r, rho


Om       = 0.29
rho_crit = 27.75e10
rhom     = Om*rho_crit
colors   = ['black','grey','yellow','orange','blue','green','magenta','red','darkred']
#Rv = array([8.0,12.0,24.0,36.0,48.0,55.0,60.0,70.0,80.0])
R_min    = 1e-3
R_max    = 100.
n_val    = 10000
RR       = logspace(log10(R_min), log10(R_max), n_val , endpoint=True)
Rv       = RR[(20.0 < RR) & (90 > RR)] #need to go through this mess to match volume filling Rvs 
lr       = len(Rv)

# Hamaus 2014 DM
rvmock   = array([8.7,11.860490045941805,16.322358346094944,21.43491577335375,26.547473200612558,34.26278713629402,46.53292496171517,69.4])
dcmock   = array([-0.9543269230769231,-0.9479768786127167,-0.9190751445086706,-0.8439306358381503,-0.7341040462427746,-0.6271676300578035,-0.5520231213872833,-0.49132947976878616])
sdcmock  = array([0.03365692,0.02730688,0.02003514,0.01220064,0.01054405,0.00937763,0.01837312,0.05622948])
rsmock   = array([0.7503607503607505,0.7956678700361012,0.8635379061371842,0.940072202166065,1.0166064981949459,1.132129963898917,1.3155234657039712,1.6577617328519856])*rvmock
srsmock  = array([0.0034632,0.00259428,0.00466345,0.00847047,0.00996869,0.01120644,0.03991019,0.23438511])*rvmock
dcinterp = CubicSpline(rvmock,dcmock)
rsinterp = CubicSpline(rvmock,rsmock)
dclinfit = poly1d(polyfit(rvmock,dcmock,1,w=1./sdcmock)) 
rslinfit = poly1d(polyfit(rvmock,rsmock,1,w=1./srsmock))
dcqadfit = poly1d(polyfit(rvmock,dcmock,2,w=1./sdcmock))
rsqadfit = poly1d(polyfit(rvmock,rsmock,2,w=1./srsmock))
print("Hamaus dc(0)="+str(dclinfit(0)))
print("Hamaus dc(0)="+str(dcqadfit(0)))

#Sutter 2014 
#DM full
rvmockS  = array([10.,15.,20.,25.,30.,35.])  
dcmockS  = array([-0.80782,-0.77907,-0.78161,-0.78414,-0.76469,-0.72326])
sdcmockS = array([0.02074,0.02805,0.02024,0.0156,0.01183,0.00796])
rsmockS  = array([0.76872,0.80900,0.86351,0.939697,1.02464,1.12773])*rvmockS
srsmockS = array([0.00257,0.00507,0.00854,0.004603,0.01376,0.02379])*rvmockS
dcinterpS = CubicSpline(rvmockS,dcmockS)
rsinterpS = CubicSpline(rvmockS,rsmockS)
dclinfitS = poly1d(polyfit(rvmockS,dcmockS,1,w=1./sdcmockS))
rslinfitS = poly1d(polyfit(rvmockS,rsmockS,1,w=1./srsmockS))
print("Sutter DM dc(0)="+str(dclinfitS(0)))

#HOD sparse
rvmockSH  = array([25.,30.,35.,40.,45.,50.,55.,60.,65.,70.,75.])
dcmockSH  = array([-0.82389,-0.87378,-0.82135,-0.78330,-0.75116,-0.77400,-0.69958,-0.67336,-0.69873,-0.65729,-0.57696])
sdcmockSH = array([0.04023, 0.03888, 0.02588, 0.03649, 0.0433 , 0.03523, 0.06005,0.05039, 0.05562, 0.04218, 0.02647])
rsmockSH  = array([0.75924,0.76517,0.77583,0.79834,0.80427,0.81848,0.83507,0.86469,0.83270,0.85284,0.83981])*rvmockSH
srsmockSH = array([0.00886,0.00623,0.00453,0.00709,0.01302,0.00916,0.021,0.0183,0.02053,0.01501,0.01053])*rvmockSH
dcinterpSH = CubicSpline(rvmockSH,dcmockSH)
rsinterpSH = CubicSpline(rvmockSH,rsmockSH)
dclinfitSH = poly1d(polyfit(rvmockSH,dcmockSH,1,w=1./sdcmockSH))
rslinfitSH = poly1d(polyfit(rvmockSH,rsmockSH,1,w=1./srsmockSH))
print("Sutter HOD sparse dc(0)="+str(dclinfitSH(0)))

#HOD dense
rvmockSHd  = array([10.,15.,20.,25.,30.,35.,40.,45.])
dcmockSHd  = array([-0.83034,-0.86412,-0.81029,-0.76358,-0.72797,-0.65092,-0.57810,-0.56966])
sdcmockSHd = array([0.01794,0.03483,0.04855,0.04194,0.03694,0.02955,0.02744,0.02955])
rsmockSHd  = array([0.69762,0.73563,0.75226,0.76651,0.79026,0.79976,0.81758,0.80570])*rvmockSHd
srsmockSHd = array([1e-3,1e-3,0.00763,0.01518,0.00756,0.00851,0.0117,0.01217])*rvmockSHd
dcinterpSHd = CubicSpline(rvmockSHd,dcmockSHd)
rsinterpSHd = CubicSpline(rvmockSHd,dcmockSHd)
dclinfitSHd = poly1d(polyfit(rvmockSHd,dcmockSHd,1,w=1./sdcmockSHd))
rslinfitSHd = poly1d(polyfit(rvmockSHd,rsmockSHd,1,w=1./srsmockSHd))
print("Sutter HOD dense dc(0)="+str(dclinfitSHd(0)))

'''
print('dc(0) Hamaus '+str(dclinfit(0)))
print('dc(0) Sutter '+str(dclinfitS(0)))
print('dc(0) Sutter HOD sparse '+str(dclinfitSH(0)))
print('dc(0) Sutter HOD dense '+str(dclinfitSHd(0)))
'''
fig1 = plt.figure(1)
plt.plot(rvmock,dcmock,'r.',label=r'Hamaus+2014')
plt.errorbar(rvmock,dcmock,xerr=0.0,yerr=sdcmock,fmt='r.')
plt.plot(Rv,dclinfit(Rv),'r:',label=r'DM full linear fit')
plt.plot(Rv,dcqadfit(Rv),'r-.',label=r'DM full quad fit')
#plt.plot(Rv,dcinterpS(Rv),'r--',label=r'DM full cubic spline')
plt.plot(rvmockSH,dcmockSH,'b.',label=r'Sutter+2014 sparse')
plt.errorbar(rvmockSH,dcmockSH,xerr=0.0,yerr=sdcmockSH,fmt='b.')
plt.plot(Rv,dclinfitSH(Rv),'b:',label=r'HOD sparse linear fit')
#plt.plot(Rv,dcinterpSH(Rv),'b--',label=r'HOD sparse cubic spline')
plt.plot(rvmockSHd,dcmockSHd,'m.',label=r'Sutter+2014 dense')
plt.errorbar(rvmockSHd,dcmockSHd,xerr=0.0,yerr=sdcmockSHd,fmt='m.')
plt.plot(Rv,dclinfitSHd(Rv),'m:',label=r'HOD dense linear fit')
#plt.plot(Rv,dcinterpSHd(Rv),'m--',label=r'HOD dense cubic spline')
plt.legend()
plt.ylabel('dc')
fig2 = plt.figure(2)
plt.plot(rvmock,rsmock,'r.',label=r'Hamaus+2014')
plt.errorbar(rvmock,rsmock,xerr=0.0,yerr=srsmock,fmt='r.')
plt.plot(Rv,rslinfit(Rv),'r:',label=r'DM full linear fit')
plt.plot(Rv,rsqadfit(Rv),'r-.',label=r'DM full quad fit')
#plt.plot(Rv,rsinterpS(Rv),'r--',label=r'DM full cubic spline')
plt.plot(rvmockSH,rsmockSH,'b.',label=r'Sutter+2014 sparse')
plt.errorbar(rvmockSH,rsmockSH,xerr=0.0,yerr=srsmockSH,fmt='b.')
plt.plot(Rv,rslinfitSH(Rv),'b:',label=r'HOD sparse linear fit')
#plt.plot(Rv,rsinterpSH(Rv),'b--',label=r'HOD sparse cubic spline')
plt.plot(rvmockSHd,rsmockSHd,'m.',label=r'Sutter+2014 dense')
plt.errorbar(rvmockSHd,rsmockSHd,xerr=0.0,yerr=srsmockSHd,fmt='m.')
plt.plot(Rv,rslinfitSHd(Rv),'m:',label=r'HOD dense linear fit')
#plt.plot(Rv,rsinterpSHd(Rv),'m--',label=r'HOD dense cubic spline')
plt.legend()
plt.ylabel('rs')
plt.show()
plt.close(fig1)
plt.close(fig2)

alpha = [0.0,0.5] #DM, gal -0.19
beta  = [1.86,4.16] #DM, gal

plt.figure()
for i in range(lr):
	print(Rv[i])
	#Hamaus DM
	#r, rho   = dprofile_Hamaus(dcinterp(Rv[i]),rsinterp(Rv[i]),Rv[i])
	r, rho   = dprofile_Hamaus(dclinfit(Rv[i]),rsqadfit(Rv[i]),Rv[i])
	norm     = trapz(rho,r)

	#Sutter DM sparse
	#r, rhoS  = dprofile_Hamaus(dcinterpS(Rv[i]),rsinterpS(Rv[i]),Rv[i]) 
	r, rhoS  = dprofile_Hamaus(dclinfitS(Rv[i]),rslinfitS(Rv[i]),Rv[i]) 
	normS    = trapz(rhoS,r)

	#Sutter HOD sparse
	#r, rhoSH = dprofile_Hamaus(dcinterpSH(Rv[i]),rsinterpSH(Rv[i]),Rv[i])
	r, rhoSH = dprofile_Hamaus(dclinfitSH(Rv[i]),rslinfitSH(Rv[i]),Rv[i])
	normSH   = trapz(rhoSH,r)

	#Sutter HOD dense
	r, rhoSHd = dprofile_Hamaus(dclinfitSHd(Rv[i]),rslinfitSHd(Rv[i]),Rv[i])
	normSHd   = trapz(rhoSHd,r)
	
	#CDF
	fcdf = open('Rv_CDF_'+str(round(Rv[i],2))+'_HamDM_SutDM_SutH_RicDM_RicGal.txt','w')
	CDF, CDFS, CDFSH, CDFSHd, CDF2, CDFg  = empty(len(r)), empty(len(r)), empty(len(r)), empty(len(r)), empty(len(r)), empty(len(r))
	
	#if Rv[i] < 24:
	for ir in range(len(r)):
		CDF[ir]   = trapz(rho[0:ir+1],r[0:ir+1])/norm
		CDFS[ir]  = trapz(rhoS[0:ir+1],r[0:ir+1])/normS
			#CDFSH[ir] = nan
		CDFSH[ir] = trapz(rhoSH[0:ir+1],r[0:ir+1])/normSH
		CDFSHd[ir] = trapz(rhoSHd[0:ir+1],r[0:ir+1])/normSHd
			
		#Ricciardelli et al 2014
		CDF2[ir] = CDFR(r[ir],Rv[i],alpha[0],beta[0])
		CDFg[ir] = CDFR(r[ir],Rv[i],alpha[1],beta[1])
		fcdf.write(str(r[ir]/Rv[i])+"\t"+str(CDF[ir])+"\t"+str(CDFS[ir])+"\t"+str(CDFSH[ir])+"\t"+str(CDF2[ir])+"\t"+str(CDFg[ir])+"\n")
	fcdf.close()

'''
	if Rv[i] == 48. or Rv[i] == 60. or Rv[i] == 70.:
		plt.plot(r/Rv[i],CDF,color=colors[i],marker='',ls='-',label=r'$R_v$ = '+str(int(Rv[i]))+' Mpc/$h$')
		plt.plot(r/Rv[i],CDFSH,color=colors[i],marker='',ls='--')#,label=r'$R_v$ = '+str(int(Rv[i]))+' Mpc/$h$')
		#plt.plot(r/Rv[i],CDFSHd,color=colors[i],marker='',ls='-.')
		#plt.plot(r/Rv[i],CDF2,color=colors[i],marker='',ls='--')
	if Rv[i] == 70.:
		plt.plot(r/Rv[i],CDFg,color='black',marker='',ls=':')	

#for legend
plt.plot([0,1],[-1,-1],'k-',label=r'DM Hamaus+ 2014')
#plt.plot([0,1],[-1,-1],'k-.',label=r'Sutter dense') #Sutter+2014
plt.plot([0,1],[-1,-1],'k--',label=r'Galaxies Sutter+ 2014') #Galaxies Sutter+2014 #Ricciardelli+ 2014 DM')
plt.plot([0,1],[-1,-1],'k:',label=r'Galaxies Ricciardelli+ 2014')
plt.xlim(0.0,1.0)
plt.ylim(0,1)
plt.legend(prop={'size':12},frameon=False)
plt.xlabel(r'$r/R_v$',fontsize=14)
plt.ylabel(r'CDF($\rho$)',fontsize=14)
plt.savefig( "CDF_rho_ofR_z=0.pdf",bbox_inches='tight')
plt.show()
'''
