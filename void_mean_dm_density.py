### compute the mean density (and thus Omega_m) in voids of various radius following the profile of Hamaus et al PRL 2014

from numpy import*
from scipy.interpolate import CubicSpline
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
colors=['black','yellow','orange','green','blue','magenta','red']
Rv = array([12.0,24.0,36.0,48.0,55.0,60.0,70.0])
lr = len(Rv)

# Hamaus 2014 DM
rvmock = array([11.860490045941805,16.322358346094944,21.43491577335375,26.547473200612558,34.26278713629402,46.53292496171517,69.4])
dcmock = array([-0.9479768786127167,-0.9190751445086706,-0.8439306358381503,-0.7341040462427746,-0.6271676300578035,-0.5520231213872833,-0.49132947976878616])
rsmock = array([0.7956678700361012,0.8635379061371842,0.940072202166065,1.0166064981949459,1.132129963898917,1.3155234657039712,1.6577617328519856])*rvmock
dcinterp = CubicSpline(rvmock,dcmock)
rsinterp = CubicSpline(rvmock,rsmock)

#Sutter 2014 
#DM
rvmockS = array([10.,15.,20.,25.,30.,35.])  
dcmockS = array([-0.80782,-0.77907,-0.78161,-0.78414,-0.76469,-0.72326])
rsmockS = array([0.76872,0.80900,0.86351,0.939697,1.02464,1.12773])*rvmockS
dcinterpS = CubicSpline(rvmockS,dcmockS)
rsinterpS = CubicSpline(rvmockS,rsmockS)

#HOD
rvmockSH = array([25.,30.,35.,40.,45.,50.,55.,60.,65.,70.,75.])
dcmockSH = array([-0.82389,-0.87378,-0.82135,-0.78330,-0.75116,-0.77400,-0.69958,-0.67336,-0.69873,-0.65729,-0.57696])
rsmockSH = array([0.75924,0.76517,0.77583,0.79834,0.80427,0.81848,0.83507,0.86469,0.83270,0.85284,0.83981])*rvmockSH
dcinterpSH = CubicSpline(rvmockSH,dcmockSH)
rsinterpSH = CubicSpline(rvmockSH,rsmockSH)

'''
fout = open("Omv_of_Rv.dat",'w')
fout.write("#Rv, Rv/f, dmv (mean delta_mv to Rv/f),  Omv (mean Om to Rv/f), Odc (center), Om_at_R\n")
'''

alpha = [-0.19,0.5] #DM, gal
beta  = [1.86,4.16] #DM, gal

plt.figure()
#plt.plot([0,2],[1,1],'k-')
for i in range(lr):
	print(Rv[i])
	#Hamaus DM
	r, rho   = dprofile_Hamaus(dcinterp(Rv[i]),rsinterp(Rv[i]),Rv[i])
	norm     = trapz(rho,r)

	#Sutter DM
	r, rhoS  = dprofile_Hamaus(dcinterpS(Rv[i]),rsinterpS(Rv[i]),Rv[i]) 
	normS    = trapz(rhoS,r)

	#Sutter HOD
	r, rhoSH = dprofile_Hamaus(dcinterpSH(Rv[i]),rsinterpSH(Rv[i]),Rv[i])
	normSH   = trapz(rhoSH,r)

	#CDF
	fcdf = open('Rv_CDF_'+str(int(Rv[i]))+'_HamDM_SutDM_SutH_RicDM_RicGal.txt','w')
	CDF, CDFS, CDFSH, CDF2, CDFg  = empty(len(r)), empty(len(r)), empty(len(r)), empty(len(r)), empty(len(r))
	
	if Rv[i] < 24:
		for ir in range(len(r)):
			CDF[ir]   = trapz(rho[0:ir+1],r[0:ir+1])/norm
			CDFS[ir]  = trapz(rhoS[0:ir+1],r[0:ir+1])/normS
			CDFSH[ir] = nan
			#Ricciardelli et al 2014
			CDF2[ir] = CDFR(r[ir],Rv[i],alpha[0],beta[0])
			CDFg[ir] = CDFR(r[ir],Rv[i],alpha[1],beta[1])
			fcdf.write(str(r[ir]/Rv[i])+"\t"+str(CDF[ir])+"\t"+str(CDFS[ir])+"\t"+str(CDFSH[ir])+"\t"+str(CDF2[ir])+"\t"+str(CDFg[ir])+"\n")
	elif 24 <= Rv[i] <= 36:
		for ir in range(len(r)):
			CDF[ir]   = trapz(rho[0:ir+1],r[0:ir+1])/norm
			CDFS[ir]  = trapz(rhoS[0:ir+1],r[0:ir+1])/normS
			CDFSH[ir] = trapz(rhoSH[0:ir+1],r[0:ir+1])/normSH
			#Ricciardelli et al 2014
			CDF2[ir] = CDFR(r[ir],Rv[i],alpha[0],beta[0])
			CDFg[ir] = CDFR(r[ir],Rv[i],alpha[1],beta[1])
			fcdf.write(str(r[ir]/Rv[i])+"\t"+str(CDF[ir])+"\t"+str(CDFS[ir])+"\t"+str(CDFSH[ir])+"\t"+str(CDF2[ir])+"\t"+str(CDFg[ir])+"\n")
	else:	
		for ir in range(len(r)):
			CDF[ir]   = trapz(rho[0:ir+1],r[0:ir+1])/norm
			CDFS[ir]  = nan
			CDFSH[ir] = trapz(rhoSH[0:ir+1],r[0:ir+1])/normSH
			#Ricciardelli et al 2014
			CDF2[ir] = CDFR(r[ir],Rv[i],alpha[0],beta[0])
			CDFg[ir] = CDFR(r[ir],Rv[i],alpha[1],beta[1])
			fcdf.write(str(r[ir]/Rv[i])+"\t"+str(CDF[ir])+"\t"+str(CDFS[ir])+"\t"+str(CDFSH[ir])+"\t"+str(CDF2[ir])+"\t"+str(CDFg[ir])+"\n")
	fcdf.close()

	#plt.plot(r/Rv[i],rho/rhom,color=colors[i],marker='',ls='--',label=r''+str(Rv[i]))
	if Rv[i] >= 48.:
		plt.plot(r/Rv[i],CDF,color=colors[i],marker='',ls='-',label=r'$R_v$ = '+str(int(Rv[i]))+' Mpc/$h$')
		plt.plot(r/Rv[i],CDFSH,color=colors[i],marker='',ls='--')
	plt.plot(r/Rv[i],CDFS,color=colors[i],marker='',ls='-.')
	plt.plot(r/Rv[i],CDF2,color=colors[i],marker='',ls='--')
	plt.plot(r/Rv[i],CDFg,color=colors[i],marker='',ls=':')	

'''
	#integrate over void radius, or a fraction of void radius f satisfying Rv/f > 12.0 to get mean rho in that region
	f = 1.
	while Rv[i]/f >= Rv[0]:
		filter   = r<=Rv[i]/f
		rrange   = r[filter]
		rhorange = rho[filter]
		mean_rho = trapz(rhorange,rrange)/(max(rrange)-min(rrange))
		Omv      = mean_rho/rho_crit
		Odc      = (dcinterp(Rv[i]) +1.)*rhom/rho_crit
		rhoatR   = rhom * (dc * ((1. - (max(rrange)/rs)**a)/(1. + (max(rrange)/Rv[i])**b)) + 1.)
		OmatR    = rhoatR/rho_crit
		fout.write(str(Rv[i])+"\t"+str(Rv[i]/f)+"\t"+str(mean_rho/rhom-1.)+"\t"+str(Omv)+"\t"+str(Odc)+"\t"+str(OmatR)+"\n")
		plt.plot(1./f,mean_rho/rhom-1.,color=colors[i],ls='',marker='.')

		f = f+1.

		del rrange, rhorange, mean_rho, Omv

	del dc, rs, a, b, r, rho
'''
#for legend
plt.plot([0,1],[-1,-1],'k-',label=r'Hamaus+ 2014')
plt.plot([0,1],[-1,-1],'k-.',label=r'Sutter+ 2014')
plt.plot([0,1],[-1,-1],'k--',label=r'Sutter+ 2014 HOD') #Ricciardelli+ 2014 DM')
plt.plot([0,1],[-1,-1],'k:',label=r'Ricciardelli+ 2014 gal')
plt.xlim(0.0,1.0)
plt.ylim(0,1)
plt.legend(prop={'size':12},frameon=False)
plt.xlabel(r'$r/R_v$',fontsize=14)
plt.ylabel(r'CDF($\rho$)',fontsize=14)
#plt.savefig( "CDF_rho_ofR_z=0.pdf",bbox_inches='tight')
plt.show()
