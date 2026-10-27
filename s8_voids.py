#COMPUTE sigma8 in void following Gottlober et al. 2003 eqs (6)-(8) 

import numpy as np
from scipy.special import hyp2f1
from astropy import cosmology
from scipy.integrate import quad
import matplotlib.pyplot as plt


#fiducial flat cosmo
a = 1./1.3 #0.7692307692307692 #sigma_8 is specified for 0. We will rescale the CAMB power spectrum linearly to 0.3 
h = 0.6831
cosmo = cosmology.Planck18.clone(H0=100.*h, Om0=0.29)
cosmo.meta['sigma8']=0.678 #0.79 at z=0, 0.678 at z=0.3, 0.598 at z=0.55
Om = cosmo.Om0
zi = 1000.
ai = 1./(1.+zi)
Df = lambda x : x*hyp2f1(1./3.,1.,11./6.,x**3.*(Om-1)/Om) #Silveira & waga 1994 corrected -> check how this compares with Komatsu's code
Omv = [0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.11,0.12,0.13,0.14,0.15,0.29]
OL  = 0.71
dc  = 1.62

#void cosmo
fout=open("Omv_s8v_Gottlober03_z0p3.txt",'w')
for om in Omv:
	f = lambda x : (1.+om*(1./x-1.)+OL*(x**2.-1.))**(-0.5)
	f3 = lambda x : (1.+om*(1./x-1.)+OL*(x**2.-1.))**(-1.5)
	Dv = lambda x : 5.*om/(2.*x*f(x))*quad(f3,0,x)[0]

	print(Dv(a)*dc/Dv(1))
	
	s8v = cosmo.meta['sigma8']*Df(ai)/Df(a)*Dv(a)/Dv(ai)
	fout.write(str(om)+"\t"+str(s8v)+"\n")
fout.close()
