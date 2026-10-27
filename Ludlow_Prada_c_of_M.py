from numpy import*
import matplotlib.pyplot as plt
from scipy.integrate import quad

#params
z  = 0.55
a  = 1./(1.+z)
OL = 0.71
Om = 0.29

x  = (OL/Om)**(1./3.)*a

#get M and sigma(M,z)
dndM   = loadtxt("M_sigma_dndM_Tinker_Om0p29_OL0p71_s80p79_z0.txt").T
M      = dndM[0]
#sigma0 = dndM[1]

#get x, D and sigma
prefac     = lambda t : 5./2. * (Om/OL)**(1./3.) * sqrt(1.+t**3.)/t**1.5         
integrand  = lambda t : t**1.5/(1.+t**3.)**1.5
D          = lambda t : prefac(t) * quad(integrand,0,t)[0]
y          = 1e12/M
sigma      = D(x) * 16.9*y**0.41 / ( 1.+1.102*y**0.2+6.22*y**0.333 ) #sufficiently accurate and can go to M=1e3 
#sigma      = sigma0 * D(x)

#plt.figure(1)
#plt.loglog() 
#plt.plot(M,sigma,'k--',label='from mass function')
#plt.plot(M,sigmaPrada,'b--',label='Prada')
#plt.legend(loc=0)

#Ludlow16
dsc   = 1.686
nu    = dsc/sigma
coL   = 3.395*(1.+z)**(-0.215)
btaL  = 0.307*(1.+z)**0.54
g1    = 0.628*(1.+z)**(-0.047)
g2    = 0.317*(1.+z)**(-0.893)
nu0   = (4.135 - 0.564/a - 0.21/a**2. + 0.557/a**3 - 0.00348/a**4)/D(x)
cofnu = coL * (nu/nu0)**(-g1) * ( 1. + (nu/nu0)**(1./btaL) )**(-btaL*(g2-g1))

fL = open("M_c_Ludlow16_z"+str(z)+".txt",'w')
for i in range(len(M)):
        fL.write(str(M[i])+"\t"+str(cofnu[i])+"\n")
fL.close()

#Prada 2012:
#B0 and B1
c0    = 3.681
c1    = 5.033
afa   = 6.948
x0    = 0.424
s0inv = 1.047
s1inv = 1.646
bta   = 7.386
x1    = 0.526

cmin        = c0 + (c1-c0)*( 1./pi*arctan(afa*(x-x0)) + 0.5 )
sinvmin     = s0inv + (s1inv-s0inv)*( 1./pi*arctan(bta*(x-x1)) + 0.5 ) 
cmin1393    = c0 + (c1-c0)*( 1./pi*arctan(afa*(1.393-x0)) + 0.5 )
sinvmin1393 = s0inv + (s1inv-s0inv)*( 1./pi*arctan(bta*(1.393-x1)) + 0.5 )

B0 = cmin/cmin1393
B1 = sinvmin/sinvmin1393

#sprime and bigC
A = 2.881
b = 1.257
c = 1.022
d = 0.06

sprime = B1 * sigma
bigC   = A * ( (sprime/b)**c + 1. ) * exp(d/sprime**2.)

#concentration
cofM = B0*bigC

f=open("M_c_Prada12_z"+str(z)+".txt",'w')
for i in range(len(M)):
	f.write(str(M[i])+"\t"+str(cofM[i])+"\n")
f.close()

plt.figure(2)
plt.plot(log10(M),log10(cofM),'k--',label='Prada+ 2012')
plt.plot(log10(M),log10(cofnu),'b--',label='Ludlow+ 2016')
plt.legend(loc=0)
plt.show()
