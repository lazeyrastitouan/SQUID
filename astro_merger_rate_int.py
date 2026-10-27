from numpy import*
import argparse
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument('--zmin',nargs=1,type=float, default=0.15, help='zmin')
parser.add_argument('--zmax',nargs=1,type=float, default=0.43, help='zmax')
args     = parser.parse_args()

def integral(x,k):
	return (1+x)**(k+1)/(k+1)

zmin = args.zmin
zmax = args.zmax
z    = arange(zmin,zmax,0.001)

k    = 4.4
errk = array([1.9,-2.0])
#error(1+z)^k is (1+z)^k log(k) * errk

R02 = 32.4
errR02 = array([18.5,-12.2])

R0 = R02/(1.+0.2)**k  #at z=0.2 
errR0 = R0*sqrt((1.2**k*log(k)*errk/1.2**k)**2. + (errR02/R02)**2.)
R = R0*(integral(zmax,k)-integral(zmin,k))
print(R)

plt.plot(z,R0*(1+z)**k)
plt.show()
