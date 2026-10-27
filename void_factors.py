from numpy import*


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
