import numpy as np
import scipy.integrate as integ


#Here the model is defined. This is the normalized Hubble function 

def Hubble(z,H0,om):
        omm = om*pow(1.+z,3)
        ol0 = 1 - om 
        return H0*pow(omm + ol0 ,0.5)


def Hub_inv(z,H0,om):
        return 1./Hubble(z,H0,om)


# This is the comoving distance
def com_dis(z,args):
        res = integ.quad(Hub_inv, 0., z, args=args)
        return 3e5*res[0]


def mu_the(z,H0,om):
        mu = np.zeros(len(z))
        args = [H0,om]
        for i in range(len(z)):
                lsd = (1+z[i])*com_dis(z[i],args)
                mu[i] = 5*np.log10(lsd) + 25.
        return mu
