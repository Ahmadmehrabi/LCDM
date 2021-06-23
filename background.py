import numpy as np
import scipy.integrate as integ


#Here the model is defined. This is the normalized Hubble function 

def Hubble(z,args):
        H0 = args[0]
        om0 = args[1]
        omm = om0*pow(1.+z,3)
        ol0 = 1 - om0 

        return H0*pow(omm + ol0 ,0.5)


def Hub_inv(z,args):
        return 1./Hubble(z,args)


# This is the comoving distance
def com_dis(z,args):
        res = integ.quad(Hub_inv, 0., z, args=args)
        return 3e3*res[0]


def lom_dis(z,args):
        return (1+z)*com_dis(z,args)
