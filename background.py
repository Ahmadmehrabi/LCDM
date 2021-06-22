import numpy as np
import scipy.integrate as integ


#Here the model is defined. This is the normalized Hubble function 

def normalized_hubble(z,args):
        om0 = args[0]
        omm = om0*pow(1.+z,3)
        ol0 = 1 - om0 

        return pow(omm + ol0 ,0.5)


def Ez_inv(z,args):
        return 1./normalized_hubble(z,args)


# This is the comoving distance
def com_dis(z,args):
        res = integ.quad(Ez_inv, 0., z, args=args)
        return res[0]


def lom_dist(z,args):
        return (1+z)*com_dis(z,args)
