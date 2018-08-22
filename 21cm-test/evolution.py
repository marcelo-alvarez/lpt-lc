import numpy as np

omegab = 0.042
omegam = 0.30
h      = 0.70
YHe    = 0.25
T0     = 23.   # in mK

def growth(z):
    x  = 1+z
    x2 = x**2
    x3 = x**3

    omegal = 1 - omegam
    
    omega = omegam*x3/ (omegam*x3 + omegal)
    lambd = omegal   / (omegam*x3 + omegal)

    g  = 2.5*omega / (omega**( 4./7.)-lambd +(1.+omega /2.)*(1+lambd /70.))
    g0 = 2.5*omegam/ (omegam**(4./7.)-omegal+(1.+omegam/2.)*(1+omegal/70.))

    D = (g/x)/g0

    return D

def dtbofz(z):
    dtb = T0 * ( (omegab*h**2/0.02) *
                 np.sqrt((0.15/omegam/h**2)*((1+z)/10.)))
    return dtb

