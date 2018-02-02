# fermi_function.py
#
# Code to calculate the beta spectrum with the ab initio approach

import numpy as np
import scipy.integrate as spint

from cmath import sin, sqrt, pi, exp

# Complex gamma function

p = [676.5203681218851
        ,-1259.1392167224028
        ,771.32342877765313
        ,-176.61502916214059
        ,12.507343278686905
        ,-0.13857109526572012
        ,9.9843695780195716e-6
        ,1.5056327351493116e-7
        ]

EPS = 1e-07
def drop_imag(z):
    if abs(z.imag) <= EPS:
        z = z.real
    return z

def gamma(z):
    z = complex(z)
    if z.real < 0.5:
        y = pi / (sin(pi*z) * gamma(1-z))  # Reflection formula
    else:
        z -= 1
        x = 0.99999999999980993
        for (i, pval) in enumerate(p):
            x += pval / (z+i+1)
        t = z + len(p) - 0.5
        y = sqrt(2*pi) * t**(z+0.5) * exp(-t) * x
    return drop_imag(y)

# Fermi function
def Fermi_fcn(Enu,Zf,Af):
    '''
    Computes the Fermi function.

    Parameters
    ----------
    Enu : float
            Kinetic energy in KeV
    Zf : int
            Charge of the parent nucleus
    Af : int
            Atomic number of the parent nucleus 
    Returns
    -------
    F : float
            Result of Fermi function
    '''

    m_e = 511  # units of KeV 
    W = Enu / m_e + 1

    a = 0.0072973525664  # the dimensionless fine structure constant
    g = np.sqrt(1 - alpha*Zf**2)

    r0 = 1.25e-15  # units of meters
    R = r0*Af**(1/3)  # nuclear radius in meters

    p = np.sqrt(Enu**2 - m_e**2)  # momentum in KeV

    z = complex(g, a*Zf*W/p)

    F = 2*(g+1)*((2*p*R)**(2*(g-1)))*np.exp(math.pi*a*Zf*W/p)*(abs(gamma(z))**2)/(gamma((2*g+1)**2))

    return F

Fermi_fcn = np.vectorize(Fermi_fcn)

def Cfb(Enu):
    return 1.0

Cfb = np.vectorize(Cfb)

# S_f calculates the neutrino spectrum for a given fission product f and
# branch b. Normalized to 1, so Int(S_f dE) = 1.
def S_f_fcn_unnormalized(Enu,Zf,Af,E0fb):
    m_e = 511
    E = E0fb-Enu
    if(E<m_e or E>E0fb):
        return 0.0
    p = np.sqrt(E**2-m_e**2)
    return p*E*(E-E0fb)**2*Fermi_fcn(E,Zf,Af)*Cfb(E)
S_f_fcn_unnormalized = np.vectorize(S_f_fcn_unnormalized)

def S_f_fcn_unnormalized_beta(Enu,Zf,Af,E0fb):
    m_e = 511
    if(Enu<m_e or Enu>E0fb):
        return 0.0
    p = np.sqrt(Enu**2-m_e**2)
    return p*Enu*(Enu-E0fb)**2*Fermi_fcn(Enu,Zf,Af)*Cfb(Enu)
S_f_fcn_unnormalized_beta = np.vectorize(S_f_fcn_unnormalized_beta)

def S_f(Enu,Zf,Af,E0fb):
    S_arr = S_f_fcn_unnormalized(Enu,Zf,Af,E0fb)
    norm = spint.quad(lambda Enu_: \
            S_f_fcn_unnormalized(Enu_,Zf,Af,E0fb),\
            0,E0fb)[0]
    return S_arr/norm

def S_f_beta(Enu,Zf,Af,E0fb):
    S_arr = S_f_fcn_unnormalized_beta(Enu,Zf,Af,E0fb)
    norm = spint.quad(lambda Enu_: \
            S_f_fcn_unnormalized(Enu_,Zf,Af,E0fb),\
            0,E0fb)[0]
    return S_arr/norm
