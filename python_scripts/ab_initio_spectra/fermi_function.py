# Preliminary code towards calculating the spectrum with the ab initio approach
# This needs to be completed and interfaced with the code to calculate the ab initio spectrum

# Fermi function
def Fermi_fcn(Enu,Zf,Af):
        return 1.0
Fermi_fcn = np.vectorize(Fermi_fcn)

def Cfb(Enu):
        return 1.0
Cfb = np.vectorize(Cfb)

# RIGHT NOW THIS IS THE BETA SPECTRUM, NOT THE NEUTRINO SPECTRUM

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
S_f_fcn_unnormalized = np.vectorize(S_f_fcn_unnormalized)

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
