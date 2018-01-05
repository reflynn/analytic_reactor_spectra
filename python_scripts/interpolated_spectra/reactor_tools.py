# ReactorTools.py
#
# Some tools for calculating the neutrino rate from a nuclear reactor.
#
# Adam Anderson
# 14 April 2016
# adama@fnal.gov
#
# Note: Convention on units:
#   --all masses are in kg
#   --all energies are in keV

import numpy as np
import scipy.interpolate as interp

import ROOT

def nuFlux(power, distance):
	'''
	Computes the total flux per fission of reactor antineutrinos at a given
	distance from the core, assuming a point-like flux, and nominal neutrino production

	Parameters
	----------
	power : float
		Reactor power in MW
	distance : float
		Distance in cm from reactor core at which flux is to be calculated

	Returns
	-------
	flux : float
		The reactor neutrino flux in fissions/s/cm^2 
	'''
        # Assume 200 MeV per fission
	flux = power/200.0/1.602176565e-19 / (4*np.pi * distance**2.)
	return flux

# Setup for the huber spectra
huber_setup_complete = False
huber_spl_loc = "data_files/interpolations/"
spl_U235 = 0
spl_Pu239 = 0
spl_Pu241 = 0
def spl_U235_eval(): return
def spl_Pu239_eval(): return
def spl_Pu241_eval(): return

def Huber_setup(file_U235=huber_spl_loc+'U235-anti-neutrino-flux-250keV.dat',
                file_Pu239=huber_spl_loc+'Pu239-anti-neutrino-flux-250keV.dat',
                file_Pu241=huber_spl_loc+'Pu241-anti-neutrino-flux-250keV.dat'):
        global huber_setup_complete
        global interp_min, interp_max
        global spl_U235, spl_Pu239, spl_Pu241
        global spl_U235_eval, spl_Pu239_eval, spl_Pu241_eval
        
        # U235
        enU235, specU235 = np.loadtxt(file_U235,usecols=(0,1),unpack=True)
        spl_U235_eval = interp.interp1d(x=enU235, y=specU235,
                                   bounds_error=False,
                                   fill_value=(specU235[0],specU235[-1]))
        
        # Pu239
        enPu239, specPu239 = np.loadtxt(file_Pu239,usecols=(0,1),unpack=True)
        spl_Pu239_eval = interp.interp1d(x=enPu239, y=specPu239,
                                   bounds_error=False,
                                   fill_value=(specPu239[0],specPu239[-1]))

        # Pu241
        enPu241, specPu241 = np.loadtxt(file_Pu241,usecols=(0,1),unpack=True)
        spl_Pu241_eval = interp.interp1d(x=enPu241, y=specPu241,
                                   bounds_error=False,
                                   fill_value=(specPu241[0],specPu241[-1]))

        huber_setup_complete = True
        return

def dRdEnu_U235(Enu):
        # check global setup
        global huber_setup_complete
        global spl_U235_eval
        if(not huber_setup_complete):
                Huber_setup()

	if type(Enu) == float:
		Enu = np.asarray([Enu])
	else:
		Enu = np.asarray(Enu)
        # input is in keV, huber spline expects MeV
        # huber spline gives results in 1/Mev/fission, we want 1/keV/fission
        spec = spl_U235_eval(Enu*1e-3)*1e-3
        spec[Enu<2.e3] = spl_U235_eval(2.0)*1e-3
        spec[spec<0] = 0
        return spec

def dRdEnu_Pu239(Enu):
        # check global setup
        global huber_setup_complete
        global spl_Pu239_eval
        if(not huber_setup_complete):
                Huber_setup()

	if type(Enu) == float:
		Enu = np.asarray([Enu])
	else:
		Enu = np.asarray(Enu)
        # input is in keV, huber spline expects MeV
        # huber spline gives results in 1/Mev/fission, we want 1/keV/fission
        spec = spl_Pu239_eval(Enu*1e-3)*1e-3
        spec[Enu<2.e3] = spl_Pu239_eval(2.0)*1e-3
        spec[spec<0] = 0
        return spec

def dRdEnu_Pu241(Enu):
        # check global setup
        global huber_setup_complete
        global spl_Pu241_eval
        if(not huber_setup_complete):
                Huber_setup()

	if type(Enu) == float:
		Enu = np.asarray([Enu])
	else:
		Enu = np.asarray(Enu)
        # input is in keV, huber spline expects MeV
        # huber spline gives results in 1/Mev/fission, we want 1/keV/fission
        spec = spl_Pu241_eval(Enu*1e-3)*1e-3
        spec[Enu<2.e3] = spl_Pu241_eval(2.0)*1e-3
        spec[spec<0] = 0
        return spec

# The fit from Mueller is used for U-238
def dRdEnu_U238(Enu):
	'''
	Reactor anti neutrino spectrum from U238 (see arXiv:1101.2663v3), per
	fission

	Parameters
	----------
	Enu : array
		Neutrino energy in keV

	Returns
	-------
	spectrum : array
		Spectrum [nu / keV / fission]
	'''
	if type(Enu) == float:
		Enu = np.asarray([Enu])
	else:
		Enu = np.asarray(Enu)
	EnuMeV = Enu / 1.e3
	spectrum = 1e-3 * np.exp((4.833e-1) + (1.927e-1)*EnuMeV - (1.283e-1)*EnuMeV**2.0 - \
						(6.762e-3)*EnuMeV**3.0 + (2.233e-3)*EnuMeV**4.0 - (1.536e-4)*EnuMeV**5.0)
	#spectrum[EnuMeV<1.0] = 1e-3 * np.exp((4.833e-1) + (1.927e-1)*1.0 - (1.283e-1)*1.0**2.0 - \
	#					(6.762e-3)*1.0**3.0 + (2.233e-3)*1.0**4.0 - (1.536e-4)*1.0**5.0)
	return spectrum
