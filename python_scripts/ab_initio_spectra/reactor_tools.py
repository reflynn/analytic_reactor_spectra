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
#
# CURRENTLY THESE METHODS RETURN 1.0e-9, AND SERVE ONLY AS PLACEHOLDERS

import numpy as np
import ab_initio_spectra.fermi_function as fermi_function

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
	flux = power/200.0/1.602176565e-19 / (4*np.pi * distance**2.)
	return flux

def dRdEnu_U235(Enu):
	'''
	Reactor anti neutrino spectrum from U235 (see arXiv:1101.2663v3), per
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
	return 1.0e-6+0.0*Enu


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
	return 1.0e-6+0.0*Enu


def dRdEnu_Pu239(Enu):
	'''
	Reactor anti neutrino spectrum from Pu239 (see arXiv:1101.2663v3), per
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
	return 1.0e-6+0.0*Enu

def dRdEnu_Pu241(Enu):
	'''
	Reactor anti neutrino spectrum from Pu239 (see arXiv:1101.2663v3), per
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
	return 1.0e-6+0.0*Enu
