# tritium.py
#
# Test case: calculating the spectra for Tritium
#
# Rionna Flynn
# 26 January 2018
# reflynn@mit.edu
#

import sys
import argparse

import numpy as np
import matplotlib.pyplot as plt

import fermi_function as ff
import reactor_tools as rt

Z = 1 # Atomic number
A = 3 # Mass number
E0 = 18.575 # Endpoint energy from arXiv:1409.0920

def S_f(Enu):
    return ff.S_f(Enu, Z, A, E0)

x = np.arange(0, E0, 0.01)
y = S_f(x)

fig = plt.figure()
plt.loglog(x, y, '-', label="Tritium, Ab Initio")
plt.xlabel("Neutrino Energy (keV)")
plt.ylabel("Neutrino Spectrum (neutrinos/keV/decay)")
plt.title("Tritium Neutrino Spectra")
plt.legend(loc=0, prop={'size':10})
plt.savefig("tritium_spectra.png")
fig.clf()
plt.close()
