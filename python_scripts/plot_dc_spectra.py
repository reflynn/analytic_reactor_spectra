#!/usr/bin/python

import sys
import argparse

import numpy as np
import matplotlib.pyplot as plt

import interpolated_spectra.reactor_tools as interp_rt
import ab_initio_spectra.reactor_tools as ab_initio_rt

def str2bool(v):
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

def parse_args(args):
    parser = argparse.ArgumentParser(description='Get settings')
    parser.add_argument('--lower-bound', type=float, default=1.0e2,
                        help='Lower bound for plot in keV')
    parser.add_argument('--upper-bound', type=float, default=8.e3,
                        help='Upper bound for plot in keV')
    parser.add_argument('--interpolated', type=str2bool, default=True,
                        help='Whether interpolated spectra should be included')
    parser.add_argument('--ab-initio', type=str2bool, default=True,
                        help='Whether ab initio spectra should be included')
    return parser.parse_args()

if __name__ == '__main__':
    arg_ns = parse_args(sys.argv[1:])

    energies = np.logspace(np.log10(arg_ns.lower_bound),
                           np.log10(arg_ns.upper_bound))

    fig = plt.figure()
    if(arg_ns.interpolated):
        plt.loglog(energies, interp_rt.dRdEnu_U235(energies), '-', label="U235, Interpolated")
        plt.loglog(energies, interp_rt.dRdEnu_U238(energies), '-', label="U238, Interpolated")
        plt.loglog(energies, interp_rt.dRdEnu_Pu239(energies), '-', label="Pu239, Interpolated")
        plt.loglog(energies, interp_rt.dRdEnu_Pu241(energies), '-', label="Pu241, Interpolated")
    if(arg_ns.ab_initio):
        plt.loglog(energies, ab_initio_rt.dRdEnu_U235(energies), '--', label="U235, Ab Initio")
        plt.loglog(energies, ab_initio_rt.dRdEnu_U238(energies), '--', label="U238, Ab Initio")
        plt.loglog(energies, ab_initio_rt.dRdEnu_Pu239(energies), '--', label="Pu239, Ab Initio")
        plt.loglog(energies, ab_initio_rt.dRdEnu_Pu241(energies), '--', label="Pu241, Ab Initio")

    plt.xlabel("Neutrino Energy (keV)")
    plt.ylabel("Neutrino Spectrum (neutrinos/keV/fission)")
    plt.title("Reactor Neutrino Spectra")
    plt.legend(loc=0,prop={'size':10})
    plt.savefig("results/spectra.png")
    fig.clf()
    plt.close()
