#!/usr/bin/python

import sys
import argparse

import numpy as np
import matplotlib.pyplot as plt

import interpolated_spectra.reactor_tools as interp_rt
import ab_initio_spectra.reactor_tools as ab_initio_rt
import cenns_rate.CNSexperiment as CNSexperiment

def str2bool(v):
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

def parse_args(args):
    parser = argparse.ArgumentParser(description='Get settings')
    parser.add_argument('--lower-bound', type=float, default=0.01,
                        help='Lower bound for plot in keV')
    parser.add_argument('--upper-bound', type=float, default=10.0,
                        help='Upper bound for plot in keV')
    parser.add_argument('--interpolated', type=str2bool, default=True,
                        help='Whether interpolated spectra should be included')
    parser.add_argument('--ab-initio', type=str2bool, default=True,
                        help='Whether ab initio spectra should be included')
    return parser.parse_args()

if __name__ == '__main__':
    arg_ns = parse_args(sys.argv[1:])

    # Use average fuel fractions
    U235frac = 0.556;
    U238frac = 0.071;
    Pu239frac = 0.326;
    Pu241frac = 0.047;

    # Settings for CNSExperiment Objects
    detMass_ = 1
    bkg_ = np.zeros(100)
    T_bkg_ = np.linspace(0,100,100)
    time_ = 1*24*60*60 # output differential rate in units of per day
    # DC Reactor has power 8.54 GW and distance 400 m (400e2 cm)
    flux_ = interp_rt.nuFlux(8.54e3,400e2)

    spectra = list()
    labels = list()

    if(arg_ns.interpolated):
        def interp_dRdEnu(Enu):
            return U235frac*interp_rt.dRdEnu_U235(Enu)+\
                U238frac*interp_rt.dRdEnu_U238(Enu)+\
                Pu239frac*interp_rt.dRdEnu_Pu239(Enu)+\
                Pu241frac*interp_rt.dRdEnu_Pu241(Enu)
        spectra.append(interp_dRdEnu)
        labels.append("Interpolated")
    
    if(arg_ns.ab_initio):
        def ab_initio_dRdEnu(Enu):
            return U235frac*ab_initio_rt.dRdEnu_U235(Enu)+\
                U238frac*ab_initio_rt.dRdEnu_U238(Enu)+\
                Pu239frac*ab_initio_rt.dRdEnu_Pu239(Enu)+\
                Pu241frac*ab_initio_rt.dRdEnu_Pu241(Enu)
        spectra.append(ab_initio_dRdEnu)
        labels.append("Ab Initio")

    energies = np.logspace(np.log10(arg_ns.lower_bound),
                           np.log10(arg_ns.upper_bound))

    fig = plt.figure()
    for i in range(len(spectra)):
        # Germanium
        dRdEnu_=spectra[i]
        
        NGe_ = 40.63
        ZGe_ = 32
        myExptGe = CNSexperiment.CNSexperiment(N=NGe_, Z=ZGe_, dRdEnu=dRdEnu_,
                                               T_bg=T_bkg_, dRdT_bg=bkg_, detMass=detMass_,
                                               time=time_, nuFlux=flux_)
        plt.loglog(energies*1000.0, myExptGe.dRdT_CNS(energies),
                   label='Ge, %s'%labels[i])

        NZn_ = 35.38
        ZZn_ = 30
        myExptZn = CNSexperiment.CNSexperiment(N=NZn_, Z=ZZn_, dRdEnu=dRdEnu_,
                                               T_bg=T_bkg_, dRdT_bg=bkg_, detMass=detMass_,
                                               time=time_, nuFlux=flux_)
        plt.loglog(energies*1000.0, myExptZn.dRdT_CNS(energies),
                   label='Zn, %s'%labels[i])

        NSi_ = 14.085
        ZSi_ = 14
        myExptSi = CNSexperiment.CNSexperiment(N=NSi_, Z=ZSi_, dRdEnu=dRdEnu_,
                                               T_bg=T_bkg_, dRdT_bg=bkg_, detMass=detMass_,
                                               time=time_, nuFlux=flux_)
        plt.loglog(energies*1000.0, myExptSi.dRdT_CNS(energies),
                   label='Si, %s'%labels[i])

    ymin, ymax = plt.gca().get_ylim()
    plt.plot((50,50),(ymin,ymax),'k-')
    plt.plot((100,100),(ymin,ymax),'k-')
    plt.plot((200,200),(ymin,ymax),'k-')
    plt.plot((1000,1000),(ymin,ymax),'k-')
    plt.legend(loc=0,prop={'size':10})
    plt.xlabel('Recoil Energy [eV]')
    plt.ylabel('Differential Event Rate (evts/eV/kg/day)')
    plt.title('Differential Event Rate at Double Chooz')
    plt.savefig('results/cenns_rate.png')
    fig.clf()
    plt.close()
