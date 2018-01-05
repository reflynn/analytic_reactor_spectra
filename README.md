# analytic_reactor_spectra
Neutrino spectra at a nuclear reactor, calculated using an analytic formula for the spectra of each isotope

# Overview

Previously, in "Coherent Neutrino Scattering with Low Temperature Bolometers at Chooz Reactor Complex" (arXiv:1612.09035 [physics.ins-det]) We used the reactor neutrino spectra from Mueller 2011 (Mueller T A, et al. 2011 Phys. Rev. C 83(5)) for U-238, and we used the spectrum from Huber 2011 (Huber P 2011 Phys. Rev. C 84(2) 024617) for U-235, Pu-239, and Pu-241. These spectra were calculated using a combination two approaches. The first was an "ab initio" approach, where the spectrum for one isotope (eg U-235) is calculated as the sum of the spectra from each of resulting fission products. The second method used electron flux from the high flux Institut Laue-Langevin (ILL) reactor in  order to compute the spectrum from data.

The combination of the "ab initio" approach with ILL data gave a precise spectrum. However, the result in the above papers were only valid in the range ~2 to 8 MeV.The lower bound was chosen because in most antineutrino expermients, the detection channel is ν¯e + p → e+ + n, which has a threshold of 1.8 MeV. However, coherent elastic neutrino-nucleus scattering (CENNS) has lower thresholds (eg ~0.3 MeV for Si). Previously, we assumed the spectrum was constant below 2 MeV. Here we would like to use the "ab initio" approach in order to calculate the spectrum below 2 MeV.

The Mueller 2011 paper walks through the ab initio calculation in detail. We have the MURE (MCNP Utility for Reactor Evolution) code used by the paper, and we will use it to calculate the spectra.

# Usage

- python python_scripts/plot_spectra.py

Creates a plot with the interpolated neutrino spectrum from the Mueller and Huber papers, and with the analytically calculated ab initio spectrum. --lower-bound and --upper-bound give the bounds for the plot in keV. --interpolated and --ab-initio take a boolean argument to specify whether the spectrum should be plotted.

- python python_scripts/plot_dc_cenns_rate.py 

Creates a plot of the resulting cenns rate. lower and upper bound are in keV.

# Code Organization

#### data_files:

Any data files used in the calculation. Currently, the text files with the spectra from the Mueller and Huber papers are stored here. The results from the MURE code will eventually be stored here.

#### results:

Resulting plots will be stored here

#### python_scripts:

Contains the code used to calculate the spectra, and the resulting differential CENNS cross sections.

- ab_initio_spectra: Code to analytically calculate the spectra using the results from the MURE code.
- interpolated_spectra: Code to generate the spectra from the Muellera and Huber papers
- cns_rate: Code to calculate the CENNS differential and total rate
- plot_results.py: Create plots of the spectra and differential rate, and print the total rate.
