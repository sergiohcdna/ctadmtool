import matplotlib.pyplot as plt
import numpy as np
from ebltable.tau_from_model import OptDepth

from auxfunctions import dminterpol , auxMan

import argparse

###########################################
#       This script generate files
#       with dm flux from annihilations
#       or decay of dm particles.
#       The format of the file is the same
#       as required by gammalib
###########################################
#
#       Parameters:
#           - mass
#           - channel
#           - energy range
#           - cross--section or lifetime
#           - process
#           - jfactor
###########################################
#
#           Sergio, September 2020
###########################################

if __name__ == '__main__' :

    options = argparse.ArgumentParser( description='Script to \
        generate files with dm flux for annihilation \
        or decay of DM particles' )

    dmmodel = options.add_argument_group( 'DM Model' ,\
        'Information about DM candidate' )

    dmmodel.add_argument( '--gammafile' ,\
        help='File with gamma production Table from PPPC4 project' ,\
        type=str ,\
        required=True ,\
        metavar='AtProduction_gammas.dat' )
    dmmodel.add_argument( '--mass' , \
        help='DM-candidate mass (in GeV)' ,\
        type=float ,\
        required=True ,\
        metavar='100.0' )
    dmmodel.add_argument( '--estart' ,\
        help='Energy to start to compute the spectra (in GeV)' ,\
        type=float ,\
        default=50.0 ,\
        metavar='50.0' )
    dmmodel.add_argument( '--sigmav' ,\
        help='Annihilation cross-section(cm**3/s)\
              or lifetime (s) according to the process' ,\
        type=float ,\
        required=True ,\
        metavar='1.e-28' )
    dmmodel.add_argument( '--jfactor' ,\
        help='Astrophysical factor according to the process,\
              (in GeV**2/cm**5) or (in GeV/cm**3)' ,\
        type=float ,\
        required=True ,\
        metavar='1.e+19' )
    dmmodel.add_argument( '--process' ,\
        help='DM process to produce gamma-rays' ,\
        type=str ,\
        default='annihilation' ,\
        metavar='annihilation' ,\
        choices=[ 'annihilation' , 'decay' ] )
    dmmodel.add_argument( '--channel' ,\
        help='Process channel to produce gamma-rays' ,\
        type=str ,\
        default='Tau' ,\
        metavar='Tau' )
    dmmodel.add_argument( '--z' ,\
        help='redshift z of the target' ,\
        type=float ,\
        default=0.0 ,\
        metavar='2.0' )
    dmmodel.add_argument( '--eblmodel' ,\
        help='EBL Model. Please check EBLTable doc to check\
            what models are available' ,\
        type=str ,\
        default='franceschini' ,\
        metavar='franceschini' )
    dmmodel.add_argument( '--npoints' ,\
        help='Number of points to compute the spectra' ,\
        type=int ,\
        default=100 ,\
        metavar='100' )
    dmmodel.add_argument( '--dfile' ,\
        help='Name of file to save data' ,\
        type=str ,\
        required=True ,\
        metavar='DMFlux.txt' )
    dmmodel.add_argument( '--outpath' ,\
        help='Output directory to save file with data' ,\
        type=str ,\
        default='./' ,\
        metavar='/path/to/directory/' )
    dmmodel.add_argument( '--plot' ,\
        help='Plot dmSpectrum' ,\
        action='store_true' )

    args = options.parse_args()

    #   Load the DM interpolator
    ewdm_int = dminterpol.pppcDM_EWinterp( args.gammafile , args.channel )

    #   Define the array with x values
    #   Remember: x = energy / DMmass
    x_start = args.estart / args.mass
    xvalues = np.logspace(  np.log10( x_start ) , 0 , args.npoints )

    thisppfactor = 0

    #   Check the process and compute the pp-factor
    if args.process == 'annihilation' :

        thisppfactor = dminterpol.ppfactor_ann( args.sigmav , args.mass )

    elif args.process == 'decay' :

        thisppfactor = dminterpol.ppfactor_dec( args.sigmav , args.mass )

    else :

        print( 'Unknown process. I assume you mean annihilation' )
        thisppfactor = dminterpol.ppfactor_ann( args.sigmav , args.mass )

    if args.z > 1.e-3 :

        #   Create Optical-Depth instance
        tau = OptDepth.readmodel( model='franceschini' )

        #   Computing EBL attenuation for E_start up to DM mass
        atten   = np.exp( -1. * tau.opt_depth( args.z , xvalues * args.mass * 1.e-3 ) )

    else :

        #   Redshift is below minimum value available
        #   Then, attenuation is equal to indentity
        atten = 1.0

    #   Computing the gamma-ray flux from dm
    flux  = ewdm_int( args.mass , xvalues )
    flux_atten = flux[ 0 ] * atten * thisppfactor * args.jfactor

    #   Converting to MeV
    #   By default, gammalib use MeV as energy unit
    eng_MeV  = xvalues * args.mass * 1.e+3
    flux_MeV = flux_atten * 1.e-3

    #   Save to file
    data    = np.array( ( eng_MeV , flux_MeV ) ).transpose()
    outfile = auxMan.createname( args.outpath , args.dfile )
    header  = ' \tE[MeV]\tFlux[1/MeV/cm**2/s]'

    np.savetxt( outfile , data , fmt='%.6e' , delimiter=' ' )

    if args.plot :

        fig , ax = plt.subplots( figsize=( 9 , 6 ) )

        ax.plot( xvalues * args.mass , flux , color=( 0.82 , 0.1 , 0.26 ) , \
            label='DM flux without EBL' , lw=2 )
        ax.plot( xvalues * args.mass , flux_atten , color=( 0.57 , 0.36 , 0.51 ) , \
            label='DM flux with EBL' , lw=2 )

        ax.set_xlim( args.estart , 10 * args.mass )
        ax.set_xlabel( 'E(GeV)' )
        ax.set_ylabel( '$\\frac{d\\Phi_{\\gamma}}{dE} (ph cm^{-2} s^{-1} GeV^{-1})$' )
        ax.set_xscale( 'log' )
        ax.set_yscale( 'log' )
        ax.legend( loc='best' , prop={ 'size': 10 } )
        pfile = 'DMflux.png'
        plt.savefig( pfile )

