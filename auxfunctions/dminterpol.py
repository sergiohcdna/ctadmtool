########################################################
#       Generate a interpolation function to
#       compute dm spectra based on
#       Cirelli PPPC4-DM
#       This should work with any version of
#       gamma-ray production tables.
#       But, I put a specific file in the same
#       directory.      :)
########################################################
#
#       shkdna, september, 2020
#
########################################################

import numpy as np
from scipy.interpolate import RegularGridInterpolator

import sys

def pppcDM_EWinterp( channel ) :

    #   Load data from PPPC4 tables
    #   The file is in the same
    #   directory as the location of
    #   this script

    #   The function create a grid with
    #   values of dN/dE (no dN/dLog10x)
    #   And the interpolation function
    #   is created using :
    #       - DM candidate mass (in GeV)
    #       - ratio of gamma-ray energy and mass

    #   The tables used in this sript
    #   correspond to the data with
    #   electroweak corrections :)
    fname = 'AtProduction_gammas.dat'
    data  = np.genfromtxt( fname , names=True , dtype=None )

    #   Tuple with column names
    #   from PPPC4 tables
    #   This include also the channels
    dchannels = data.dtype.names

    #   check if channel is in dchannels
    #   if not, then exit
    if channel not in dchannels :
        print( 'Channel not available,\
            in PPPC4-DM tables.\n\
            The available channels are:' )
        print( dchannels )
        sys.exit( 'Unknown channel' )

    #   Get values of DM mass and
    #   x_steps to compute the 
    #   number of gamma-rays
    masses  = np.unique( data[ 'mDM' ] )

    xvalues = np.unique( data[ 'Log10x' ] )
    xvalues = np.power( 10 , xvalues )

    #   This is needed to create the
    #   interpolating function
    dmgrid  = np.zeros( ( masses.size , xvalues.size ) )

    #   Filling the grid
    for m_index in range( masses.size ) :

        mass    = masses[ m_index ]

        #   Get the indices where data[ 'mDM']
        #   is equal to mass
        indices = np.where( data[ 'mDM'] == mass )

        #   Loop to extract dm_flux
        for index in indices :
            for x_index in range( xvalues.size ) :

                xvalue = xvalues[ x_index ]

                #   Getting the flux
                flux = data[ channel ][ index ][ x_index ]
                flux = flux / xvalue / mass / np.log( 10 )

                #   filling the grid
                dmgrid[ m_index ][ x_index ] = flux

    interpolator = RegularGridInterpolator( [ masses , xvalues ] ,\
        dmgrid , bounds_error=False , fill_value=0 )

    return interpolator

def ppfactor_ann( sigmav , dmmass ) :

    #   This function compute the
    #   Particle Physics term to
    #   get the gamma-ray flux
    #   from annihilations
    ppfactor = sigmav / 8 / np.pi / np.power( dmmass , 2 )

    return ppfactor

def ppfactor_dec( lifetime , dmmass ) :

    #   This function compute the
    #   Particle Physics term to
    #   get the gamma-ray flux
    #   from decay
    ppfactor = 1 / 4 / np.pi / dmmass / lifetime

    return ppfactor


