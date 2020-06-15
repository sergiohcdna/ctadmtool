#! /usr/bin/env python

# ====================================================================
# This script perform the DM limits using the Cirelli et al. spectrum.
# ====================================================================

import ctools
import gammalib
import cscripts

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import SymLogNorm
import math
import numpy as np
import scipy

import glob
import os
import sys
import csv
import time
import argparse

from auxfunctions import *

#=====================#
# Routine entry point #
#=====================#

def ctamap( obs , emin , emax , \
    coordsys , ra , dec , rad , \
    caldb , irf , name , plotfile ) :
    skymap = ctools.ctskymap() 
    skymap[ 'inobs' ] = obs
    skymap[ 'emin' ] = emin
    skymap[ 'emax' ] = emax
    skymap[ 'nxpix' ] = 100
    skymap[ 'nypix' ] = 100
    skymap[ 'binsz' ] = 0.02
    skymap[ 'proj' ] = 'CAR'
    skymap[ 'coordsys' ] = coordsys
    skymap[ 'xref' ] = ra
    skymap[ 'yref' ] = dec
    skymap[ 'bkgsubtract' ] = 'IRF'
    skymap[ 'caldb' ] = caldb
    skymap[ 'irf' ] = irf
    skymap[ 'outmap' ] = name

    skymap.execute()

    skymap.skymap().smooth( 'GAUSSIAN', 0.02 )

    ax = plt.subplot()
    plt.imshow( skymap.skymap().array() , \
        origin='lower' , \
        extent=[ ra + rad , \
            ra - rad , \
            dec - rad , \
            dec + rad ] , \
        norm=SymLogNorm( 1 , base=10 ) )
    xlabel = '' ; ylabel = '' ;
    if coordsys == 'GAL' :
        xlabel = 'Longitude (deg)'
        ylabel = 'Latitude (deg)'
    elif coordsys == 'CEL' :
        xlabel = 'R.A. (deg)'
        ylabel = 'DEC (deg)'
    else :
        print( 'Unknown coordinate System. \
            I will assume Celestial Coordinate System' )
        xlabel = 'R.A. (deg)'
        ylabel = 'DEC (deg)'
    ax.set_xlabel( xlabel )
    ax.set_ylabel( ylabel )
    cbar = plt.colorbar()
    cbar.set_label( 'Counts' )

    if len( plotfile ) > 0 :
        plt.savefig( plotfile )
    else :
        plt.show()
    plt.close()

def getONOFFCounts( onoffobs , models , ebounds , with_Components ) :

    bkgclasses = [ 'GCTAModelAeffBackground' , \
        'GCTAModelBackground' , 'GCTAModelCubeBackground' , \
        'GCTAModelIrfBackground' , 'GCTAModelRadialAcceptance' ]

    names = [ source.name() for source in models if source.classname() not in bkgclasses ]
    typelist = [ ( 'Total' , float ) , ( 'Signal' , float ) , ( 'TotalOff' , float ) ]
    for name in names :
        typelist.append( (  name , float ) )

    n_modelcounts = np.zeros( ebounds.size() , dtype=typelist )

    counts       = onoffobs.on_spec().counts_spectrum()
    background   = onoffobs.model_background( models ).counts_spectrum()
    alpha        = onoffobs.on_spec().backscal_spectrum()
    modelcounts  = background.copy()
    modelcounts *= alpha
    modelcounts += onoffobs.model_gamma( models ).counts_spectrum()

    coff         = onoffobs.off_spec().counts_spectrum()
    # counts_off   = onoffobs.off_spec().counts_spectrum()
    for i in range( ebounds.size() ) :
        n_modelcounts[ 'Total' ][ i ]    = counts[ i ]
        n_modelcounts[ 'Signal' ][ i ]   = modelcounts[ i ]
        n_modelcounts[ 'TotalOff' ][ i ] = coff[ i ]

    if with_Components :
        for source in models :
            if source.classname() not in bkgclasses :
                modelContainer = gammalib.GModels()
                modelContainer.append( source )
                modelsource = onoffobs.model_gamma( modelContainer )
                modelsource = modelsource.counts_spectrum()
                for i in range( ebounds.size() ) :
                    n_modelcounts[ source.name() ][ i ] = modelsource[ i ]

    return n_modelcounts

def getONOFFCountsAfterLike( onoffobs , models , ebounds , with_Components ) :

    bkgclasses = [ 'GCTAModelAeffBackground' , \
        'GCTAModelBackground' , 'GCTAModelCubeBackground' , \
        'GCTAModelIrfBackground' , 'GCTAModelRadialAcceptance' ]

    names = [ source.name() for source in models if source.classname() not in bkgclasses ]
    typelist = [ ( 'Total' , float ) , ( 'Model' , float ) , ( 'Background' , float ) , ( 'TotalOff' , float ) ]
    for name in names :
        typelist.append( (  name , float ) )

    n_modelcounts = np.zeros( ebounds.size() , dtype=typelist )

    counts       = onoffobs.on_spec().counts_spectrum()

    background   = onoffobs.model_background( models ).counts_spectrum()
    alpha        = onoffobs.on_spec().backscal_spectrum()
    modelcounts  = background.copy()
    modelcounts *= alpha
    modelcounts += onoffobs.model_gamma( models ).counts_spectrum()
    bkg          = background.copy()
    backscal     = onoffobs.on_spec().backscal_spectrum()
    bkg         *= backscal

    coff         = onoffobs.off_spec().counts_spectrum()

    # counts_off   = onoffobs.off_spec().counts_spectrum()
    for i in range( ebounds.size() ) :
        n_modelcounts[ 'Total' ][ i ]      = counts[ i ]
        n_modelcounts[ 'Model' ][ i ]      = modelcounts[ i ]
        n_modelcounts[ 'Background' ][ i ] = bkg[ i ]
        n_modelcounts[ 'TotalOff' ][ i ]   = coff[ i ]

    if with_Components :
        for source in models :
            if source.classname() not in bkgclasses :
                modelContainer = gammalib.GModels()
                modelContainer.append( source )
                modelsource = onoffobs.model_gamma( modelContainer )
                modelsource = modelsource.counts_spectrum()
                for i in range( ebounds.size() ) :
                    n_modelcounts[ source.name() ][ i ] = modelsource[ i ]

    return n_modelcounts

def getEmptyCountsArray( enumbins , srcnames ) :
    typelist = [ ( 'Total' , float ) , ( 'TotalSigmaUp' , float ) , ( 'TotalSigmaDn' , float ) , ( 'Signal' , float ) , ( 'SignalSigmaUp' , float ) , ( 'SignalSigmaDn' , float ) , ( 'TotalOff' , float ) , ( 'TotalOffSigmaUp' , float ) , ( 'TotalOffSigmaDn' , float ) ]
    for name in srcnames :
        typelist.append( (  name , float ) )
        typelist.append( (  name + 'SigmaUp' , float ) )
        typelist.append( (  name + 'SigmaDn' , float ) )

    modelcounts = np.zeros( enumbins , dtype=typelist )
    return modelcounts

def getEmptyCountsLikeArray( enumbins , srcnames ) :
    typelist = [ ( 'Total' , float ) , ( 'TotalSigmaUp' , float ) , ( 'TotalSigmaDn' , float ) , ( 'Model' , float ) , ( 'ModelSigmaUp' , float ) , ( 'ModelSigmaDn' , float ) , ( 'Background' , float ) , ( 'BackgroundSigmaUp' , float ) , ( 'BackgroundSigmaDn' , float ) , ( 'TotalOff' , float ) , ( 'TotalOffSigmaUp' , float ) , ( 'TotalOffSigmaDn' , float ) ]
    for name in srcnames :
        typelist.append( (  name , float ) )
        typelist.append( (  name + 'SigmaUp' , float ) )
        typelist.append( (  name + 'SigmaDn' , float ) )

    modelcounts = np.zeros( enumbins , dtype=typelist )
    return modelcounts

def getMeanCountsFromList( MeanModelCounts , listaObsContainer , confidence=0.68 ) :

    for name in listaObsContainer[ 0 ].dtype.names :
        for eindex in range( args.enumbins ) :
            temparray = np.zeros( args.nsims )
            for rindex in range( args.nsims ) :
                temparray[ rindex ] = listaObsContainer[ rindex ][ name ][ eindex ]
            mean = np.mean( temparray )
            ndata = len( temparray )
            stdDev = scipy.stats.sem( temparray )
            stdData = stdDev *  scipy.stats.t.ppf( ( 1 + confidence ) / 2. , ndata - 1 )
            MeanModelCounts[ name ][ eindex ] = mean
            MeanModelCounts[ name + 'SigmaUp' ][ eindex ] = mean + stdData
            MeanModelCounts[ name + 'SigmaDn' ][ eindex ] = mean - stdData

def writeMeanCountsArray( MeanModelCounts , file ) :
    fmt = '#'
    for index in range( len( MeanModelCounts.dtype.names ) ) :
        name = MeanModelCounts.dtype.names[ index ]
        fmt += '{:s}'.format( name )
        if index != ( len( MeanModelCounts.dtype.names ) - 1 ) :
            fmt += '\t'
    
    fmt += '\n'

    with open( file , 'w' ) as f :
        f.write( fmt )
        for eindex in range( MeanModelCounts.size ) :
            dummystr = ''
            for index in range( len( MeanModelCounts.dtype.names ) ) :
                name      = MeanModelCounts.dtype.names[ index ]
                dummystr += '{:.5e}'.format( MeanModelCounts[ name ][ eindex ] )
                if index != ( len( MeanModelCounts.dtype.names ) - 1 ) :
                    dummystr += '\t'
                else :
                    dummystr += '\n'
            f.write( dummystr )

def writeParametersArray( parameters , file ) :

    fmt = '#'

    for nindex in range( len( parameters.dtype.names ) ) :
        name = parameters.dtype.names[ nindex ]
        for parindex in range( len( parameters[ name ].dtype.names ) ) :
            par = parameters[ name ].dtype.names[ parindex ]
            fmt += '{:s}{:s}'.format( name , par )
            if parindex != ( len( parameters[ name ].dtype.names ) - 1 ) :
                fmt += '\t'
        if nindex != ( len( parameters.dtype.names ) - 1 ) :
            fmt += '\t'

    fmt += '\n'

    with open( file , 'w' ) as f :
        f.write( fmt )
        for sim in range( parameters.size ) :
            dummystr = ''
            for nindex in range( len( parameters.dtype.names ) ) :
                name = parameters.dtype.names[ nindex ]
                for parindex in range( len( parameters[ name ].dtype.names ) ) :
                    par = parameters[ name ].dtype.names[ parindex ]
                    dummystr += '{:.5e}'.format( parameters[ name ][ par ][ sim ] )
                    if parindex != ( len( parameters[ name ].dtype.names ) - 1 ) :
                        dummystr += '\t'
                if nindex != ( len( parameters.dtype.names ) - 1 ) :
                    dummystr += '\t'
                else :
                    dummystr += '\n'
            f.write( dummystr )

def writeEnergyBins( energybins , file ) :
    fmt = '#'
    for index in range( len( energybins.dtype.names ) ) :
        name = energybins.dtype.names[ index ]
        fmt += '{:s}'.format( name )
        if index != ( len( energybins.dtype.names ) - 1 ) :
            fmt += '\t'
    
    fmt += '\n'

    with open( file , 'w' ) as f :
        f.write( fmt )
        for eindex in range( energybins.size ) :
            dummystr = ''
            for index in range( len( energybins.dtype.names ) ) :
                name      = energybins.dtype.names[ index ]
                dummystr += '{:.5e}'.format( energybins[ name ][ eindex ] )
                if index != ( len( energybins.dtype.names ) - 1 ) :
                    dummystr += '\t'
                else :
                    dummystr += '\n'
            f.write( dummystr )

def plotParametersDist( parameters , models ) :
    for name in parameters.dtype.names :
        for par in parameters[ name ].dtype.names :
            unit = models[ name ].spectral()[ par ].unit()
            myCTAfuncs.plot_onedDist( parameters[ name ][ par ] , 10 , 0 ,\
                ( 0.44 , 0.16 , 0.39 ) , \
                '{:s}, {:s} ({:s})'.format( name , par , unit ) , \
                parameters.size , \
                auxMan.createname( args.outpath , '{:s}{:s}Distribution.png'.format( name , par ) ) )

def select_onoff_obs( obs , emin , emax ) :
    """
    Select an energy interval from one CTA On/Off observation

    Parameters
    ----------
    obs : `~gammalib.GCTAOnOffObservation`
        Minimum energy
    emin : `~gammalib.GEnergy()`
        Minimum energy
    emax : `~gammalib.GEnergy()`
        Maximum energy

    Returns
    -------
    obs : `~gammalib.GCTAOnOffObservation`
        CTA On/Off observation
    """
    # Select energy bins in etrue and ereco. All etrue energy bins are
    # selected. A 0.1% margin is added for reconstructed energies to
    # accomodate for rounding errors.

    etrue     = obs.rmf().etrue()
    ereco     = gammalib.GEbounds()
    itrue     = [ i for i in range( obs.rmf().etrue().size() ) ]
    ireco     = []
    for i in range( obs.rmf().emeasured().size() ) :
        ereco_bin_min = obs.rmf().emeasured().emin( i )
        ereco_bin_max = obs.rmf().emeasured().emax( i )
        if ereco_bin_min * 1.001 >= emin and ereco_bin_max * 0.999 <= emax:
            ereco.append(ereco_bin_min, ereco_bin_max)
            ireco.append(i)

    # Extract PHA
    pha_on  = gammalib.GPha( ereco )
    pha_off = gammalib.GPha( ereco )
    pha_on.exposure( obs.on_spec().exposure() )
    pha_off.exposure( obs.off_spec().exposure() )
    for idst, isrc in enumerate( ireco ) :
        # On
        pha_on[ idst ] = obs.on_spec()[ isrc ]
        pha_on.areascal( idst , obs.on_spec().areascal( isrc ) )
        pha_on.backscal( idst , obs.on_spec().backscal( isrc ) )
        # Off
        pha_off[ idst ] = obs.off_spec()[ isrc ]
        pha_off.areascal( idst , obs.off_spec().areascal( isrc ) )
        pha_off.backscal( idst , obs.off_spec().backscal( isrc ) )

    # Extract BACKRESP
    pha_backresp = obs.off_spec()[ 'BACKRESP' ]
    backresp     = []
    for idst, isrc in enumerate( ireco ) :
        backresp.append( pha_backresp[ isrc ] )
    pha_off.append( 'BACKRESP' , backresp )

    # Extract ARF
    arf = gammalib.GArf( etrue )
    for idst, isrc in enumerate( itrue ) :
        arf[idst] = obs.arf()[ isrc ]

    # Extract RMF
    rmf = gammalib.GRmf( etrue , ereco)
    for idst_true , isrc_true in enumerate( itrue ) :
        for idst_reco , isrc_reco in enumerate( ireco ) :
            rmf[ idst_true , idst_reco ] = obs.rmf()[ isrc_true, isrc_reco ]

    # Set On/Off observations
    obsid      = obs.id()
    statistic  = obs.statistic()
    instrument = obs.instrument()
    obs = gammalib.GCTAOnOffObservation( pha_on , pha_off , arf , rmf )
    obs.id( obsid )
    obs.statistic( statistic )
    obs.instrument(instrument)

    # Return observation
    return obs

def ONOFFAnalysisObservationContainer( meanobservation , caldb , irf , onoffmodelname , outpath , \
    nameSource , emin , emax , enumbins , algorithm  , ROIradius , models , gname , is_onoff , \
    sensData , srcSens , additional , addLabel , prefix ) :


    like = ctools.ctlike( meanobservation )
    like[ 'caldb' ]         = caldb
    like[ 'irf' ]           = irf
    like[ 'inmodel' ]       = onoffmodelname
    likeModelName           = auxMan.createname( outpath , '{:s}LikeOutModel.xml'.format( prefix ) )
    like[ 'outmodel' ]      = likeModelName
    like[ 'like_accuracy' ] = 1.e-4
    like[ 'max_iter' ]      = 100
    like[ 'nthreads' ]      = 2
    like.execute()

    sspec = mycsspec.csspec( like.obs() )

    specname = auxMan.createname( outpath , '{:s}Spectrum.fits'.format( prefix ) )

    sspec[ 'inmodel' ]  = onoffmodelname
    sspec[ 'srcname' ]  = nameSource
    sspec[ 'caldb' ]    = caldb
    sspec[ 'irf' ]      = irf
    sspec[ 'outfile' ]  = specname
    sspec[ 'method' ]   = 'SLICE'
    sspec[ 'ebinalg' ]  = 'LOG'
    sspec[ 'emin' ]     = emin
    sspec[ 'emax' ]     = emax
    sspec[ 'enumbins' ] = enumbins
    sspec[ 'nthreads' ] = 2

    sspec.execute()

    resname = auxMan.createname( outpath , '{:s}ResidualMap.fits'.format( prefix ) )

    respec                 = cscripts.csresspec( like.obs() )
    respec[ 'algorithm' ]  = algorithm
    respec[ 'mask' ]       = True
    respec[ 'ra' ]         = models[ nameSource ].spatial()[ 'RA' ].value()
    respec[ 'dec' ]        = models[ nameSource ].spatial()[ 'DEC' ].value()
    respec[ 'rad' ]        = ROIradius
    respec[ 'components' ] = True
    respec[ 'outfile' ]    = resname

    respec.execute()

    compname = auxMan.createname( outpath , '{:s}ComponentCounts.png'.format( prefix ) )
    myCTAfuncs.plot_residuals( resname , compname , 0 )

    thisenergies , dataon , dataoff , comps , label = myCTAfuncs.get_residuals( resname , 0 )

    myCTAfuncs.plot_spectrum( specname , auxMan.createname( outpath , '{:s}Spectrum{:s}.png'.format( prefix , gname ) ) ,\
        models , nameSource , dataon , dataoff , comps , is_onoff ,\
        1.e-9 , 1.e-15 , emax , emin , sensData , srcSens , additional , addLabel )

    os.remove( specname )
    os.remove( resname )
    for i in range( enumbins ) :
        rem = specname[ : -5 ] + 'ScanDataForBin{:d}.txt'.format( i + 1 )
        os.remove( rem )


if __name__ == '__main__':

    options = argparse.ArgumentParser( description='This script compute \
        the DM limits using Cirelli et al. Spectrum' )
    source = options.add_argument_group( 'Source' , \
        'All the relevant information about the source is passed in the inmodel\
        argument, where a XML file is passed describing parameters for the\
        observation, as position of the source, name, and spatial and spectral\
        models for the source. It is assumed that the file provides information about\
        the background model used in the observation. For more details, the user\
        is referred to the gammalib and ctools documentations :)' )
    source.add_argument( '--inmodel' ,\
        help='File with xml model' ,\
        type=str ,\
        required=True ,
        metavar='PerseusModel.xml' )
    source.add_argument( '--gname' ,\
        help='Name of the simulation' ,\
        type=str ,\
        required=True ,
        metavar='Perseus' )
    source.add_argument( '--coordsys' ,\
        help='Coordinate system used to perform the simulation. \
            Options are: [ GAL , CEL ] ( Galactic , Celestial )' ,\
        type=str ,\
        default='CEL' ,\
        metavar='CEL' ,\
        choices=[ 'GAL' , 'CEL' ] )
    source.add_argument( '--ROIra' ,\
        help='RA for center of Region of Interest (according to Coord. Sys.)' ,\
        type=float ,\
        required=True ,\
        metavar='47.0' )
    source.add_argument( '--ROIdec' ,\
        help='DEC forcenter of Region of Interestn (according to Coord. Sys.)' ,\
        type=float ,\
        required=True ,\
        metavar='39.0' )
    source.add_argument( '--ROIradius' ,\
        help='Radius of the Region of Interest (degrees)' ,\
        type=float ,\
        required=True ,\
        metavar='10.0' )
    
    instrument = options.add_argument_group( 'CTA-Intrument options' , \
        'Information about CTA IRF, observation time, etc.' )
    instrument.add_argument( '--pntra' ,\
        help='RA for pointing direction (according to Coord. Sys.)' ,\
        type=float ,\
        required=True ,\
        metavar='47.0' )
    instrument.add_argument( '--pntdec' ,\
        help='DEC for pointing direction (according to Coord. Sys.)' ,\
        type=float ,\
        required=True ,\
        metavar='39.0' )
    instrument.add_argument( '--pntradius' ,\
        help='Radius of the observation (degrees)' ,\
        type=float ,\
        required=True ,\
        metavar='10.0' )
    instrument.add_argument( '--caldb' ,\
        help='Database production for the IRF file. \
            Options are: [ prod2 , prod3b-v1 , prod3b-v2 ] ' ,\
        type=str ,\
        default='prod3b-v1' ,\
        metavar='prod3b-v1' ,\
        choices=[ 'prod2' , 'prod3b-v1' , 'prod3b-v2' ] )
    instrument.add_argument( '--irf' ,\
        help='Instrument Response Function for CTA' ,\
        type=str ,\
        default='North_z20_average_5h' ,\
        metavar='North_z20_average_50h' )
    instrument.add_argument( '--hours' ,\
        help='Time for simulation of observation (in hours)' ,\
        type=float ,\
        default=5.0 ,\
        metavar='50.0' )
    instrument.add_argument( '--enumbins',\
        help='Number of energy bins' ,\
        type=int ,\
        default=10 ,\
        metavar='20' )
    instrument.add_argument( '--emin',\
        help='Minimum energy for events (TeV)' ,\
        type=float ,\
        default=0.05 ,\
        metavar='0.01' )
    instrument.add_argument( '--emax',\
        help='Maximum energy for events (TeV)' ,\
        type=float ,\
        default=10.0 ,\
        metavar='100.0' )
    instrument.add_argument( '--id' ,\
        help='A number identifier for events in simulation. ' ,\
        type=str ,\
        default='0' ,\
        metavar='0' )
    instrument.add_argument( '--nsims' ,\
        help='Number of simulations' ,\
        type=int ,\
        required=True ,\
        metavar='10' )
    instrument.add_argument( '--outpath' , \
        help='Path to save files' ,\
        type=str ,\
        required=False ,\
        default='./' ,\
        metavar='path/to/file/' )
    instrument.add_argument( '--is_onoff' ,\
        help='Boolean to indicate ON/OFF observation type' ,\
        action='store_true' )
    instrument.add_argument( '--algorithm' ,\
        help='Method used to plot error for residual plot' ,\
        type=str ,\
        default='SIGNIFICANCE' ,\
        metavar='SIGNIFICANCE' )
    instrument.add_argument( '--with_sensitivity' ,\
        help='Boolean to indicate if spectrum plots need to show sensitivity curves' ,\
        action='store_true' )
    instrument.add_argument( '--sensFile' ,\
        help='File with data for sensitivity' ,\
        metavar='CrabSens.txt' )
    instrument.add_argument( '--srcSens' ,\
        help='Name of Source used to compute the sensitivity' ,\
        type=str ,\
        default='Crab' ,\
        metavar='Crab' )
    instrument.add_argument( '--addData' ,\
        help='Additional data to plot in spectrum plot' ,\
        type=str ,\
        metavar='data.txt' )


    args = options.parse_args()

    MeVtoTeV = 1.e-6 ;

    if args.with_sensitivity :
        if len( args.sensFile ) > 0 :
            sensData = np.genfromtxt( args.sensFile , names=True , delimiter=',' )
        else :
            print( 'Sorry, the file with data for sensitivity is not provided.' )
            print( 'No sensitivity curve wil be plotted' )
    else :
        if len( args.sensFile ) > 0 :
            print( 'Sorry, it is not required to plot sensitivity curve' )
            print( 'No sensitivity curve wil be plotted' )
        else :
            print( 'No sensitivity curve wil be plotted' )

    if len( args.addData ) > 0 :
        try :
            additionalData       = np.genfromtxt( args.addData , names=True , delimiter='  ' )
        except ValueError :
            print( 'Wrong delimiter' )
            additionalData       = np.zeros( 0 )
    else :
        print( 'No additional Data will be plotted' )
        additionalData           = np.zeros( 0 )

    bkgclasses = [ 'GCTAModelAeffBackground' , 'GCTAModelBackground' , 'GCTAModelCubeBackground' , 'GCTAModelIrfBackground' , 'GCTAModelRadialAcceptance' ]

    auxMan.checkDir( args.outpath )

    models = gammalib.GModels( args.inmodel )
    dir = gammalib.GSkyDir()
    dir.radec_deg( args.pntra , args.pntdec )
    duration = args.hours * 3600.

    Etypes     = np.dtype( [ ( 'Emean' , float ) , ( 'Emin' , float ) , ( 'Emax' , float ) ] )
    energybins = np.zeros( args.enumbins , dtype=Etypes )

    nameSource = [ source.name() for source in models if source.classname() not in bkgclasses ]
    partypes   = [ ( source.name() , [ ( par.name() , float ) for par in source.spectral() ] ) for source in models if source.classname() not in bkgclasses ]
    stypes     = np.dtype( partypes )
    parameters = np.zeros( args.nsims , dtype=stypes )


    MeanModelCounts     = getEmptyCountsArray( args.enumbins , nameSource )
    MeanModelCountsLike = getEmptyCountsLikeArray( args.enumbins , nameSource )
    listaObsContainer   = []
    listaLikeContainer  = []


    flux   = []
    ulflux = []
    eflux  = np.empty( shape=( 0 , 0 ) )
    euflux = np.empty( shape=( 0 , 0 ) )


    #################################################################
    #####                       The Simuations
    #################################################################
    for sim in range( args.nsims ) :

        print( '**********************************' )
        print( 'Observation: {:d}'.format( sim + 1 ) )
        print( '**********************************' )

        cntname = auxMan.createname( args.outpath , \
            'CountCubeObsSim{:d}.fits'.format( sim ) )


        obslist = gammalib.GObservations()
        type = ''
        eventlist=''

        obs_binned = myCTAfuncs.single_obs( type=type , \
            eventfile=eventlist ,\
            pntdir=dir ,\
            tstart=0.0 ,\
            duration=duration ,\
            deadc=0.95 ,\
            emin=args.emin ,\
            emax=args.emax ,\
            rad=args.pntradius ,\
            irf=args.irf ,\
            caldb=args.caldb ,\
            id=args.id, \
            instrument="CTA" )
        obslist.append( obs_binned )


        obssim                = ctools.ctobssim( obslist )
        obssim[ 'inmodel' ]   = args.inmodel
        obssim[ 'outevents' ] = cntname
        obssim[ 'nthreads' ]  = 2
        obssim[ 'seed' ]      = int( time.time() )
        obssim[ 'edisp' ]     = False

        obssim.execute()

        obsname = auxMan.createname( args.outpath , \
            'Obs{:d}.xml'.format( sim ) )

        obsxml  = gammalib.GXml()
        obslist = obsxml.append( 'observation_list title="observation library"' )
        obs     = obslist.append( 'observation name="%s" id="%s" \
            instrument="CTA"' % ( args.gname , args.id + str( sim ) ) )
        obs.append( 'parameter name="EventList" file="%s"' % cntname )
        obsxml.save( obsname )

        if sim == 0 :
            for source in models :
                if source.classname() not in bkgclasses :
                    mapname = auxMan.createname( args.outpath , \
                        'Source{:s}CTAMAP{:d}.png'.format( source.name() , sim ) )
                    ctamap( cntname , args.emin , args.emax , args.coordsys , \
                        source.spatial()[ 'RA' ].value() , \
                        source.spatial()[ 'DEC' ].value() , \
                        2.0 , args.caldb , args.irf , \
                        source.name() , mapname )
            mapname = auxMan.createname( args.outpath , \
                'Global{:s}CTAMAP{:d}.png'.format( args.gname , sim ) ) 
            ctamap( cntname , args.emin , args.emax , args.coordsys ,\
                args.ROIra , args.ROIdec , args.ROIradius , \
                args.caldb , args.irf , args.gname , mapname )

        onoffgen = cscripts.csphagen()

        onoffname = auxMan.createname( args.outpath , \
            'ONOFFObs{:d}.xml'.format( sim ) )
        onoffmodelname = auxMan.createname( args.outpath , \
            'ONOFFModel{:d}.xml'.format( sim ) )
        prefix = os.path.join( args.outpath , \
            'onoff{:d}'.format( sim ) )

        onoffgen[ 'inobs' ]      = obsname
        onoffgen[ 'inmodel' ]    = args.inmodel
        onoffgen[ 'srcname' ]    = nameSource[ 0 ]
        onoffgen[ 'caldb' ]      = args.caldb
        onoffgen[ 'irf' ]        = args.irf
        onoffgen[ 'outobs' ]     = onoffname
        onoffgen[ 'outmodel' ]   = onoffmodelname
        onoffgen[ 'ebinalg' ]    = 'LOG'
        onoffgen[ 'emin' ]       = args.emin
        onoffgen[ 'emax' ]       = args.emax
        onoffgen[ 'enumbins' ]   = args.enumbins
        onoffgen[ 'coordsys' ]   = args.coordsys
        onoffgen[ 'ra' ]         = models[ nameSource[ 0 ] ].spatial()[ 'RA' ].value()
        onoffgen[ 'dec' ]        = models[ nameSource[ 0 ] ].spatial()[ 'DEC' ].value()
        onoffgen[ 'rad' ]        = args.ROIradius
        onoffgen[ 'bkgmethod' ]  = 'REFLECTED'
        onoffgen[ 'bkgregskip' ] = 0
        onoffgen[ 'etruemin' ]   = args.emin
        onoffgen[ 'etruemax' ]   = args.emax
        onoffgen[ 'nthreads' ]   = 2
        onoffgen[ 'stack' ]      = False
        onoffgen[ 'prefix' ]     = prefix

        onoffgen.execute()

        onoffobs = onoffgen.obs()[ 0 ]
        ebounds  = onoffobs.on_spec().ebounds()

        modelcounts = getONOFFCounts( onoffobs , models , ebounds , 1 )

        listaObsContainer.append( modelcounts )

        like = ctools.ctlike()

        like[ 'inobs' ]    = onoffname
        like[ 'caldb' ]    = args.caldb
        like[ 'irf' ]      = args.irf
        like[ 'inmodel' ]  = onoffmodelname
        likeModelName      = auxMan.createname( args.outpath , 'LikeOutModel{:d}.xml'.format( sim ) )
        like[ 'outmodel' ] = likeModelName
        like[ 'like_accuracy' ] = 1.e-4
        like[ 'max_iter' ] = 100
        like.execute()

        likeonoffobs      = like.obs()[ 0 ]
        thismodel         = like.obs().models()
        modelcountsLike   = getONOFFCountsAfterLike( likeonoffobs , thismodel , ebounds , 1 )
        listaLikeContainer.append(  modelcountsLike )

        for source in thismodel :
            if source.classname() not in bkgclasses :
                for par in source.spectral() :
                    parameters[ source.name() ][ par.name() ][ sim ] = thismodel[ source.name() ][ par.name() ].value()

        if sim == 0 :

            meanObs    = select_onoff_obs( onoffobs , ebounds.emin( 0 ) , ebounds.emax( ebounds.size() - 1 ) )
            # SigmaUpObs = select_onoff_obs( onoffobs , ebounds.emin( 0 ) , ebounds.emax( ebounds.size() - 1 ) )
            # SigmaDnObs = select_onoff_obs( onoffobs , ebounds.emin( 0 ) , ebounds.emax( ebounds.size() - 1 ) )
        else :
            os.remove( cntname )
            os.remove( obsname )
            os.remove( onoffname )
            os.remove( onoffmodelname )
            os.remove( likeModelName )
            os.remove( prefix + '_arf.fits' )
            os.remove( prefix + '_off.reg' )
            os.remove( prefix + '_on.reg' )
            os.remove( prefix + '_pha_off.fits' )
            os.remove( prefix + '_pha_on.fits' )
            os.remove( prefix + '_rmf.fits' )

        #os.system( 'say -v Kyoko Simulation {:d} Done'.format( sim + 1 ) )

        auxMan.delete_last_lines( 3 )


    ebounds = meanObs.on_spec().ebounds()
    for eindex in range( ebounds.size() ) :
        meanObs.on_spec().counts_spectrum()[ eindex ]     = MeanModelCounts[ 'Total' ][ eindex ]
        meanObs.off_spec().counts_spectrum()[ eindex ]    = MeanModelCounts[ 'TotalOff' ][ eindex ]
        # SigmaUpObs.on_spec().counts_spectrum()[ eindex ]  = MeanModelCounts[ 'TotalSigmaUp' ][ eindex ]
        # SigmaUpObs.off_spec().counts_spectrum()[ eindex ] = MeanModelCounts[ 'TotalOffSigmaUp' ][ eindex ]
        # SigmaDnObs.on_spec().counts_spectrum()[ eindex ]  = MeanModelCounts[ 'TotalSigmaDn' ][ eindex ]
        # SigmaDnObs.off_spec().counts_spectrum()[ eindex ] = MeanModelCounts[ 'TotalOffSigmaDn' ][ eindex ]

    meanobservation    = gammalib.GObservations()
    # sigmaupobservation = gammalib.GObservations()
    # sigmadnobservation = gammalib.GObservations()

    meanobservation.append( meanObs )
    # sigmaupobservation.append( SigmaUpObs )
    # sigmadnobservation.append( SigmaDnObs )

    for i in range( ebounds.size() ) :
        energybins[ 'Emean' ][ i ] = ebounds.emean( i ).TeV()
        energybins[ 'Emin' ][ i ]  = ebounds.emin( i ).TeV()
        energybins[ 'Emax' ][ i ]  = ebounds.emax( i ).TeV()

    onoffmodelname = auxMan.createname( args.outpath , \
        'ONOFFModel0.xml' )

    addLabel = 'CTA Gamma Prop. Paper'
    ONOFFAnalysisObservationContainer( meanobservation , args.caldb , args.irf , onoffmodelname , \
        args.outpath , nameSource[ 0 ] , args.emin , args.emax , args.enumbins , \
        args.algorithm , args.ROIradius , models , args.gname , args.is_onoff , \
        sensData , args.srcSens ,\
        additionalData , addLabel , 'meanObs' )

    # ONOFFAnalysisObservationContainer( sigmaupobservation , args.caldb , args.irf , onoffmodelname , \
    #   args.outpath , nameSource[ 0 ] , args.emin , args.emax , args.enumbins , \
    #   args.algorithm , args.ROIradius , models , args.gname , args.is_onoff , \
    #   sensData , args.srcSens , 'SigmaUpObs' )

    # ONOFFAnalysisObservationContainer( sigmadnobservation , args.caldb , args.irf , onoffmodelname , \
    #   args.outpath , nameSource[ 0 ] , args.emin , args.emax , args.enumbins , \
    #   args.algorithm , args.ROIradius , models , args.gname , args.is_onoff , \
    #   sensData , args.srcSens , 'SigmaDnObs' )

    writeEnergyBins( energybins ,\
        auxMan.createname( args.outpath , 'AnalysisEnergyBins.txt' ) )


    print( '**************************************************' )

    getMeanCountsFromList( MeanModelCounts , listaObsContainer )
    getMeanCountsFromList( MeanModelCountsLike , listaLikeContainer )

    writeMeanCountsArray( MeanModelCounts , \
        auxMan.createname( args.outpath , '{:s}MeanCountsObsContainer.txt'.format( args.gname ) ) )

    writeMeanCountsArray( MeanModelCountsLike , \
        auxMan.createname( args.outpath , '{:s}MeanCountsObsLikeContainer.txt'.format( args.gname ) ) )

    writeParametersArray( parameters ,\
        auxMan.createname( args.outpath , '{:s}ParametersModels.txt'.format( args.gname ) ) )

    plotParametersDist( parameters , models )

    #os.system( 'say -v Kyoko This is the END' )


















































