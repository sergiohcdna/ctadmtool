#! /usr/bin/env python

# ====================================================================
# This script perform the DM limits using the Cirelli et al. spectrum.
# ====================================================================

import ctools 
from gammalib import *
import math

import glob
import os
import sys
from sys import argv
import os.path
import csv
import time
import argparse
from myCTAfuncs import *
    
#=====================#
# Routine entry point #
#=====================#

def checkPathandExit( thispath ) :
    if not os.path.exists( thispath ) :
        print( '\t\tError. You are trying to pass a file or directory that does not exists.' )
        sys.exit( 1 )

def checkDir( thispath ) :
    if os.path.exists( thispath ) :
        if not os.path.isdir( thispath ) :
            print( 'Sorry, you are trying to pass any other file-type as the outpath directory' )
            sys.exit( 1 )
        else :
            print( 'Good, this a directory. Nothing to do' )
    else :
        print( 'It seems like the path does not exists. The scripts create the directory for you' )
        os.makedirs( thispath )

if __name__ == '__main__':

    doSim = False
    options = argparse.ArgumentParser( description='This script compute the DM limits using Cirelli et al. Spectrum' )
    source = options.add_argument_group( 'Source' , 'Information about the source where the limits are computed' )
    source.add_argument( '--name' ,\
        help='Name of the source' ,\
        type=str ,\
        required=True ,
        metavar='Perseus' )
    source.add_argument( '--longitude' ,\
        help='Position of the source: Longitud in degrees' ,\
        type=float ,\
        required=True ,\
        metavar='145.5' )
    source.add_argument( '--latitude' ,\
        help='Position of the source: Latitude in degrees' ,\
        type=float ,\
        required=True ,\
        metavar='-23.45' ) 
    source.add_argument( '--radius' ,\
        help='Angular extension of the source (in degrees)' ,\
        type=float ,\
        required=True ,\
        metavar='2.0')
    
    dmmodel = options.add_argument_group( 'Dark Matter (DM) Model' ,\
        'Information about the spectra and spatial model of DM' )
    dmmodel.add_argument( '--mass' ,\
        help='Mass of the DM particle (in TeV)' ,\
        type=float ,\
        required=True ,\
        metavar=100.0 )
    dmmodel.add_argument( '--sigmav' ,\
        help='Value for annihilation cross-section (in cm**3/s) or decay lifetime (in s)' ,\
        type=float ,\
        required=True ,\
        metavar='1.e-28' )
    dmmodel.add_argument( '--process' ,\
        help='DM annihilation or DM decay' ,\
        type=str ,\
        required=True ,\
        metavar='annihilation' ,\
        choices=[ 'annihilation' , 'decay' ] )
    dmmodel.add_argument( '--channel' ,\
        help='DM channel in annihilation or decay process. Possible options are: - 5 : mumu - 8 : tautau- 11 : bb- 15 : ww ' ,\
        type=int ,\
        required=True ,\
        metavar='11' ,\
        choices= [ 5 , 8 , 11 , 15 ] )

    instrument = options.add_argument_group( 'CTA-Intrument related options' , 'Information about CTA IRF, observation time, etc.' )
    instrument.add_argument( '--prod' ,\
        help='Database production for the IRF file. Options are: [ prod2 , prod3b-v1 ] ' ,\
        type=str ,\
        required=True ,\
        metavar='prod3b-v1' ,\
        choices=[ 'prod2' , 'prod3b-v1' ] )
    instrument.add_argument( '--irf' ,\
        help='Instrument Response Function for CTA' ,\
        type=str ,\
        required=True ,\
        metavar='South_50h' )
    instrument.add_argument( '--hours' ,\
        help='Time for simulation of observation (in hours)' ,\
        type=float ,\
        required=True ,\
        metavar='50.0' )
    instrument.add_argument( '--emin',\
        help='Minimum energy for events (in TeV)' ,\
        type=float ,\
        required=True ,\
        metavar='0.1' )
    instrument.add_argument( '--emax',\
        help='Maximum energy for events (in TeV)' ,\
        type=float ,\
        required=True ,\
        metavar='100.0' )
    instrument.add_argument( '--id' ,\
        help='A number identifier for events in simulation. ' ,\
        type=str ,\
        required=True ,\
        metavar='0' )
    instrument.add_argument( '--nsims' ,\
        help='Number of simulations' ,\
        type=int ,\
        required=True ,\
        metavar='50' )
    instrument.add_argument( '--xmlmodel' ,\
        help='File with model for CTA observation' ,\
        type=str ,\
        required=True ,\
        metavar='xmlmodel.xml' )
    instrument.add_argument( '--outpath' , \
        help='Path to save files' ,\
        type=str ,\
        required=False ,\
        default='./' ,\
        metavar='path/to/file/' )


    args = options.parse_args()

    checkPathandExit( args.xmlmodel )
    checkDir( args.outpath )
        
   # Create observation container
    models = GModels( args.xmlmodel )
    dir = GSkyDir()
    dir.lb_deg( args.longitude , args.latitude )

    duration = args.hours * 3600.0
    
    type = ""
    eventlist= ""
    
    print( '#======================================#\
        \n#          ENTER THE LOOP!!            #\
        \n#======================================#\
        \n# channel_number = [11,5,8]  \
        \n# Channels: bb, mumu, tautau\
        \n# channel_number = [11,5,12,8,15]  \
        \n# Channels: bb, mumu, tt, tautau, ww ' )

    channel_number_dic = { '5' : 'mumu' , '8':'tautau' , '11':'bb' , '15':'ww' }       # Channels: tautau

    for sims in range( args.nsims ) :
        thisrunID = int( args.id ) + sims
        print ( '\n%%%%%%%%%%%%%%%%%' )
        print ( 'Mass: ' , args.mass , ' TeV' )
        print ( 'Cross Section: ' , args.sigmav , ' cm3/s' )
        print ( 'Channel: %s' % channel_number_dic[ str( args.channel ) ] ) 
        print ( 'Simulation: %s' % str( thisrunID ) )
        print ( '%%%%%%%%%%%%%%%%%\n' )
    
        print ( 'Copy model...' )
        time.sleep(1)
        
        obslist = GObservations()
        
        obs_binned = single_obs( type=type , \
            eventfile=eventlist ,\
            pntdir=dir ,\
            tstart=0.0 ,\
            duration=duration ,\
            deadc=0.95 ,\
            emin=args.emin ,\
            emax=args.emax ,\
            rad=args.radius ,\
            irf=args.irf ,\
            caldb=args.prod ,\
            id=args.id, \
            instrument="CTA" )
        obslist.append( obs_binned )
        
        obslist.models( models )
        theoretical_spectrum = models[ args.name ].spectral()
        
        thiseref = 0.0 ; thisemax = 0.0 ;
        if args.process == 'decay' :
            thiseref = args.mass * 0.9 / 2.
            thisemax = args.mass / 2.
        elif args.process == 'annihilation' :
            thiseref = args.mass * 0.9 
            thisemax = args.mass 
        else :
            print( 'Unknown process. I will assume DM annihilation' )
            thiseref = args.mass * 0.9
            thisemax = args.mass
        thisemin = 0.01
        e_ref   = GEnergy( thiseref , "TeV" )
        e_min   = GEnergy( thisemin , "TeV" )
        e_max   = GEnergy( thisemax , "TeV" )
        print( 'Reference Energu: %.2f TeV\nMaximum Energy: %.2f TeV\nMinimum Energy: %.2f TeV' % ( thiseref , thisemax , thisemin ) )
        
        th_flux_ref = theoretical_spectrum.eval( e_ref , GTime() )
        print( 'Theoretical Spectrum at Eref %.4f TeV: %.4e' % ( thiseref , th_flux_ref ) ) ;
        try:
            print( "Doing the simulation" )
            cube_name = args.name + "_" + args.process + "_" + str( args.channel ) + "_" + str( args.mass ) + "TeV_" + str( thisrunID ) + "_" + str( args.hours ) + "h_" + str( args.sigmav ) + "_cntcube.fits"
            cube_name = os.path.join( args.outpath , cube_name )
            print( 'Cube saved with name %s' % cube_name )
            obssim = sim( obslist ,\
                log=False ,\
                debug=False ,\
                chatter=4 ,\
                edisp=False ,\
                seed=int( time.time() ) ,\
                outfile=cube_name )
            print ( 'Number of CUBE events: ' ,  obssim[ 0 ].events().size() , 'Total photons: ' , obssim[ 0 ].events().size() )
            # Optimize forDM analysis
            print ( '\nDM fitting...' )
            like = ctools.ctlike()
            like[ 'inobs' ] = cube_name
                
            like[ 'expcube' ] = "NONE"
            like[ 'psfcube' ] = "NONE"
            like[ 'bkgcube' ] = "NONE"
            like[ 'caldb' ] = args.prod
            like[ 'irf' ] = args.irf
            like[ 'inmodel' ] = args.xmlmodel
            
            like.run()

            models_fit = like.obs().models();
            ts = 0
            for m in models_fit :

                if( m.tscalc()) :
                    print ( "\n" , m.name() , "\t TS: " , m.ts() )
                    ts = m.ts()
                else :
                    print ( "\n", m.name() )

                for i in m :
                    print ( i )

            print ( '\nDM calc upper limits... ' )
            limit = ctools.ctulimit()
            limit[ 'inobs' ] = cube_name
            
            limit[ 'srcname' ].string( models[ args.name ].name() )
            limit[ 'inmodel' ] = args.xmlmodel
            limit[ 'expcube' ] = 'NONE'
            limit[ 'psfcube' ] = 'NONE'
            limit[ 'bkgcube' ] = 'NONE'
            limit[ 'eref' ].real( thiseref )
            limit[ 'emin' ].real( thisemin )
            limit[ 'emax' ].real( thisemax )
            limit[ 'caldb' ] = args.prod
            limit[ 'irf' ] = args.irf
            limit[ 'statistic' ] = 'CHI2'
            # limit[ 'debug' ] = True
            limit.run()

            print ( '\nRESULTS:' )
            cs = 0
            flux    = limit.diff_ulimit()
            scale   = flux / th_flux_ref
            if args.process == "decay" :
                cs = args.sigmav / scale
            elif args.process == "annihilation" :
                cs = args.sigmav * scale
            else :
                print( "Sorry, unknown process. I will assume that you refer to annihilation process" )
                print( "But, All is not missing, you can get the scale factor and compute your parameter of interest" )
                cs = sigma_v * scale
            print( '\t\tFlux: %.5e\n\t\tScale: %.5e\n\t\tPOI: %.5e\n' % ( flux , scale , cs ) )
            file_name = args.name + '_' + args.process + '_' + str( args.channel ) + '_' + str( args.mass ) + 'TeV_' + str( thisrunID ) + '_' + str( args.hours ) + 'h_' + str( args.sigmav ) + '.txt'
            file_name = os.path.join( args.outpath , file_name )
            with open( file_name , 'a' ) as outfile:
                outfile.write('%d\t%d\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\n' % ( thisrunID , args.channel , obslist[0].events().number() , cs , scale , ts , args.sigmav ) )
            # Free memory
            del limit
            del like
            del obssim
            del obs_binned
            del obslist
            del cube_name
        except RuntimeError:
            pass
    print( 'ok!' )

