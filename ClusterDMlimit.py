#! /usr/bin/env python

# ====================================================================
# This script perform the DM limits using the Cirelli et al. spectrum.
# ====================================================================

import ctools
import gammalib
import cscripts
import math

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

if __name__ == '__main__':

    doSim = False
    options = argparse.ArgumentParser( description='This script compute \
                            the DM limits using Cirelli et al. Spectrum' )
    source = options.add_argument_group( 'Source' , \
        'Information about the source where the limits are computed' )
    source.add_argument( '--name' ,\
        help='Name of the source' ,\
        type=str ,\
        required=True ,
        metavar='Perseus' )
    source.add_argument( '--ra' ,\
        help='Position of the source: RA in degrees' ,\
        type=float ,\
        required=True ,\
        metavar='47.0' )
    source.add_argument( '--dec' ,\
        help='Position of the source: Dec in degrees' ,\
        type=float ,\
        required=True ,\
        metavar='39.0' )
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
        help='Value for annihilation cross-section (in cm**3/s) \
              or decay lifetime (in s)' ,\
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
        help='DM channel in annihilation or decay process. \
            Possible options are: - 5 : mumu - 8 : tautau- 11 : bb- 15 : ww ' ,\
        type=int ,\
        required=True ,\
        metavar='11' ,\
        choices= [ 5 , 8 , 11 , 15 ] )

    instrument = options.add_argument_group( 'CTA-Intrument related options' , \
    'Information about CTA IRF, observation time, etc.' )
    instrument.add_argument( '--prod' ,\
        help='Database production for the IRF file. \
        Options are: [ prod2 , prod3b-v1 , prod3b-v2 ] ' ,\
        type=str ,\
        default='prod3b-v2' ,\
        metavar='prod3b-v1' ,\
        choices=[ 'prod2' , 'prod3b-v1' , 'prod3b-v2' ] )
    instrument.add_argument( '--irf' ,\
        help='Instrument Response Function for CTA' ,\
        type=str ,\
        default='North_z20_average_5h' ,\
        metavar='North_z20_average_5h' )
    instrument.add_argument( '--hours' ,\
        help='Time for simulation of observation (in hours)' ,\
        type=float ,\
        default=50.0 ,\
        metavar='50.0' )
    instrument.add_argument( '--emin',\
        help='Minimum energy for events (in TeV)' ,\
        type=float ,\
        default=0.01 ,\
        metavar='0.1' )
    instrument.add_argument( '--emax',\
        help='Maximum energy for events (in TeV)' ,\
        type=float ,\
        default=100.0 ,\
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

    auxMan.checkPathandExit( args.xmlmodel )
    auxMan.checkDir( args.outpath )

    #   Create model container
    models = GModels( args.xmlmodel )

    #   Setting the direction of the pointing
    #   This is needed to create the obs container
    dir = GSkyDir()
    dir.radec_deg( args.ra , args.dec )

    #   Setting the duration of the obs container
    duration = args.hours * 3600.0

    #   I think, this can be avoided if we need to create new
    #   sims every time we want to compute the upper-limits
    #   BAH
    type = ""
    eventlist= ""

    #   Now, this is to print useful information about the simulation
    print( '#======================================#\
        \n#          Gral. Information:!\
        \n#======================================#\
        \n# channel_number = [11,5,8]  \
        \n# Channels: bb, mumu, tautau\
        \n# channel_number = [11,5,12,8,15]  \
        \n# Channels: bb, mumu, tt, tautau, ww ' )

    channel_number_dic = { '5' : 'mumu' , '8':'tautau' , '11':'bb' , '15':'ww' }
    # Channels: tautau

    for sims in range( args.nsims ) :
        thisrunID = int( args.id ) + sims
        print ( '\n%%%%%%%%%%%%%%%%%' )
        print ( 'Mass: ' , args.mass , ' TeV' )
        print ( 'Cross Section: ' , args.sigmav , ' cm3/s' )
        print ( 'Channel: %s' % channel_number_dic[ str( args.channel ) ] )
        print ( 'Simulation: %s' % str( thisrunID ) )
        print ( '%%%%%%%%%%%%%%%%%\n' )

        print ( 'Copy model...' )
        time.sleep( 1 )

        #   Observation container
        obslist = gammalib.GObservations()

        #   Setting observation parameters
        obs_binned = myCTAfuncs.single_obs( type=type , \
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
        #   Append obs to the container
        obslist.append( obs_binned )

        #   I think, this is not needed
        #   Setting model to the obs container
        obslist.models( models )
        #   Get spectrum reference
        theoretical_spectrum = models[ args.name ].spectral()

        #   Minimum and maximum energy according to the DM process
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

        #   Create GEnergy objects
        e_ref   = gammalib.GEnergy( thiseref , "TeV" )
        e_min   = gammalib.GEnergy( thisemin , "TeV" )
        e_max   = gammalib.GEnergy( thisemax , "TeV" )
        print( 'Reference Energu: %.2f TeV\nMaximum Energy: \
                %.2f TeV\nMinimum Energy: %.2f TeV' \
                % ( thiseref , thisemax , thisemin ) )

        #   Computing theoretical flux
        th_flux_ref = theoretical_spectrum.eval( e_ref , gammalib.GTime() )
        print( 'Theoretical Spectrum at Eref %.4f TeV: %.4e' \
            % ( thiseref , th_flux_ref ) ) ;
        try:
            print( "Doing the simulation" )

            #   Cube name
            cube_name = args.name + "_" + args.process \
                + "_" + str( args.channel ) + "_" \
                + str( args.mass ) + "TeV_" + str( thisrunID ) \
                + "_" + str( args.hours ) + "h_" \
                + str( args.sigmav ) + "_cntcube.fits"
            cube_name = auxMan.createname( args.outpath , cube_name )

            #   Doing the simulation
            obssim = myCTAfuncs.sim( obslist ,\
                log=False ,\
                debug=False ,\
                chatter=4 ,\
                edisp=False ,\
                seed=int( time.time() ) ,\
                nbins=10 ,\
                binsz=0.02 ,\
                npix=100 ,\
                outfile=cube_name )

            print ( 'Number of Cube events: ' ,  obssim[ 0 ].events().size() , \
                'Total photons: ' , obssim[ 0 ].events().size() )

            #   Creating OnOff observation library XML

            obsname =  auxMan.createname( args.outpath , \
                'Obs{:d}.xml'.format( thisrunID ) )

            obsxml  = gammalib.GXml()
            obslist = obsxml.append( 'observation_list title="observation library"' )
            obs     = obslist.append( 'observation name="%s" id="%s" \
                instrument="CTA"' % ( args.name , str( thisrunID ) ) )
            obs.append( 'parameter name="EventList" file="%s"' % cube_name )
            obsxml.save( obsname )

            #   csphagen to compute all relevant parameters for OnOff obs.

            onoffgen       = cscripts.csphagen()

            onfoffname     = auxMan.createname( args.outpath , \
                'ONOFFObs{:d}.xml'.format( thisrunID ) )
            onoffmodelname = auxMan.createname( args.outpath , \
                'ONOFFModel{:d}.xml'.format( thisrunID ) )
            prefix         = os.path.join( args.outpath , \
                'onoff{:d}'.format( thisrunID ) )

            onoffgen[ 'inobs' ]      = obsname
            onoffgen[ 'inmodel' ]    = args.xmlmodel
            onoffgen[ 'srcname' ]    = args.name
            onoffgen[ 'caldb' ]      = args.prod
            onoffgen[ 'irf' ]        = args.irf
            onoffgen[ 'outobs' ]     = onoffname
            onoffgen[ 'outmodel' ]   = onoffmodelname
            onoffgen[ 'ebinalg' ]    = 'LOG'
            onoffgen[ 'emin' ]       = thisemin
            onoffgen[ 'emax' ]       = thisemax
            onoffgen[ 'enumbins' ]   = 10
            onoffgen[ 'coordsys' ]   = 'CEL'
            onoffgen[ 'ra' ]         = args.ra
            onoffgen[ 'dec' ]        = args.dec
            onoffgen[ 'rad' ]        = args.radius
            onoffgen[ 'bkgmethod' ]  = 'REFLECTED'
            onoffgen[ 'bkgregskip' ] = 0
            onoffgen[ 'etruemin' ]   = thisemin
            onoffgen[ 'etruemax' ]   = thisemax
            onoffgen[ 'nthreads' ]   = 8
            onoffgen[ 'stack' ]      = False
            onoffgen[ 'prefix' ]     = prefix

            onoffgen.execute()

            # Optimize forDM analysis
            print ( '\nDM fitting...' )
            like                    = ctools.ctlike()
            like[ 'inobs' ]         = onoffname
            like[ 'expcube' ]       = "NONE"
            like[ 'psfcube' ]       = "NONE"
            like[ 'bkgcube' ]       = "NONE"
            like[ 'caldb' ]         = args.prod
            like[ 'irf' ]           = args.irf
            like[ 'inmodel' ]       = onoffmodelname
            likeModelName           = auxMan.createname( args.outpath , \
                'LikeOutModel{:d}.xml'.format( sim ) )
            like[ 'outmodel' ]      = likeModelName
            like[ 'like_accuracy' ] = 1.e-4
            like[ 'max_iter' ]      = 100
            like[ 'nthreads' ]      = 8
            like.execute()

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
            limit                = ctools.ctulimit()
            limit[ 'inobs' ]     = onoffname
            limit[ 'srcname' ]   = args.name
            limit[ 'inmodel' ]   = onoffmodelname
            limit[ 'expcube' ]   = 'NONE'
            limit[ 'psfcube' ]   = 'NONE'
            limit[ 'bkgcube' ]   = 'NONE'
            limit[ 'eref' ]      = thiseref
            limit[ 'emin' ]      = thisemin
            limit[ 'emax' ]      = thisemax
            limit[ 'caldb' ]     = args.prod
            limit[ 'irf' ]       = args.irf
            limit[ 'statistic' ] = 'DEFAULT'
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
                print( "Sorry, unknown process. \
                    I will assume that you refer to annihilation process" )
                print( "But, All is not missing, \
                    you can get the scale factor and \
                    compute your parameter of interest" )
                cs = sigma_v * scale
            print( '\t\tFlux: %.5e\n\t\tScale: %.5e\n\t\tPOI: %.5e\n' \
                % ( flux , scale , cs ) )
            file_name = args.name + '_' + args.process + '_' \
                + str( args.channel ) + '_' \
                + str( args.mass ) + 'TeV_' + str( thisrunID ) \
                + '_' + str( args.hours ) + 'h_' + str( args.sigmav ) + '.txt'
            file_name = os.path.join( args.outpath , file_name )
            with open( file_name , 'a' ) as outfile:
                outfile.write('%d\t%d\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\n' \
                    % ( thisrunID , args.channel ,\
                        obslist[0].events().number() ,\
                        cs , scale , \
                        ts , args.sigmav ) )
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

