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
    source.add_argument( '--longitude' ,\
        help='Position of the source: Dec in degrees' ,\
        type=float ,\
        required=True ,\
        metavar='-13.22' )
    source.add_argument( '--latitude' ,\
        help='Position of the source: RA in degrees' ,\
        type=float ,\
        required=True ,\
        metavar='150.57' )
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
        default='North_z20_5h' ,\
        metavar='North_z20_5h' )
    instrument.add_argument( '--coordsys' ,\
        help='Coordinate system selected' ,\
        type=str ,\
        default='GAL' ,\
        metavar='GAL' ,\
        choices=[ 'GAL' , 'CEL' ] )
    instrument.add_argument( '--pntlon' ,\
        help='Pointing direction (Longitude in degrees)' ,\
        type=float ,\
        required=True ,\
        metavar='-15.00' )
    instrument.add_argument( '--pntlat' ,\
        help='Ponting direction (Latitude in degrees)' ,\
        type=float ,\
        required=True ,\
        metavar='153.00' )
    instrument.add_argument( '--pntrad' ,\
        help='Radius of obs ROI' ,\
        type=float ,\
        required=True ,\
        metavar='4.5' )
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
    instrument.add_argument( '--fitsfile' ,\
        help='Name fo fits file' ,\
        type=str ,\
        required=True ,\
        metavar='myfits.fits' )


    args = options.parse_args()

    auxMan.checkPathandExit( args.xmlmodel )
    auxMan.checkDir( args.outpath )

    #   Create model container
    models = gammalib.GModels( args.xmlmodel )

    #   Setting the direction of the pointing
    #   This is needed to create the obs container
    dir = gammalib.GSkyDir()

    #   Check the coordsys argument and set
    #   the source position properly
    if args.coordsys == 'GAL' :
        dir.lb_deg( args.longitude , args.latitude )
    elif args.coordsys == 'CEL' :
        dir.radec_deg( args.longitude , args.latitude )
    else :
        print( 'Coordsystem unknown. Using default value' )
        dir.lb_deg( args.longitude , args.latitude )

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

    #   Now, creating a fits file to save the results
    #   All the information related to DM process,
    #   mass, and channel used to compute the spectra;
    #   is added to the header
    #   The number of parameters to save is 5:
    #       - run ID for the repetition
    #       - Number of events in the obs container
    #       - Upper-limit obtained for cross-section
    #         or lifetime
    #       - Actual scale factor obtained from the
    #         ratio of the differential flux to the
    #         theoretical flux
    #       - TS obtained during the likelihood calc.

    #   Creating the fits file
    dmfits     = gammalib.GFits()

    #   Creating table with rows equal to 
    #   the number of sims
    table      = gammalib.GFitsBinTable( args.nsims )

    #   Creating columns for every parameter to save
    col_runID  = gammalib.GFitsTableShortCol( 'RunID' , args.nsims )
    col_events = gammalib.GFitsTableDoubleCol( 'ObsEvents'  , args.nsims )

    if args.process == "decay" :

        col_cs = gammalib.GFitsTableDoubleCol( 'LogLifetime' , args.nsims )
    elif args.process == "annihilation" :

        col_cs = gammalib.GFitsTableDoubleCol( 'LogCrossSection' , args.nsims )
    else :

        print( 'Unknown process' )
        print( 'I will assume annihilation, but you can get\n\
            the scale factor and compute the parameter of interest' )

        col_cs = gammalib.GFitsTableDoubleCol( 'LogCrossSection' , args.nsims )

    col_scale  = gammalib.GFitsTableDoubleCol( 'ScaleFactor' , args.nsims )
    col_ts     = gammalib.GFitsTableDoubleCol( 'TS' , args.nsims )

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
            thiseref = args.mass * 0.5 / 2.
            thisemax = args.mass / 2.
        elif args.process == 'annihilation' :
            thiseref = args.mass * 0.5
            thisemax = args.mass
        else :
            print( 'Unknown process. I will assume DM annihilation' )
            thiseref = args.mass * 0.5
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

            onoffname      = auxMan.createname( args.outpath , \
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
            onoffgen[ 'coordsys' ]   = args.coordsys
            if args.coordsys == 'GAL' :
                onoffgen[ 'glon' ]   = args.pntlon
                onoffgen[ 'glat' ]   = args.pntlat
            elif args.coordsys == 'CEL' :
                onoffgen[ 'ra' ]     = args.pntlon
                onoffgen[ 'dec' ]    = args.pntlat
            else :
                print( 'Unknow coordsystem. I will use default value' )
                onoffgen[ 'glon' ]   = args.pntlon
                onoffgen[ 'glat' ]   = args.pntlat
            onoffgen[ 'rad' ]        = args.pntrad
            onoffgen[ 'bkgmethod' ]  = 'REFLECTED'
            onoffgen[ 'bkgregskip' ] = 0
            onoffgen[ 'bkgregmin' ]  = 1
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
                'LikeOutModel{:d}.xml'.format( thisrunID ) )
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

            #   Then, set the different results
            #   for every parameter
            col_runID[ sims ]  = thisrunID
            col_events[ sims ] = obslist[ 0 ].events().number()
            col_cs[ sims ]     = math.log10( cs )
            col_scale[ sims ]  = math.log10( scale )
            col_ts[ sims ]     = ts

            #file_name = args.name + '_' + args.process + '_' \
            #    + str( args.channel ) + '_' \
            #    + str( args.mass ) + 'TeV_' + str( thisrunID ) \
            #    + '_' + str( args.hours ) + 'h_' + str( args.sigmav ) + '.txt'
            #file_name = os.path.join( args.outpath , file_name )
            #with open( file_name , 'a' ) as outfile:
            #    outfile.write('%d\t%d\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\n' \
            #        % ( thisrunID , args.channel ,\
            #            obslist[0].events().number() ,\
            #            cs , scale , \
            #            ts , args.sigmav ) )

            # Free memory
            del limit
            del like
            del obssim
            del obs_binned
            del obslist
            del cube_name
        except RuntimeError:
            pass

    #   Above, end of loop for simulations

    #   Append columns to the table
    table.append( col_runID )
    table.append( col_events )
    table.append( col_cs )
    table.append( col_scale )
    table.append( col_ts )

    #   Now, create cards to save information
    #   about the DM model
    table.card( 'OBJECT' , args.name , 'Name of target' )

    if args.process == 'decay' :

        table.card( 'DMPROCESS' , args.process , 'Process... Dah!' )
        table.card( 'REFCROSSSECTION' , math.log( args.sigmav ) ,\
            'log10(sigmav/[cm**3/s])' , \
            'Logarithmic Annihilation cross-section' )

    elif args.process == 'annihilation' :

        table.card( 'DMPROCESS' , args.process , 'Process... Dah!' )
        table.card( 'REFLIFETIME' , math.log( args.sigmav ) ,\
            'log10(tau/[s])' , \
            'Logarithmic Decay Lifetime' )

    else :

        print( 'Unknown DM process' )
        print( 'Assuming annihilation ... (?)' )
        table.card( 'DMPROCESS' , 'Annihilation' , 'The actual process is ?' )

    table.card( 'CHANNEL' , args.channel , \
        'DM channel to generate gamma-rays' )
    table.card( 'DMMASS' , args.mass , 'TeV' , 'Mass of DM Candidate' )
    table.card( 'HOURS' , args.hours , 'h' , 'Observation Time' )

    #   Then, append table to fits file
    dmfits.append( table )

    #   And save file
    fits_file = auxMan.createname( args.outpath , args.fitsfile )
    dmfits.saveto( fits_file , True )

    print( 'ok!' )

