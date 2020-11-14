# ====================================================================
# This script perform the DM limits using the Cirelli et al. spectrum.
# ====================================================================

import ctools
import gammalib
import cscripts

import numpy as np

import os
import sys
import time
import argparse

from auxfunctions import *

#=====================#
# Routine entry point #
#=====================#

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
    source.add_argument( '--ROIra' ,\
        help='Right Ascension of ROI (in degrees)'  ,\
        type=float ,\
        required=True ,\
        metavar='47.0' )
    source.add_argument( '--ROIdec' ,\
        help='Declination of ROI (in degrees)' ,\
        type=float ,\
        required=True ,\
        metavar='39.0' )
    source.add_argument( '--ROIradius' ,\
        help='Radius of ROI (in degrees)' ,\
        type=float ,\
        required=True ,\
        metavar='10.0' )

    instrument = options.add_argument_group( 'CTA-Intrument options' , \
        'Information about CTA IRF, observation time, etc.' )
    instrument.add_argument( '--pntra' ,\
        help='Right Ascension of pointing (in degrees)' ,\
        type=float ,\
        required=True ,\
        metavar='47.0' )
    instrument.add_argument( '--pntdec' ,\
        help='Declination of pointing (in degrees)' ,\
        type=float ,\
        required=True ,\
        metavar='39.0' )
    instrument.add_argument( '--pntradius' ,\
        help='Radius of the observation (in degrees)' ,\
        type=float ,\
        required=True ,\
        metavar='10.0' )
    instrument.add_argument( '--caldb' ,\
        help='Database production for the IRF file. \
            Options are: [ prod2 , prod3b-v1 , prod3b-v2 ] ' ,\
        type=str ,\
        default='prod3b-v2' ,\
        metavar='prod3b-v2' ,\
        choices=[ 'prod2' , 'prod3b-v1' , 'prod3b-v2' ] )
    instrument.add_argument( '--irf' ,\
        help='Instrument Response Function for CTA' ,\
        type=str ,\
        default='North_z20_5h' ,\
        metavar='North_z20_50h' )
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
        default=0.03 ,\
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
    instrument.add_argument( '--pts' ,\
        help='Value of TS to compute Exclusion limits (TS < pts).' ,\
        type=float ,\
        default='25.0' ,\
        metavar='25.0' )
    instrument.add_argument( '--nsims' ,\
        help='Number of simulations' ,\
        type=int ,\
        required=True ,\
        metavar='10' )
    instrument.add_argument( '--nthreads' ,\
        help='Number of threads to speed the calculations' ,\
        type=int ,\
        metavar='10' ,\
        default=1 )
    instrument.add_argument( '--outpath' , \
        help='Path to save files' ,\
        type=str ,\
        default='./' ,\
        metavar='path/to/file/' )
    instrument.add_argument( '--is_onoff' ,\
        help='Indicate ON/OFF observation type' ,\
        action='store_true' )

    #   Create two parsers to manage the DM process selection
    subparsers = options.add_subparsers( help='Selecting DM processess \
        to produce gamma-rays' , dest='dmprocess' )
    subparsers.required = True

    #   Parser for annihilation
    anna       = subparsers.add_parser( 'anna' , help='DM annihilation Process' )
    anna.add_argument( '--sigmav' ,\
        help='Annihilation cross-section (cm**3/s)' ,\
        type=float ,\
        metavar='1.e-28' )
    anna.add_argument( '--jfactor' ,\
        help='Astrophysical J factor (in GeV**2/cm**5)' ,\
        type=float ,\
        metavar='1.e+20' )
    anna.add_argument( '--channel' ,\
        help='Annihilation channel (According to PPPC4DM Tables)' ,\
        type=int ,\
        metavar='11' ,\
        choices=[ 5 , 8 , 11 , 15 ] )
    anna.add_argument( '--mass' ,\
        help='DM mass (in TeV)' ,\
        type=float ,\
        metavar='100.0' )

    #   Parser for decay
    decs       = subparsers.add_parser( 'decs' , help='DM decay Process' )
    decs.add_argument( '--lifetime' ,\
        help='DM Decay Lifetime (s)' ,\
        type=float ,\
        metavar='1.e+35' )
    decs.add_argument( '--dfactor' ,\
        help='Astrophysical D factor (in GeV/cm**2)' ,\
        type=float ,\
        metavar='1.e+20' )
    decs.add_argument( '--channel' ,\
        help='Annihilation channel (According to PPPC4DM Tables)' ,\
        type=int ,\
        metavar='11' ,\
        choices=[ 5 , 8 , 11 , 15 ] )
    decs.add_argument( '--mass' ,\
        help='DM mass (in TeV)' ,\
        type=float ,\
        metavar='100.0' )


    #   Parsing arguments
    args = options.parse_args()

    MeVtoTeV = 1.e-6 ;

    #   This list put all the CTA-Background models
    bkgclasses = [ 'GCTAModelAeffBackground' , 'GCTAModelBackground' , \
        'GCTAModelCubeBackground' , 'GCTAModelIrfBackground' , \
        'GCTAModelRadialAcceptance' ]

    #   Check if xmlfile and outpath exist
    auxMan.checkDir( args.outpath )
    auxMan.checkPathandExit( args.inmodel )

    is_anna_or_decs = False
    process         = ''

    if args.dmprocess.lower() == 'anna' :
        is_anna_or_decs = True
        process         = 'Annihilation'

    elif args.dmprocess.lower() == 'decs' :

        is_anna_or_decs = False
        process         = 'Decay'

    #   Just printing useful information
    print( '************************************************' )
    print( '**    Mass:            {:.2e} TeV'.format( args.mass ) )
    print( '**    Process:         {:s}'.format( process ) )
    print( '**    Channel:         {:d}'.format( args.channel ) )

    if is_anna_or_decs :

        print( '**    Cross-Section:   {:.2e} cm**3/s'.format( args.sigmav ) )
        print( '**    J Factor:        {:.2e} GeV**2/cm**5'.format( args.jfactor ) )

    else :

        print( '**    Lifetime:        {:.2e} s'.format( args.lifetime ) )
        print( '**    D Factor:        {:.2e} GeV/cm**2'.format( args.dfactor ) )

    print( '************************************************' )

    #   Then, set a GModels container
    models = gammalib.GModels( args.inmodel )
    nameSources = [ source.name() \
        for source in models if source.classname() not in bkgclasses ]

    if args.gname in nameSources :

        print( '\n**    {:s} found in the models\n'.format( args.gname ) )

    else :

        print( '\n\tERROR: {:s} not found in Model container.\
            \n\tAvailbale models are:' )
        for name in nameSources :

            print( '\t - {:s}'.format( name ) )

        sys.exit()

    #   Extracting ts flag (calculation of TS)
    ts_flag = models[ args.gname ].tscalc()

    #   Extracting free parameters
    freepars = [ par.name() for par in models[ args.gname ].spectral() \
        if par.is_free() ]

    #   Set the GSkydir according to the coordinate system especified
    #   in the input parameters
    pdir = gammalib.GSkyDir()
    pdir.radec_deg( args.pntra , args.pntdec )


    #   Set duration of observation
    duration = args.hours * 3600.

    #   Set GEnergy objects and compute theoretical spectrum
    thiseref = 0.0

    if is_anna_or_decs :

        thiseref = args.mass / 2.0

    else :

        thiseref = args.mass / 4.0

    e_ref = gammalib.GEnergy( thiseref  , 'TeV' )
    e_min = gammalib.GEnergy( args.emin , 'TeV' )
    e_max = gammalib.GEnergy( args.emax , 'TeV' )

    #   Get expected dmflux at reference energy
    #   for all sources in Models
    theo_flux = 0.0

    for nameSource in nameSources :

        spec = models[ nameSource ].spectral()

        flux = spec.eval( e_ref , gammalib.GTime() )
        theo_flux += flux

    print( '\t* dF/dE(@ {:.2f} TeV): {:.3e} ph/(MeV cm**2 s)'.format( thiseref , theo_flux ) )

    #   Creating th fits file to save all the relevant results
    dmfits = gammalib.GFits()

    #   Creating table with number of rows equal to the number of sims
    table  = gammalib.GFitsBinTable( args.nsims )

    #   Creating columns for every parameter to save
    col_runID   = gammalib.GFitsTableShortCol( 'RunID' , args.nsims )
    col_events  = gammalib.GFitsTableDoubleCol( 'CubeEvents' , args.nsims )
    # col_oevents = gammalib.GFitsTableDoubleCol( 'ObsEvents' , args.nsims )
    col_ts      = gammalib.GFitsTableDoubleCol( 'TS' , args.nsims )
    col_scale   = gammalib.GFitsTableDoubleCol( 'ScaleFactor' , args.nsims )

    if is_anna_or_decs :

        col_cs  = gammalib.GFitsTableDoubleCol( 'LogCrossSection' , args.nsims )

    else :

        col_cs  = gammalib.GFitsTableDoubleCol( 'LogLifetime' , args.nsims )


    #   List of GFitsColumns
    colpars = []

    for freepar in freepars :

        col_par = gammalib.GFitsTableDoubleCol( freepar , args.nsims , 2 )
        colpars.append( col_par )

    #   At this moment I am not using this
    #   But, the idea is to try to use as an optional
    #   input parameter
    cubetype  = ''
    eventlist = ''


    #################################################################
    #####                       The Simuations
    #################################################################
    for nsim in range( args.nsims ) :

        #   ID to save in fits file
        thisrunID = args.id + '{:d}'.format( nsim )
        print( '\t**********************' )
        print( '\t*    Observation: {:d}'.format( nsim + 1 ) )
        print( '\t**********************' )

        #   GObservations Container
        obslist    = gammalib.GObservations()

        #   Setting observation parameters
        obs_binned = myCTAfuncs.single_obs( type=cubetype , eventfile=eventlist ,\
            pntdir=pdir , tstart=0.0 , duration=duration , deadc=0.95 ,\
            emin=args.emin , emax=args.emax , rad=args.pntradius ,\
            irf=args.irf , caldb=args.caldb , id=thisrunID , instrument='CTA' )

        #   Append obs_binned to obslist container
        obslist.append( obs_binned )

        #   And set the models to obslist container
        obslist.models( models )

        #   Name of the fits file to save events
        cubename = '{:s}_{:s}_{:d}_{:.2f}TeV_{:.2e}_Sim{:s}_{:d}min_cntcube.fits'.format( \
            args.gname , args.dmprocess , args.channel , args.mass ,\
            args.sigmav , thisrunID , int( duration ) )
        cubename = auxMan.createname( args.outpath , cubename )

        #   And now, perform the simulation using ctobssim
        #   and ctbin
        obssim = myCTAfuncs.sim( obslist , log=False , debug=False , \
            chatter=4 , edisp=False , seed=int( time.time() ) , \
            nbins=args.enumbins , binsz=0.02 , npix=200 , outfile=cubename )

        #   This is to save in the fits file
        events = obssim[ 0 ].events().number()

        #   Now, creating the observation definition XML file
        obsname = auxMan.createname( args.outpath , \
            'DMSimObs{:s}.xml'.format( thisrunID ) )

        obsxml  = gammalib.GXml()
        lists   = obsxml.append( 'observation_list title="observation library"' )
        obsss   = lists.append( 'observation name ="{:s}" id="{:s}"\
            instrument="CTA"'.format( args.gname , thisrunID ) )
        obsss.append( 'parameter name="EventList" file="{:s}"'.format( cubename ) )
        obsxml.save( obsname )

        #   Now, if the user request to generate an OnOff observation
        #   we use csphagen to handle this
        if args.is_onoff :

            onoffname      = auxMan.createname( args.outpath , \
                '{:s}OnOffObs.xml'.format( thisrunID ) )
            onoffmodelname = auxMan.createname( args.outpath ,\
                '{:s}OnOffModel.xml'.format( thisrunID ) )
            prefix         = os.path.join( args.outpath ,\
                '{:s}OnOff'.format( thisrunID ) )

            onoffgen       = cscripts.csphagen()

            onoffgen[ 'inobs' ]      = obsname
            onoffgen[ 'inmodel' ]    = args.inmodel
            onoffgen[ 'srcname' ]    = args.gname
            onoffgen[ 'caldb' ]      = args.caldb
            onoffgen[ 'irf' ]        = args.irf
            onoffgen[ 'outobs' ]     = onoffname
            onoffgen[ 'outmodel' ]   = onoffmodelname
            onoffgen[ 'ebinalg' ]    = 'LOG'
            onoffgen[ 'emin' ]       = args.emin
            onoffgen[ 'emax' ]       = args.emax
            onoffgen[ 'enumbins' ]   = args.enumbins
            onoffgen[ 'coordsys' ]   = 'CEL'
            onoffgen[ 'ra' ]         = args.ROIra
            onoffgen[ 'dec' ]        = args.ROIdec
            onoffgen[ 'rad' ]        = args.ROIradius
            onoffgen[ 'bkgmethod' ]  = 'REFLECTED'
            onoffgen[ 'bkgregskip' ] = 1
            onoffgen[ 'bkgregmin' ]  = 2
            onoffgen[ 'etruemin' ]   = args.emin
            onoffgen[ 'etruemax' ]   = args.emax
            onoffgen[ 'nthreads' ]   = args.nthreads
            onoffgen[ 'stack' ]      = False
            onoffgen[ 'prefix' ]     = prefix

            onoffgen.execute()

            #   Number of observed events.
            nobserved = onoffgen.obs()[ 0 ].nobserved()

        #   Now, Optimization using MLE
        print( '\t**    DM fitting...' )

        likeModelName = auxMan.createname( args.outpath ,\
            '{:s}LikeOutModel.xml'.format( thisrunID ) )

        like = ctools.ctlike()

        if args.is_onoff :

            like[ 'inobs' ]     = onoffname
            like[ 'inmodel' ]   = onoffmodelname

        else :

            like[ 'inobs' ]     = obsname
            like[ 'inmodel' ]   = args.inmodel

        like[ 'caldb' ]         = args.caldb
        like[ 'irf' ]           = args.irf
        like[ 'outmodel' ]      = likeModelName
        like[ 'like_accuracy' ] = 1.e-4
        like[ 'max_iter' ]      = 100
        like[ 'nthreads' ]      = args.nthreads
        like[ 'debug' ]         = False

        like.execute()

        #   Likelihood including DM Model
        #   (Alternative Hipotesis)
        like_alth = like.opt().value()

        #   Get the best fit value
        spec_fit  = like.obs().models()[ args.gname ].spectral()
        pars      = []

        for freepar in freepars :

            listpar = [ freepar , spec_fit[ freepar ].value() ,\
                spec_fit[ freepar ].error() ]
            pars.append( listpar )

        #Perform calculation of TS

        TS = 0.0

        if not ts_flag :

            #   Now, removing the source model to
            #   compute the likelihood for Null hypothesis
            like.obs().models().remove( args.gname )
            like.run()

            #   Likelihood of null hypothesis
            like_null = like.opt().value()

            #   And, calculation of TS
            TS        = - 2.0 * ( like_alth - like_null )

        else :

            TS = like.obs().models()[ args.gname ].ts()


        #   Now, performin calculation of UL if TS < 25

        #  Initializing Scale factor and UL
        scale        = 0.0
        parameter_UL = 0.0

        if TS < args.pts :

            print( '\t**    UL(95%C.L.) calculation' )

            limit = ctools.ctulimit()

            if args.is_onoff :

                limit[ 'inobs' ]   = onoffname
                limit[ 'inmodel' ] = onoffmodelname

            else :

                limit[ 'inobs' ]   = obsname
                limit[ 'inmodel' ] = args.inmodel

            limit[ 'srcname' ]     = args.gname
            limit[ 'eref' ]        = thiseref
            limit[ 'emin' ]        = args.emin
            limit[ 'emax' ]        = args.emax
            limit[ 'caldb' ]       = args.caldb
            limit[ 'irf' ]         = args.irf
            limit[ 'statistic' ]   = 'DEFAULT'

            limit.run()

            #   Get differential UL-flux
            diff_flux = limit.diff_ulimit()
            scale     = diff_flux / theo_flux

            if is_anna_or_decs :

                parameter_UL = args.sigmav * scale

            else :

                parameter_UL = args.lifetime / scale

        #   Saving to fits file
        col_runID[ nsim ]  = nsim + 1
        col_events[ nsim ] = events
        col_ts[ nsim ]     = TS

        #   If TS < 25, then save log(scale) and log(UL)
        #   else, save 0 for scale and UL
        if TS < args.pts :

            col_scale[ nsim ]  = np.log10( scale )
            col_cs[ nsim ]     = np.log10( parameter_UL )

        else :

            col_scale[ nsim ]  = 0.0
            col_cs[ nsim ]     = 0.0

        for ipar in range( len( pars ) ) :

            colpars[ ipar ][ nsim , 0 ] = pars[ ipar ][ 1 ]
            colpars[ ipar ][ nsim , 1 ] = pars[ ipar ][ 2 ]

        if nsim == 0 :

            print( 'Saving some files' )

        else :

            os.remove( obsname )
            os.remove( cubename )
            del obssim
            del obs_binned
            del obslist
            del obsxml

            if args.is_onoff :

                os.remove( onoffname )
                os.remove( onoffmodelname )

                os.remove( prefix + '_arf.fits' )
                os.remove( prefix + '_off.reg' )
                os.remove( prefix + '_on.reg' )
                os.remove( prefix + '_pha_off.fits' )
                os.remove( prefix + '_pha_on.fits' )
                os.remove( prefix + '_rmf.fits' )
                del onoffgen

            os.remove( likeModelName )
            del like

    print( '**    Finishing simulations' )

    #   Append columns to the table
    table.append( col_runID )
    table.append( col_events )
    table.append( col_ts )
    table.append( col_scale )
    table.append( col_cs )

    for ipar in range( len( colpars ) ) :

        table.append( colpars[ ipar ] )

    #   Cards to save information about the DM model
    table.card( 'OBJECT' , args.gname , 'Name of target' )

    #   Card for the DM process
    if is_anna_or_decs :

        table.card( 'DMPROCESS' , args.dmprocess , 'X(+) + X(-) -> gammas' )
        table.card( 'REF_CROSS_SECTION' , np.log10( args.sigmav ) ,\
            'log10(sigmav/[cm**3/s])' )

    else :

        table.card( 'DMPROCESS' , args.dmprocess , 'X -> gammas' )
        table.card( 'REF_LIFETIME' , np.log10( args.lifetime ) ,\
            'log10( lifetime/[s])' )

    table.card( 'CHANNEL' , args.channel , \
        'DM channel to generate gamma-rays' )
    table.card( 'DM_MASS' , args.mass , 'Mass of DM candidate (TeV)' )
    table.card( 'HOURS' , args.hours , 'Observation time (h)' )

    #   Then append to fits file
    dmfits.append( table )

    #   And save file
    fits_file = '{:s}_{:s}_{:s}_{:d}_Results.fits'.format( args.id , \
        args.gname ,args.dmprocess , args.channel )
    fits_file = auxMan.createname( args.outpath , fits_file )

    dmfits.saveto( fits_file , True )

    print( 'This is the END!...' )


