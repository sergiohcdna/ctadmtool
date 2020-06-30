#! /usr/bin/env python

# ====================================================================
# This script perform the DM limits using the Cirelli et al. spectrum.
# ====================================================================

import ctools 
from gammalib import *
import math

import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import colors
import numpy as np
import scipy.stats

import glob
import os
import sys
from sys import argv
import os.path
import csv
import time

def mean_confidence_interval( data , confidence=0.95 ) :
    n = len( data )
    m , se = np.mean( data ), scipy.stats.sem( data )
    h = se * scipy.stats.t.ppf( ( 1 + confidence ) / 2. , n - 1 )
    return m , m-h , m+h

# ============================== #
# Setup a single CTA observation #
#             T = 500h           #         
# ============================== #
def single_obs(type, eventfile, pntdir, tstart=0.0, duration=1800.0, deadc=0.95, \
                emin=0.03, emax=80.0, rad=5.0, \
                irf="South_50h", caldb="prod2", id="000001", instrument="CTA"):
    """
    Returns a single CTA observation.

    Parameters:
     pntdir   - Pointing direction [GSkyDir]
    Keywords:
     tstart     - Start time [seconds] (default: 0.0)
     duration   - Duration of observation [seconds] (default: 1800.0)
     deadc      - Deadtime correction factor (default: 0.95)
     emin       - Minimum event energy [TeV] (default: 0.1)
     emax       - Maximum event energy [TeV] (default: 100.0)
     rad        - ROI radius used for analysis [deg] (default: 5.0)
     irf        - Instrument response function (default: cta_dummy_irf)
     caldb      - Calibration database path (default: "dummy")
     id         - Run identifier (default: "000000")
     instrument - Intrument (default: "CTA")
    """
    # Allocate CTA observation
    obs_cta = GCTAObservation()

    # Set calibration database
    db = GCaldb()
    if (os.path.isdir(caldb)):
        db.rootdir(caldb)
    else:
        db.open("cta", caldb)

    # Set pointing direction
    pnt = GCTAPointing()
    pnt.dir(pntdir)
    obs_cta.pointing(pnt)

    # Set ROI
    roi     = GCTARoi()
    instdir = GCTAInstDir()
    instdir.dir(pntdir)
    roi.centre(instdir)
    roi.radius(rad)

    # Set GTI
    gti = GGti()
    gti.append(GTime(tstart), GTime(tstart+duration))

    # Set energy boundaries
    ebounds = GEbounds(GEnergy(emin, "TeV"), \
                                GEnergy(emax, "TeV"))

    # Allocate event cube
    if(type == "cube"):
        events = GCTAEventCube(eventfile);
        obs_cta.events(events);
    elif(type=='list'): 
        events = GCTAEventList(eventfile);
        obs_cta.events(events);
    else:
        events = GCTAEventList()
        events.roi(roi)
        events.gti(gti)
        events.ebounds(ebounds)
        obs_cta.events(events)

    # Set instrument response
    obs_cta.response(irf, db)

    # Set ontime, livetime, and deadtime correction factor
    obs_cta.ontime(duration)
    obs_cta.livetime(duration*deadc)
    obs_cta.deadc(deadc)
    obs_cta.id(id)

    # Return CTA observation
    return obs_cta


# ===================== #
# Simulate observations #
# ===================== #
def sim(obs, log=False, debug=False, chatter=2, edisp=False, seed=0, nbins=31,
        binsz=0.02, npix=50, proj="CAR", coord="CEL", outfile=""):
    """
    Simulate events for all observations in the container.

    Parameters:
     obs   - Observation container
    Keywords:
     log   - Create log file(s)
     debug - Create console dump?
     edisp - Apply energy dispersion?
     seed  - Seed value for simulations (default: 0)
     nbins - Number of energy bins (default: 0=unbinned)
     binsz - Pixel size for binned simulation (deg/pixel)
     npix  - Number of pixels in X and Y for binned simulation
    """

    # Allocate ctobssim application and set parameters
    sim = ctools.ctobssim(obs)
    sim["seed"] = seed
    sim["edisp"] = edisp
    sim["outevents"] = outfile
    sim[ 'nthreads' ] = 2

    # Optionally open the log file
    if log:
        sim.logFileOpen()

    # Optionally switch-on debugging model
    if debug:
        sim["debug"] = True

    # Set chatter level
    sim["chatter"] = chatter

    # Run ctobssim application. This will loop over all observations in the
    # container and simulation the events for each observation. Note that
    # events are not added together, they still apply to each observation
    # separately.
    sim.run()

    # Binned option?
    if nbins > 0:

        # Determine common energy boundaries for observations
        emin = None
        emax = None
        for run in sim.obs():
            run_emin = run.events().ebounds().emin().TeV()
            run_emax = run.events().ebounds().emax().TeV()
            if emin == None:
                emin = run_emin
            elif run_emin > emin:
                emin = run_emin
            if emax == None:
                emax = run_emax
            elif run_emax > emax:
                emax = run_emax


        # Allocate ctbin application and set parameters
        bin = ctools.ctbin(sim.obs())
        bin["ebinalg"] = "LOG"
        bin["emin"] = emin
        bin["emax"] = emax
        bin["enumbins"] = nbins
        bin["usepnt"] = True # Use pointing for map centre
        bin["nxpix"] = npix
        bin["nypix"] = npix
        bin["binsz"] = binsz
        bin["coordsys"] = coord
        bin["proj"] = proj
        bin[ 'nthreads' ] = 2

        # Optionally open the log file
        if log:
            bin.logFileOpen()

        # Optionally switch-on debugging model
        if debug:
            bin["debug"] = True

        # Set chatter level
        bin["chatter"] = chatter

        # Run ctbin application. This will loop over all observations in
        # the container and bin the events in counts maps
        bin.run()

        # Make a deep copy of the observation that will be returned
        # (the ctbin object will go out of scope one the function is
        # left)
        obs = bin.obs().copy()

    else:

        # Make a deep copy of the observation that will be returned
        # (the ctobssim object will go out of scope one the function is
        # left)
        obs = sim.obs().copy()

    # Optionally save file
    if (len(outfile) > 0):
        sim.save()

    # Delete the simulation
    del sim

    # Return observation container
    return obs

def log_likelihood_ONOFF( counts_on , source_on , bkg_on , counts_off , bkg_off ) :
    #Following convention from gammalib and ctools, the term for factorial of the
    #observed events is removed from the calculation of the likelihood
    like_on  = counts_on * np.log( source_on + bkg_on ) - ( source_on + bkg_on )
    like_off = counts_off * np.log( bkg_off ) - bkg_off
    return -( like_on + like_off )

def log_likelihood_ON( counts_on , source_on , bkg_on ) :
    #Following convention from gammalib and ctools, the term for factorial of the
    #observed events is removed from the calculation of the likelihood
    like_on  = counts_on * np.log( source_on + bkg_on ) - ( source_on + bkg_on )
    return -( like_on )

def ts_ONOFF( counts_on , source_on , bkg_on , counts_off , bkg_off ) :
    like_signal = log_likelihood_ONOFF( counts_on , source_on , bkg_on , counts_off , bkg_off )
    like_bkg    = log_likelihood_ONOFF( counts_on , 0 , bkg_on , counts_off , bkg_off )
    ts = -2.0 * ( like_signal - like_bkg )
    return ts

def ts_ON( counts_on , source_on , bkg_on ) :
    like_signal = log_likelihood_ON( counts_on , source_on , bkg_on )
    like_bkg    = log_likelihood_ON( counts_on , 0 , bkg_on )
    ts = -2.0 * ( like_signal - like_bkg )
    return ts

def plot_spectrum( filename , plotfile , model , source , con , coff , comp , \
        is_onoff , yup_lim , ydn_lim , xup_lim , xdn_lim , sensData , srcSens , \
        additional , addLabel , p_ts ) :
    """
    Plot spectrum

    Parameters
    ----------
    filename : str
        Name of spectrum FITS file
    plotfile : str
        Plot file name
    """
    # Read spectrum file  
    fits     = GFits( filename )
    table    = fits.table( 1 )
    thisfile = fits.filename().path() + fits.filename().file()
    thisfile = thisfile[ : -5 ] + 'ScanDataForBin'
    c_energy = table[ 'Energy' ]
    c_ed     = table[ 'ed_Energy' ]
    c_eu     = table[ 'eu_Energy' ]
    c_flux   = table[ 'Flux' ]
    c_eflux  = table[ 'e_Flux' ]
    c_ts     = table[ 'TS' ]
    c_upper  = table[ 'UpperLimit' ]

    # Initialise arrays to be filled
    energies    = []
    flux        = []
    ed_engs     = []
    eu_engs     = []
    e_flux      = []
    ul_energies = []
    ul_ed_engs  = []
    ul_eu_engs  = []
    ul_flux     = []
    thflux      = []
    allengs     = []
    likelihoods = []

    # Loop over rows of the file
    nrows = table.nrows()
    spectral = model[ source ].spectral()
    like = 0

    fig , ax = plt.subplots()
    plt.xscale( 'log' )
    plt.yscale( 'log' )

    for row in range( nrows ) :

        file = thisfile + '{:d}.txt'.format( row + 1 ) 

        thisdata = np.genfromtxt( file , comments='#' , delimiter='\t' , names=True )
        thisdata = np.array( ( thisdata[ 'Normalization' ] , thisdata[ 'TS' ] ) )
        xdata = np.zeros( ( 1 , 2 ) )
        xdata[ 0 ][ 0 ] = c_energy.real( row ) - c_ed.real( row )
        xdata[ 0 ][ 1 ] = c_energy.real( row ) + c_eu.real( row )
        yvalue = thisdata[ 0 ]
        points = thisdata[ 1 ]
        X , Y = np.meshgrid( xdata , yvalue )
        Z = np.zeros( ( 2 , len( yvalue ) ) )
        for i in range( len( yvalue ) ) :
            Z[ 0 ][ i ] = points[ i ]
            Z[ 1 ][ i ] = points[ i ]

        # Get Test Statistic, flux and flux error
        ts    = c_ts.real( row )
        flx   = c_flux.real( row )
        e_flx = c_eflux.real( row )

        allengs.append( c_energy.real( row ) )
        thiseng = GEnergy( c_energy.real( row ) , "TeV" )
        thf = spectral.eval( thiseng , GTime() )
        engd = c_energy.real( row ) * 1.e+6 # energy from TeV to ergs
        thflux.append( engd**2 * thf * MeV2erg )

        # If Test Statistic is larger than 25 and flux error is smaller than
        # flux then append flux plots ...
        if ts > p_ts and e_flx < flx:
            energies.append( c_energy.real( row ) )
            flux.append( c_flux.real( row) )
            ed_engs.append( c_ed.real( row) )
            eu_engs.append( c_eu.real( row) )
            e_flux.append( c_eflux.real( row) )

        # ... otherwise append upper limit
        else:
            ul_energies.append( c_energy.real( row ) )
            ul_flux.append( c_upper.real( row ) )
            ul_ed_engs.append( c_ed.real( row ) )
            ul_eu_engs.append( c_eu.real( row ) )

        caxes = ax.pcolormesh( X , Y , Z.T , vmin=-4.0 , vmax=0 , cmap='Blues' )
        if row == 0 :
            ax.figure.colorbar( caxes , ax=ax , label='TS Difference' )

        ###Section to compute likelihood for events
        #the calculation is over the energy bins

        # if is_onoff :
        #     l = log_likelihood_ONOFF( con[ row ][ 0 ] , comp[ row ][ 0 ] , \
        #         comp[ row ][ 1 ] , coff[ row ][ 0 ] , coff[ row ][ 1 ] )
        #     l0 = log_likelihood_ONOFF( con[ row ][ 0 ] , 0 , \
        #         comp[ row ][ 1 ] , coff[ row ][ 0 ] , coff[ row ][ 1 ] )
        #     likelihoods.append( l )
        #     myts = ts_ON( con[ row ][ 0 ] , comp[ row ][ 0 ] , comp[ row ][ 1 ] )
        #     like += l
        # else :
        #     l = log_likelihood_ON( con[ row ][ 0 ] , comp[ row ][ 0 ] , comp[ row ][ 1 ] ) 
        #     l0 = log_likelihood_ON( con[ row ][ 0 ] , 0 , comp[ row ][ 1 ] ) 
        #     likelihoods.append( l )
        #     myts = ts_ON( con[ row ][ 0 ] , comp[ row ][ 0 ] , comp[ row ][ 1 ] )
        #     like += l



    # Set upper limit errors
    yerr = [ 0.6 * x for x in ul_flux ]

    # Plot the spectrum 
    # plt.grid()
    plt.errorbar( energies , flux , yerr=e_flux , xerr=[ ed_engs , eu_engs ] , fmt='ro' , label='Flux Points' )
    plt.errorbar( ul_energies , ul_flux , xerr=[ ul_ed_engs , ul_eu_engs ] , yerr=yerr , uplims=True , fmt='ro' , label='UL (C.L.=0.95)' )
    plt.plot( allengs , thflux , color='black' , linewidth=2 , label='Model' )
    plt.xlabel( 'Energy (TeV)' )
    plt.ylabel( r'E$^2$ $\times$ dN/dE (erg cm$^{-2}$ s$^{-1}$)' )
    plt.ylim( ydn_lim , yup_lim )
    plt.xlim( xdn_lim , xup_lim )

    plt.plot( np.power( 10 , ( sensData[ 'loge' ] ) ) , sensData[ 'sensitivity' ] , color=( 0.57 , 0.36 , 0.51 ) , linewidth=2 , label=srcSens )

    if additional.size == 0 :
        print( 'No additional Data in the plot' )
    else :
        plt.plot( additional[ 'e' ] , additional[ 'e2dnde' ] / 0.624 , color=( 0.29 , 0.33 , 0.13 ) , linewidth=2 , label=addLabel )

    plt.legend( loc='best' , prop={'size': 12} )

    # lev_exp = np.arange( np.floor( np.log10( 5.e+7 ) - 1 ) , \
    #     np.ceil( np.log10( 1.e+2 ) - 1 ) )
    # levels  = np.power( 10 , lev_exp )
    # cs      = ax.contourf( allengs , thflux , likelihoods , levels , norm=colors.LogNorm() )
    # cbar    = fig.colorbar( cs )



    # Show figure
    if len( plotfile ) > 0:
        plt.savefig( plotfile )
    else:
        plt.show()

    plt.close()

    # Return
    return

# ======================= #
# Set residual label text #
# ======================= #
def set_residual_label(algorithm):
    """
    Set residual label text for a given algorithm.

    Parameters
    ----------
    algorithm : str
        Algorithm

    Returns
    -------
    label : str
        Label text
    """
    # Set algorithm dependent residual label text
    if algorithm == 'SUB':
        label = 'Residuals (counts)'
    elif algorithm == 'SUBDIV':
        label = 'Residuals (fraction)'
    elif algorithm == 'SUBDIVSQRT' or algorithm == 'SIGNIFICANCE':
        label = r'Residuals ($\sigma$)'
    else:
        label = ''

    # Return label
    return label


# ================================= #
# Fill counts/model/residual arrays #
# ================================= #
def fill_cmr(row, counts, model, resid, e_counts, e_resid, c_counts, c_model,
             c_resid, algorithm):
    """
    Helper function that fills in the counts, model, residuals lists and
    calculate relative errors for given table row

    Parameters
    ----------
    row : int
        Table row index
    counts : list of float
        Counts
    model : list of float
        Model
    resid : list of float
        Residual
    e_counts : list of float
        Counts error
    e_resid : list of float
        Residual error
    c_counts : `~gammalib.GFitsTableCol'
        Counts FITS table column
    c_model : `~gammalib.GFitsTableCol'
        Model FITS table column
    c_resid : `~gammalib.GFitsTableCol'
        Residual FITS table column
    algorithm : str
        Algorithm

    Returns
    -------
    counts, model, resid, e_counts, e_resid : tuple of float
        Counts, model and residual arrays
    """
    # Extract values
    counts.append(c_counts.real(row))
    model.append(c_model.real(row))
    resid.append(c_resid.real(row))

    # Calculate count error
    err = math.sqrt(c_counts.real(row))
    if err == 1:
        err = 0.99 # This prevents visualization problem in matplotlib
    e_counts.append(err)

    # Calculate residual error
    if algorithm == 'SUB':
        e_resid.append(err)
    elif algorithm == 'SUBDIV':
        if c_model.real(row) > 0.:
            e_resid.append(err / c_model.real(row))
        else:
            e_resid.append(0.)
    elif algorithm == 'SUBDIVSQRT':
        if c_model.real(row) > 0.:
            e_resid.append(err / math.sqrt(c_model.real(row)))
        else:
            e_resid.append(0.)
    elif algorithm == 'SIGNIFICANCE':
        e_resid.append(1.0)
    else:
        e_resid.append(0.0)

    # Return tuple
    return counts, model, resid, e_counts, e_resid


# ====================== #
# Plot residual spectrum #
# ====================== #
def plot_residuals( filename , plotfile , hdu ) :
    """
    Plot spectrum

    Parameters
    ----------
    filename : str
        Name of spectrum FITS file
    plotfile : str
        Plot file name
    hdu : int
        Number of observation to plot
    """
    # Read spectrum file    
    fits     = GFits(filename)
    table    = fits.table(1 + hdu)
    c_emin   = table['Emin']
    c_emax   = table['Emax']
    c_counts = table['Counts']
    c_model  = table['Model']
    c_resid  = table['Residuals']
    try:
        c_counts_off = table['Counts_Off']
        c_model_off  = table['Model_Off']
        c_resid_off  = table['Residuals_Off']
        is_onoff     = True
    except:
        is_onoff = False

    # Initialise arrays to be filled
    ebounds      = []
    ed_engs      = []
    eu_engs      = []
    em_engs      = []
    counts       = []
    e_counts     = []
    model        = []
    resid        = []
    e_resid      = []
    counts_off   = []
    e_counts_off = []
    model_off    = []
    resid_off    = []
    e_resid_off  = []

    # Residual algorithm
    algorithm = table.header()['ALGORITHM'].string()

    # Loop over rows of the file
    nrows = table.nrows()

    # Add first energy boundary
    ebounds.append(c_emin.real(0))
    for row in range(nrows):

        # Boundaries
        ebounds.append(c_emax.real(row))

        # Geometrical mean energy and errors
        emean = math.sqrt(c_emin.real(row) * c_emax.real(row))
        em_engs.append(emean)
        ed_engs.append(emean - c_emin.real(row))
        eu_engs.append(c_emax.real(row) - emean)

        # counts, model, residuals
        counts, model, resid, e_counts, e_resid = fill_cmr(row, counts, model,
                                                           resid, e_counts,
                                                           e_resid, c_counts,
                                                           c_model, c_resid,
                                                           algorithm)
        # On/off
        if is_onoff:
            counts_off, model_off, resid_off, e_counts_off, e_resid_off = fill_cmr(
                row, counts_off,
                model_off, resid_off,
                e_counts_off,
                e_resid_off,
                c_counts_off,
                c_model_off, c_resid_off,
                algorithm)

    # Add model value to be compatible with plt.step
    model = [model[0]] + model
    if is_onoff:
        model_off = [model_off[0]] + model_off

    # Initialise figure
    axarr = []
    if is_onoff:
        f   = plt.figure(figsize=(12,7))
        ax1 = f.add_subplot(221)
        ax2 = f.add_subplot(222)
        ax3 = f.add_subplot(223)
        ax4 = f.add_subplot(224)
        axarr.append(ax1)
        axarr.append(ax2)
        axarr.append(ax3)
        axarr.append(ax4)
    else:
        f   = plt.figure()
        ax1 = f.add_subplot(211)
        ax2 = f.add_subplot(212)
        axarr.append(ax1)
        axarr.append(ax2)
        f.subplots_adjust(hspace=0)
    axarr[0].set_xscale('log')
    axarr[1].set_xscale('log')
    axarr[0].set_yscale('log')

    # Counts and model
    axarr[0].errorbar(em_engs, counts, yerr=e_counts, xerr=[ed_engs, eu_engs],
                      fmt='ko', capsize=0, linewidth=2, zorder=2, label='Data')
    axarr[0].step(ebounds, model, color='0.5', linewidth=2, zorder=1,
                  label='Model')
    axarr[0].set_ylabel('Counts')

    # Residuals
    axarr[1].errorbar(em_engs, resid, yerr=e_resid, xerr=[ed_engs, eu_engs],
                      fmt='ko', capsize=0, linewidth=2, zorder=2)
    axarr[1].axhline(0, color='0.5', linestyle='--')
    axarr[1].set_xlabel('Energy (TeV)')
    axarr[1].set_ylabel(set_residual_label(algorithm))

    # On/Off
    if is_onoff:

        # Set titles and axes labels
        axarr[0].set_title('ON')
        axarr[2].set_title('OFF')
        axarr[2].set_ylabel('Counts')
        axarr[2].set_xscale('log')
        axarr[3].set_xscale('log')
        axarr[2].set_yscale('log')
        axarr[3].set_ylabel(set_residual_label(algorithm))

        # Counts and model
        axarr[2].errorbar(em_engs, counts_off, yerr=e_counts_off,
                          xerr=[ed_engs, eu_engs],
                          fmt='ko', capsize=0, linewidth=2, zorder=2,
                          label='Data')
        axarr[2].step(ebounds, model_off, color='0.5', linewidth=2, zorder=1,
                      label='Model')
        axarr[2].set_xlabel('Energy (TeV)')

        # Residuals
        axarr[3].errorbar(em_engs, resid_off, yerr=e_resid_off,
                          xerr=[ed_engs, eu_engs],
                          fmt='ko', capsize=0, linewidth=2, zorder=2)
        axarr[3].axhline(0, color='0.5', linestyle='--')
        axarr[3].set_xlabel('Energy (TeV)')

    # Add spectra of individual components. Skip all standard columns.
    skiplist = ['Counts', 'Model', 'Residuals', 'Counts_Off', 'Model_Off',
                'Residuals_Off', 'Emin', 'Emax']
    for s in range(table.ncols()):
        if table[s].name() not in skiplist:
            component = []
            for row in range(nrows):
                component.append(table[s].real(row))
            component = [component[0]] + component
            axarr[0].step(ebounds, component, zorder=0, label=table[s].name())

    # Add legend
    axarr[0].legend(loc='best', prop={'size': 10})

    # Show figure
    if len(plotfile) > 0:
        plt.savefig(plotfile)
    else:
        plt.show()

    plt.close()
    # Return
    return

def npfill( row , c_counts , c_model , c_resid , algorithm ) :
    """
    Helper function that fills in the counts, model, residuals lists and
    calculate relative errors for given table row

    Parameters
    ----------
    row : int
        Table row index
    counts : list of float
        Counts
    model : list of float
        Model
    resid : list of float
        Residual
    e_counts : list of float
        Counts error
    e_resid : list of float
        Residual error
    c_counts : `~gammalib.GFitsTableCol'
        Counts FITS table column
    c_model : `~gammalib.GFitsTableCol'
        Model FITS table column
    c_resid : `~gammalib.GFitsTableCol'
        Residual FITS table column
    algorithm : str
        Algorithm

    Returns
    -------
    counts, model, resid, e_counts, e_resid : tuple of float
        Counts, model and residual arrays
    """
    # Extract values
    counts = c_counts.real( row ) 
    model  = c_model.real( row ) 
    resid  = c_resid.real( row ) 

    # Calculate count error
    e_counts = np.sqrt( c_counts.real( row ) )
    if e_counts == 1:
        e_counts = 0.99 # This prevents visualization problem in matplotlib

    # Calculate residual error
    if algorithm == 'SUB':
        e_resid =  e_counts 
    elif algorithm == 'SUBDIV':
        if c_model.real( row ) > 0.:
            e_resid = e_counts / c_model.real( row )
        else:
            e_resid = np.append( 0. )
    elif algorithm == 'SUBDIVSQRT':
        if c_model.real( row ) > 0.:
            e_resid = e_counts / math.sqrt( c_model.real( row ) )
        else:
            e_resid = 0.
    elif algorithm == 'SIGNIFICANCE':
        e_resid = 1.0
    else:
        e_resid = e_resid , 0.0 

    # Return tuple
    return counts, model, resid, e_counts, e_resid


def get_residuals( filename , hdu ) :

    fits     = GFits(filename)
    table    = fits.table(1 + hdu)
    c_emin   = table['Emin']
    c_emax   = table['Emax']
    c_counts = table['Counts']
    c_model  = table['Model']
    c_resid  = table['Residuals']
    try:
        c_counts_off = table['Counts_Off']
        c_model_off  = table['Model_Off']
        c_resid_off  = table['Residuals_Off']
        is_onoff     = True
    except:
        is_onoff = False

    # Residual algorithm
    algorithm = table.header()['ALGORITHM'].string()

    # Loop over rows of the file
    nrows = table.nrows()
    ncols = table.ncols()

    skiplist = ['Counts', 'Model', 'Residuals', 'Counts_Off', 'Model_Off',
        'Residuals_Off', 'Emin', 'Emax']

    lcomps  = [ table[ s ].name() for s in range( ncols ) if table[ s ].name() not in skiplist ]
    ncomps = len( lcomps )

    energybins = np.zeros( shape=( nrows , 3 ) )
    oncounts   = np.zeros( shape=( nrows , 5 ) )
    offcounts  = np.zeros( shape=( nrows , 5 ) )
    components = np.zeros( shape=( nrows , ncomps ) )


    # Add first energy boundary
    # ebounds = np.append( ebounds , c_emin.real( 0 ) )
    for row in range(nrows):

        # Geometrical mean energy and errors
        emean = math.sqrt( c_emin.real( row ) * c_emax.real( row ) )
        energybins[ row , 0 ] = emean - c_emin.real( row )
        energybins[ row , 1 ] = emean
        energybins[ row , 2 ] = c_emax.real( row ) - emean 

        # counts, model, residuals
        counts, model, resid, e_counts, e_resid = npfill( row, c_counts, c_model , \
            c_resid , algorithm )
        oncounts[ row , 0 ] = counts
        oncounts[ row , 1 ] = model
        oncounts[ row , 2 ] = resid
        oncounts[ row , 3 ] = e_counts
        oncounts[ row , 4 ] = e_resid
        # On/off
        if is_onoff:
            counts_off, model_off, resid_off, e_counts_off, e_resid_off = npfill(
                row, c_counts_off , c_model_off , c_resid_off , algorithm )
            offcounts[ row , 0 ] = counts_off
            offcounts[ row , 1 ] = model_off
            offcounts[ row , 2 ] = resid_off
            offcounts[ row , 3 ] = e_counts_off
            offcounts[ row , 4 ] = e_resid_off

    c = 0
    for comp in lcomps :
        for row in range( nrows ) :
            components[ row , c ] = table[ comp ].real( row )
        c += 1

    if is_onoff :
        return energybins , oncounts , offcounts , components , lcomps
    else :
        return energybins , oncounts , components , lcomps

def get_mean_residuals( data ) :

    if np.ndim( data ) != 3 :
        print( 'Incorrect number of dimensions for data' )
        sys.exit( 1 )

    ss = data.shape
    nbins = ss[ 0 ]
    columns = ss[ 1 ]

    newdata = np.empty( shape=( nbins , columns , 5 ) )

    for column in range( columns ) :
        thisarray = np.empty( shape=( nbins , 5 ) )
        for nbin in range( nbins ) :
            c , c1 , c2 = mean_confidence_interval( 
                data[ nbin , column , : ] , 0.68 )
            d , d1 , d2 = mean_confidence_interval( 
                data[ nbin , column , : ] , 0.95 )
            thisarray[ nbin , : ] = np.array( [ d1 , c1 , c , c2 , d2 ] )
        newdata[ : , column , : ] = thisarray[ : , : ]

    return newdata 



def plot_Mean_residuals( con , coff , comp , energybins , plotfile , is_onoff , \
    algorithm , label_models ) :


    axarr = []

    if is_onoff:
        f   = plt.figure( figsize=( 16 , 9 ) )
        ax1 = f.add_subplot( 221 )
        ax2 = f.add_subplot( 222 )
        ax3 = f.add_subplot( 223 )
        ax4 = f.add_subplot( 224 )
        axarr.append( ax1 )
        axarr.append( ax2 )
        axarr.append( ax3 )
        axarr.append( ax4 )
    else:
        f   = plt.figure( figsize=( 9 , 5 ) )
        ax1 = f.add_subplot( 211 )
        ax2 = f.add_subplot( 212 )
        axarr.append( ax1 )
        axarr.append( ax2 )
        f.subplots_adjust( hspace=0 )
    axarr[ 0 ].set_xscale( 'log' )
    axarr[ 1 ].set_xscale( 'log' )
    axarr[ 0 ].set_yscale( 'log' )

    # Counts and model
    axarr[ 0 ].errorbar( energybins[ : , 1 ] , \
        con[ : , 0 , 2 ] , \
        yerr=con[ : , 3 , 2 ] , \
        xerr=[ energybins[ : , 0 ] , energybins[ : , 2 ] ] ,\
        fmt='ko' , capsize=0 , linewidth=2 , zorder=2 , label='Data' )
    axarr[ 0 ].plot( energybins[ : , 1 ] , con[ : , 1 , 2 ] , \
        color='green' , linewidth=2 , zorder=1 , label='Model' )
    axarr[ 0 ].fill_between( energybins[ : , 1 ] ,\
        con[ : , 0 , 0 ] , con[ : , 0 , 4 ] , \
        facecolor='gray' , alpha=0.95 , interpolate=True )
    axarr[ 0 ].fill_between( energybins[ : , 1 ] ,\
        con[ : , 0 , 1 ] , con[ : , 0 , 3 ] , \
        facecolor='gray' , alpha=0.75 , interpolate=True )
    axarr[ 0 ].fill_between( energybins[ : , 1 ] ,\
        con[ : , 1 , 0 ] , con[ : , 1 , 4 ] , \
        facecolor='green' , alpha=0.95 , interpolate=True )
    axarr[ 0 ].fill_between( energybins[ : , 1 ] ,\
        con[ : , 1 , 1 ] , con[ : , 1 , 3 ] , \
        facecolor='green' , alpha=0.75 , interpolate=True )
    axarr[ 0 ].set_ylabel( 'Counts' )

    # Residuals
    axarr[ 1 ].errorbar( energybins[ : , 1 ] , \
        con[ : , 2 , 2 ] , yerr=con[ : , 4 , 2 ] , \
        xerr=[ energybins[ : , 0 ] , energybins[ : , 2 ] ] ,\
        fmt='ko' , capsize=0 , linewidth=2 , zorder=2 )
    axarr[ 1 ].fill_between( energybins[ : , 1 ] ,\
        con[ : , 2 , 0 ] , con[ : , 2 , 4 ] , \
        facecolor='gray' , alpha=0.95 )
    axarr[ 1 ].fill_between( energybins[ : , 1 ] ,\
        con[ : , 2 , 1 ] , con[ : , 2 , 3 ] , \
        facecolor='gray' , alpha=0.75 )
    axarr[ 1 ].axhline( 0 , color='0.5' , linestyle='--' )
    axarr[ 1 ].set_xlabel( 'Energy (TeV)' )
    axarr[ 1 ].set_ylabel( set_residual_label( algorithm ) )

    # On/Off
    if is_onoff:

        # Set titles and axes labels
        axarr[0].set_title( 'ON' )
        axarr[2].set_title( 'OFF' )
        axarr[2].set_ylabel( 'Counts' )
        axarr[2].set_xscale( 'log' )
        axarr[3].set_xscale( 'log' )
        axarr[2].set_yscale( 'log' )
        axarr[3].set_ylabel( set_residual_label( algorithm ) )

        # Counts and model
        axarr[ 2 ].errorbar( energybins[ : , 1 ] , \
            coff[ : , 0 , 2 ] , \
            yerr=coff[ : , 3 , 2 ] , \
            xerr=[ energybins[ : , 0 ] , energybins[ : , 2 ] ] ,\
            fmt='ko' , capsize=0 , linewidth=2 , zorder=2 , label='Data' )
        axarr[ 2 ].plot( energybins[ : , 1 ] , coff[ : , 1 , 2 ] , \
            color='green' , linewidth=2 , zorder=1 , label='Model' )
        axarr[ 2 ].fill_between( energybins[ : , 1 ] ,\
            coff[ : , 0 , 0 ] , coff[ : , 0 , 4 ] , \
            facecolor='gray' , alpha=0.95 , interpolate=True )
        axarr[ 2 ].fill_between( energybins[ : , 1 ] ,\
            coff[ : , 0 , 1 ] , coff[ : , 0 , 3 ] , \
            facecolor='gray' , alpha=0.75 , interpolate=True )
        axarr[ 2 ].fill_between( energybins[ : , 1 ] ,\
            coff[ : , 1 , 0 ] , coff[ : , 1 , 4 ] , \
            facecolor='green' , alpha=0.95 , interpolate=True )
        axarr[ 2 ].fill_between( energybins[ : , 1 ] ,\
            coff[ : , 1 , 1 ] , coff[ : , 1 , 3 ] , \
            facecolor='green' , alpha=0.75 , interpolate=True )
        axarr[2].set_xlabel( 'Energy (TeV)' )

        # Residuals
        axarr[ 3 ].errorbar( energybins[ : , 1 ] , \
            coff[ : , 2 , 2 ] , yerr=coff[ : , 4 , 2 ] , \
            xerr=[ energybins[ : , 0 ] , energybins[ : , 2 ] ] ,\
            fmt='ko' , capsize=0 , linewidth=2 , zorder=2 )
        axarr[ 3 ].fill_between( energybins[ : , 1 ] ,\
            coff[ : , 2 , 0 ] , coff[ : , 2 , 4 ] , \
            facecolor='gray' , alpha=0.95 , interpolate=True )
        axarr[ 3 ].fill_between( energybins[ : , 1 ] ,\
            coff[ : , 2 , 1 ] , coff[ : , 2 , 3 ] , \
            facecolor='gray' , alpha=0.75 , interpolate=True )
        axarr[ 3 ].axhline( 0 , color='0.5' , linestyle='--' )
        axarr[ 3 ].set_xlabel( 'Energy (TeV)' )

    for s in range( comp.shape[ 1 ] ) :
        axarr[ 0 ].step( energybins[ : , 1 ] , comp[ : , s ] , \
                zorder=0 , label=label_models[ s ] )

    # Add legend
    axarr[ 0 ].legend( loc='best' , prop={ 'size': 10 } )

    # Show figure
    if len( plotfile ) > 0:
        plt.savefig( plotfile )
    else:
        plt.show()

    # Return
    return


def plot_onedDist( data , nbins , is_log , color , xlabel , nsims , plotfile ) :

    f = plt.figure( figsize=( 9 , 5 ) )
    ax = f.add_subplot( 111 )

    mean , c , d = mean_confidence_interval( data , confidence=0.68 )

    ax.hist( data , bins=nbins , log=is_log , color=color )
    ax.axvline( mean , color='black' , linewidth=2.0 )
    ax.axvspan( c , d , color='black' , alpha=0.65 )
    if len( xlabel ) > 0 :
        ax.set_xlabel( xlabel )
    ax.set_ylabel( 'Number of simulations' )
    ax.set_ylim( 0 , nsims * 0.75 )
    if len( plotfile ) > 0 :
        plt.savefig( plotfile )
    else :
        plt.show()

    return


def get_spectrum_data( filename ) :
    """
    Getting spectrum points in order to save for later calculations.
    Useful when it's needed to compute the mean and standar deviation
    for several simulation

    Parameters
    ----------
    filename : str
        Name of spectrum FITS file
    """
    # Read spectrum file    
    fits     = GFits( filename )
    table    = fits.table( 1 )
    c_energy = table[ 'Energy' ]
    c_ed     = table[ 'ed_Energy' ]
    c_eu     = table[ 'eu_Energy' ]
    c_flux   = table[ 'Flux' ]
    c_eflux  = table[ 'e_Flux' ]
    c_ts     = table[ 'TS' ]
    c_upper  = table[ 'UpperLimit' ]

    # Initialise arrays to be filled

    energies    = np.empty( shape=( 0 ) )
    flux        = np.empty( shape=( 0 ) )
    ed_engs     = np.empty( shape=( 0 ) )
    eu_engs     = np.empty( shape=( 0 ) )
    e_flux      = np.empty( shape=( 0 ) )
    ul_energies = np.empty( shape=( 0 ) )
    ul_ed_engs  = np.empty( shape=( 0 ) )
    ul_eu_engs  = np.empty( shape=( 0 ) )
    ul_flux     = np.empty( shape=( 0 ) )

    # Loop over rows of the file
    nrows = table.nrows()
    ncols = table.ncols()

    for row in range(nrows) :

        # Get Test Statistic, flux and flux error
        ts    = c_ts.real( row )
        flx   = c_flux.real( row )
        e_flx = c_eflux.real( row )

        # If Test Statistic is larger than 9 and flux error is smaller than
        # flux then append flux plots ...
        if ts > 9.0 and e_flx < flx :
            energies = np.append( energies , c_energy.real( row ) )
            flux     = np.append( flux , c_flux.real( row ) )
            ed_engs  = np.append( ed_engs , c_ed.real( row ) )
            eu_engs  = np.append( eu_engs , c_eu.real( row ) )
            e_flux   = np.append( e_flux , c_eflux.real( row ) )

        # ... otherwise append upper limit
        else :
            ul_energies = np.append( ul_energies , c_energy.real( row ) )
            ul_flux     = np.append( ul_flux , c_upper.real( row ) )
            ul_ed_engs  = np.append( ul_ed_engs , c_ed.real( row ) )
            ul_eu_engs  = np.append( ul_eu_engs , c_eu.real( row ) )

    energyFlux = np.array( [ energies , ed_engs , eu_engs ] )
    dataFlux   = np.array( [ flux , e_flux ] )
    ulenegFlux = np.array( [ ul_energies , ul_ed_engs , ul_eu_engs ] )
    uldataFlux = np.array( [ ul_flux , 0.6 * ul_flux ] )

    # Return
    return energyFlux , dataFlux , ulenegFlux , uldataFlux


def get_mean_spectrum( data ) :

    if np.ndim( data ) != 3 :
        print( 'Incorrect number of dimensions for data' )
        sys.exit( 1 )

    ss = data.shape
    nbins = ss[ 2 ]
    columns = ss[ 1 ]

    newdata = np.empty( shape=( columns , nbins , 5 ) )

    for column in range( columns ) :
        thisarray = np.empty( shape=( nbins , 5 ) )
        for nbin in range( nbins ) :
            c , c1 , c2 = mean_confidence_interval( 
                data[ : , column , nbin ] , 0.68 )
            d , d1 , d2 = mean_confidence_interval( 
                data[ : , column , nbin ] , 0.95 )
            thisarray[ nbin , : ] = np.array( [ d1 , c1 , c , c2 , d2 ] )
        newdata[ column , : , : ] = thisarray[ : , : ]

    return newdata 


# def plot_spectrum(filename, plotfile):
#     """
#     Plot spectrum

#     Parameters
#     ----------
#     filename : str
#         Name of spectrum FITS file
#     plotfile : str
#         Plot file name
#     """
#     # Read spectrum file    
#     fits     = GFits(filename)
#     table    = fits.table(1)
#     c_energy = table['Energy']
#     c_ed     = table['ed_Energy']
#     c_eu     = table['eu_Energy']
#     c_flux   = table['Flux']
#     c_eflux  = table['e_Flux']
#     c_ts     = table['TS']
#     c_upper  = table['UpperLimit']

#     # Initialise arrays to be filled
#     energies    = []
#     flux        = []
#     ed_engs     = []
#     eu_engs     = []
#     e_flux      = []
#     ul_energies = []
#     ul_ed_engs  = []
#     ul_eu_engs  = []
#     ul_flux     = []

#     # Loop over rows of the file
#     nrows = table.nrows()
#     for row in range(nrows):

#         # Get Test Statistic, flux and flux error
#         ts    = c_ts.real(row)
#         flx   = c_flux.real(row)
#         e_flx = c_eflux.real(row)

#         # If Test Statistic is larger than 9 and flux error is smaller than
#         # flux then append flux plots ...
#         if ts > 9.0 and e_flx < flx:
#             energies.append(c_energy.real(row))
#             flux.append(c_flux.real(row))
#             ed_engs.append(c_ed.real(row))
#             eu_engs.append(c_eu.real(row))
#             e_flux.append(c_eflux.real(row))

#         # ... otherwise append upper limit
#         else:
#             ul_energies.append(c_energy.real(row))
#             ul_flux.append(c_upper.real(row))
#             ul_ed_engs.append(c_ed.real(row))
#             ul_eu_engs.append(c_eu.real(row))

#     # Set upper limit errors
#     yerr = [0.6 * x for x in ul_flux]

#     # Plot the spectrum 
#     plt.figure()
#     plt.xscale( 'log' )
#     plt.yscale( 'log' )
#     plt.grid()
#     plt.errorbar(energies, flux, yerr=e_flux, xerr=[ed_engs, eu_engs],
#                  fmt='ro')
#     plt.errorbar(ul_energies, ul_flux, xerr=[ul_ed_engs, ul_eu_engs],
#                  yerr=yerr, uplims=True, fmt='ro')
#     plt.xlabel('Energy (TeV)')
#     plt.ylabel(r'E$^2$ $\times$ dN/dE (erg cm$^{-2}$ s$^{-1}$)')
#     plt.ylim( 1.e-13 , 1.e-9 )

#     # Show figure
#     if len(plotfile) > 0:
#         plt.savefig(plotfile)
#     else:
#         plt.show()

#     # Return
#     return
















