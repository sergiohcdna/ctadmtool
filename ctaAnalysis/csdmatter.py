#! /usr/bin/env python
# ==========================================================================
# Perform Analysis for DM models
#
# Copyright (C) 2020 Sergio, Judit, Miguel
# Modified 2021 Sergio
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# ==========================================================================
import sys
import gammalib
import ctools
from cscripts import mputils

import os
import math

#====================================================================#
#                                                                    #
#   First, add a path to PFILES environment path to be able to load  #
#   parameters file. This change only survives during the execution  #
#   time of the application.                                         #
#                                                                    #
#====================================================================#

BASEDIR = os.path.abspath( os.path.dirname( __file__ ) )
pfiles  = os.path.join( BASEDIR , 'pfiles' )
os.environ[ 'PFILES' ] += os.pathsep + pfiles

channels = {0:'eL', 1:'eR', 2:'e', 3:'MuL', 4:'MuR',
    5:'Mu', 6:'TauL', 7:'TauR', 8:'Tau', 9:'q',
    10:'c', 11:'b', 12:'t', 13:'WL', 14:'WT',
    15:'W', 16:'ZL', 17:'ZT', 18:'Z', 19:'g',
    20:'Gamma', 21:'h', 22:'Nue', 23:'NuMu', 24:'NuTau',
    25:'Ve', 26:'VMu', 27:'VTau'}

# =============== #
# csdmatter class #
# =============== #

class csdmatter( ctools.csobservation ) :
    """
    This class perform a search for DM signals on a CTA-Observation.
    At this moment, please use for observations without any extra source
    """

    def __init__( self , *argv ) :
        """
        Default Constructor
        """

        #   Initialise application by calling the appropiate class constructor
        self._init_csobservation( self.__class__.__name__ , ctools.__version__ , argv )

        #   Initialise data members
        self._fits        = None
        self._binned_mode = False
        self._onoff_mode  = False
        self._nthreads    = 0
        self._masses      = None

        #   Return
        return

    def __del__( self ) :
        """
        Destructor
        """

        #   Return
        return

    def __getstate__( self ) :
        """
        Extend ctools.csobservation getstate method to include some members
        """

        #   Set pickled dictionary
        state = { 'base'        : ctools.csobservation.__getstate__( self ) ,
                  'fits'        : self._fits ,
                  'binned_mode' : self._binned_mode ,
                  'onoff_mode'  : self._onoff_mode ,
                 # 'dmass'       : self._dmass ,
                 # 'sigmav'      : self._sigmav ,
                  'nthreads'    : self._nthreads ,
                  'masses'      : self._masses }

        #   Return dictionary
        return state

    def __setstate__( self , state ) :
        """
        Extend ctools.csobservation getstate method to include some members
        """

        #   Set dictionary
        ctools.csobservation.__setstate__( self , state[ 'base' ] )

        self._fits        = state[ 'fits' ]
        self._binned_mode = state[ 'binned_mode' ]
        self._onoff_mode  = state[ 'onoff_mode' ]
        # self._dmass       = state[ 'dmass' ]
        # self._sigmav      = state[ 'sigmav' ]
        self._nthreads    = state[ 'nthreads' ]
        self._masses      = state[ 'masses' ]

        #   Return
        return

    def _get_parameters( self ) :
        """
        Get parameters and setup the observation
        """

        #   Set observation if not done before
        if self.obs().is_empty() :
            self._require_inobs( 'csdmatter::get_parameters' )
            self.obs( self._get_observations() )

        #   Set Obs statistic
        self._set_obs_statistic( gammalib.toupper( self[ 'statistic' ].string() ) )

        #   Set Models
        #if self.obs().models().is_empty() :
        #    self.obs().models( self[ 'inmodel' ].filename() )

        #   Query source name
        self[ 'srcname' ].string()

        #   Collect number of unbinned, binned and OnOff obs
        #   in observation container
        n_unbinned = 0
        n_binned   = 0
        n_onoff    = 0

        for obs in self.obs() :
            if obs.classname() == 'GCTAObservation' :
                if obs.eventtype() == 'CountsCube' :
                    n_binned += 1
                else :
                    n_unbinned += 1
            elif obs.classname() == 'GCTAOnOffObservation' :
                n_onoff += 1
        n_cta = n_unbinned + n_binned + n_onoff
        n_other = self.obs().size() - n_cta

        #   Query other parameters
        self[ 'edisp' ].boolean()
        self[ 'calc_ulim' ].boolean()
        self[ 'calc_ts' ].boolean()
        self[ 'fix_bkg' ].boolean()
        self[ 'fix_srcs' ].boolean()

        #   Query all dark-matter related parameters
        self[ 'mmin' ].real()
        self[ 'mmax' ].real()
        self[ 'mnumpoints' ].integer()
        self[ 'process' ].string()
        self[ 'channel' ].string()
        self[ 'ewcorrections' ].boolean()
        self[ 'logsigmav' ].real()
        self[ 'logastfactor' ].real()
        self[ 'redshift' ].real()
        self[ 'eblmodel' ].string()
        self[ 'emin' ].real()
        self[ 'emax' ].real()
        self[ 'modtype' ].string()
        self['dmspecfits'].filename()

        #   Query parameters according to the Source Model type
        if self[ 'modtype' ].string() == 'PointSource' :

            self[ 'ra' ].real()
            self[ 'dec' ].real()

        if self[ 'modtype' ].string() == 'DiffuseSource' :

            self[ 'map_fits' ].filename()


        #   Set mass points
        self._mlogspace()

        #   Read ahead output parameters
        if self._read_ahead() :
            self[ 'outfile' ].filename()


        #   Write into logger
        self._log_parameters( gammalib.TERSE )

        #   Set number of processes for multiprocessing
        self._nthreads = mputils.nthreads( self )

        self._log_header1( gammalib.TERSE , 'DM analysis' )
        self._log_value( gammalib.TERSE , 'Unbinned observations' , n_unbinned )
        self._log_value( gammalib.TERSE , 'Binned observations' , n_binned )
        self._log_value( gammalib.TERSE , 'OnOff Observations' , n_onoff )
        self._log_value( gammalib.TERSE , 'NonCTA Observations' , n_other )

        if n_other == 0 :

            if n_unbinned == 0 and n_binned != 0 and n_onoff == 0 :
                self._binned_mode = True

            elif n_unbinned == 0 and n_binned == 0 and n_onoff != 0 :
                self._onoff_mode = True

            elif n_unbinned == 0 and n_binned != 0 and n_onoff != 0 :
                msg = 'Mixing of binned and OnOff Observations'
                raise RuntimeError( msg )

            elif n_unbinned != 0 and ( n_binned != 0 or n_onoff != 0 ) :
                msg = 'Mixing of different CTA Observations'
                raise RuntimeError( msg )

        else :

            msg = 'csdmatter only supports CTA-observations'
            raise RuntimeError( msg )

        return

    def _mlogspace( self ) :
        """
        Generate vector with points separated logarithmically.

        Return
        ------
        GVector with n mass points, including endpoint
        """

        mmin    = self[ 'mmin' ].real()
        mmax    = self[ 'mmax' ].real()
        mlogmin = math.log10( mmin )
        mlogmax = math.log10( mmax )
        mpoints = self[ 'mnumpoints' ].integer()

        #   initialize an empty list
        # masses  = gammalib.GVector( mpoints )
        masses  = []

        #   Compute width
        width = ( mlogmax - mlogmin ) / ( mpoints - 1 )

        for index in range( mpoints ) :

            # masses[ index ] = 10**( mlogmin + index * width )
            masses.append(math.pow(10., mlogmin + index * width))

        self._masses = masses

        #   Return
        return

    def _get_channel_key(self):

        #   channel from input parameters
        channel = self['channel'].string()
        ckey    = 0

        for key, value in channels.items():
             if channel.upper() == value.upper():
                 ckey = key

        return ckey

    def _gen_model_anna(self, i) :
        """
        in-fly Creation of GModel

        Return
        ------
        GModel for a dark matter annihilation
        """

        #   Create xml element for the spectral type

        #   These are the min and max values that can take
        #   the prefactor parameter
        #   Note: I am not sure why I need to rise so much
        #   the maxval to specify the valid range od the
        #   parameter
        #   I hope this issue change when implementing the
        #   GModelSpectralDMMmodel class
        minval  = 1.0e-40
        maxval  = 1.0e+40

        #   Model type
        modtype  = self['modtype'].string()

        #   recover filename created with _gen_dmfile_anna
        srcname = self['srcname'].string()

        #   Create spectral container using GModelSpectralTable
        #   and setup the mass, channel and normalization of the model
        ch_number = self._get_channel_key()
        dmmass    = self._masses[i]
        jfactor   = math.pow(10., self['logastfactor'].real())
        sigmav    = math.pow(10., self['logsigmav'].real())
        fluxnorm  = sigmav * jfactor / 8. / gammalib.pi / dmmass / dmmass
        fluxnorm *= 1.e-3
        ffile     = self['dmspecfits'].filename()

        dmspec    = gammalib.GModelSpectralTable()
        dmspec.load(ffile)
        dmspec['Mass'].value(dmmass)
        dmspec['Channel'].value(ch_number)
        dmspec['Channel'].scale(1)
        dmspec['Normalization'].value(fluxnorm)
        dmspec['Normalization'].range(minval, maxval)

        #   Mass and Channel parameters should be fixed
        #   But, just to be sure
        dmspec['Mass'].fix()
        dmspec['Channel'].fix()
        dmspec['Normalization'].free()

        #   Creating Spatial container
        if self[ 'modtype' ].string() == 'PointSource' :

            ra     = self['ra'].real()
            dec    = self['dec'].real()
            dmspat = gammalib.GModelSpatialPointSource(ra, dec)

        elif self[ 'modtype' ].string() == 'DiffuseSource' :

            mfile  = self['map_fits'].filename()
            dmspat = gammalib.GModelSpatialDiffuseMap(mfile)

        #   Then generate GModel from spatial and spectral part
        #   This must avoid to create a lot of XML Templates
        #   to specify DM models :P
        dmmodel = gammalib.GModelSky(dmspat, dmspec)
        dmmodel.name(srcname)
        dmmodel.tscalc(True)

        #   Return
        return dmmodel

    # def _gen_bkgmodel(self) :
    #     """
    #     Create bkg model. This model assume that the bkg type
    #     is modeled by the IRF provided by the user.

    #     Return
    #     ------
    #     Gmodel for bkg
    #     """

    #     #   spectral correction
    #     genergy = gammalib.GEnergy(1, 'TeV')
    #     spectral = gammalib.GModelSpectralPlaw(1, 0, genergy)
        
    #     # create background model
    #     bkgmodel = gammalib.GCTAModelIrfBackground(spectral)
    #     bkgmodel.name('Background')
    #     bkgmodel.instruments('CTA')
        
    #     return bkgmodel

    def _fit_mass_point(self, i) :
        """
        Fit Model to DATA in the observation for a specific
        value of dark matter mass

        Return
        ------
        Result , dictionary with relevant fit results
        """

        #   Get value of mass (already in GeV)
        dmmass = self._masses[i]

        self._log_header2(gammalib.EXPLICIT, 'Mass point ' + str(i + 1))

        #   Set reference energy and energy range for calculations
        geref = gammalib.GEnergy(dmmass / 2., 'GeV')
        gemin = gammalib.GEnergy(self['emin'].real(), 'GeV')
        gemax = gammalib.GEnergy(dmmass, 'GeV')

        #   Create file with flux according to process
        #   Well, at this moment, just annihilation :P
        self._log_header1( gammalib.TERSE , 'Compute DM model' )

        if self['process'].string() == 'ANNA' :

            thisdmmodel = self._gen_model_anna(i)

        # elif self['process'].string() == 'DEC' :

        #     print( 'Not Implemented' )
        #     sys.exit()

        #   Then create GModel containers for source and bkg
        # thisbkgmodel = self._gen_bkgmodel()

        #   GModels source and append dm and bkg models
        # mymodels = gammalib.GModels()
        # mymodels.append(thisdmmodel)
        # mymodels.append(thisbkgmodel)

        #   Show mymodels in logfile, just to check that everything is Ok!
        self._log_string(gammalib.EXPLICIT , str(self.obs().models()))

        self._log_header1(gammalib.TERSE, 'Set or replace by Dark matter model')

        for model in self.obs().models() :

            if model.classname() != 'GCTAModelIrfBackground' :

                self.obs().models().remove(model.name())

        self.obs().models().append(thisdmmodel)
        obssim = self.obs().copy()

        #   Now, all the analysis is the same as in csspec script

        #   Get expected dmflux between emin and emax
        #   for the source of interest
        # srcname  = self['srcname'].string()
        # srcmodel = obssim.models()[srcname]
        # srcspec  = srcmodel.spectral()
        theoflux  = thisdmmodel.spectral().flux(gemin, gemax)
        # theoflux *= thisdmmodel.spectral()['Normalization'].value()
        sigmav    = math.pow(10., self['logsigmav'].real())

        #   Header
        self._log_header1(gammalib.TERSE, 'Fitting DM Model')

        #   So, at this moment interesting results to save are:
        #       - Min and Max Energy to compute integrated fluxes
        #       - Mass of dark matter candidate
        #       - Differential flux obtained in the fit
        #       - Error
        #       - TS
        #       - Upper-limit on the differential flux
        #       - Reference value of sigmav
        #       - UL computed of sigmav
        #       - Scale factor computed to obtain the UL on sigmav
        #   This may be change when including Spectral class
        #   for DM annihilation
        result = {'e_min'     : gemin.TeV(),
                  'e_max'     : gemax.TeV(),
                  'mass'      : dmmass,
                  'flux'      : 0.0,
                  'flux_err'  : 0.0,
                  'logL'      : 0.0,
                  'TS'        : 0.0,
                  'ulimit'    : 0.0,
                  'sigma_ref' : sigmav,
                  'sigma_lim' : 0.0 ,
                  'sc_factor' : 0.0}

        #   Header for ctlike instance :)
        self._log_header3(gammalib.EXPLICIT, 'Performing likelihood fit')

        #   Maximum likelihood fit via ctlike
        like             = ctools.ctlike(obssim)
        like['edisp']    = self['edisp'].boolean()
        like['nthreads'] = 1

        #   Chatter
        if self._logVerbose() and self._logDebug() :

            like['debug'] = True

        like.run()

        #   Extract fit results
        model    = like.obs().models()[self['srcname'].string()]
        spectrum = model.spectral()
        logL0    = like.obs().logL()

        result['logL'] = logL0

        #   Write models results
        self._log_string(gammalib.EXPLICIT, str(like.obs().models()))

        #   Continue only if logL0 is different from zero
        if logL0 != 0.0 :

            #   Extract TS value
            result['TS'] = model.ts()


            #   Calculation of upper-limit via ctulimit
            ulimit_value = -1.0

            if self['calc_ulim'].boolean() :

                #   Print to log
                self._log_header3(gammalib.EXPLICIT,
                                   'Computing Upper Limit')

                #   Instance for ctulimit
                ulimit             = ctools.ctulimit(like.obs())
                ulimit['srcname']  = self['srcname'].string()
                ulimit['eref']     = geref.TeV()
                ulimit['emin']     = gemin.TeV()
                ulimit['emax']     = gemax.TeV()

                #   Set chatter
                if self._logVerbose() and self._logDebug() :

                    ulimit['debug'] = True

                ulimit.run()
                ulimit_value = ulimit.diff_ulimit()
                #    Catching exceptions
                # try :

                #     ulimit.run()
                #     ulimit_value = ulimit.diff_ulimit()

                # except :

                #     self._log_string(gammalib.EXPLICIT, 'UL Calculation failed')
                #     ulimit_value = -1.0

                #   Compute quantities related to ulimit

                if ulimit_value > 0.0 :

                    flimit_value          = ulimit.flux_ulimit()
                    result[ 'ulimit' ]    = flimit_value
                    scfactor              = flimit_value / theoflux
                    result[ 'sc_factor' ] = scfactor
                    result[ 'sigma_lim' ] = scfactor * sigmav

                #   Get flux and error
                fitted_flux = spectrum.eval(geref)
                parvalue    = spectrum[0].value()

                if parvalue != 0.0 :

                    rel_error = spectrum[0].error() / parvalue
                    e_flux    = fitted_flux * rel_error

                else :

                    e_flux    = 0.0

                # If a cube, then compute corresponding weight
                if model.spatial().classname() == 'GModelSpatialDiffuseCube' :

                    dir          = gammalib.GSkyDir()
                    model.spatial().set_mc_cone( dir , 180 )
                    norm         = model.spatial().spectrum().eval( geref )
                    fitted_flux *= norm
                    e_flux      *= norm

                #   Convert to nuFnu
                eref2              = geref.MeV() * geref.MeV()
                result['flux']     = fitted_flux * eref2 * gammalib.MeV2erg
                result['flux_err'] = e_flux      * eref2 * gammalib.MeV2erg

                #   Logging
                value = '%e +/- %e' % (fitted_flux, e_flux)
                svmsg = ''

                if self[ 'calc_ulim' ].boolean() and result[ 'ulimit' ] > 0.0 :

                    value += ' [< %e]' % (result['ulimit'])
                    svmsg += ' [%e]' % (result['sc_factor'])

                value += ' 1/cm**2/s'

                if self['calc_ts'].boolean() and result['TS'] > 0.0:

                    value += ' (TS = %.3f)' % (result['TS'])

                self._log_value(gammalib.TERSE, 'Flux', value)

                if len( svmsg ) > 0 :

                    self._log_value(gammalib.TERSE, 'ScaleFactor', svmsg)

        #   If logL0 == 0, then failed :(
        #   but, this does not raise any error
        else :

            value = 'Likelihood is zero. Something is weird. Check model'
            self._log_value(gammalib.TERSE, 'Warning: ', value)

        #   Return
        return result

    def _fit_mass_points(self) :
        """
        Fit for GVector masses

        Return
        ------
        results: dictionary with result for every mass point
        """

        self._log_header1(gammalib.TERSE, 'Fitting models for different masses')
        self._log_string(gammalib.TERSE, str(self._masses))

        # Initialise results
        results = []

        #   I need to check how to use multiprocessing
        #   I hope, this is not caused by a bug in the code
        if self._nthreads > 1 :

            #   force to set nthreads to one

            # self._nthreads = 1
            # Compute for mass points
            args        = [(self, '_fit_mass_point', i)
                            for i in range(len(self._masses))]
            poolresults = mputils.process(self._nthreads, mputils.mpfunc, args)

            # Construct results
            for i in range(len(self._masses)) :

                results.append(poolresults[i][0])
                self._log_string(gammalib.TERSE,poolresults[i][1]['log'],False)

        # Otherwise, loop over energy bins
        else:
            for i in range(len(self._masses)) :

                # Fit energy bin
                result = self._fit_mass_point(i)

                # Append results
                results.append(result)

        # Return results
        return results


    def _create_fits(self, results) :
        """
        Create fits file

        Parameters
        ----------
        result: Dictionary with results obtained from fit
        """

        #   Create columns (><'! Now, added for n mass points')
        nrows = len(self._masses)

        e_min        = gammalib.GFitsTableDoubleCol('MinEnergy', nrows)
        e_max        = gammalib.GFitsTableDoubleCol('MaxEnergy', nrows)
        mass         = gammalib.GFitsTableDoubleCol('Mass', nrows)
        flux         = gammalib.GFitsTableDoubleCol('Flux', nrows)
        flux_err     = gammalib.GFitsTableDoubleCol('EFlux', nrows)
        slogl        = gammalib.GFitsTableDoubleCol('LogL', nrows)
        TSvalues     = gammalib.GFitsTableDoubleCol('TS', nrows)
        ulim_values  = gammalib.GFitsTableDoubleCol('UpperLimit', nrows)
        sigma_lim    = gammalib.GFitsTableDoubleCol('ULCrossSection', nrows)
        sigma_ref    = gammalib.GFitsTableDoubleCol('RefCrossSection', nrows)
        sc_factor    = gammalib.GFitsTableDoubleCol('ScaleFactor', nrows)

        e_min.unit('TeV')
        e_max.unit('TeV')
        mass.unit('TeV')
        flux.unit('erg/cm2/s')
        flux_err.unit('erg/cm2/s')
        ulim_values.unit('1/cm2/s')
        sigma_lim.unit('cm3/s')
        sigma_ref.unit('cm3/s')

        #   Fill fits
        for i, result in enumerate(results) :

            e_min[i]       = result['e_min']
            e_max[i]       = result['e_max']
            mass[i]        = result['mass']
            flux[i]        = result['flux']
            flux_err[i]    = result['flux_err']
            slogl[i]       = result['logL']
            TSvalues[i]    = result['TS']
            ulim_values[i] = result['ulimit']
            sigma_lim[i]   = result['sigma_lim']
            sigma_ref[i]   = result['sigma_ref']
            sc_factor[i]   = result['sc_factor']

        #   Initialise FITS Table with extension "DMATTER"
        table = gammalib.GFitsBinTable(nrows)
        table.extname('DMATTER')

        #   Add Header for compatibility with gammalib.GMWLSpectrum
        table.card('INSTRUME', 'CTA', 'Name of Instrument')
        table.card('TELESCOP', 'CTA', 'Name of Telescope')

        #   Append filled columns to fits table
        table.append(e_min)
        table.append(e_max)
        table.append(mass)
        table.append(flux)
        table.append(flux_err)
        table.append(slogl)
        table.append(TSvalues)
        table.append(ulim_values)
        table.append(sigma_lim)
        table.append(sigma_ref)
        table.append(sc_factor)

        #   Create the FITS file now
        self._fits = gammalib.GFits()
        self._fits.append(table)

        #   Return
        return

    def run(self) :
        """
        Run the script
        """

        #   Screen logging
        if self._logDebug() :
            self._log.cout(True)

        #   Getting parameters
        self._get_parameters()
        srcname = self['srcname'].string()

        #   Write input observation container into logger
        self._log_observations(gammalib.NORMAL, self.obs(), 'Input observation')

        #   Adjust model parameters dependent on input user parameters
        # self._adjust_models()

        #   Fit model
        results = self._fit_mass_points()

        #   Create FITS file
        self._create_fits(results)

        #   Publishing...?
        if self['publish'].boolean() :
            self.publish()

        #   Return
        return

    def save(self) :
        """
        Save dmatter analysis results
        """

        #   Write header
        self._log_header1(gammalib.TERSE, 'Save dmatter analysis results')

        #   Continue only if FITS file is valid
        if self._fits != None :

            #   Get outmap parameter
            outfile = self['outfile'].filename()

            #   Log file name
            self._log_value(gammalib.NORMAL, 'dmatter file', outfile.url())

            #   Save results
            self._fits.saveto(outfile, self['clobber'].boolean())

        #   Return
        return

    def publish(self, name='') :
        """
        Publish results
        """

        #   Write header
        self._log_header1(gammalib.TERSE, 'Publishing Results')

        #   Checking fits
        if self._fits != None :

            if not name :
                user_name = self._name()
            else :
                user_name = name

        #   Log file
        self._log_value(gammalib.NORMAL, 'DM anna name', user_name)

        #   Publishin
        self._fits.publish('DM Anna', user_name)

        #   Return
        return

    def dmatter_fits(self) :
        """
        Return fits with results
        """

        #   Return
        return self._fits

    # def models( self ) :
    #     """
    #     Set models
    #     """

    #     #   Copy models
    #     self.obs().models( models.clone() )

    #     #   Return
    #     return

# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Create instance of application
    app = csdmatter(sys.argv)

    # Execute application
    app.execute()

