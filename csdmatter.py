#! /usr/bin/env python
# ==========================================================================
# Perform Analysis for DM models
#
# Copyright (C) 2020 Sergio, Judit, Miguel
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
                  'nthreads'    : self._nthreads }

        #   Return dictionary
        return state

    def _setstate__( self ) :
        """
        Extend ctools.csobservation getstate method to include some members
        """

        #   Set dictionary
        ctools.csobservation.__setstate__( self , state[ 'base' ] )

        self._binned_mode = state[ 'binned_mode' ]
        self._onoff_mode  = state[ 'onoff_mode' ]
        # self._dmass       = state[ 'dmass' ]
        # self._sigmav      = state[ 'sigmav' ]
        self._nthreads    = state[ 'nthreads' ]

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
        if self.obs().models().is_empty() :
            self.obs().models( self[ 'inmodel' ].filename() )

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

        #   Set energy bounds
        self._set_ebounds()

        #   Query other parameters
        self[ 'edisp' ].boolean()
        self[ 'calc_ulim' ].boolean()
        self[ 'calc_ts' ].boolean()
        self[ 'fix_bkg' ].boolean()
        self[ 'fix_srcs' ].boolean()
        # self[ 'dmass' ].real()
        # self[ 'sigmav' ].real()

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

    def _set_ebounds( self ) :
        """
        Set energy boundaries
        """

        #   Taking energy boundaries from observation
        #   in binned or onoff modes
        if self._binned_mode or self._onoff_mode :

            #   GEbounds instance
            self._ebounds = gammalib.GEbounds()

            #   Energy boundaries
            ebounds = self._create_ebounds()

            #   Extract energy boundaries from first observation
            if self._binned_mode :

                cube_ebounds = self.obs()[ 0 ].events().ebounds()

            else :

                cube_ebounds = self.obs()[ 0 ].rmf().emeasured()

            #   Collect all overlapping, user energy bins
            for ibin in range( ebounds.size() ) :

                #   Extract boundaries
                emin = ebounds.emin( ibin ).TeV() * 0.999
                emax = ebounds.emax( ibin ).TeV() * 1.001

                #   Overlapping bins (counter)
                nbins = 0

                #   Search in the first cube
                emin_value = -1.0
                for k in range( cube_ebounds.size() ) :

                    Emin = cube_ebounds.emin( k ).TeV()
                    Emax = cube_ebounds.emax( k ).TeV()

                    if Emin >= emin and Emax <= emax :

                        emin_value = Emin

                        break

                if emin_value < 0.0 :

                    continue

                for k in range( cube_ebounds.size() ) :

                    Emin = cube_ebounds.emin( k ).TeV()
                    Emax = cube_ebounds.emax( k ).TeV()

                    if Emin >= emin and Emax <= emax :

                        emax_value = Emax
                        nbins     += 1

                #   Append overlapping energy bin
                if nbins > 0 :

                    self._ebounds.append( gammalib.GEnergy( emin_value , 'TeV' ) ,
                                          gammalib.GEnergy( emax_value , 'TeV' ) )

            #   Exception if there are not overlapping bins
            if len( self._ebounds ) == 0 :
                msg = 'Energy range ['+str(cube_ebounds.emin()) + \
                      ', '+str(cube_ebounds.emax())+'] of counts ' + \
                      'cube does not overlap with specified energy ' + \
                      'range [' + \
                      str(ebounds.emin())+', '+str(ebounds.emax())+'].' + \
                      ' Specify overlapping energy range.'
                raise RuntimeError( msg )

        else :

            self._ebounds = self._create_ebounds()

        return

    def _adjust_models( self ) :
        """
        Fixing parameters for bkg and other src in the model
        """

        #   Header
        self._log_header1( gammalib.TERSE , 'Adjust model parameters' )

        #   Loop over models
        for model in self.obs().models() :

            #   set tscalc to false
            model.tscalc( False )

            #   Model name
            self._log_header3( gammalib.EXPLICIT , model.name() )

            #   Deal with the source of interest
            if model.name == self[ 'srcname' ].string() :

                #   At this moment, only check that norm parameter
                #   is free. Later, when including a DMSpectrumclass
                #   then, more operations will be included

                normpar = model.spectral()[ 0 ]
                if normpar.is_fixed() :

                    self._log_string( gammalib.EXPLICIT ,
                                      'Now,' + normpar.name() + ' is free' )

                normpar.free()

                #   Compute TS?
                if self[ 'calc_ts' ].boolean() :

                    model.tscalc( True )

            #   Deal with bkg models
            #   Fixing parameters
            elif self[ 'fix_bkg' ].boolean() and not model.classname() == 'GModelSky' :

                for par in model :

                    if par.is_free() :

                        self._log_string( gammalib.EXPLICIT,
                                          'Now, ' + par.name() + ' is fixed' )

                    par.fix()

            #   Deal with other source
            #   Also, fixing parameters :)
            elif self[ 'fix_srcs' ].boolean() and model.classname() == 'GModelSky' :

                for par in model :

                    if par.is_free() :

                        self._log_string( gammalib.EXPLICIT,
                                          'Now, ' + par.name() + ' is fixed' )

                    par.fix()

        return

    def _fit_model( self ) :
        """
        Fit Model to DATA in the observation

        Return
        ------
        Result , dictionary with relevant fit results
        """

        #   Set reference energy for calculations
        eref = gammalib.GEnergy( self[ 'dmass' ].real() / 2.0 , 'TeV' )

        #   Get expected dmflux at reference energy
        #   for the source of interest
        srcmodel = self.obs().models()[ self[ 'srcname' ].string() ]
        srcspec  = srcmodel.spectral()
        theoflux = srcspec.eval( eref )

        #   Header
        self._log_header1( gammalib.TERSE , 'Fitting DM Model' )

        #   So, at this moment interesting results to save are:
        #       - Reference Energy
        #       - Differential flux obtained in the fit
        #       - Error
        #       - TS
        #       - Upper-limit on the differential flux
        #       - Reference value of sigmav
        #       - UL computed of sigmav
        #       - Scale factor computed to obtain the UL on sigmav
        #   This may be change when including Spectral class
        #   for DM annihilation
        result = { 'energy'    : eref.TeV() ,
                   'flux'      : 0.0 ,
                   'flux_err'  : 0.0 ,
                   'TS'        : 0.0 ,
                   'ulimit'    : 0.0 ,
                   # 'sigma_ref' : self[ 'sigmav' ] ,
                   # 'sigma_lim' : 0.0 ,
                   'sc_factor' : 0.0 }

        #   Header for ctlike instance :)
        self._log_header3( gammalib.EXPLICIT , 'ctlike instance' )

        #   Maximum likelihood fit via ctlike
        like               = ctools.ctlike( self.obs() )
        like[ 'edisp' ]    = self[ 'edisp' ].boolean()
        like[ 'nthreads' ] = 1

        #   Chatter
        if self._logVerbose() and self._logDebug() :

            like[ 'debug' ] = True

        like.run()

        #   Extract fit results
        model    = like.obs().models()[ self[ 'srcname' ].string() ]
        spectrum = model.spectral()
        logL0    = like.obs().logL()

        #   Write models results
        self._log_string( gammalib.EXPLICIT , str( like.obs().models() ) )

        #   Continue only if logL0 is different from zero
        if logL0 != 0.0 :

            #   Extract TS value
            result[ 'TS' ] = model.ts()

            #   Calculation of upper-limit via ctulimit
            ulimit_value = -1.0

            if self[ 'calc_ulim' ].boolean() :

                #   Print to log
                self._log_header3( gammalib.EXPLICIT ,
                                   'Computing Upper Limit' )

                #   Instance for ctulimit
                ulimit              = ctools.ctulimit( like.obs() )
                ulimit[ 'srcname' ] = self[ 'srcname' ].string()
                ulimit[ 'eref' ]    = eref.TeV()

                #   Set chatter
                if self._logVerbose() and self._logDebug() :

                    ulimit[ 'debug' ] = True

                #    Catching exceptions

                try :

                    ulimit.run()
                    ulimit_value = ulimit.diff_ulimit()

                except :

                    self._log_string( gammalib.EXPLICIT , 'UL Calculation failed :(' )
                    ulimit_value = -1.0

                #   Compute quantities related to ulimit
                if ulimit_value > 0.0 :

                    result[ 'ulimit' ]    = ulimit_value * eref.MeV() * \
                                            eref.MeV() * gammalib.MeV2erg
                    result[ 'sc_factor' ] = ulimit_value / theoflux
                    # result[ 'sigma_lim' ] = ulimit_value / theoflux * \
                    #                         self[ 'sigmav' ]

                #   Get flux and error
                fitted_flux = spectrum.eval( eref )
                parvalue    = spectrum[ 0 ].value()

                if parvalue != 0.0 :

                    rel_error = spectrum[ 0 ].error() / parvalue
                    e_flux    = fitted_flux * rel_error

                else :

                    e_flux    = 0.0

                # If a cube, then compute corresponding weight
                if model.spatial().classname() == 'GModelSpatialDiffuseCube' :

                    dir          = gammalib.GSkyDir()
                    model.spatial().set_mc_cone( dir , 180 )
                    norm         = model.spatial().spectrum().eval( eref )
                    fitted_flux *= norm
                    e_flux      *= norm

                #   Convert to nuFnu
                eref2                = eref.MeV() * eref.MeV()
                result[ 'flux' ]     = fitted_flux * eref2 * gammalib.MeV2erg
                result[ 'flux_err' ] = e_flux      * eref2 * gammalib.MeV2erg

                #   Logging
                value = '%e +/- %e' % ( result[ 'flux' ] , result[ 'flux_err' ] )
                svmsg = ''

                if self[ 'calc_ulim' ].boolean() and result[ 'ulimit' ] > 0.0 :

                    value += ' [< %e]' % ( result[ 'ulimit' ] )
                    svmsg += ' [%e]' % ( result[ 'sc_factor' ] )

                value += ' erg/cm**2/s'

                if self[ 'calc_ts' ].boolean() and result[ 'TS' ] > 0.0:

                    value += ' (TS = %.3f)' % ( result[ 'TS' ] )

                self._log_value( gammalib.TERSE , 'Flux' , value )

                if len( svmsg ) > 0 :

                    self._log_value( gammalib.TERSE , 'ScaleFactor' , svmsg )

        #   If logL0 == 0, then failed :(
        #   but, this does not raise any error
        else :

            value = 'Likelihood is zero. Something is weird. Check model'
            self._log_value( gammalib.TERSE , 'Warning: ' , value )

        #   Return
        #   At this moment, only save for individual mass and channel
        #   Later, include loop to compute over several masses and channels
        return result

    def _create_fits( self , result ) :
        """
        Create fits file

        Parameters
        ----------
        result: Dictionary with results obtained from fit
        """

        #   Create columns (just one row ><'!)
        #   This will change when adding loop over different masses
        #   and/or channels
        #   By now, hard-coding the number of rows, 1
        nrows = 1
        energy       = gammalib.GFitsTableDoubleCol( 'RefEnergy' , nrows )
        flux         = gammalib.GFitsTableDoubleCol( 'Flux' , nrows )
        flux_err     = gammalib.GFitsTableDoubleCol( 'EFlux', nrows )
        TSvalues     = gammalib.GFitsTableDoubleCol( 'TS' , nrows )
        ulim_values  = gammalib.GFitsTableDoubleCol( 'UpperLimit' , nrows )
        # sigma_lim    = gammalib.GFitsTableDoubleCol( 'ULCrossSection' , nrows )
        # sigma_ref    = gammalib.GFitsTableDoubleCol( 'RefCrossSection' , nrows )
        sc_factor    = gammalib.GFitsTableDoubleCol( 'ScaleFactor' , nrows )

        energy.unit( 'TeV' )
        flux.unit( 'erg/cm2/s' )
        flux_err.unit( 'erg/cm2/s' )
        ulim_values.unit( 'erg/cm2/s' )
        # sigma_lim.unit( 'cm3/s' )
        # sigma_ref.unit( 'cm3/s' )

        #   Fill fits
        energy[ 0 ]      = result[ 'energy' ]
        flux[ 0 ]        = result[ 'flux' ]
        flux_err[ 0 ]    = result[ 'flux_err' ]
        TSvalues[ 0 ]    = result[ 'TS' ]
        ulim_values[ 0 ] = result[ 'ulimit' ]
        # sigma_lim[ 0 ]   = result[ 'sigma_lim' ]
        # sigma_ref[ 0 ]   = result[ 'sigma_ref' ]
        sc_factor[ 0 ]   = result[ 'sc_factor' ]

        #   Initialise FITS Table with extension "DMATTER"
        table = gammalib.GFitsBinTable( nrows )
        table.extname( 'DMATTER' )

        #   Add Header for compatibility with gammalib.GMWLSpectrum
        table.card( 'INSTRUME' , 'CTA' , 'Name of Instrument' )
        table.card( 'TELESCOP' , 'CTA' , 'Name of Telescope' )

        #   Append filled columns to fits table
        table.append( energy )
        table.append( flux )
        table.append( flux_err )
        table.append( TSvalues )
        table.append( ulim_values )
        # table.append( sigma_lim )
        # table.append( sigma_ref )
        table.append( sc_factor )

        #   Create the FITS file now
        self._fits = gammalib.GFits()
        self._fits.append( table )

        #   Return
        return

    def run( self ) :
        """
        Run the script
        """

        #   Screen logging
        if self._logDebug() :
            self._log.cout( True )

        #   Getting parameters
        self._get_parameters()

        #   Write input observation container into logger
        self._log_observations( gammalib.NORMAL , self.obs() , 'Input observation' )

        #   Adjust model parameters dependent on input user parameters
        self._adjust_models()

        #   Fit model
        results = self._fit_model()

        #   Create FITS file
        self._create_fits( results )

        #   Publishing...?
        if self[ 'publish' ].boolean() :
            self.publish()

        #   Return
        return

    def save(self) :
        """
        Save dmatter analysis results
        """

        #   Write header
        self._log_header1( gammalib.TERSE , 'Save dmatter analysis results' )

        #   Continue only if FITS file is valid
        if self._fits != None :

            #   Get outmap parameter
            outfile = self[ 'outfile' ].filename()

            #   Log file name
            self._log_value( gammalib.NORMAL, 'dmatter file' , outfile.url() )

            #   Save results
            self._fits.saveto( outfile , self['clobber'].boolean() )

        #   Return
        return

    def publish( self , name='' ) :
        """
        Publish results
        """

        #   Write header
        self._log_header1( gammalib.TERSE , 'Publishing Results' )

        #   Checking fits
        if self._fits != None :

            if not name :
                user_name = self._name()
            else :
                user_name = name

        #   Log file
        self._log_value( gammalib.NORMAL , 'DM anna name' , user_name )

        #   Publishin
        self._fits.publish( 'DM Anna' , user_name )

        #   Return
        return

    def dmatter_fits( self ) :
        """
        Return fits with results
        """

        #   Return
        return self._fits

    def models( self ) :
        """
        Set models
        """

        #   Copy models
        self.obs().models( models.clone() )

        #   Return
        return

# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Create instance of application
    app = csdmatter(sys.argv)

    # Execute application
    app.execute()

