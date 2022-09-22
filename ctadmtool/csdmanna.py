#! /usr/bin/env python
# ==========================================================================
# Perform Analysis for DM models: Annihilation
#
# Copyright (C) 2022 Sergio, Judit, Miguel
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
import gammalib
import ctools
from cscripts import mputils

from ctadmtool.dmspectrum.dmspectra import dmspectrum

import math
import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import minimize,root_scalar,fsolve

import sys
import os
#====================================================================#
#                                                                    #
#   First, add a path to PFILES environment path to be able to load  #
#   parameters file. This change only survives during the execution  #
#   time of the application.                                         #
#                                                                    #
#====================================================================#

BASEDIR = os.path.abspath(os.path.dirname(__file__))
pfiles  = os.path.join(BASEDIR, 'pfiles')
os.environ['PFILES'] += os.pathsep + pfiles

# This is the dictionary of channels available in PPPC4DMID tables
# for the data with electroweak corrections
channels = {0:'eL', 1:'eR', 2:'e', 3:'MuL', 4:'MuR',
    5:'Mu', 6:'TauL', 7:'TauR', 8:'Tau', 9:'q',
    10:'c', 11:'b', 12:'t', 13:'WL', 14:'WT',
    15:'W', 16:'ZL', 17:'ZT', 18:'Z', 19:'g',
    20:'Gamma', 21:'h', 22:'Nue', 23:'NuMu', 24:'NuTau',
    25:'Ve', 26:'VMu', 27:'VTau'}

# This is the dictionary of channels available in PPPC4DMID tables
# for the data without electroweak corrections
channelsNoEw = {0:'e', 1:'Mu', 2:'Tau', 3:'q',
    4:'c', 5:'b', 6:'t', 7:'W', 8:'Z', 9:'g'}

# ================ #
#  csdmanna class  #
# ================ #

class csdmanna(ctools.csobservation):
    """
        This class performs a search for gamma-ray signals
        induced by annihilation of DM particles in a target
        The analysis use the TS (or profile LLR) to set
        an upper-limit to annihilation cross-section
        if there is no detection of the DM signal.
        Extra sources (or components to the total emission)
        in the region of interest are passed via the inmodel
        parameter
    """

    # Initialize class
    def __init__(self,*argv):
        """
        Default Constructor
        """

        #   Initialise application by calling the appropiate class constructor
        self._init_csobservation(self.__class__.__name__,
            ctools.__version__,argv)

        #   Initialise data members
        self._fits        = None
        self._binned_mode = False
        self._onoff_mode  = False
        self._nthreads    = 0
        self._masses      = []
        self._jfactor     = 0.0
        self._sigmavref   = 0.0
        self._sigmavUL    = 0.0
        self._deltats     = []

        #   Return
        return

    def __del__(self):
        """
        Destructor
        """

        #   Return
        return

    def __getstate__(self):
        """
        Extend ctools.csobservation getstate method to include some members
        """

        #   Set pickled dictionary
        state = {'base'       : ctools.csobservation.__getstate__(self),
                'fits'        : self._fits,
                'binned_mode' : self._binned_mode,
                'onoff_mode'  : self._onoff_mode,
                'nthreads'    : self._nthreads,
                'masses'      : self._masses,
                'jfactor'     : self._jfactor,
                'sigmavref'   : self._sigmavref,
                'sigmavUL'    : self._sigmavUL,
                'deltats'     : self._deltats}

        #   Return dictionary
        return state

    def __setstate__(self, state):
        """
        Extend ctools.csobservation getstate method to include some members
        """

        #   Set dictionary
        ctools.csobservation.__setstate__(self, state['base'])

        self._fits        = state['fits']
        self._binned_mode = state['binned_mode']
        self._onoff_mode  = state['onoff_mode']
        self._nthreads    = state['nthreads']
        self._masses      = state['masses']
        self._jfactor     = state['jfactor']
        self._sigmavref   = state['sigmavref']
        self._sigmavUL    = state['sigmavUL']
        self._deltats     = state['deltats']

        #   Return
        return

    def _get_parameters(self):
        """
        Get parameters and setup the observation
        """

        #   Set observation if not done before
        if self.obs().is_empty():
            self._require_inobs('csdmatter::get_parameters')
            self.obs(self._get_observations())

        #   Set Obs statistic
        self._set_obs_statistic(gammalib.toupper(self['statistic'].string()))

        #   Set models if there are not in the observation
        #   This is the usual behavior for cscripts
        #   The DM is still computed during execution time,
        #   but only appended to the models container
        #   passed by the user
        if self.obs().models().is_empty():
            self.obs().models(self['inmodel'].filename())

        #   Query source name
        self['srcname'].string()

        #   Collect number of unbinned, binned and OnOff obs
        #   in observation container
        n_unbinned = 0
        n_binned   = 0
        n_onoff    = 0

        for obs in self.obs():
            if obs.classname() == 'GCTAObservation':
                if obs.eventtype() == 'CountsCube':
                    n_binned += 1
                else :
                    n_unbinned += 1
            elif obs.classname() == 'GCTAOnOffObservation':
                n_onoff += 1
        n_cta   = n_unbinned + n_binned + n_onoff
        n_other = self.obs().size() - n_cta

        #   Query other parameters
        self['edisp'].boolean()
        self['calc_ulim'].boolean()
        # self['fix_bkg'].boolean()
        # self['fix_srcs'].boolean()

        #   Query all dark-matter related parameters
        #   This class is specific for Annihilation of DM
        self['mmin'].real()
        self['mmax'].real()
        self['mnumpoints'].integer()
        self['channel'].string()
        self['ewcorrections'].boolean()
        self['eblmodel'].string()
        self['redshift'].real()
        self['emin'].real()
        self['emax'].real()
        self['modtype'].string()
        self['logsigmav'].real()
        self['logjfactor'].real()
        self['interp_method'].string()
        self['npoints'].integer()

        #   Save values of jfactor and sigmav to future use
        self._jfactor   = math.pow(10.,self['logjfactor'].real())
        self._sigmavref = math.pow(10.,self['logsigmav'].real())

        #   Query parameters according to the Source Model type
        if self['modtype'].string() == 'PointSource':
            self['ra'].real()
            self['dec'].real()

        if self['modtype'].string() == 'DiffuseSource':
            self['map_fits'].filename()

        #   Set mass points
        self._masses = self._mlogspace()


        #   Read ahead output parameters
        if self._read_ahead() :
            self['outfile'].filename()


        #   Write into logger
        self._log_parameters(gammalib.TERSE)

        #   Set number of processes for multiprocessing
        self._nthreads = mputils.nthreads(self)

        self._log_header1(gammalib.TERSE,'DM analysis')
        self._log_value(  gammalib.TERSE,'Unbinned observations',n_unbinned)
        self._log_value(  gammalib.TERSE,'Binned observations'  ,n_binned)
        self._log_value(  gammalib.TERSE,'OnOff Observations'   ,n_onoff)
        self._log_value(  gammalib.TERSE,'NonCTA Observations'  ,n_other)

        if n_other == 0:
            if n_unbinned == 0 and n_binned != 0 and n_onoff == 0:
                self._binned_mode = True
            elif n_unbinned == 0 and n_binned == 0 and n_onoff != 0:
                self._onoff_mode = True
            elif n_unbinned == 0 and n_binned != 0 and n_onoff != 0:
                msg = 'Mixing of binned and OnOff Observations'
                raise RuntimeError(msg)
            elif n_unbinned != 0 and (n_binned != 0 or n_onoff != 0):
                msg = 'Mixing of different CTA Observations'
                raise RuntimeError(msg)
        else:
            msg = 'csdmatter only supports CTA-observations'
            raise RuntimeError(msg)

        return

    def _mlogspace(self):
        """
        Generate list with mass points separated logarithmically.

        Return
        ------
        List with n mass points, including endpoint
        """

        mmin    = self['mmin'].real()
        mmax    = self['mmax'].real()
        mlogmin = np.log10(mmin)
        mlogmax = np.log10(mmax)
        mpoints = self['mnumpoints'].integer()

        masses  = np.logspace(mlogmin,mlogmax,mpoints)

        thisarray = np.zeros(masses.size)
        for i, mass in enumerate(masses):
            m            = float(math.ceil(mass))
            thisarray[i] = m

        #   Return
        return thisarray.tolist()

    def _gen_dmmodel(self,i):
        """
        Creation of DM model for annihilation during execution time

        Return
        ------
        GModel for a dark matter annihilation
        The GModel container used is the GModelSpectralFunc.
        This is the class used to read spectra from File,
        but I do not create any file (either temporal or another type)
        The normalization is computed using (anna):
            N_0 = J * sigma_v / (8*pi*mass**2)
        """

        #   These are the min and max values that can take
        #   the prefactor parameter
        #   Note: I am not sure why I need to rise so much
        #   the maxval to specify the valid range of the
        #   parameter. Probably because the likelihood profile
        #   is extremely flat
        #   NOTE: In fact, the Profile Likelihood is flat
        #   but only in the range of small values for normalization
        self._log_header3(gammalib.EXPLICIT, 'New DM-spectral model')
        minval  = 0.0
        maxval  = 1.0e+20

        #   Number of energy points used to compute the gamma-ray flux
        epoints = 200

        #   Model type
        modtype = self['modtype'].string()

        #   recover filename created with _gen_dmfile_anna
        srcname = self['srcname'].string()

        #   Create spectral container using GModelSpectralTable
        #   and setup the mass, channel and normalization of the model
        channel   = self['channel'].string()
        dmmass    = self._masses[i]
        jfactor   = self._jfactor
        sigmav    = self._sigmavref
        fluxnorm  = sigmav*jfactor/(8.*gammalib.pi*dmmass*dmmass)
        fluxnorm *= 1.e-3
        #scale     = 1.0
        #fluxnorm /= scale

        msg = 'New range around mass {0}'.format(dmmass)
        self._log_header3(gammalib.EXPLICIT, msg)

        #=================================================================
        #   The photon spectra is computed using th PPPC4DMID tables
        #   The energy range is from the min energy indicated by the
        #   user in the parameters. Maximum energy is computed
        #   according to the mass of the DM candidate. This is to
        #   avoid a lot of rows with zeros, that can lead to bad
        #   calculations of the integral flux
        #=================================================================

        #   Create interpolation function for DM spectra
        emin     = self['emin'].real()
        thisemax = 0.9*dmmass
        hasew    = self['ewcorrections'].boolean()
        eblmodel = self['eblmodel'].string()
        redshift = self['redshift'].real()

        dminterp = dmspectrum(dmmass,emin,thisemax,channel,
            redshift,process='anna',eblmod=eblmodel,
            has_EW=hasew,epoints=epoints)

        #   Get spectra induced by DM, and the energy values
        dnde     = dminterp.spectra()
        dnde    *= fluxnorm
        energies = dminterp.energy

        #   Create ModelSpectra container
        dmspec = gammalib.GModelSpectralFunc()
        for i,energy in enumerate(energies):
            eng    = gammalib.GEnergy(energy,'GeV')
            dphide = dnde[i]
            dmspec.append(eng,dphide)

        #   Tunning the DM-spectral component
        dmspec['Normalization'].scale(1.0)
        dmspec['Normalization'].factor_value(1)
        dmspec['Normalization'].factor_min(0)
        dmspec['Normalization'].factor_max(1.e+20)

        #   Create the spatial component of the emission model
        if self['modtype'].string() == 'PointSource':
            ra     = self['ra'].real()
            dec    = self['dec'].real()
            dmspat = gammalib.GModelSpatialPointSource(ra,dec)
        elif self['modtype'].string() == 'DiffuseSource':
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

    def _get_fittedpars(self,gmodels):
        """
        Get Value parameters from GModel Container.
        Save only the free parameters in the container

        Return
        ------
        Dictionary with values and errors for parameters
        """

        result = {}
        for gmodel in gmodels:
            gname = gmodel.name()

            if gmodel.classname() == 'GModelSky':
                spectral   = gmodel.spectral()
                spatial    = gmodel.spatial()
                temporal   = gmodel.temporal()
                components = [spectral,spatial,temporal]
                cnames     = ['Spec','Spat','Temp']
                for i,component in enumerate(components):
                    for par in component:
                        if par.is_free():
                            entry = '-'.join([gname,cnames[i],par.name()])
                            result[entry] = [par.value(),par.error()]
            else:
                spectral = gmodel.spectral()
                for par in spectral:
                    if par.is_free():
                        entry = '-'.join([gname,'Spec',par.name()])
                        result[entry] = [par.value(),par.error()]

        return result

    def _get_profilets(self,like):
        """
        Compute TS profile for a given Likelihood instance

        Return
        ------
        Array with TS values
        """
        src = self['srcname'].string()

        #   Make a deep copy of the observation in the ctlike container
        obs = like.obs().copy()

        #   Get value of the fitted normalization
        fnorm = like.obs().models()[src].spectral()['Normalization'].factor_value()

        #   Get number of points used to create the spectra
        points = self['npoints'].integer()

        #   Create array for normalization values
        if fnorm == 0.0:
            norms = np.logspace(0,5,points)
            norms = np.insert(norms,0,0.0)
        else:
            logfnorm = np.log10(fnorm)
            norms    = np.logspace(logfnorm-4,logfnorm+4,points+1)

        #   Fix the parameters of the DM Model
        for par in obs.models()[src].spectral():
            if par.is_free():
                par.fix()

        #   Get TS Values
        tsvals = np.zeros(points+1)
        for i,norm in enumerate(norms):
            obs.models()[src].spectral()['Normalization'].factor_value(norm)
            fit = ctools.ctlike(obs)
            fit.process()

            tsval = fit.obs().models()[src].ts()
            tsvals[i] = tsval

        return norms,tsvals

    def _get_parbracket(self,norms,tsvals,deltats):
        """
        Get Initial Parameter bracket to start root finding

        Return
        ------
        Bracket
        """
        bracket = []
        y       = -tsvals - deltats

        for i in range(1,len(y)):
            dot = y[i-1]*y[i]
            if dot < 0:
                bracket.append(norms[i-1])
                bracket.append(norms[i])
                break;

        return bracket

    def _fit_mass_point(self,i):
        """
        Fit Model to DATA in the observation for a specific
        value of dark matter mass

        Return
        ------
        Result, dictionary with relevant fit results
        """

        #   Get value of mass (already in GeV)
        dmmass  = self._masses[i]
        srcname = self['srcname'].string()

        header = 'Fit the DM model for annihilation'
        self._log_header1(gammalib.TERSE, header)
        msg = 'Mass point {0}: {1}'.format(i+1,dmmass)
        self._log_header2(gammalib.EXPLICIT, msg)

        # Get the Sky DM-model for the target
        self._log_header3(gammalib.EXPLICIT,'Create the DM model')
        msg = '{:.3e} (cm**3/s)'.format(self._sigmavref)
        self._log_value(gammalib.TERSE,'Ref. Cross-section',msg)
        msg = '{:.3e} (GeV**2/cm**5/s)'.format(self._jfactor)
        self._log_value(gammalib.TERSE,'J factor',msg)

        thisdmmodel = self._gen_dmmodel(i)

        # Set reference energy and energy range for calculations
        # according to the DM process
        # the 0.9 factor is to avoid problems at the end of the spectrum
        thiseref = dmmass/2.
        thisemax = 0.9*dmmass
        thisemin = self['emin'].real()

        geref = gammalib.GEnergy(thiseref,'GeV')
        gemin = gammalib.GEnergy(thisemin,'GeV')
        gemax = gammalib.GEnergy(thisemax,'GeV')

        #   Get expected dmflux between emin and emax
        #   for the source of interest
        theoflux  = thisdmmodel.spectral().flux(gemin,gemax)
        if theoflux == 0.0 :
            theoflux = -1.0
        msg = '{:.3e} (ph/cm**2/s)'.format(theoflux)
        self._log_value(gammalib.TERSE,'Photon Flux',msg)

        #   Making a deep copy of the observation
        obssim = self.obs().copy()

        self._log_header2(gammalib.TERSE,'Insert the DM model')
        obssim.models().insert(0,thisdmmodel)

        #   Also activate TS calculation for all the GModelSky containers
        for model in obssim.models():
            if model.classname() == 'GModelSky':
                model.tscalc(True)

        #   Show mymodels in logfile, just to check that everything is Ok!
        self._log_header2(gammalib.TERSE,'Emission Models')
        self._log_string(gammalib.EXPLICIT,str(obssim.models()))

        #   Header
        self._log_header2(gammalib.TERSE,'Fitting DM Model')

        #   So, at this moment interesting results to save are:
        #       - Min and Max Energy to compute integrated fluxes
        #       - Mass of dark matter candidate
        #       - Differential flux obtained in the fit
        #       - Error
        #       - logL
        #       - TS
        #       - Upper-limit on the differential flux
        #       - Scale factor computed to obtain the UL on sigmav
        #   And according to the process, either:
        #       - Reference value of sigmav
        #       - UL computed of sigmav
        #   or:
        #       - Reference value of lifetime
        #       - UL computed of lifetime
        result = {'e_min'      : gemin.TeV(),
                  'e_max'      : gemax.TeV(),
                  'mass'       : dmmass,
                  'flux'       : 0.0,
                  'flux_err'   : 0.0,
                  'e2flux'     : 0.0,
                  'e2flux_err' : 0.0,
                  'logL'       : 0.0,
                  'TS'         : 0.0,
                  'ulimit'     : 0.0,
                  'sc_factor'  : 0.0,
                  'sigmav_ref' : self._sigmavref,
                  'sigmav_lim' : 0.0}

        #   Header for ctlike instance :)
        self._log_header3(gammalib.EXPLICIT,'Performing likelihood fit')

        #   Maximum likelihood fit via ctlike
        #   I increased some values for the accuracy
        #   and the acceptance decrease
        #   This is because, for some observations, the
        #   minimization and search of UL via Profile LL
        #   requires to find change until the 7th or greater
        #   digit. This also increase the number of iterations
        #   and the time spend to find a solution
        like                  = ctools.ctlike(obssim)
        like['edisp']         = self['edisp'].boolean()
        like['like_accuracy'] = 1.e-10
        like['accept_dec']    = 10
        like['statistic']     = self['statistic'].string()
        like['nthreads']      = 1

        #   Chatter
        if self._logVerbose() and self._logDebug():
            like['debug'] = True

        like.process()

        #   Extract best fit model
        fitted_model   = like.obs().models()[srcname]
        fitted_spectra = fitted_model.spectral()
        logL0          = like.obs().logL()
        fittednorm     = fitted_spectra['Normalization'].factor_value()

        # Fill some results
        result['logL'] = logL0
        par_results    = self._get_fittedpars(like.obs().models())
        for k,v in par_results.items():
            result[k] = v

        #   Write models results
        self._log_header2(gammalib.TERSE,'Results from fitting: Optimization')
        self._log_string( gammalib.EXPLICIT,str(like.opt()))

        self._log_header2(gammalib.TERSE,'Results from fitting')
        self._log_string( gammalib.EXPLICIT,str(like.obs().models()))

        #   Continue only if logL0 is different from zero
        if logL0 != 0.0:
            #   Extract TS of the best fit
            ts           = fitted_model.ts()
            result['TS'] = ts

            #   Get Covariance Matrix
            covariance = like.obs().function().covariance()

            #   Convert Covariance Matrix to numpy array
            nrows  = covariance.rows()
            ncols  = covariance.columns()
            np_cov = np.zeros((nrows,ncols))

            for i in range(ncols):
                for j in range(nrows):
                    val = covariance[j,i]
                    if val != 0:
                        np_cov[j,i] = val

            self._log_header3(gammalib.TERSE,'Compute Correlation Matrix')
            #   Compute Correlation Matrix
            #   Remove all zeros from covariance matrix
            np_cov = np_cov[:,~np.all(np_cov == 0,axis=0)]
            np_cov = np_cov[  ~np.all(np_cov == 0,axis=1)]

            #   Create diagonal matrix
            d = np.sqrt(np.diag(np_cov))

            dmat = np.zeros(np_cov.shape)
            np.fill_diagonal(dmat,d)

            # Get inverse matrix
            dmatinv = np.linalg.inv(dmat)

            # Get Correlation matrix
            rho = dmatinv@np_cov@dmatinv.T

            result['correlation'] = rho

            parnames = [k for k in par_results.keys()]
            result['parnames'] = parnames

            self._log_value(gammalib.TERSE,'TS',ts)
            self._log_value(gammalib.TERSE,'Fitted Norm',fittednorm)
            self._log_string(gammalib.EXPLICIT,str(rho))

            if self['calc_ulim'].boolean():
                self._log_header2(gammalib.TERSE,'UL Calculation')

                #   I start with the calculation of Profile TS
                #   This is because for some cases
                #   the UL calculation using LLR requires to
                #   find increments only visible in the 7th
                #   or greater digit, reaching the precision
                #   in the computer

                msg = 'Calculation of Profile TS'
                self._log_header3(gammalib.TERSE,msg)
                norms,tsvals    = self._get_profilets(like)
                result['dts']   = tsvals.tolist()
                result['norms'] = norms.tolist()
                self._deltats.append(tsvals.tolist())

                # self._log_header3(gammalib.TERSE,'TS Scan')
                # self._log_string(gammalib.EXPLICIT,str(tsvals.tolist()))
                #   Now, start calculation of UL
                #   Because, we want to compute the 95% CL
                #   for the UL to the flux, we search
                #   the value of normalization (factor_value)
                #   that give us an increment of 3.841
                #   The values for the normalization are hard coded
                #   So,I need to compute again, hahaha
                #   Also, because the normalization parameter
                #   is bounded to be always positive
                #   Some times, the best fit parameter
                #   is found to be zero (correctly, should be negative)
                #   but, then the calculation in this case is slightly
                #   different.
                points = self['npoints'].integer()

                #   Now, I get the interpolation function
                kind    = self['interp_method'].string()
                ulnorm  = 0.0

                if fittednorm != 0.0:
                    thisf = interp1d(np.log10(norms),-tsvals,kind=kind)

                    #   Get the minimum for -TS profile
                    rmin = minimize(thisf,np.log10(fittednorm),method='Nelder-Mead',tol=1.e-6)
                    self._log_value(gammalib.TERSE,'Minimum found: Norm',np.power(10,rmin.x[0]))

                    #   The following are the TS at the minimum
                    #   and the value of the TS after the increment
                    #   to set a 95% CL Upper limit to the flux
                    tsmin = thisf(rmin.x[0])
                    tsnew = tsmin + 3.841

                    self._log_value(gammalib.TERSE,'TS (at min)',tsmin)
                    #   But, I need to redefine the interp1d function
                    #   to search the value of norm where the increase occurs
                    modf = interp1d(np.log10(norms),-tsvals-tsnew,kind=kind)

                    #   Get initial bracket to start root finding
                    #   Now, I will start from the minimum and increase
                    #   in two orders of magnitud (only 2 in log(norms)
                    bracket = self._get_parbracket(np.log10(norms),tsvals,tsnew)
                    self._log_value(gammalib.TERSE,
                        'Initial Parameter Range',str(bracket))

                    #   Root finding using the brent Method
                    sol  = root_scalar(modf,bracket=bracket,method='brentq')
                    # sol  = fsolve(lambda x:thisf(x)-tsnew,)
                    ulnorm = np.power(10,sol.root)
                    #root = fsolve(lambda x:thisf(x)-tsnew,3,xtol=1.e-6)

                else:
                    thisf = interp1d(norms,-tsvals,kind=kind)

                    #   Evaluated at zero, cuz the min is there
                    tsmin = thisf(0)

                    tsnew   = tsmin+3.841
                    modf    = interp1d(norms,-tsvals-tsnew,kind=kind)
                    bracket = self._get_parbracket(norms,tsvals,tsnew)
                    sol     = root_scalar(modf,bracket=bracket,method='brentq')
                    ulnorm  = sol.root

                #   Compute UL to cross-section
                #   and other related quantities
                dmmodel   = like.obs().models()[srcname].copy()
                dmspectra = dmmodel.spectral()

                dmspectra['Normalization'].factor_value(ulnorm)
                fluxUL               = dmspectra.flux(gemin,gemax)
                scale                = fluxUL/theoflux
                self._sigmavUL       = scale*self._sigmavref
                result['sigmav_lim'] = self._sigmavUL
                result['sc_factor']  = scale
                result['ulimit']     = fluxUL

                #   Get diff-flux and error
                fitted_flux = fitted_spectra.eval(geref)
                parvalue    = fitted_spectra['Normalization'].value()

                if parvalue != 0.0:
                    rel_error = fitted_spectra['Normalization'].error()/parvalue
                    e_flux    = fitted_flux*rel_error
                else:
                    e_flux    = 0.0

                #   Save fitted flux
                result['flux']     = fitted_flux/gammalib.MeV2erg
                result['flux_err'] = e_flux     /gammalib.MeV2erg

                #   Convert to nuFnu
                eref2                = geref.MeV()*geref.MeV()
                result['e2flux']     = fitted_flux*eref2*gammalib.MeV2erg
                result['e2flux_err'] = e_flux     *eref2*gammalib.MeV2erg

                #   Logging
                value = '%e +/- %e' % (fitted_flux,e_flux)
                svmsg = ''

                if self['calc_ulim'].boolean() and result['ulimit'] > 0.0:
                    value += ' [< %e]' % (result['ulimit'])
                    svmsg += ' [%e]'   % (result['sc_factor'])
                value += ' 1/cm**2/s'
                if result['TS'] > 0.0:
                    value += ' (TS = %.3f)' % (result['TS'])
                self._log_value(gammalib.TERSE, 'Flux', value)
                if len(svmsg) > 0:
                    self._log_value(gammalib.TERSE,'ScaleFactor',svmsg)
                    sigmalim = result['sigmav_lim']
                    msg  = 'UL on anna cross-section is: %.3e' % (sigmalim)
                    msg += ' cm**3/s'
                    self._log_value(gammalib.TERSE,'Cross-Section',msg)

        #   If logL0 == 0, then failed :(
        #   but, this does not raise any error
        else:
            value = 'Likelihood is zero. Something is weird. Check model'
            self._log_value(gammalib.TERSE,'Warning: ',value)

        #   Return
        return result

    def _fit_mass_points(self):
        """
        Fit for list of masses

        Return
        ------
        results: dictionary with result for every mass point
        """
        self._log_header1(gammalib.TERSE,'Fitting models for different masses')
        self._log_string(gammalib.TERSE,str(self._masses))

        # Initialise results
        results = []

        #   Now, running multiprocessing
        if self._nthreads > 1:
            # Compute for mass points
            args=[(self,'_fit_mass_point',i) for i in range(len(self._masses))]
            poolresults = mputils.process(self._nthreads,mputils.mpfunc,args)

            # Construct results
            for i in range(len(self._masses)):
                results.append(poolresults[i][0])
                self._log_string(gammalib.TERSE,poolresults[i][1]['log'],False)

        # Otherwise, loop over energy bins
        else:
            for i in range(len(self._masses)):
                # Fit energy bin
                result = self._fit_mass_point(i)

                # Append results
                results.append(result)

        # Return results
        return results

    def _create_fits(self,results):
        """
        Create fits file

        Parameters
        ----------
        result: Dictionary with results obtained from fit
        """

        #   Create columns (><'! Now, added for n mass points)
        nrows = len(self._masses)
        ncols = len(results[0]['dts'])

        e_min        = gammalib.GFitsTableDoubleCol('MinEnergy'    ,nrows)
        e_max        = gammalib.GFitsTableDoubleCol('MaxEnergy'    ,nrows)
        mass         = gammalib.GFitsTableDoubleCol('Mass'         ,nrows)
        flux         = gammalib.GFitsTableDoubleCol('Flux'         ,nrows)
        flux_err     = gammalib.GFitsTableDoubleCol('ErrFlux'      ,nrows)
        e2flux       = gammalib.GFitsTableDoubleCol('E2Flux'       ,nrows)
        e2flux_err   = gammalib.GFitsTableDoubleCol('E2ErrFlux'    ,nrows)
        slogl        = gammalib.GFitsTableDoubleCol('LogL'         ,nrows)
        TSvalues     = gammalib.GFitsTableDoubleCol('TS'           ,nrows)
        ulim_values  = gammalib.GFitsTableDoubleCol('UpperLimit'   ,nrows)
        sc_factor    = gammalib.GFitsTableDoubleCol('ScaleFactor'  ,nrows)
        tsscan       = gammalib.GFitsTableDoubleCol('Delta TS'     ,nrows,ncols)
        norms        = gammalib.GFitsTableDoubleCol('log10NormScan',nrows,ncols)

        #   Add units for the relevant columns
        e_min.unit('TeV')
        e_max.unit('TeV')
        mass.unit('TeV')
        flux.unit('1/erg/cm2/s')
        flux_err.unit('1/erg/cm2/s')
        e2flux.unit('erg/cm2/s')
        e2flux_err.unit('erg/cm2/s')
        ulim_values.unit('1/cm2/s')

        #   Creating columns for parameter of interest
        #   according to the process:
        #   Anna  --> CrossSection
        paroi_lim = gammalib.GFitsTableDoubleCol('ULCrossSection' ,nrows)
        paroi_ref = gammalib.GFitsTableDoubleCol('RefCrossSection',nrows)
        paroi_lim.unit('cm3/s')
        paroi_ref.unit('cm3/s')

        #   Fill fits
        for i, result in enumerate(results):
            e_min[i]       = result['e_min']
            e_max[i]       = result['e_max']
            mass[i]        = result['mass']
            flux[i]        = result['flux']
            flux_err[i]    = result['flux_err']
            slogl[i]       = result['logL']
            TSvalues[i]    = result['TS']
            ulim_values[i] = result['ulimit']
            sc_factor[i]   = result['sc_factor']
            paroi_lim[i]   = result['sigmav_lim']
            paroi_ref[i]   = result['sigmav_ref']

            # Special Case: Add Column for Delta TS:
            for npoint in range(ncols):
                tsscan[i,npoint] = result['dts'][npoint]
                norms[i,npoint]  = result['norms'][npoint]

        #   Initialise FITS Table with extension "DMATTER"
        table = gammalib.GFitsBinTable(nrows)
        table.extname('DMATTER')

        #   Add Header for compatibility with gammalib.GMWLSpectrum
        table.card('INSTRUME','CTA','Name of Instrument')
        table.card('TELESCOP','CTA','Name of Telescope')

        #   Append filled columns to fits table
        table.append(e_min)
        table.append(e_max)
        table.append(mass)
        table.append(flux)
        table.append(flux_err)
        table.append(e2flux)
        table.append(e2flux_err)
        table.append(slogl)
        table.append(TSvalues)
        table.append(ulim_values)
        table.append(sc_factor)
        table.append(paroi_lim)
        table.append(paroi_ref)
        table.append(tsscan)
        table.append(norms)

        #   Create table for covariance matrix
        thisresult = results[0]
        parnames   = thisresult['parnames']
        size       = len(parnames)
        covmat     = gammalib.GFitsBinTable()
        par        = gammalib.GFitsTableStringCol('Parameters',1,100,size)

        #   Fill columns
        for i in range(size):
            par[0,i] = parnames[i]

        covmat.append(par)

        for i,result in enumerate(results):
            covar = result['correlation']
            colnm = 'Covariance (Mass {})'.format(i)
            cov   = gammalib.GFitsTableDoubleCol(colnm,1,size*size)

            for nrow in range(size):
                for ncol in range(size):
                    cov[0,nrow*size+ncol] = covar[nrow,ncol]

            #   Set dimension of covmatrix
            cov.dim([size,size])

            #   Append columns to covariance table
            covmat.append(cov)

        #   Set Extension name
        covmat.extname("Correlation Matrix")


        #   Create table for the model parameters:
        #   I think, I am hard-coding several things
        #   But, the 2 here is the size of the parameter list
        #   Value and error
        parvalues = gammalib.GFitsBinTable(nrows)
        for name in parnames:
            thiscol = gammalib.GFitsTableDoubleCol(name,nrows,2)
            for i,result in enumerate(results):
                # print(result[name])
                thiscol[i,0] = result[name][0]
                thiscol[i,1] = result[name][1]
            parvalues.append(thiscol)

        #   Set extension name
        parvalues.extname('Fitted parameters')

        #   Create the FITS file now
        self._fits = gammalib.GFits()
        self._fits.append(table)
        self._fits.append(parvalues)
        self._fits.append(covmat)

        #   Return
        return

    def process(self):
        """
        Run the script
        """

        #   Screen logging
        if self._logDebug():
            self._log.cout(True)

        #   Getting parameters
        self._get_parameters()
        srcname = self['srcname'].string()

        #   Write input observation container into logger
        self._log_observations(gammalib.NORMAL,self.obs(),'Input observation')

        #   Fit model
        results = self._fit_mass_points()

        #   Create FITS file
        self._create_fits(results)

        #   Publishing...?
        if self['publish'].boolean():
            self.publish()

        #   Return
        return

    def save(self):
        """
        Save dmatter analysis results
        """
        #   Write header
        self._log_header1(gammalib.TERSE,'Save dmatter analysis results')

        #   Continue only if FITS file is valid
        if self._fits != None:
            #   Get outmap parameter
            outfile = self['outfile'].filename()

            #   Log file name
            self._log_value(gammalib.NORMAL,'dmatter file',outfile.url())

            #   Save results
            self._fits.saveto(outfile,self['clobber'].boolean())

        #   Return
        return

    def publish(self,name=''):
        """
        Publish results
        """

        #   Write header
        self._log_header1(gammalib.TERSE,'Publishing Results')

        #   Checking fits
        if self._fits != None:
            if not name:
                user_name = self._name()
            else :
                user_name = name

        #   Log file
        self._log_value(gammalib.NORMAL,'DM anna name',user_name)

        #   Publishing
        self._fits.publish('DM Anna',user_name)

        #   Return
        return

    def dmatter_fits(self):
        """
        Return fits with results
        """
        #   Return
        return self._fits

    def jfactor(self):
        """
        Return jfactor value
        """
        #   Return
        return self._jfactor

    def sigmav_ref(self):
        """
        Return Reference Cross-section value
        """
        #   Return
        return self._sigmavref

    def sigmav_ul(self):
        """
        Return Reference Cross-section value
        """
        #   Return
        return self._sigmavUL

    def masses(self):
        """
        Return list of masses used to compute hte ULs
        """
        #   Return
        return self._masses

    def deltaTS(self):
        """
        Return lista of profile TS
        """

        #   Return
        return self._deltats

# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Create instance of application
    app = csdmatter(sys.argv)

    # Execute application
    app.execute()
