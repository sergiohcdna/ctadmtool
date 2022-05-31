import gammalib
import math
import numpy as np
from ctadmtool.dmspectrum.dmspectra import dmspectrum
from ctadmtool.tools.misc import ValidValue , ValidString
from tqdm import tqdm

import warnings

ALLOWED_FERMIONS = ('Majorana', 'Dirac')

ALLOWED_CHANNELS = ('eL', 'eR', 'e',
    'MuL', 'MuR', 'Mu', 'TauL', 'TauR', 'Tau',
    'q', 'c', 'b', 't',
    'WL', 'WT', 'W', 'ZL', 'ZT', 'Z', 'g', 'Gamma', 'h',
    'Nue', 'NuMu', 'NuTau',
    'Ve', 'VMu', 'VTau')

ALLOWED_CHANNELSNOEW = ('e','Mu','Tau','q','c','b','t','W','Z','g')

@ValidValue("_dfactor", min_val=1.e-40)
@ValidValue('_lifetime', min_val=1.e-40)
@ValidValue("_jfactor", min_val=1.e-40)
@ValidValue('_sigmav', min_val=1.e-40)
@ValidString('_delta', empty_allowed=False, options=ALLOWED_FERMIONS)
@ValidValue('_mmin', min_val=10.0)
@ValidValue('_mmax', max_val=1.e+5)
@ValidString('_srcname', empty_allowed=False)
class dmtable() :
    """
    Class to compute the flux generated
    from annihilation of dark matter
    particles.
    dmflux is a derived class from dmspectrum.
    """

    #   Init
    def __init__(self, srcname, mmin, mmax, mpoints, dminterp, delta='Majorana',
        sigmav=3.6e-26, jfactor=1.e+19, lifetime=1.e+30, dfactor=1.e+19) :
        """
        Initialize the dmflux_anna class

        Parameters:
        ----------
            srcname  : Name of the target or family targets
            mmin     : Min Mass of dark matter candidate
            mmax     : Max Mass of dark matter candidate
            mpoints  : Number of mass points to create the Fits table
            dminterp : dmspectrum class instance (I avoid to write a lot
                       of code I already have)
            delta    : Parameter to describe if dark matter candidate
                       is a Majorana (delta=2) fermion or a
                       Dirac (delta=4) fermion
            sigmav   : Annihilation cross-section (in cm**3/s)
            jfactor  : Astrophysical factor in (GeV**2/cm**5)
            lifetime : Decay lifetime (in s)
            dfactor  : Astrophysical factor in (GeV/cm**2)
        """
        #   And, I check that mmin < mmax, if not, then reverse the order
        if mmin > mmax :
            msg = ('\nI found that Minimum mass {0} '.format(mmin) +
                'is greater than Maximum mass {0}.\n'.format(mmax) +
                'Changing the order...')
            warnings.warn(msg, RuntimeWarning)
            m_min = mmax
            m_max = mmin
        else :
            m_min = mmin
            m_max = mmax

        #   Initialize parameters of dmflux_ana class
        self._srcname  = srcname
        self._sigmav   = sigmav
        self._jfactor  = jfactor
        self._lifetime = lifetime
        self._dfactor  = dfactor
        self._delta    = delta
        self._mmin     = m_min
        self._mmax     = m_max
        self._mpoints  = mpoints

        if not isinstance(dminterp, dmspectrum) :
            msg = 'dminterp must be an instance of dmspectrum class'
            raise TypeError(msg)
        else :
            self._dminterp = dminterp

        self._masses  = self._marray(m_min, m_max, mpoints)

        if dminterp.hasEW :
            self._allowed_channels = ALLOWED_CHANNELS
        else :
            self._allowed_channels = ALLOWED_CHANNELSNOEW

        self._model = None
        #   Return
        return

    @property
    def sigmav(self) :
        """
        Return value of the annihilation cross-section
        used to compute the flux
        """
        #   Return
        return self._sigmav

    @sigmav.setter
    def sigmav(self, sigmav) :
        """
        Set the value of Annihilation cross-section (in cm**3/s)
        used to compute the flux

        Parameters
        ----------
            sigmav : Annihilation cross-section (cm**3/s)
        """
        #   Check that sigmav is greater than 1.e-35
        if sigmav < 1.e-40 :
            raise ValueError(('\nValue of annihilation cross-section ' +
                ' must be greater than 1.e-40.\n' +
                'This is just to avoid possible round errors'))

        #   Set sigmav
        self._sigmav = sigmav

        #   Return
        return

    @property
    def lifetime(self) :
        """
        Return value of the decay lifetime
        used to compute the flux
        """
        #   Return
        return self._lifetime

    @lifetime.setter
    def lifetime(self, tau_chi) :
        """
        Set the value of decay lifetime (in s)
        used to compute the flux

        Parameters
        ----------
            tau_chi : Annihilation cross-section (cm**3/s)
        """
        #   Check that sigmav is greater than 1.e-35
        if tau_chi < 1.e-40 :
            raise ValueError(('\nValue of decay lifetime ' +
                ' must be greater than 1.e-40.\n' +
                'This is just to avoid possible round errors'))

        #   Set sigmav
        self._lifetime = tau_chi

        #   Return
        return

    @property
    def jfactor(self) :
        """
        Return the value of the Astrophysical factor
        used to compute the flux
        """
        #   Return
        return self._jfactor

    @jfactor.setter
    def jfactor(self, jfactor) :
        """
        Set the value of the Astrophysical factor (GeV**2/cm**5)
        to compute the dm flux

        Parameters
        ----------
            jfactor : Astrophysical factor J (GeV**2/cm**5)
        """
        if jfactor < 1.e-40 :
            raise ValueError('\nValue of jfactor must be greater than 1.e-40.')

        #   Set the jfactor
        self._jfactor = jfactor

        #   Return
        return

    @property
    def dfactor(self) :
        """
        Return the value of the Astrophysical factor
        used to compute the flux
        """
        #   Return
        return self._dfactor

    @dfactor.setter
    def dfactor(self, dfactor) :
        """
        Set the value of the Astrophysical factor (GeV/cm**2)
        to compute the dm flux

        Parameters
        ----------
            dfactor : Astrophysical factor D (GeV/cm**2)
        """
        #   Set the jfactor
        self._dfactor = dfactor

        #   Return
        return

    @property
    def mmin(self) :
        """
        Return Minimum value mass (GeV) used to compute
        the dm flux
        """
        #   Return
        return self._mmin

    @mmin.setter
    def mmin(self, m_min) :
        """
        Set the value of minimum mass (GeV) used to compute
        the dm flux
        """
        #   Just check that the minimum mass is greater than
        #   10.0 GeV.
        if m_min < 10. :
            raise ValueError(('\nMinimum mass {0} GeV '.format(m_min) +
                'is below the allowed value (10GeV)'))

        #   Set minimum energy
        self._mmin = m_min

        #   Update masses
        mvalues      = self._marray(self._mmin, self._mmax, self._mpoints)
        self._masses = mvalues

        #   Return
        return

    @property
    def mmax(self) :
        """
        Return Maximum value of mass (GeV) used to compute
        the dm flux
        """
        #   Return
        return self._mmax

    @mmax.setter
    def mmax(self, m_max) :
        """
        Set the value of minimum mass (GeV) used to compute
        the dm flux
        """
        if m_max > 1.e+5 :
            raise ValueError(('\nMaximum mass {0} GeV '.format(m_max) +
                'is above the allowed value (1.e+5GeV)'))

        #   Set minimum energy
        self._mmax = m_max

        #   Update masses
        mvalues      = self._marray(self._mmin, self._mmax, self._mpoints)
        self._masses = mvalues

        #   Return
        return

    @property
    def masses(self) :
        """
        Return the values of the energy array used to compute the spectrum
        """
        #   Return
        return self._masses

    @masses.setter
    def masses(self, m_vals) :
        """
        Set the masses used to compute the spectrum
        Parameters
        ----------
            - evals   : tuple with:
                - mmin    : Minimum mass (GeV)
                - mmax    : Maximum mass (GeV)
                - mpoints : Number of points to create the array
        """

        mmin, mmax, mpoints = m_vals

        #   Check if emin and emax are valid
        if mmin < 10.0 :

            raise ValueError(('Mass {0} '.format(mmin) +
                'is lower than the allowed value 10.0'))

        if mmax > 1.e+5 :

            raise ValueError(('Mass {0} '.format(mmax) +
                'is greater than the allowed value 1.e+5'))

        #   Create energy array
        mvalues = self._marray(mmin, mmax, mpoints)

        self._masses = mvalues

        #   Return
        return

    @staticmethod
    def _marray(mmin, mmax, mpoints) :
        """
        Create list of masses to generate the fits table.
        The calculation is based in the number of points
        The masses are computed assuming logarithmic distance
        """
        logmmin = np.log10(mmin)
        logmmax = np.log10(mmax)
        width   = (logmmax - logmmin)/(mpoints-1)
        masses  = []

        for index in range(mpoints) :
            masses.append(math.pow(10., logmmin+index*width))

        #   Return
        return masses

    @property
    def delta(self) :
        """
        Return what kind of dark matter particle is
        used to compute the dm flux
        """
        #   Return
        return self._delta

    @delta.setter
    def delta(self, delta) :
        """
        Set the value of delta to describe what kind of
        dark matter particle is used to compute the
        dm flux.

        Parameters
        ----------
            delta : String, either Majorana or Dirac
        """

        #   Just to check that delta is valid
        if delta not in ALLOWED_FERMIONS :
            raise ValueError(('\nKind of Dark matter particle not ' +
                'supported.\nOptions are:{0}'.format(ALLOWED_FERMIONS)))

        #   Set minimum energy
        self._delta = delta

        #   Return
        return

    @property
    def hasEW(self) :
        """
        Return whether EW corrections are included or not
        """
        #   Return
        return self._dminterp.hasEW

    @hasEW.setter
    def hasEW(self, has_EW) :
        """
        Include EW corrections in computation of DM spectra
        """
        self._dminterp.hasEW = has_EW

        #   Update the tuple of allowed channels
        if has_EW :
            self._allowed_channels = ALLOWED_CHANNELS
        else :
            self._allowed_channels = ALLOWED_CHANNELSNOEW

        #   Return
        return

    @property
    def allowed_channels(self) :
        """
        Return tuple of allowed channels according to
        whether or not to include EW corrections in spectra
        """
        #   Return
        return self._allowed_channels

    @property
    def tablemodel(self) :
        """
        Return GModelSpectralTable
        """
        #   Return
        return self._model

    @property
    def process(self) :
        """
        Return dm process
        """
        #   Return
        return self._dminterp.process

    @process.setter
    def process(self, process_vals) :
        """
        Set annihilation (anna) or decay process in dminterp
        Also update the properties jfactor and sigmav for anna
        or dfactor and lifetime for decay
        """
        #   Extract values
        dmprocess = process_vals[0]
        astfactor = process_vals[1]
        paroi     = process_vals[2]

        #   Check that process is valid
        VALID_PROCESSES = ['anna', 'decay']
        if dmprocess not in VALID_PROCESSES :
            msg = 'Valid options are: {0}'.format(VALID_PROCESSES)
            raise ValueError(msg)

        if astfactor < 1.e-40 or paroi < 1.e-40 :
            raise ValueError('\nParameters must be greater than 1.e-40.')

        #   Update properties
        if dmprocess == 'anna' :
            self._jfactor = astfactor
            self._sigmav  = paroi
        elif dmprocess == 'decay' :
            self._dfactor  = astfactor
            self._lifetime = paroi

        self._dminterp.process = dmprocess

        #   Update

        #   Return
        return

    @property
    def elist(self) :
        """
        Return list of energy values used to compute the spectrum
        """
        #   Return
        return self._dminterp.energy

    @elist.setter
    def elist(self, evals) :
        """
        Update energy values used to compute the spectrum
        evals[0]  --> emin
        evals[1]  --> emax
        evals[2]  --> epoints
        """

        #   Check that emin and emax are ok
        #   Note, that I set the minimum to 500 MeV
        #   There is no meaning to go to lower energies
        #   In the case of CTA
        if evals[0] < 5.0e-3 or evals[1] > 1.e+5 :
            raise ValueError('\nParameters outside of range')

        #   Update properties
        self._dminterp.energy = evals

        #   Return
        return

    @staticmethod
    def _norm_anna(sigmav, mass, delta, jfactor) :
        """
        Compute normalization of the dm flux compatible with gammalib

        Parameters
        ----------
            sigmav  : Value of annihilation cross-section (cm**3/s)
            mass    : Mass of dark matter particles (GeV)
            delta   : String to indicate if dark matter is a
                      Majorana or Dirac fermion
            jfactor : Astrophysica factor for annihilation

        Return
        ------
            norm : (1/[MeV* cm^2 * s])
        """
        d = 0.
        #   Check delta
        if delta == 'Majorana' :
            d = 2.
        elif delta == 'Dirac' :
            d = 4.

        #   Compute ppfactor
        ppfactor = sigmav / (d*4.*gammalib.pi*mass*mass)
        norm     = ppfactor * jfactor

        return norm * 1.0e-3

    @staticmethod
    def _norm_decay(lifetime, mass, dfactor) :
        """
        Compute normalization of the dm flux compatible with gammalib

        Parameters
        ----------
            lifetime : Value of decay lifetime (s)
            mass     : Mass of dark matter particles (GeV)
            dfactor  : Astrophysical factor for ecay

        Return
        ------
            norm : (1/[MeV* cm^2 * s])
        """
        #   Compute ppfactor
        ppfactor = 1 / (4.*gammalib.pi*mass*lifetime)
        norm     = ppfactor * dfactor

        return norm * 1.0e-3

    def create_modeltable(self) :
        """
        Create fits table with spectrum and channels
        """
        #   Get list of channel indices
        #   First, I get the number of channels and energy points
        #   I don't want to access a private member from dmspectrum
        #   class, but I can get the number of points from the
        #   energy array
        ch_indices = [i for i in range(len(self._allowed_channels))]
        n_chs   = len(ch_indices)
        n_eng   = len(self._dminterp.energy)

        # Array with definitions of energy bins
        gemin = gammalib.GEnergy(self._dminterp.emin, 'GeV')
        gemax = gammalib.GEnergy(self._dminterp.emax, 'GeV')
        ebins = gammalib.GEbounds(n_eng, gemin, gemax)

        #   Then create the GModelPar objects for mass and channel
        #   I know, default channel is hard coded, but we don't need
        #   to select any particular channel at this moment.
        #   Select Tau channel is just for initialization
        dmmass    = gammalib.GModelPar('Mass', self._mmin, 1.0)
        dmmass.unit('GeV')
        index     = self._allowed_channels.index('Tau')
        dmchannel = gammalib.GModelPar('Channel', index, 1.0)

        #   Create the GSpectralTablePar objects
        par_mass    = gammalib.GModelSpectralTablePar(dmmass, self._masses)
        par_channel = gammalib.GModelSpectralTablePar(dmchannel, ch_indices)

        #   Create the container GSpectralTablePars and append the pars
        pars = gammalib.GModelSpectralTablePars()
        pars.append(par_mass)
        pars.append(par_channel)

        #   Get ppfactor and normalization
        #   This normalization computed here
        #   is not neccessary. You can change the normalization
        #   of the GModelSpectralTable later during simulation
        #   or analysis steps via GModelSpectralTable methods
        norm   = 0.0
        minval = 0.0
        maxval = 1.0e+20

        if self._dminterp.process == 'anna' :
            norm = self._norm_anna(self._sigmav, self._mmin,
                self._delta, self._jfactor)
        elif self._dminterp.process == 'decay' :
            norm = self._norm_decay(self._lifetime, self._mmin, self._dfactor)

        #   GNdarray to save the spectra
        spectra = gammalib.GNdarray(self._mpoints,n_chs,n_eng)

        #   filling the spectrum
        desc = 'Computing {}-spectrrum'.format(self._dminterp.process)
        for index, mass in tqdm(enumerate(self._masses),desc=desc,leave=False):
            #   Change the value of the mass
            self._dminterp.mass = mass
            for cindex, thisch in enumerate(self._allowed_channels):
                #    Modified the instance of dmspectrum
                #   to match values for every channel
                #   I don't need to change the array for energy
                #   And also, I don't need to check whether I want
                #   to include EW corrections or not
                self._dminterp.channel = thisch
                dmspec                 = self._dminterp.spectra()
                for eindex in range(n_eng):
                    spectra[index, cindex, eindex] = dmspec[eindex]

        #   Tuning the ModelSpectralTable
        #   I set the interpolation method of masses to logarithmic
        #   Mass and channel are fixed.
        #   Particularly, it's mandatory that channel parameter is fixed
        model = gammalib.GModelSpectralTable(ebins, pars, spectra)
        model.table_par('Mass').method(1)
        model.table_par('Channel').method(0)
        model['Mass'].fix()
        model['Channel'].fix()
        model['Normalization'].value(norm)
        model['Normalization'].scale(1.0)
        model['Normalization'].range(minval,maxval)

        self._model = model

        #   Return
        return

    def save(self) :
        """
        Save the DM table
        """
        process = self._dminterp.process
        ew      = int(self._dminterp.hasEW)
        name = 'DMModel{0}{1}EW{2}.fits'.format(process, self._srcname, ew)

        self._model.save(name, True)

        return

@ValidValue("_dfactor", min_val=1.e-40)
@ValidValue('_lifetime', min_val=1.e-40)
@ValidValue("_jfactor", min_val=1.e-40)
@ValidValue('_sigmav', min_val=1.e-40)
@ValidString('_delta', empty_allowed=False, options=ALLOWED_FERMIONS)
@ValidValue('_mmin', min_val=10.0)
@ValidValue('_mmax', max_val=1.e+5)
@ValidString('_srcname', empty_allowed=False)
class dmtable_ch() :
    """
    Class to compute the flux generated
    from annihilation of dark matter
    particles.
    dmflux is a derived class from dmspectrum.
    the suffix 'ch' stands for single channel,
    so the class only create table-models for
    specific channels
    """

    #   Init
    def __init__(self, srcname, mmin, mmax, mpoints, dminterp,
        channel='Tau', delta='Majorana', sigmav=3.6e-26, jfactor=1.e+19,
        lifetime=1.e+30, dfactor=1.e+19) :
        """
        Initialize the dmflux_anna class

        Parameters:
        ----------
            srcname  : Name of the target or family targets
            mmin     : Min Mass of dark matter candidate
            mmax     : Max Mass of dark matter candidate
            mpoints  : Number of mass points to create the Fits table
            dminterp : dmspectrum class instance (I avoid to write a lot
                       of code I already have)
            delta    : Parameter to describe if dark matter candidate
                       is a Majorana (delta=2) fermion or a
                       Dirac (delta=4) fermion
            sigmav   : Annihilation cross-section (in cm**3/s)
            jfactor  : Astrophysical factor in (GeV**2/cm**5)
            lifetime : Decay lifetime (in s)
            dfactor  : Astrophysical factor in (GeV/cm**2)
        """
        #   And, I check that mmin < mmax, if not, then reverse the order
        if mmin > mmax :
            msg = ('\nI found that Minimum mass {0} '.format(mmin) +
                'is greater than Maximum mass {0}.\n'.format(mmax) +
                'Changing the order...')
            warnings.warn(msg, RuntimeWarning)
            m_min = mmax
            m_max = mmin
        else :
            m_min = mmin
            m_max = mmax

        #   Initialize parameters of dmflux_ana class
        self._srcname  = srcname
        self._sigmav   = sigmav
        self._jfactor  = jfactor
        self._lifetime = lifetime
        self._dfactor  = dfactor
        self._delta    = delta
        self._mmin     = m_min
        self._mmax     = m_max
        self._mpoints  = mpoints
        self._channel  = channel

        if not isinstance(dminterp, dmspectrum) :
            msg = 'dminterp must be an instance of dmspectrum class'
            raise TypeError(msg)
        else :
            self._dminterp = dminterp

        self._masses  = self._marray(m_min, m_max, mpoints)

        if dminterp.hasEW :
            self._allowed_channels = ALLOWED_CHANNELS
        else :
            self._allowed_channels = ALLOWED_CHANNELSNOEW

        #   Check if channel is valid
        if channel not in self._allowed_channels:
            msg = ('\nChannel {0} not found in'.format(channel) +
                'allowed channels. Options are: {0}'.format(ALLOWED_FERMIONS))
            raise ValueError(msg)

        #   Update channel property of spectrum interpolator dminterp
        #   Only if the channels are different
        if dminterp.channel != channel :
            dminterp.channel = channel

        self._model = None
        #   Return
        return

    @property
    def sigmav(self) :
        """
        Return value of the annihilation cross-section
        used to compute the flux
        """
        #   Return
        return self._sigmav

    @sigmav.setter
    def sigmav(self, sigmav) :
        """
        Set the value of Annihilation cross-section (in cm**3/s)
        used to compute the flux

        Parameters
        ----------
            sigmav : Annihilation cross-section (cm**3/s)
        """
        #   Check that sigmav is greater than 1.e-35
        if sigmav < 1.e-40 :
            raise ValueError(('\nValue of annihilation cross-section ' +
                ' must be greater than 1.e-40.\n' +
                'This is just to avoid possible round errors'))

        #   Set sigmav
        self._sigmav = sigmav

        #   Return
        return

    @property
    def lifetime(self) :
        """
        Return value of the decay lifetime
        used to compute the flux
        """
        #   Return
        return self._lifetime

    @lifetime.setter
    def lifetime(self, tau_chi) :
        """
        Set the value of decay lifetime (in s)
        used to compute the flux

        Parameters
        ----------
            tau_chi : Annihilation cross-section (cm**3/s)
        """
        #   Check that sigmav is greater than 1.e-35
        if tau_chi < 1.e-40 :
            raise ValueError(('\nValue of decay lifetime ' +
                ' must be greater than 1.e-40.\n' +
                'This is just to avoid possible round errors'))

        #   Set sigmav
        self._lifetime = tau_chi

        #   Return
        return

    @property
    def jfactor(self) :
        """
        Return the value of the Astrophysical factor
        used to compute the flux
        """
        #   Return
        return self._jfactor

    @jfactor.setter
    def jfactor(self, jfactor) :
        """
        Set the value of the Astrophysical factor (GeV**2/cm**5)
        to compute the dm flux

        Parameters
        ----------
            jfactor : Astrophysical factor J (GeV**2/cm**5)
        """
        if jfactor < 1.e-40 :
            raise ValueError('\nValue of jfactor must be greater than 1.e-40.')

        #   Set the jfactor
        self._jfactor = jfactor

        #   Return
        return

    @property
    def dfactor(self) :
        """
        Return the value of the Astrophysical factor
        used to compute the flux
        """
        #   Return
        return self._dfactor

    @dfactor.setter
    def dfactor(self, dfactor) :
        """
        Set the value of the Astrophysical factor (GeV/cm**2)
        to compute the dm flux

        Parameters
        ----------
            dfactor : Astrophysical factor D (GeV/cm**2)
        """
        #   Set the jfactor
        self._dfactor = dfactor

        #   Return
        return

    @property
    def mmin(self) :
        """
        Return Minimum value mass (GeV) used to compute
        the dm flux
        """
        #   Return
        return self._mmin

    @mmin.setter
    def mmin(self, m_min) :
        """
        Set the value of minimum mass (GeV) used to compute
        the dm flux
        """
        #   Just check that the minimum mass is greater than
        #   10.0 GeV.
        if m_min < 10. :
            raise ValueError(('\nMinimum mass {0} GeV '.format(m_min) +
                'is below the allowed value (10GeV)'))

        #   Set minimum energy
        self._mmin = m_min

        #   Update masses
        mvalues      = self._marray(self._mmin, self._mmax, self._mpoints)
        self._masses = mvalues

        #   Return
        return

    @property
    def mmax(self) :
        """
        Return Maximum value of mass (GeV) used to compute
        the dm flux
        """
        #   Return
        return self._mmax

    @mmax.setter
    def mmax(self, m_max) :
        """
        Set the value of minimum mass (GeV) used to compute
        the dm flux
        """
        if m_max > 1.e+5 :
            raise ValueError(('\nMaximum mass {0} GeV '.format(m_max) +
                'is above the allowed value (1.e+5GeV)'))

        #   Set minimum energy
        self._mmax = m_max

        #   Update masses
        mvalues      = self._marray(self._mmin, self._mmax, self._mpoints)
        self._masses = mvalues

        #   Return
        return

    @property
    def masses(self) :
        """
        Return the values of the energy array used to compute the spectrum
        """
        #   Return
        return self._masses

    @masses.setter
    def masses(self, m_vals) :
        """
        Set the masses used to compute the spectrum
        Parameters
        ----------
            - evals   : tuple with:
                - mmin    : Minimum mass (GeV)
                - mmax    : Maximum mass (GeV)
                - mpoints : Number of points to create the array
        """

        mmin, mmax, mpoints = m_vals

        #   Check if emin and emax are valid
        if mmin < 10.0 :

            raise ValueError(('Mass {0} '.format(mmin) +
                'is lower than the allowed value 10.0'))

        if mmax > 1.e+5 :

            raise ValueError(('Mass {0} '.format(mmax) +
                'is greater than the allowed value 1.e+5'))

        #   Create energy array
        mvalues = self._marray(mmin, mmax, mpoints)

        self._masses = mvalues

        #   Return
        return

    @staticmethod
    def _marray(mmin, mmax, mpoints) :
        """
        Create list of masses to generate the fits table.
        The calculation is based in the number of points
        The masses are computed assuming logarithmic distance
        """
        logmmin = np.log10(mmin)
        logmmax = np.log10(mmax)
        width   = (logmmax - logmmin)/(mpoints-1)
        masses  = []

        for index in range(mpoints) :
            masses.append(math.pow(10., logmmin+index*width))

        #   Return
        return masses

    @property
    def delta(self) :
        """
        Return what kind of dark matter particle is
        used to compute the dm flux
        """
        #   Return
        return self._delta

    @delta.setter
    def delta(self, delta) :
        """
        Set the value of delta to describe what kind of
        dark matter particle is used to compute the
        dm flux.

        Parameters
        ----------
            delta : String, either Majorana or Dirac
        """

        #   Just to check that delta is valid
        if delta not in ALLOWED_FERMIONS :
            raise ValueError(('\nKind of Dark matter particle not ' +
                'supported.\nOptions are:{0}'.format(ALLOWED_FERMIONS)))

        #   Set minimum energy
        self._delta = delta

        #   Return
        return

    @property
    def hasEW(self) :
        """
        Return whether EW corrections are included or not
        """
        #   Return
        return self._dminterp.hasEW

    @hasEW.setter
    def hasEW(self, has_EW) :
        """
        Include EW corrections in computation of DM spectra
        """
        self._dminterp.hasEW = has_EW

        #   Update the tuple of allowed channels
        if has_EW :
            self._allowed_channels = ALLOWED_CHANNELS
        else :
            self._allowed_channels = ALLOWED_CHANNELSNOEW

        #   Return
        return

    @property
    def allowed_channels(self) :
        """
        Return tuple of allowed channels according to
        whether or not to include EW corrections in spectra
        """
        #   Return
        return self._allowed_channels

    @property
    def channel(self) :
        '''
        Return channel used to compute the gamma-ray flux
        '''
        #   Return
        return self._channel

    @channel.setter
    def channel(self, ch) :
        '''
        Set channel used to compute the dmspectrum.
        Also updates the channel parameter of the 
        spectrum interpolator dminterp
        If channel is not valid, raise value error
        '''
        #   Check if channel is valid
        if ch not in self._allowed_channels :
            msg = ('\nChannel {0} not found in'.format(channel) +
                'allowed channels. Options are: {0}'.format(ALLOWED_FERMIONS))
            raise ValueError(msg)

        #   Set channel
        self._channel = ch

        #   Update dminterp instance
        self._dminterp.channel = ch

        #   Return
        return

    @property
    def tablemodel(self) :
        """
        Return GModelSpectralTable
        """
        #   Return
        return self._model

    @property
    def process(self) :
        """
        Return dm process
        """
        #   Return
        return self._dminterp.process

    @process.setter
    def process(self, process_vals) :
        """
        Set annihilation (anna) or decay process in dminterp
        Also update the properties jfactor and sigmav for anna
        or dfactor and lifetime for decay
        """
        #   Extract values
        dmprocess = process_vals[0]
        astfactor = process_vals[1]
        paroi     = process_vals[2]

        #   Check that process is valid
        VALID_PROCESSES = ['anna', 'decay']
        if dmprocess not in VALID_PROCESSES :
            msg = 'Valid options are: {0}'.format(VALID_PROCESSES)
            raise ValueError(msg)

        if astfactor < 1.e-40 or paroi < 1.e-40 :
            raise ValueError('\nParameters must be greater than 1.e-40.')

        #   Update properties
        if dmprocess == 'anna' :
            self._jfactor = astfactor
            self._sigmav  = paroi
        elif dmprocess == 'decay' :
            self._dfactor  = astfactor
            self._lifetime = paroi

        self._dminterp.process = dmprocess

        #   Update

        #   Return
        return

    @property
    def elist(self) :
        """
        Return list of energy values used to compute the spectrum
        """
        #   Return
        return self._dminterp.energy

    @elist.setter
    def elist(self, evals) :
        """
        Update energy values used to compute the spectrum
        evals[0]  --> emin
        evals[1]  --> emax
        evals[2]  --> epoints
        """

        #   Check that emin and emax are ok
        #   I set the minimum to 500 MeV
        if evals[0] < 5.0e-3 or evals[1] > 1.e+5 :
            raise ValueError('\nParameters outside of range')

        #   Update properties
        self._dminterp.energy = evals

        #   Return
        return

    @staticmethod
    def _norm_anna(sigmav, mass, delta, jfactor) :
        """
        Compute normalization of the dm flux compatible with gammalib

        Parameters
        ----------
            sigmav  : Value of annihilation cross-section (cm**3/s)
            mass    : Mass of dark matter particles (GeV)
            delta   : String to indicate if dark matter is a
                      Majorana or Dirac fermion
            jfactor : Astrophysica factor for annihilation

        Return
        ------
            norm : (1/[MeV* cm^2 * s])
        """
        d = 0.
        #   Check delta
        if delta == 'Majorana' :
            d = 2.
        elif delta == 'Dirac' :
            d = 4.

        #   Compute ppfactor
        ppfactor = sigmav / (d*4.*gammalib.pi*mass*mass)
        norm     = ppfactor * jfactor

        return norm * 1.0e-3

    @staticmethod
    def _norm_decay(lifetime, mass, dfactor) :
        """
        Compute normalization of the dm flux compatible with gammalib

        Parameters
        ----------
            lifetime : Value of decay lifetime (s)
            mass     : Mass of dark matter particles (GeV)
            dfactor  : Astrophysical factor for ecay

        Return
        ------
            norm : (1/[MeV* cm^2 * s])
        """
        #   Compute ppfactor
        ppfactor = 1 / (4.*gammalib.pi*mass*lifetime)
        norm     = ppfactor * dfactor

        return norm * 1.0e-3

    def create_modeltable(self) :
        """
        Create fits table with spectrum for specific channel
        """
        #   Number of points in energy array
        n_eng   = len(self._dminterp.energy)

        # Array with definitions of energy bins
        # The min and max values are encapsulated in the
        # dm spectrum interpolator dminterp
        gemin = gammalib.GEnergy(self._dminterp.emin, 'GeV')
        gemax = gammalib.GEnergy(self._dminterp.emax, 'GeV')
        ebins = gammalib.GEbounds(n_eng, gemin, gemax)

        #   Then create the GModelPar objects for mass
        dmmass    = gammalib.GModelPar('Mass', self._mmin, 1.0)
        dmmass.unit('GeV')

        #   Create the GSpectralTablePar objects
        par_mass    = gammalib.GModelSpectralTablePar(dmmass, self._masses)

        #   Create the container GSpectralTablePars and append the pars
        pars = gammalib.GModelSpectralTablePars()
        pars.append(par_mass)

        #   GNdarray to save the spectra
        spectra = gammalib.GNdarray(self._mpoints,n_eng)

        #   Get ppfactor and normalization
        #   This normalization computed here
        #   is not neccessary. You can change the normalization
        #   of the GModelSpectralTable later during simulation
        #   or analysis steps via GModelSpectralTable methods
        norm   = 0.0
        minval = 0.0
        maxval = 1.0e+20

        if self._dminterp.process == 'anna' :
            norm = self._norm_anna(self._sigmav, self._mmin,
                self._delta, self._jfactor)
        elif self._dminterp.process == 'decay' :
            norm = self._norm_decay(self._lifetime, self._mmin, self._dfactor)

        #   filling the spectrum
        desc = 'Computing {}-spectrrum'.format(self._dminterp.process)
        for index, mass in tqdm(enumerate(self._masses),desc=desc,leave=False):
            #   Change the value of the mass
            self._dminterp.mass = mass
            dmspec              = self._dminterp.spectra()
            for eindex in range(n_eng):
                spectra[index, eindex] = dmspec[eindex]

        #   Tuning the ModelSpectralTable
        #   I set the interpolation method of masses to logarithmic
        #   Mass is a fixed parameter
        model = gammalib.GModelSpectralTable(ebins, pars, spectra)
        model.table_par('Mass').method(1)
        model['Mass'].scale(1.)
        model['Mass'].fix()
        model['Normalization'].value(norm)
        model['Normalization'].scale(1.0)
        model['Normalization'].range(minval,maxval)

        self._model = model

        #   Return
        return

    def save(self) :
        """
        Save the DM table
        """
        process = self._dminterp.process
        ew      = int(self._dminterp.hasEW)
        src     = self._srcname
        ch      = self._channel
        name = 'DMModel{0}{1}EW{2}Ch{3}.fits'.format(process, src, ew, ch)

        self._model.save(name, True)

        return
