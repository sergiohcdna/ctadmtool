import numpy as np
from scipy.interpolate import interp2d
from ebltable.tau_from_model import OptDepth
from ctaAnalysis.tools.misc import ValidnpArray , ValidString , ValidValue

import os

#   List of allowed channels for annihilation of DM
#   from PPPC4DMID tables
ALLOWED_CHANNELS = ('eL', 'eR', 'e',
    'MuL', 'MuR', 'Mu', 'TauL', 'TauR', 'Tau',
    'q', 'c', 'b', 't',
    'WL', 'WT', 'W', 'ZL', 'ZT', 'Z', 'g', 'Gamma', 'h',
    'Nue', 'NuMu', 'NuTau',
    'Ve', 'VMu', 'VTau')

ALLOWED_CHANNELSNOEW = ('e','Mu','Tau','q','c','b','t','W','Z','g')

ALLOWED_EBLMODELS = ('franceschini', 'franceschini2017',
    'kneiske', 'finke', 'dominguez', 'dominguez-upper',
    'dominguez-lower','inuoe', 'inuoe-low-pop3',
    'inuoe-up-pop3','gilmore', 'gilmore-fixed')

ALLOWED_PROCESSES = ('anna', 'decay')

#   1. For mass, the allowed values are taken from PPPC4DMID project
#      (5GeV , 100TeV)
#   2. For channels, the allowed values are listed in
#      ALLOWED_CHANNELS tuple
#   3. For redshift z, I check for positive values
#      from ebl-table project 
#      (https://github.com/me-manu/ebltable/blob/master/ebltable/tau_from_model.py)
@ValidString("_eblmodel", empty_allowed=False, options=ALLOWED_EBLMODELS)
@ValidString("_process", empty_allowed=False, options=ALLOWED_PROCESSES)
@ValidValue("_z", min_val=0)
@ValidString("_channel", empty_allowed=False, options=ALLOWED_CHANNELS)
@ValidValue("_mass", min_val=5, max_val=1.e+5)
@ValidValue("_emin", min_val=5.e-9)
@ValidValue("_emax", max_val=1.e+5)
# @ValidnpArray("_energy", min_val=5.e-9, max_val=1.e+5)
class dmspectrum() :
    """
    Implementation of a class to compute photon spectra
    due to annihilations or decay of dark matter particles.
    The calculation is based tables from PPPC4DMID project:
        http://www.marcocirelli.net/PPPC4DMID.html
    You can select between spectrum wit or without EW
    corrections.
    Theres is a local copy of files with PPPC4DMID tables
    """

    def __init__(self, dm_mass, emin, emax, channel, z, process='anna',\
        eblmod='franceschini2017', has_EW=True, epoints=100) :
        """
        Initiate dark matter class

        Parameters
        ----------

        dm_mass : Mass of dark matter particle in GeV
                  Because PPPC4DMID uses in GeV
        emin    : Minimum energy to compute spectra (GeV)
        emax    : Maximum energy to compute spectra (GeV)
        channel : Annihilation/Decay channel
        z       : Redshift of source
                  This is to be able to include
                  EBL attenuation
        process : Annihilation (anna) or Decay (decay) of
                  dark matter particles
        eblmod  : EBL model used to compute attenuation
        epoints : Number of points in energy spectrum
        """

        self._emin     = emin
        self._emax     = emax
        self._channel  = channel
        self._z        = z
        self._mass     = dm_mass
        self._process  = process
        self._eblmodel = eblmod
        self._ew       = has_EW
        self._epoints  = epoints
        self._energy   = self._earray(emin, emax, epoints)

        #   Return
        return

    @classmethod
    def dminterp2d(cls, dm_channel, z) :
        """
        Initialize class with some default parameters

        Parameters
        ----------
            - dm_channel
            - z
        """

        m      = 100.0
        emin   = 30.0
        emax   = 1.0e+4
        dmspec = cls(m, emin, emax, dm_channel, z)

        #   Return
        return dmspec

    @property
    def mass(self) :
        """
        Return value of mass of the dark matter candidate
        """
        #   Return
        return self._mass

    @mass.setter
    def mass(self, dm_mass) :
        """
        Set the value of mass of the dark matter candidate

        Parameters
        ----------
            dm_mass : Mass (in GeV) [Valid range: 5GeV-100TeV]
        """

        #   Check if mdm has a valid value
        #   The valid range is between 5GeV and 100 TeV.
        #   If not, then raise Value Error
        if not (5 <= dm_mass <= 1.e+5) :

            raise ValueError(('\nMass of DM particle ' +
                'with value {0} '.format(dm_mass) +
                'is out of range: [5,1.e+5] GeV'))

        #   Set value of mass
        self._mass = dm_mass

        #   Return
        return

    @property
    def emin(self) :
        """
        Return the values of min energy used to compute the spectrum
        """

        #   Return
        return self._emin

    @emin.setter
    def emin(self, e_min) :
        """
        Set the minimum energy to compute the spectrum
        Parameters
        ----------
            - e_min : Energy (GeV)

        Additionally, check if the values is valid (>5.e-9)
        """
        if e_min < 5.e-9 :
            msg = 'Energy {0} GeV is lower than the value allowed'.format(e_min)
            raise ValueError(msg)

        # Set energy
        self._emin = e_min

        #   Return
        return

    @property
    def emax(self) :
        """
        Return the values of min energy used to compute the spectrum
        """

        #   Return
        return self._emax

    @emax.setter
    def emax(self, e_max) :
        """
        Set the minimum energy to compute the spectrum
        Parameters
        ----------
            - e_max : Energy (GeV)

        Additionally, check if the values is valid (>5.e-9)
        """
        if e_max > 1.0e+5 :
            msg = 'Energy {0} GeV is greater than value allowed'.format(e_max)
            raise ValueError(msg)

        # Set energy
        self._emax = e_max

    @property
    def energy(self) :
        """
        Return the values of the energy array used to compute the spectrum
        """

        #   Return
        return self._energy

    @energy.setter
    def energy(self, e_vals) :
        """
        Set the energy used to compute the spectrum
        Parameters
        ----------
            - evals   : tuple with:
                - emin    : Minimum energy (GeV)
                - emax    : Maximum energy (GeV)
                - epoints : Number of points to create the array
        """

        emin, emax, epoints = e_vals

        #   Check if emin and emax are valid
        if emin < 5.e-9 :

            raise ValueError(('Energy {0} '.format(emin) +
                'is lower than the allowed value 5.e-9'))

        if emax > 1.e+5 :

            raise ValueError(('Energy {0} '.format(emax) +
                'is greater than the allowed value 1.e+5'))

        #   Create energy array
        energies = self._earray(emin, emax, epoints)

        self._energy = energies

        #   Return
        return

    @staticmethod
    def _earray(emin, emax, epoints) :
        """
        Create energy array to compute the spectra.
        The calculation is based in the number of points
        Return an np.array instance. The energies are
        computed assuming logarithmic distance
        """
        logemin  = np.log10(emin)
        logemax  = np.log10(emax)
        energies = np.logspace(logemin, logemax, epoints)

        #   Return
        return energies

    @property
    def channel(self) :
        """
        Return the channel used to compute the spectra.
        """
        #   Return
        return self._channel

    @channel.setter
    def channel(self, dm_channel) :
        """
        Set the channel used to compute the spectra

        Parameters
        ----------
            dm_channel: Channel

        Options:
        'eL' , 'eR' , 'e' ,
        'MuL' , 'MuR' , 'Mu' ,
        'TauL' , 'TauR' , 'Tau' ,
        'q' , 'c' , 'b' , 't' ,
        'WL' , 'WT' , 'W' ,
        'ZL' , 'ZT' , 'Z' ,
        'g' , 'Gamma' , 'h' ,
        'Nue' , 'NuMu' , 'NuTau' ,
        'Ve' , 'VMu' , 'VTau'
        """

        #   Check if channel is valid
        #   according to the dm process
        #   if not, raise value error
        if self._ew :
            if dm_channel not in ALLOWED_CHANNELS :

                raise ValueError(('\nChannel is not valid\n' +
                    'Valid options are {0}'.format(ALLOWED_CHANNELS)))

        else :
            if dm_channel not in ALLOWED_CHANNELSNOEW :

                raise ValueError(('\nChannel is not valid\n' +
                    'Valid options are {0}'.format(ALLOWED_CHANNELSNOEW)))
        #   Set channel
        self._channel = dm_channel

        #   Return
        return

    @property
    def z(self) :
        """
        Return the value of the redshift
        """

        #   Return
        return self._z

    @z.setter
    def z(self, src_z) :
        """
        Set the value of source's redshift
        """

        #   Check if redshift is positive
        #   If not, raise value error
        if src_z < 0 :

            raise ValueError(('\nRedshift with value {0} '.format(src_z) +
                'is not allowed. Redshift must be positive.'))

        #   Set the redshift
        self._z = src_z

        #   Return
        return

    @property
    def process(self) :
        """
        Return the process used to compute the DM spectrum
        """

        #   Return
        return self._process

    @process.setter
    def process(self, dm_process) :
        """
        Set the process used to compute the spectrum

        Parameters
        ----------
            dm_process : Process

        Options:
            - anna (Annihilation)
            - dec  (Decay)
        """

        #   Check if process is valid
        #   If not, raise value error
        if dm_process not in ALLOWED_PROCESSES :

            raise ValueError(('\nProcess {0} not valid\n'.format(dm_process) +
                'allowed options are: {0}'.format(ALLOWED_PROCESSES)))

        self._process = dm_process

        #   Return
        return

    @property
    def ebl_model(self) :
        """
        Return the EBL model used to compute attenuation
        """

        #   Return
        return self._eblmodel

    @ebl_model.setter
    def ebl_model(self, eblmodel) :
        """
        Set the EBL model used to compute attenuation

        Parameters
        ----------
            eblmodel : EBL model

        Options (taken from ebl-table project):
            'franceschini' , 'franceschini2017' ,
            'kneiske' , 'finke' ,
            'dominguez' , 'dominguez-upper' , 'dominguez-lower' ,
            'inuoe' , 'inuoe-low-pop3' , 'inuoe-up-pop3' ,
            'gilmore' , 'gilmore-fixed'
        """

        #   Check if eblmodel is valid
        #   If not, raise ValueError
        if eblmodel not in ALLOWED_EBLMODELS :

            raise ValueError(('\nEBL model {0} is '.format(eblmodel) +
                'not valid.\nOptions are: {0}'.format(ALLOWED_EBLMODELS)))

        #   Set ebl_model
        self._eblmodel = eblmodel

        #   Return
        return

    @property
    def hasEW(self) :
        """
        Return if the spectrum is computed using 
        EW corrections (hasEW=True) or not (hasEW=False)
        """

        #   Return
        return self._ew

    @hasEW.setter
    def hasEW(self, has_EW) :
        """
        Set if the spectrum is computed using
        EW corrections or not
        """

        #   Set hasEW
        self._ew = has_EW

        #   Return
        return

    @staticmethod
    def _dminterp(dm_channel, has_EW) :
        """
        Create DM interpolating function using tables
        from PPPC4DMID project.

        Interpolation is computed using interp2d from scipy.
        By default, only using linear interpolation

        Parameters
        ----------
            dm_channel : Channel
            has_EW     : using EW corrections (True)
                         or not (False)
        """

        #   Name of file to read data
        BASEDIR = os.path.abspath(os.path.join(os.path.dirname(__file__),'..'))
        data_path = os.path.join(BASEDIR, 'data/')

        fname = 'AtProduction'

        if not has_EW :

            fname += 'NoEW'

        fname += '_gammas.dat'
        fname  = os.path.join(data_path, fname)

        # print( 'Reading data from {:s}'.format( fname ) )

        #   Loading data
        data = np.genfromtxt(fname, names=True)

        #   Get unique values of masses and logxvals
        masses   = np.unique(data['mDM'])
        logxvals = np.unique(data['Log10x'])

        #   Matirx to extract data
        dndlogx  = np.zeros((logxvals.size, masses.size))

        #   Filling dndlogx matrix
        for mindex, mass in enumerate(masses) :

            mass    = int(mass)
            indices = np.where(data['mDM'] == mass)
            phis    = data[dm_channel][indices]

            for index , phi in enumerate(phis) :

                logxval                = logxvals[index]
                dndlogx[index][mindex] = phi

        #   Interpolating function using interp2d
        #   By default, I am only using linear interpolation
        dminterp = interp2d(masses, logxvals, dndlogx,
            kind='linear', fill_value=0.0)

        #   Return
        return dminterp

    @staticmethod
    def _decayinterp(dm_channel, has_EW) :
        """
        Create DM interpolating function using tables
        from PPPC4DMID project.

        Interpolation is computed using interp2d from scipy.
        By default, only using linear interpolation

        Parameters
        ----------
            dm_channel : Channel
            has_EW     : using EW corrections (True)
                         or not (False)
        """

        #   Name of file to read data
        BASEDIR = os.path.abspath(os.path.join(os.path.dirname(__file__),'..'))
        data_path = os.path.join(BASEDIR, 'data/')

        fname = 'AtProduction'

        if not has_EW :

            fname += 'NoEW'

        fname += '_gammas.dat'
        fname  = os.path.join(data_path, fname)

        # print( 'Reading data from {:s}'.format( fname ) )

        #   Loading data
        data = np.genfromtxt(fname, names=True)

        #   Get unique values of masses and logxvals
        masses   = np.unique(data['mDM'])
        logxvals = np.unique(data['Log10x'])

        #   Matirx to extract data
        dndlogx  = np.zeros((logxvals.size, masses.size))

        #   Filling dndlogx matrix
        for mindex, mass in enumerate(masses) :

            mass    = int(mass)
            indices = np.where(data['mDM'] == mass)
            phis    = data[dm_channel][indices]

            for index , phi in enumerate(phis) :

                logxval                = logxvals[index]
                dndlogx[index][mindex] = phi

        #   Now, the logxvals are going from 0 to 0.5
        #   I need to apply the rule x --> x/2
        logxvals -= np.log10(2.)

        #   Interpolating function using interp2d
        #   By default, I am only using linear interpolation
        dminterp = interp2d(masses, logxvals, dndlogx,
            kind='linear', fill_value=0.0)

        #   Return
        return dminterp

    @staticmethod
    def _ebl_atten(z, eblmodel, egev) :
        """
        Compute EBL attenuation, using ebl-table project
        By default, using kx=ky=2 as interpolation order

        Parameters
        ----------
            z        : Redshift
            eblmodel : EBL Model
            egev     : Energy (in GeV)

        Return
        ------
            atten    : attenuation
        """

        #   Check if redshift is greater than 1.e-3
        #   if not, then atten = 1

        if z >= 1.e-3 :

            #   Initiate instance of OptDepth
            tau   = OptDepth.readmodel(eblmodel)
            atten = np.exp(-1. * tau.opt_depth(z, egev*1.e-3))

        else :

            atten = 1.0

        #   Return
        return atten

    def spectra(self) :
        """
        Compute the dm_spectrum
        """

        #   Compute attenuation
        atten = self._ebl_atten(self._z, self._eblmodel, self._energy)

        if self._process == 'anna' :

            #   Get Interpolator
            dm_interp  = self._dminterp(self._channel, self._ew)

            #   Compute number of photons at energy self._energy
            xval    = self._energy / self._mass
            dndlogx = dm_interp(self._mass, np.log10(xval))
            dndlogx = dndlogx.flatten()
            dnde    = dndlogx / self._energy / np.log(10)

        elif self._process == 'decay' :

            #   Get interpolator
            dm_interp = self._decayinterp(self._channel, self._ew)

            #   Compute number of photons at energy self._energy
            xval    = self._energy / self._mass
            dndlogx = dm_interp(self._mass, np.log10(xval))
            dndlogx = dndlogx.flatten()
            dnde    = dndlogx / self._energy / np.log(10)

        #   Return spectrum attenuated by EBL
        dnde_atten = dnde * atten


        #   Return
        return dnde_atten

