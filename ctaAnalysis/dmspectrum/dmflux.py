import numpy as np
from ctaAnalysis.dmspectrum.dmspectra import dmspectrum
from ctaAnalysis.tools.misc import ValidValue , ValidString

import warnings

ALLOWED_FERMIONS = ( 'Majorana' , 'Dirac' )

@ValidString( '_delta' , empty_allowed=False , options=ALLOWED_FERMIONS )
@ValidValue( '_emin' , min_val=5.e-9 )
@ValidValue( '_emax' , max_val=1.e+5 )
@ValidValue( '_sigmav' , min_val=1.e-35 )
# @ValidValue( "_jfactor" , min_val=1.e+5 )
class dmflux_anna( dmspectrum ) :
    """
    Class to compute the flux generated
    from annihilation of dark matter
    particles.
    dmflux is a derived class from dmspectrum.
    """

    #   Init
    def __init__( self , sigmav , jfactor , dm_mass , emin , emax , channel ,\
        z , delta='Majorana' , npoints=100 , eblmod='franceschini2017' , has_EW=True ) :
        """
        Initialize the dmflux_anna class

        Parameters:
        ----------
            sigmav  : Annihilation cross-section (in cm**3/s)
            jfactor : Astrophysical factor in (GeV**2/cm**5)
            dm_mass : Mass of dark matter candidate
            emin    : Minimum energy to compute the flux (in GeV)
            emax    : Maximum energy to compute the flux (in GeV)
            channel : Annihilation channel
            z       : Redshift to which photons are emitted
            delta   : Parameter to describe if dark matter candidate
                      is a Majorana (delta=2) fermion or a
                      Dirac (delta=4) fermion
            npoints : Number of points to compute the flux
            eblmod  : EBL model used to compute attenuation
            has_EW  : Boolean to indicate if EW corrections
                      are taken into account or not
        """

        #   Additionally to check that emin and emax have valid values
        #   I check that emin < emax
        if emin > emax :

            msg = ( '\nI found that Minimum energy {0} '.format( emin ) +
                'is greater than Maximum energy {0}.\n'.format( emax ) +
                'Changing the order...' )

            warnings.warn( msg , RuntimeWarning )
            e_min = emax
            e_max = emin

        else :

            e_min = emin
            e_max = emax

        #   Array to store energy values to compute the dmflux
        e_array = np.logspace( np.log10( e_min ) , np.log10( e_max ) , npoints )

        #   By default, I am using annihilation (anna)
        process = 'anna'

        #   Initialize dm_spectrum super class
        super().__init__( dm_mass , e_array , channel , z , \
            process=process , eblmod=eblmod , has_EW=has_EW )

        #   Initialize parameters of dmflux_ana class
        self._sigmav  = sigmav
        self._jfactor = jfactor
        self._emin    = e_min
        self._emax    = e_max
        self._delta   = delta

        #   Return
        return

    @property
    def sigmav( self ) :
        """
        Return value of the annihilation cross-section
        used to compute the flux
        """

        #   Return
        return self._sigmav

    @sigmav.setter
    def sigmav( self , sigmav ) :
        """
        Set the value of Annihilation cross-section (in cm**3/s)
        used to compute the flux

        Parameters
        ----------
            sigmav : Annihilation cross-section (cm**3/s)
        """

        #   Check that sigmav is greater than 1.e-35
        if sigmav < 1.e-35 :

            raise ValueError( ( '\nValue of annihilation cross-section ' +
                ' must be greater than 1.e-35.\n' +
                'This is just to avoid possible round errors' ) )

        #   Set sigmav
        self._sigmav = sigmav

        #   Return
        return

    @property
    def jfactor( self ) :
        """
        Return the value of the Astrophysical factor
        used to compute the flux
        """

        #   Return
        return self._jfactor

    @jfactor.setter
    def jfactor( self , jfactor ) :
        """
        Set the value of the Astrophysical factor (GeV**2/cm**5)
        to compute the dm flux

        Parameters
        ----------
            jfactor : Astrophysical factor J (GeV**2/cm**5)
        """

        #   Set the jfactor
        self._jfactor = jfactor

        #   Return
        return

    @property
    def emin( self ) :
        """
        Return Minimum value of energy (GeV) used to compute
        the dm flux
        """

        #   Return
        return self._emin

    @emin.setter
    def emin( self , emin ) :
        """
        Set the value of minimum energy (GeV) used to compute
        the dm flux
        """

        #   Just check that the minimum energy is greater than
        #   5.e-9.
        #   In reality, at this point, is factible to check
        #   that emin is greater than mass * 1.e-9
        #   By the way, the dminterpolator put 1.e-40
        #   in values outside the range of interpolation
        if emin < 5.e-9 :

            raise ValueError( ( '\nMinimum energy {0} GeV '.format( emin ) +
                'is below the allowed value (5.e-9GeV)') )

        #   Set minimum energy
        self._emin = emin

    @property
    def emax( self ) :
        """
        Return Maximum value of energy (GeV) used to compute
        the dm flux
        """

        #   Return
        return self._emax

    @emax.setter
    def emax( self , emax ) :
        """
        Set the value of minimum energy (GeV) used to compute
        the dm flux
        """

        #   Just check that the minimum energy is greater than
        #   1.e+5.
        #   In reality, at this point, is factible to check
        #   that emax is lower than mass
        #   By the way, the dminterpolator put 1.e-40
        #   in values outside the range of interpolation
        if emax > 1.e+5 :

            raise ValueError( ( '\nMaximum energy {0} GeV '.format( emax ) +
                'is above the allowed value (1.e+5GeV)') )

        #   Set minimum energy
        self._emax = emax

    @property
    def delta( self ) :
        """
        Return what kind of dark matter particle is
        used to compute the dm flux
        """

        #   Return
        return self._delta

    @delta.setter
    def delta( self , delta ) :
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

            raise ValueError( ( '\nKind of Dark matter particle not ' +
                'supported.\nOptions are:{0}'.format( ALLOWED_FERMIONS ) ) )

        #   Set minimum energy
        self._delta = delta

    @staticmethod
    def _ppfactor( sigmav , mass , delta ) :
        """
        Compute the Particle physics factor for the dm flux

        Parameters
        ----------
            sigmav : Value of annihilation cross-section (cm**3/s)
            mass   : Mass of dark matter particles (GeV)
            delta  : String to indicate if dark matter is a
                     Majorana or Dirac fermion

        Return
        ------
            ppfactor : (cm**3/GeV**2/s)
        """

        #   Check delta
        if delta == 'Majorana' :

            d = 2

        elif delta == 'Dirac' :
            d = 4

        #   Compute ppfactor
        ppfactor = sigmav / d / 4 / np.pi / np.power( mass , 2 )

        return ppfactor

    def flux( self ) :
        """
        Compute the DM flux
        """

        #   Get dnde
        dnde     = self.spectra()

        #   Get ppfactor
        ppfactor = self._ppfactor( self._sigmav , self._mass , self._delta )

        #   Get the flux
        dphide   = ppfactor * dnde * self._jfactor

        return dphide
