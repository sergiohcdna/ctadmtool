#=======================================#
#   Script to create gammalib-models    #
#   in runtime mode (during execution)  #
#                                       #
#       Created: Dec-2020               #
#                Sergio Hernández       #
#=======================================#

import gammalib

from ctaAnalysis.tools.misc import ValidValue , ValidString


MODEL_TYPES   = ( 'PointSource' , 'DiffuseSource' )
SPATIAL_TYPE  = ( 'PointSource' , 'DiffuseMap' )
SPECTRAL_TYPE = ( 'FileFunction' )

class dm_spectral() :

    def __init__( self , fname , minval=0.0 , maxval=1.e+3 , norm_free=True ) :

        self._file          = fname
        self._min_val       = minval
        self._max_val       = maxval
        self._is_norm_free  = norm_free

    @property
    def file( self ) :

        return self._file

    @file.setter
    def file( fname ) :

        self._file = fname

        return

    @property
    def min( self ) :

        return self._min_val

    @min.setter
    def min( min_val ) :

        self._min_val = min_val

        return

    @property
    def max( self ) :

        return self._max_val

    @max.setter
    def max( max_val ) :

        self._max_val = max_val

        return

    @property
    def norm_free( self ) :

        return self._is_norm_free

    @norm_free.setter
    def norm_free( norm_free ) :

        self._is_norm_free = norm_free

        return

    @staticmethod
    def _set_ffspectrum( file , min_val=0.0 ,
        max_val=1.e+8 , norm_free=True ) :

        #   Create GModelSpectralFunc
        spectrum = gammalib.GModelSpectralFunc( file , 1.0 )

        #   Set min and max value
        norm     = spectrum[ 'Normalization' ]
        norm.factor_min( min_val )
        norm.factor_max( max_val )

        return spectrum

    def dmspectrum( self ) :

        dmspec = self._set_ffspectrum( self._file , self._min_val ,
            self._max_val , self._is_norm_free )

        return dmspec

class dm_pointsource() :

    def __init__( self , srcra , srcdec , ra_free=False , dec_free=False ) :

        self._ra           = srcra
        self._dec          = srcdec
        self._is_ra_free   = ra_free
        self._is_dec_free  = dec_free

    @property
    def ra( self ) :

        return self._ra

    @ra.setter
    def ra( srcra ) :

        self._ra = srcra

        return

    @property
    def dec( self ) :

        return self._dec

    @dec.setter
    def dec( srcdec ) :

        self._dec = srcdec

        return

    @property
    def ra_free( self ) :

        return self._is_ra_free

    @ra_free.setter
    def ra_free( ra_free ) :

        self._is_ra_free = ra_free

        return

    @property
    def dec_free( self ) :

        return self._is_dec_free

    @dec_free.setter
    def dec_free( dec_free ) :

        self._is_dec_free

        return

    @staticmethod
    def _set_ps( ra , dec , ra_free=False , dec_free=False ) :

        #   Creating SkyDir object
        srcdir  = gammalib.GSkyDir()
        srcdir.radec_deg( ra , dec )

        #   Create PointSource model
        spatial = gammalib.GModelSpatialPointSource( srcdir )

        #   check if it's needed to free, then set to free
        if ra_free :

            spatial[ 'RA' ].free()

        if dec_free :

            spatial[ 'DEC' ].free()

        return spatial

    def spatial( self ) :

        dmspat = self._set_ps( self._ra , self._dec ,
            self._is_ra_free , self._is_dec_free )

        return dmspat

class dm_extended() :

    def __init__( self , map_fits ) :

        self._fits         = map_fits

    @property
    def fits_file( self ) :

        return self._fits

    @fits_file.setter
    def fits_file( map_fits ) :

        self._fits = map_fits

        return

    def spatial( self ) :

        dmspat = gammalib.GModelSpatialDiffuseMap( self._fits )

        return dmspat

class DMModel() :

    def __init__( self , name , spec , spat , tscalc=True ) :

        self._name     = name
        self._tscalc   = tscalc
        self._spectral = spec
        self._spatial  = spat

    @staticmethod
    def _set_model( name , spectrum , spatial , compute_ts=True ) :

        source = gammalib.GModelSky( spatial , spectrum )

        #   Set name of the source
        source.name( name )

        #   Set computation of TS
        source.tscalc( compute_ts )

        return source

    def model( self ) :

        src  = self._set_model( self._name , self._spectral ,
            self._spatial , self._tscalc )

        return src
