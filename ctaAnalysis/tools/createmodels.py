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

class dm_spectral_xml() :

    def __init__( self , spec_type , fname ,
        minval=0.0 , maxval=1.e+3 , norm_free=True ) :

        self._spectral_type = spec_type
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

        spectrum = gammalib.GXmlElement( ( 'spectrum type="FileFunction" ' +
            'file="{0}"'.format( file ) ) )
        spectrum.append( ( 'parameter scale="1.0" name="Normalization" ' +
            'min="{:.2f}" max="{:.2e}" '.format( min_val , max_val ) +
            'value="1.0" free="{0}"'.format( int( norm_free ) ) ) )

        return spectrum

    def xml_spectrum( self ) :

        xml = self._set_ffspectrum( self._file , self._min_val ,
            self._max_val , self._is_norm_free )

        return xml

class dm_pointsource_xml() :

    def __init__( self , srcra , srcdec , ra_free=False , dec_free=False ) :

        self._ra           = srcra
        self._dec          = srcdec
        self._is_ra_free   = ra_free
        self._is_dec_free  = dec_free
        self._spatial_type = 'PointSource'

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

    @property
    def spatial_type( self ) :

        return self._spatial_type

    @staticmethod
    def _set_ps( ra , dec , ra_free=False , dec_free=False ) :

        spatial = gammalib.GXmlElement( 'spatialModel type="PointSource"' )
        spatial.append( ( 'parameter name="RA" scale="1.0" ' +
            'value="{:2f}" min="-360" '.format( ra ) +
            'max="360" free="{0}"'.format( int( ra_free ) ) ) )
        spatial.append( ( 'parameter name="DEC" scale="1.0" ' +
            'value="{:2f}" min="-90" '.format( dec ) +
            'max="90" free="{0}"'.format( int( dec_free ) ) ) )

        return spatial

    def xml_spatial( self ) :

        xml = self._set_ps( self._ra , self._dec ,
            self._is_ra_free , self._is_dec_free )

        return xml

class dm_extended_xml() :

    def __init__( self , map_fits ) :

        self._fits         = map_fits
        self._spatial_type = 'DiffuseMap'

    @property
    def fits_file( self ) :

        return self._fits

    @fits_file.setter
    def fits_file( map_fits ) :

        self._fits = map_fits

        return

    @property
    def spatial_type( self ) :

        return self._spatial_type

    @staticmethod
    def _set_diffmap( map_fits ) :

        spatial = gammalib.GXmlElement( ( 'spatialModel type="DiffuseMap" ' +
            'file="{0}" normalize="1"'.format( map_fits ) ) )
        spatial.append( ( 'parameter scale="1" name="Prefactor" ' +
            'min="0.001" max="1.e+5" value="1" free="0"/' ) )

        return spatial


    def xml_spatial( self ) :

        xml = self._set_diffmap( map_fits )

        return xml

# @ValidValue( '_dec' , min_val=-90 , max_val=90 )
# @ValidValue( '_ra' , min_val=-360 , max_val=360 )
# @ValidString( '_file' , empty_allowed=False )
# @ValidString( '_modtype' , empty_allowed=False , options=MODEL_TYPES )
# @ValidString( '_name' , empty_allowed=False )
class DMModel() :

    def __init__( self , name , modtype ,
        xml_spec , xml_spat , tscalc=True ) :

        self._name         = name
        self._modtype      = modtype
        self._tscalc       = tscalc
        self._xml_spectral = xml_spec

        if self._modtype == 'PointSource' :

            if not xml_spat.spatial_type == self._modtype :

                raise AssertionError( ( 'Spatial component ' +
                    'must be of the same type as the source ' +
                    'model: {0}'.format( self._modtype ) ) )

        if self._modtype == 'DiffuseSource' :

            if not 'Diffuse' in xml_spat.spatial_type :

                raise AssertionError( ( 'Spatial component ' +
                    'must be of the same type as the source ' +
                    'model: {0}'.format( self._modtype ) ) )

        self._xml_spatial  = xml_spat

    @staticmethod
    def _set_model( name , modtype , spectrum , spatial , compute_ts=True ) :

        srcxml = gammalib.GXmlElement( ( 'source name="{0}" '.format( name ) +
            'type="{0}" '.format( modtype ) +
            'tscalc="{0}"'.format( int( compute_ts ) ) ) )
        spec   = spectrum.xml_spectrum()
        spat   = spatial.xml_spatial()

        srcxml.append( spec )
        srcxml.append( spat )

        source = gammalib.GModelSky( srcxml )

        return source

    def model( self ) :

        src  = self._set_model( self._name , self._modtype ,
            self._xml_spectral , self._xml_spatial , self._tscalc )

        return src
