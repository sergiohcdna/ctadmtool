###############################################
#   This is an incredible implementation of   #
#   Class decorators taken from:              #
#       Informit                              #
#       Mark Summerfield: Python Descriptors  #
#       https://www.informit.com/articles/    #
###############################################
import numpy as np

class GenericDescriptor :

    def __init__( self , getter , setter ) :

        self.getter = getter
        self.setter = setter

    def __get__( self , instance , owner=None ) :

        if instance is None :

            return self

        return self.getter( instance )

    def __set__(self , instance , value ) :

        return self.setter( instance , value )

#   Check if string is valid
def ValidString( attr_name , empty_allowed=True , options=None ) :

    def decorator( cls ) :

        name = "__" + attr_name

        def getter( self ) :

            return getattr( self , name )

        def setter( self , value ) :

            msg = '\n\t{0} must be a string'.format( attr_name )

            assert isinstance( value , str ) , msg

            if not empty_allowed and not value :

                raise ValueError( ( '{0} may not be empty'.format( attr_name ) ) )

            if options is not None and value not in options :

                msg = ( '\n\t{0} with value {1} '.format( attr_name , value ) +
                    'is not allowed. Valid options are:\n{0}'.format( options ) )

                raise ValueError( msg )

            setattr( self , name , value )

        setattr( cls , attr_name , GenericDescriptor( getter , setter ) )

        return cls

    return decorator

def ValidValue( attr_name , min_val=None , max_val=None ) :

    def decorator( cls ) :

        name = "__" + attr_name

        def getter( self ) :

            return getattr( self , name )

        def setter( self , value ) :

            msg = '\n\t{0} must be a float'.format( attr_name )

            assert isinstance( value , float ) , msg

            if min_val is not None and value < min_val :

                msg = ( '\n\t{0} with value {1} '.format( attr_name , value ) +
                    'is below the minimum value allowed: {0}'.format( min_val ) )

                raise ValueError( msg )

            if max_val is not None and value > max_val :

                msg = ( '\n\t{0} with value {1} '.format( attr_name , value ) +
                    'is above the maximum value allowed: {0}'.format( min_val ) )

                raise ValueError( msg )

            setattr( self , name , value )

        setattr( cls , attr_name , GenericDescriptor( getter , setter ) )

        return cls

    return decorator

def ValidnpArray( attr_name , min_val=None , max_val=None ) :

    def decorator( cls ) :

        name = "__" + attr_name

        def getter( self ) :

            return getattr( self , name )

        def setter( self , value ) :

            msg = '\n\t{0} must be a numpy ndarray'.format( attr_name )

            assert isinstance( value , np.ndarray ) , msg

            if min_val is not None and ( value < min_val ).all() :

                msg = ( '\n\t{0} has some value that  '.format( attr_name ) +
                    'is below the minimum value allowed: {0}'.format( min_val ) )

                raise ValueError( msg )

            if max_val is not None and ( value > max_val ).all() :

                msg = ( '\n\t{0} has some value that '.format( attr_name ) +
                    'is above the maximum value allowed: {0}'.format( min_val ) )

                raise ValueError( msg )

            setattr( self , name , value )

        setattr( cls , attr_name , GenericDescriptor( getter , setter ) )

        return cls

    return decorator

