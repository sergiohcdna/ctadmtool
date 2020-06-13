#!/usr/bin/python

########################################################
#       Definition of some management functions
#       for files, directories, and screen lines
########################################################
#
#       shkdna, June 12, 2020
#
########################################################

import os
import sys

#   CURSOR_UP_ONE     : Move the cursor one line above
#   ERASE_LINE        : keyword command to erase current line

CURSOR_UP_ONE = '\x1b[1A'
ERASE_LINE = '\x1b[2K'

def delete_last_lines( n=1 ) :
    '''Erase the last n lines from the terminal screen.
    Parameters:
        (int)   n   Number of lines'''
    for _ in range( n ) :
        sys.stdout.write( CURSOR_UP_ONE )
        sys.stdout.write( ERASE_LINE )

#   Wait until run the next python task.
def Wait( waiting ) :
    '''Wait x time until run the next python block code. You can use it to control
    the code flow and then print some useful messages.
    Parameters:
        (int)   waiting     Time'''
    for i in range( waiting ) :
        print( '\t%d    %s' % ( waiting - i , ( waiting - i ) * '*' ) )
        delete_last_lines( 1 )
        time.sleep( 1 )
    print( '\n' )

########################################################
#   The following functions are for the management
#   of directories and files
########################################################

def createname( mypath , name ) :
    '''Return a string with a proper name'''
    return os.path.join( mypath , name )

def checkPathandExit( thispath ) :
    '''Check if the path to file or directory'''
    if not os.path.exists( thispath ) :
        print( '\t\tError. You are trying to pass a file or directory that does not exists.' )
        sys.exit( 1 )

def checkDir( thispath ) :
    '''This function allow to check if the path belongs to a directory or any other file-type.
    The functions is useful to check if any out-put path corresponds to a valid directory'''
    if os.path.exists( thispath ) :
        if not os.path.isdir( thispath ) :
            print( 'Sorry, you are trying to pass any other file-type as the outpath directory' )
            sys.exit( 1 )
        else :
            print( 'Good, this is a directory. Nothing to do' )
    else :
        print( 'It seems like the path does not exists.\nCreating directory' )
        os.makedirs( thispath )



