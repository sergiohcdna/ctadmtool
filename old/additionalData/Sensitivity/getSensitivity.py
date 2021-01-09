#! /usr/bin/env python

import cscripts
import gammalib
import ctools

import argparse

options = argparse.ArgumentParser( description='This script compute \
	the DM limits using Cirelli et al. Spectrum' )
source = options.add_argument_group( 'Source' , \
	'All the relevant information about the source is passed in the inmodel\
	argument, where a XML file is passed describing parameters for the\
	observation, as position of the source, name, and spatial and spectral\
	models for the source. It is assumed that the file provides information about\
	the background model used in the observation. For more details, the user\
	is referred to the gammalib and ctools documentations :)' )
source.add_argument( '--inmodel' ,\
	help='File with xml model' ,\
	type=str ,\
	required=True ,
	metavar='PerseusModel.xml' )
source.add_argument( '--srcname' ,\
	help='Name of the source in the xml model' ,\
	type=str ,\
	required=True ,
	metavar='Perseus' )
source.add_argument( '--rad' ,\
	help='Radius of the Region of Interest (degrees)' ,\
	type=float ,\
	default=1.0 ,\
	metavar='10.0' )
source.add_argument( '--caldb' ,\
	help='Database production for the IRF file. \
		Options are: [ prod2 , prod3b-v1 ] ' ,\
	type=str ,\
	default='prod3b-v1' ,\
	metavar='prod3b-v1' ,\
	choices=[ 'prod2' , 'prod3b-v1' ] )
source.add_argument( '--irf' ,\
	help='Intrument Response Function for CTA' ,\
	type=str ,\
	default='North_z20_average_5h' ,\
	metavar='North_z20_average_50h' )
source.add_argument( '--hours' ,\
	help='Time for simulation of observation (in hours)' ,\
	type=float ,\
	default=5.0 ,\
	metavar='50.0' )
source.add_argument( '--enumbins' ,\
	help='Number of energy bins' ,\
	type=int ,\
	default=10 ,\
	metavar='20' )
source.add_argument( '--emin' ,\
	help='Minimum energy for events (TeV)' ,\
	type=float ,\
	default=0.05 ,\
	metavar='0.01' )
source.add_argument( '--emax' ,\
	help='Maximum energy for events (TeV)' ,\
	type=float ,\
	default=10.0 ,\
	metavar='100.0' )

args = options.parse_args()


sensitivity = cscripts.cssens()


sensitivity[ 'inmodel' ]  = args.inmodel
sensitivity[ 'srcname' ]  = args.srcname
sensitivity[ 'caldb' ]    = args.caldb
sensitivity[ 'irf' ]      = args.irf
sensitivity[ 'outfile' ]  = args.srcname + 'Sens{:.1f}h.txt'.format( args.hours )
sensitivity[ 'duration' ] = args.hours * 3600.
sensitivity[ 'rad' ]      = args.rad
sensitivity[ 'emin' ]     = args.emin
sensitivity[ 'emax' ]     = args.emax
sensitivity[ 'bins' ]     = args.enumbins
sensitivity[ 'enumbins' ] = args.enumbins
sensitivity[ 'nthreads' ] = 2

sensitivity.execute() 


