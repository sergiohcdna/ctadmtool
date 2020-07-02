/*******************************************************************************
 *                 ONOFFAnalysis.cpp - Simulation and analysis                 *
 *-----------------------------------------------------------------------------*
 *    shkdna                                                                   *
 *-----------------------------------------------------------------------------*
 *                                                                             *
 * This program is free software: you can redistribute it and/or modify it     *
 * under the terms of the GNU General Public License.                          *
 *******************************************************************************
 *  @file ONOFFAnalysis.cpp                                                    *
 *  @brief Simulation and analysis of an observation from a XML-model template *
 *  @author shkdna                                                             *
 ******************************************************************************/

/*_____________Includes_______________________________________________________*/
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>
#include <stdlib.h>
#include <unistd.h>
#include "GammaLib.hpp"
#include "ctools.hpp"

/*******************************************************************************
 * @brief Simulation and analysis of a (wobble) observation                    *
 *                                                                             *
 *  This code show how to initialize an observation container,                 *
 *  run the simulation, perform a MLE fit, and get parameters                  *
 *  to compute the recovered spectrum from the observation.                    *
 *                                                                             *
 *  This code parses some argumments from the command line                     *
 *  but, this is in a simple and dummy way.                                    *
 *  At this moments, the options are:                                          *
 *      - XML template model                                                   *
 *      - Pointing right ascension in degrees                                  *
 *      - Pointin declination in degrees                                       *
 *      - Name of the source. This must be a valid name                        *
 *        of a model described in XML template                                 *
 *      - Duration of observation in hours                                     *
 *      - Instrument response function                                         *
 *      - Calibration data base                                                *
 *      - Dead time in percentage [ 0 , 1 ]                                    *
 *      - Minimum energy in TeV                                                *
 *      - Maximum energy in TeV                                                *
 *      - Radius of region of interest in degrees                              *
 ******************************************************************************/

static void show_usage(std::string name)
{
    // Print help message
    std::cerr << "Usage: " << name << " <args> CTA-related"
              << "Arguments:\n"
              << "\tXMLfile  : (string) XML file with model\n"
              << "\tpntRA    : (double) R.A. of pointing, in degrees\n"
              << "\tpntDeC   : (double) Dec. of pointing, in degrees\n"
              << "\tsource   : (string) Name of source of interest\n"
              << "\tduration : (double) Duration of observation, in hours\n"
              << "\tirf      : (string) Intrument response function\n"
              << "\tcaldb    : (string) Calibration database\n"
              << "\tdeadc    : (double) Dead time, [0,1]\n"
              << "\temin     : (double) Minimum energy, in TeV\n"
              << "\temax     : (double) Maximum energy, in TeV\n"
              << "\trad      : (double) Radius of region of interest\n"
              << std::endl ;
}

/********************************************
 *  Function to convert string to double    *
 *  Required parameters:                    *
 *      - strnumber: (string) number        *
 *******************************************/
double strTodouble( std::string strnumber )
{
    // Declare double var 'number'
    double number ;
    // Declare stringstream var 'ssnumber' with strnumber
    std::stringstream ssnumber( strnumber ) ;

    // Conversion from string to double
    // If the conversion fails, exit
    if ( ! ( ssnumber >> number ) ) {
        std::cout << "Error converting string to double" << std::endl ;
        exit( EXIT_FAILURE ) ;
    }

    // Return double number
    return number ;
}

/********************************************
 *  Print number of observation id_obs      *
 *  Required parameters:                    *
 *      -id_obs: (int) observation number   *
 *******************************************/
void print_obs( int id_obs )
{
    // Print message
    std::cout << "******************************\n"
              << "\tObservation " << id_obs << "\n"
              << "******************************"
              << std::endl ;

}

/************************************************************
 *      Create empty CTA observation container              *
 *      Required parameters:                                *
 *          - pntdir: (GSkyDir) Direction to pointing       *
 *          - caldb : (string)  CTA-calibration database    *
 *          - id    : (string)  Identifier                  *
 *          - tstart: (double)  Start time, in seconds      *
 *          - tinter: (double)  Duration, in seconds        *
 *          - deadc : (double)  Dead time, [0,1]            *
 *          - emin  : (double)  Minimum Energy, in TeV      *
 *          - emax  : (double)  Maximum Energy, in TeV      *
 *          - rad   : (double)  Radius of ROI               *
 *          - irf   : (string)  CTA-IRF                     *
 ***********************************************************/
GCTAObservation single_obs( GSkyDir pntdir , 
                            std::string irf , std::string caldb , 
                            double tstart = 0.0 , double tinter = 1800.0 , 
                            double deadc = 0.95 ,
                            double emin = 0.01 , double emax = 100.0 , 
                            double rad = 5.0 , 
                            std::string id = "001" )
{
    
    // Declare GCTAObservation container
    GCTAObservation obs_cta ;

    // Open cta-associated calibration data base
    GCaldb db ;
    db.open( "cta" , caldb ) ;

    // CTA pointing using GSkyDir object
    GCTAPointing pnt ;
    pnt.dir( pntdir ) ;
    obs_cta.pointing( pnt ) ;

    // Setting CTA region of interest with centre at pntdir and radius rad
    GCTARoi roi ;
    GCTAInstDir instdir ;
    instdir.dir( pntdir ) ;
    roi.centre( instdir ) ;
    roi.radius( rad ) ;

    // Setting start and duration of observation
    GGti gti ;
    gti.append( GTime( tstart ) , GTime( tstart + tinter ) ) ;

    // Setting energy bounds in TeV
    GEbounds ebounds( GEnergy( emin , "TeV" ) , GEnergy( emax , "TeV" ) ) ;

    // GCTAEvents container for roi, gti, and ebounds objects
    GCTAEventList events ;
    events.roi( roi ) ;
    events.gti( gti ) ;
    events.ebounds( ebounds ) ;

    // Setting events object to obs_cta container
    obs_cta.events( events ) ;

    // Setting instrument response function irf from db object
    obs_cta.response( irf , db ) ;

    // Setting time-related options in observation
    obs_cta.ontime( tinter ) ;
    obs_cta.livetime( tinter * deadc ) ;
    obs_cta.deadc( deadc ) ;
    
    //Setting identifier for observation
    obs_cta.id( id ) ;

    // Return
    return obs_cta ;

}

int main( int argc , char* argv[] )
{

    // Check the number of arguments in the command line
    if ( ( argc != 12 ) ) {
        
        // If number of arguments is different from 6
        // print message and exit
        std::cout << "You must check the number of arguments" << std::endl ;
        show_usage( argv[ 0 ] ) ;
        
        //Return
        return 1 ;

    } else {
        
        //Print message
        std::cout << "Nice\nParsing arguments..." << std::endl ;
    }

    // Getting the arguments from command-line
    
    // XML file
    std::string xmlfile = argv[ 1 ] ;
    
    // pointing right ascension
    double pntRA        = strTodouble( argv[ 2 ] ) ;
    
    // pointing declination
    double pntDeC       = strTodouble( argv[ 3 ] ) ;
    
    // Name of the source
    std::string source  = argv[ 4 ] ;
    
    // Duration 
    double duration     = strTodouble( argv[ 5 ] ) * 3600 ;

    // Print right ascensiond and declination got from command-line
    std::cout << std::setprecision( 4 ) << "Pointing at ( RA , Dec ) =  ( " 
              << pntRA << " , " << pntDeC << " )" << std::endl ;
    
    // Declare GModels object models with xmlfile template
    GModels models( xmlfile ) ;

    // Check if models contains source
    // If not, print message and exit
    if ( ! ( models.contains( source ) ) ) {
        std::cout << "Sorry, the model with name " << source 
                  << "is not available in the xml template"
                  << std:: endl ;
        exit( EXIT_FAILURE ) ;
    }

    // Declare GSkyDir objecto for the pointing direction
    GSkyDir pntdir ;
    
    // Set Pointing to pntRA and pntDeC
    pntdir.radec_deg( pntRA , pntDeC ) ;

    // Get GModel source from models
    GModel* src_roi = models[ source ] ;
    
    // Declare GSkyDir for the source position
    GSkyDir srcdir ;
    
    // Get source RA
    double srcRA  = src_roi -> at( 0 ).value() ;
    
    // Get source declination
    double srcDeC = src_roi -> at( 1 ).value() ;
    
    // Set Direction to source position
    srcdir.radec_deg( srcRA , srcDeC ) ;

    // Print number of observation
    print_obs( 10 ) ;

    // Declare GObservations container
    GObservations obslist ;

    // CTA observation container
    GCTAObservation thisobs = single_obs( pntdir , 
                                          "North_z20_50h" , "prod3b-v2" ,
                                          0.0 , duration , 0.95 ,
                                          0.1 , 10.0 , 3.0 , "0001"  ) ;

    // Return
    return 0 ;

}

