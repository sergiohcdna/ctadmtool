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
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <ctime>
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
 *      - Output directory                                                     *
 *      - Obs-list XML file name                                               *
 *      - Number of threads used to parallelization                            *
 ******************************************************************************/

static void show_usage( std::string name )
{
    // Print help message
    std::cerr << "Usage: " << name << " <args> CTA-related"
              << "Arguments:\n"
              << "\tXMLfile  : (string) XML file with model\n"
              << "\tpntRA    : (double) R.A. of pointing, in degrees\n"
              << "\tpntDeC   : (double) Dec. of pointing, in degrees\n"
              << "\tsource   : (string) Name of source of interest\n"
              << "\tinterval : (double) Duration of observation, in hours\n"
              << "\tirf      : (string) Intrument response function\n"
              << "\tcaldb    : (string) Calibration database\n"
              << "\tdeadc    : (double) Dead time, [0,1]\n"
              << "\temin     : (double) Minimum energy, in TeV\n"
              << "\temax     : (double) Maximum energy, in TeV\n"
              << "\tenbins   : (int)    Number of energy bins\n"
              << "\trad      : (double) Radius of region of interest\n"
              << "\toutpath  : (string) Path to output directory\n"
              << "\tobsname  : (string) Name of the XML-file with all the\n"
              << "\t                    details about CTAObservation\n"
              << "\tthreads  : (int)    Number of threads used to\n"
              << "\t                    parallelization when possible\n"
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
 *  Function to convert string to double    *
 *  Required parameters:                    *
 *      - strnumber: (string) number        *
 *******************************************/
double strToint( std::string strnumber )
{
    // Declare double var 'number'
    int number ;
    // Declare stringstream var 'ssnumber' with strnumber
    std::stringstream ssnumber( strnumber ) ;

    // Conversion from string to int
    // If the conversion fails, exit
    if ( ! ( ssnumber >> number ) ) {
        std::cout << "Error converting string to double" << std::endl ;
        exit( EXIT_FAILURE ) ;
    }

    // Return double number
    return number ;
}


/************************************************************
 *      Function to check if path-to-dir exists.            *
 *      If path does not exist, create directory.           *
 *      If path-to-dir exists and it is not a directory     *
 *      then exit. But if path-to-dir exists and is         *
 *      a directory, then print message with path.          *
 *      Required parameters:                                *
 *          - dirname : (string) path to directory          *
 ***********************************************************/
void checkdir( std::string dirname )
{
    // Declare stat object to manage file
    struct stat sb ;
    
    // Check if path-to directory exists and if is a directory
    if ( stat( dirname.c_str() , &sb ) != 0 ) {
        std::cout << "Cannot access " << dirname 
                  << "\nCreating directory..."
                  << std::endl ;
        int check = mkdir( dirname.c_str() , 0777 ) ;
        if ( !check ) {
            std::cout << "Directory" << dirname << " created successfully"
                      << ":)" << std::endl ;
        } else {
            std::cout << "There was problems to create directory"
                      << ":(" << std::endl ;
            exit( EXIT_FAILURE ) ;
        }
    } else if ( sb.st_mode & S_IFDIR ) {
        std::cout << "Good!\n\tI love it when I have not to work\n"
                  << dirname << " is a directory"
                  << std::endl ;
    } else {
        std::cout << dirname << " is not a directory"
                  << "Exit\nExit\nExit" << std::endl ;
        exit( EXIT_FAILURE ) ;
    }
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

/********************************************
 *      UTC to MJD (in seconds) converter   *
 *      Required parameters:                *
 *          - gmt_time: (struct tm) Time    *
 *******************************************/

unsigned long int get_mjd( std::tm *gmt_time )
{

    //  Get year, month, day, hour, minutes and seconds
    //  from gmt_time struct
    int year = gmt_time -> tm_year ;
    int month = gmt_time -> tm_mon ;
    int day = gmt_time -> tm_mday ;
    int hour = gmt_time -> tm_hour ;
    int mins = gmt_time -> tm_min ;
    int secs = gmt_time -> tm_sec ;

    //  Compute total seconds from hours, minutes and seconds
    int tsecs = hour * 3600 + mins * 60 + secs ;

    //  MJD Calculation
    int l = 0 ;
    unsigned long int mjd = 0 ;

    if ( ( month == 1 ) || ( month == 2 ) ) {
        l = 1 ;
    }

    mjd = 14956 + day + ( int )( ( year - l ) * 365.25 )
          + ( int )( ( month + 1 + l * 12 ) * 30.60001 ) ;

    //  Converting to seconds
    mjd *= 86400 ;
    mjd += tsecs ;

    //  Return
    return mjd ;

}

/************************************************************
 *      Create empty CTA observation container              *
 *      Required parameters:                                *
 *          - pntdir : (GSkyDir) Direction to pointing      *
 *          - caldb  : (string)  CTA-calibration database   *
 *          - id     : (string)  Identifier                 *
 *          - tstart : (double)  Start time, in seconds     *
 *          - tinter : (double)  Duration, in seconds       *
 *          - deadc  : (double)  Dead time, [0,1]           *
 *          - emin   : (double)  Minimum Energy, in TeV     *
 *          - emax   : (double)  Maximum Energy, in TeV     *
 *          - enbins : (int)     Number of energy bins      *
 *          - rad    : (double)  Radius of ROI              *
 *          - irf    : (string)  CTA-IRF                    *
 ***********************************************************/
GCTAObservation single_obs( GSkyDir pntdir , 
                            std::string irf , std::string caldb , 
                            double tstart = 0.0 , double tinter = 1800.0 , 
                            double deadc = 0.95 ,
                            double emin = 0.01 , double emax = 100.0 , 
                            int enbins = 20 , double rad = 5.0 , 
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
    if ( ( argc != 16 ) ) {
        
        // If number of arguments is different from 15
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
    std::string xmlfile = argv[ 1 ] ;
    double pntRA        = strTodouble( argv[ 2 ] ) ;
    double pntDeC       = strTodouble( argv[ 3 ] ) ;
    std::string source  = argv[ 4 ] ;
    double interval     = strTodouble( argv[ 5 ] ) * 3600 ;
    std::string irf     = argv[ 6 ] ;
    std::string caldb   = argv[ 7 ] ;
    double deadc        = strTodouble( argv[ 8 ] ) ;
    double emin         = strTodouble( argv[ 9 ] ) ;
    double emax         = strTodouble( argv[ 10 ] ) ;
    int enbins          = strToint( argv[ 11 ] ) ;
    double rad          = strTodouble( argv[ 12 ] ) ;
    std::string outpath = argv[ 13 ] ;
    std::string obsname = argv[ 14 ] ;
    int threads         = strToint( argv[ 15 ] ) ;

    std::time_t  thistime     = std::time( 0 ) ;
    std:tm *gmt_time          = std::gmtime( &thistime ) ;
    unsigned long int thismjd = get_mjd( gmt_time ) ;

    // Check if output directory exists. 
    // If output path is not valid, then exit
    checkdir( outpath ) ;

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

    //  Check if distance between pointin and source centre
    //  is greater tahn radius of source ROI. The last is
    //  taking to be a 1 degree, but in the future, the
    //  user will be able to pass as an option. :)
    
    double offset = pntdir.dist_deg( srcdir ) ;

    if ( offset < 1.10 ) {
        std::cout << "Sorry, the distance between the centers"
                  << " of the source and the pointing are so"
                  << " close to get background regions"
                  << std::endl ;
        exit( EXIT_FAILURE ) ;
    } else if ( offset > 10.0 ) {
        std::cout << "Sorry, the distance between the centers"
                  << " of the source and the pointing are so"
                  << " far to compute background regions"
                  << std::endl ;
    }

    // Print number of observation
    print_obs( 10 ) ;

    // Declare GObservations container
    GObservations obslist ;

    // CTA observation container
    GCTAObservation myobs = single_obs( pntdir , irf , caldb , 0.0 , interval ,
                                        deadc , emin , emax , enbins , 
                                        rad , "0001"  ) ;

    // Append GCTAObservation myobs to observation_list
    obslist.append( myobs ) ;
    // Save XML-definition file
    obsname = outpath + "/" + obsname ;
    obslist.save( obsname ) ;

    /****************************************************************
     *      In Order to pass all the arguments to ctobssim tool     *
     *      and to avoid asking to the user for all the arguments   *
     *      then, I am using the contruction method to pass an      *
     *      array with all the options.                             *
     *      As usual, the options are passed via:                   *
     *          option=val_option                                   *
     *      For detailed information, see ctobssim help.            *
     ***************************************************************/
    std::string thislist = "inobs=" + obsname ;
    std::string inmodel  = "inmodel=" + xmlfile ;
    std::string scaldb   = "caldb=" + caldb ;
    std::string sirf     = "irf=" + irf ;
    std::string edisp    = "edisp=no" ;
    std::string chatter  = "chatter=4" ;
    std::string outevent = "outevents=" + outpath 
                           + "/" + source + "Events.fits" ;
    std::string seed     = "seed=" + std::to_string( ( int ) std::time( 0 ) ) ;
    std::string sra      = "ra=" + std::to_string( pntRA ) ;
    std::string sdec     = "dec=" + std::to_string( pntDeC ) ;
    std::string srad     = "rad=" + std::to_string( rad ) ;
    std::string tmin     = "tmin=" + std::to_string( thismjd ) ;
    std::string tmax     = "tmax=" + std::to_string( thismjd + interval ) ;
    std::string semin    = "emin=" + std::to_string( emin ) ;
    std::string semax    = "emax=" + std::to_string( emax ) ;
    std::string sdeadc   = "deadc=" + std::to_string( deadc ) ;
    std::string nthreads = "nthreads=" + std::to_string( threads ) ;
    std::string logfile  = "logfile=" + outpath + "/" 
                           + source + "ctobssim.log" ;
    
    //  Array of arguments with ctobssim options
    char* thisargv[]     = { &thislist[ 0 ] ,
                             &inmodel[ 0 ] ,
                             &scaldb[ 0 ] ,
                             &sirf[ 0 ] ,
                             &edisp[ 0 ] ,
                             &outevent[ 0 ] ,
                             &seed[ 0 ] ,
                             &sra[ 0 ] ,
                             &sdec[ 0 ] ,
                             &srad[ 0 ] ,
                             &tmin[ 0 ] ,
                             &tmax[ 0 ] ,
                             &semin[ 0 ] ,
                             &semax[ 0 ] ,
                             &sdeadc[ 0 ] ,
                             &nthreads[ 0 ] ,
                             &chatter[ 0 ] , 
                             &logfile[ 0 ] } ;
    
    int charsize         = ( sizeof thisargv) / ( sizeof thisargv[ 0 ] ) ;
    
    //  Initializaing ctobssim tool
    ctobssim sim( charsize , thisargv ) ;
    
    //  Run and save of ctobssim events
    sim.run() ;
    sim.save() ;

    const GObservation* simobs = sim.obs()[ 0 ] ;

    //  The next section is based in csphagen script
    //  But, I need to test some "new classes" and
    //  I want to avoid compile all the source code
    //  to see that I forget that I was using a pointer
    //  :: :)
    std::cout << "Calculating background regions" << std::endl ;

    //  I will assume that there is only one background region
    //  near to the source center that must be skipped.
    //  Then, the number of background regions to start is 1 + Nskip
    int Nskip  = 1 ;
    int Nstart = 1 + Nskip ;
    int Nlim   = 1 + 2 * Nskip ;

    //  Angular separation between reflected regions
    double alpha = 1.05 * 2.0 * rad / offset ;

    // Return
    return 0 ;

}


