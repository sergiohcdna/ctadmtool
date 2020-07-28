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
              << "\tsrcrad   : (double) Radius of SRC-region of interest\n"
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

    //  Setting energy bounds in TeV
    GEbounds ebounds( enbins , GEnergy( emin , "TeV" ) , GEnergy( emax , "TeV" ) ) ;

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

/****************************************************************
 *      Create ctobssim object and run simulation               *
 *      Required parameters:                                    *
 *          - obsname   : (string)  XML obs-definition file     *
 *          - xmlmodel  : (string)  XML source model            *
 *          - caldb     : (string)  CTA-calibration database    *
 *          - irf       : (string)  CTA-IRF                     *
 *          - outname   : (string)  Output Event fits           *
 *          - seed      : (int)     Seed for random gen         *
 *          - pntRA     : (double)  pointing RA in degrees      *
 *          - pntDec    : (double)  poiniting Dec in degrees    *
 *          - pntrad    : (double)  radius of ROI               *
 *          - tmin      : (ul int)  obs-Start time (MJD)        *
 *          - tmax      : (ul int)  obs-End time (MJD)          *
 *          - emin      : (double)  Minimum Energy (TeV)        *
 *          - emax      : (double)  Maximum Energy (TeV)        *
 *          - deadc     : (double)  Deadtime (0-1)              *
 *          - nthreads  : (int)     Number of threads (par)     *
 *          - logfile   : (string)  Ouput logfile               *
 ***************************************************************/
ctobssim simobs( std::string obsname , std::string xmlmodel , std::string caldb ,
                 std::string irf ,     std::string outname ,  int seed ,
                 double pntRA ,        double pntDec ,        double pntrad ,
                 unsigned long int tmin ,            
                 unsigned long int tmax ,             
                 double emin ,         double emax ,         
                 double deadc ,        int nthreads ,
                 std::string logfile )
{
    /****************************************************************
     *      In Order to pass all the arguments to ctobssim tool     *
     *      and to avoid asking to the user for all the arguments   *
     *      then, I am using the contruction method to pass an      *
     *      array with all the options.                             *
     *      As usual, the options are passed via:                   *
     *          option=val_option                                   *
     *      For detailed information, see ctobssim help.            *
     ***************************************************************/

     std::string s_inobs     = "inobs=" + obsname ;
     std::string s_inmodel   = "inmodel=" + xmlmodel ;
     std::string s_caldb     = "caldb=" + caldb ;
     std::string s_irf       = "irf=" + irf ;
     std::string s_edisp     = "edisp=no" ;
     std::string s_chatter   = "chatter=4" ;
     std::string s_outevents = "outevents=" + outname ;
     std::string s_seed      = "seed=" + std::to_string( seed ) ;
     std::string s_ra        = "ra=" + std::to_string( pntRA ) ;
     std::string s_dec       = "dec=" + std::to_string( pntDec ) ;
     std::string s_rad       = "rad=" + std::to_string( pntrad ) ;
     std::string s_tmin      = "tmin=" + std::to_string( tmin ) ;
     std::string s_tmax      = "tmax=" + std::to_string( tmax ) ;
     std::string s_emin      = "emin=" + std::to_string( emin ) ;
     std::string s_emax      = "emax=" + std::to_string( emax ) ;
     std::string s_deadc     = "deadc=" + std::to_string( deadc ) ;
     std::string s_nthreads  = "nthreads=" + std::to_string( nthreads ) ;
     std::string s_logfile   = "logfile=" + logfile ;

     //  Array of arguments with ctobssim options
     char* s_args [] = { &s_inobs[ 0 ] , &s_inmodel[ 0 ] , &s_caldb[ 0 ] ,
                         &s_irf[ 0 ] , &s_edisp[ 0 ] , &s_outevents[ 0 ] ,
                         &s_seed[ 0 ] , &s_ra[ 0 ] , &s_dec[ 0 ] ,
                         &s_rad[ 0 ] , &s_tmin[ 0 ] , &s_tmax[ 0 ] ,
                         &s_emin[ 0 ] , &s_emax[ 0 ] , &s_deadc[ 0 ] ,
                         &s_nthreads[ 0 ] , &s_chatter[ 0 ] , &s_logfile[ 0 ]
                       } ;
    int s_sargs      = ( sizeof s_args ) / ( sizeof s_args[ 0 ] ) ;

    //  Initializaing ctobssim tool
    ctobssim sim( s_sargs , s_args ) ;

    //  Run and save of ctobssim events
    sim.execute() ;

    //  Return ctobssim object
    return sim ;

}

/************************************************************************
 *      Create ctlike object and run likelihood calculation             *
 *      Required parameters:                                            *
 *          - obsxml    : (string)          XML obs-definition file     *
 *          - obslit    : (GObservations)   Obs-list object             *
 *          - model     : (string)          XML model                   *
 *          - caldb     : (string)          CTA-calibration database    *
 *          - irf       : (string)          CTA-IRF                     *
 *          - outmodel  : (string)          XML out-model file          *
 *          - accuracy  : (double)          Accuraccy for calculations  *
 *          - max_iters : (int)             Maximum number of iters     *
 *          - fix_spat  : (bool)            Fix spatial model?          *
 *          - nthreads  : (int)             Number of threads (par)     *
 *          - logfile   : (string)          Ouput logfile               *
 ***********************************************************************/
ctlike obslike( std::string obsxml , GObservations obslist , std::string model ,
                std::string caldb , std::string irf , std::string outmodel ,
                double accuracy , int max_iters , bool fix_spat ,
                int nthreads , std::string logfile )
{
    //  So, now it's time to compute the likelihood
    std::string sfix_spat = "no" ;
    if ( fix_spat ) {
        sfix_spat = "yes" ;
    }

    std::string l_inobs    = "inobs=" + obsxml ;
    std::string l_inmodel  = "inmodel=" + model ;
    std::string l_caldb    = "caldb=" + caldb ;
    std::string l_irf      = "irf=" + irf ;
    std::string l_outmodel = "outmodel=" + outmodel ;
    std::string l_accuracy = "like_accuracy=" + std::to_string( accuracy ) ;
    std::string l_iters    = "max_iter=" + std::to_string( max_iters ) ;
    std::string l_fixspat_ = "fix_spat_for_ts=" + sfix_spat ;
    std::string l_nthreads = "nthreads=" + std::to_string( nthreads ) ;
    std::string l_logfile  = "logfile=" + logfile ;

    //  Array of arguments with ctobssim options
    char* l_args [] = { &l_inobs[ 0 ] , &l_inmodel[ 0 ] , &l_caldb[ 0 ] ,
                        &l_irf[ 0 ] , &l_outmodel[ 0 ] ,
                        &l_accuracy[ 0 ] , &l_iters[ 0 ] , 
                        &l_nthreads[ 0 ] , &l_logfile[ 0 ] } ;
    int l_sargs     = ( sizeof l_args ) / ( sizeof l_args[ 0 ] ) ;

    //  Initializaing ctobssim tool
    ctlike like( l_sargs , l_args ) ;

    //  At this point, I don't know why I must to set the observation
    //  if I already pass an observation-list file, :)
    like.obs( obslist ) ;

    //  Execute...
    like.execute() ;

    //  Return the object after the fit
    return like ;

}

GSkyRegions get_onregions( GSkyDir srcdir, double srcrad , 
                           bool save=false , std::string fname="" )
{
    //  Container for On-region
    GSkyRegions on ;

    //  Setting circle region
    GSkyRegionCircle reg( srcdir , srcrad ) ;

    on.append( reg ) ;

    if( save and ( fname.size() > 0 ) ) {
        on.save( fname )
    } else{
        std::cout << "Sky (on) region not saved" << std::endl ;
    }

    //  Return skyregion
    return on ;
}

GSkyRegions get_offregions( GSkyDir pnt , GskyDir src , double srcrad ,
                            int Nskip = 1 , int Nmin=2 
                            bool save=false , std::string fname="" )
{
    //  Container for Off-regions
    GSkyRegions off ;

    //  Setting conditions to compute number of circle off regions
    int Nstart    = 1 + Nskip ;
    int Nlim      = 1 + 2 * Nskip ;

    //  Distance between pointing and center of source
    double offset = pntdir.dist_deg( srcdir ) ;

    //  Compute position angle between sky directions in degrees
    double posang = pnt.posang_deg( src ) ;

    //  Angular separation between reflected regions
    double alpha  = 1.05 * 2.0 * srcrad / offset ;

    //Calculation of the number of background regions
    int N = ( int ) ( 2.0 * gammalib::pi / alpha ) ;

    //  If number of bkg regions is less than Nmin, then exit
    if ( N < ( Nmin + Nlim ) ) {
        std::cout << "\tThe number of reflected regions "
                  << "for background estimation is "
                  << "smaller than the minimum required."
                  << "\n\t\t\tAborting..."
                  << std::endl ;
        exit( EXIT_FAILURE ) ;
    }

    std::cout << "Number of background regions "
              << "computed: " << N - Nlim
              << std::endl ;

    //  Angular separation between regions
    double a        = 360.0 / N ;
    double dphi_max = 360.0 - a * ( 1 + Nskip ) ;
    double dphi     = a ;

    //  Creating CIRCLE--off-regions
    while ( dphi < dphi_max ) {
        GSkyDir ctrdir( pntdir ) ;
        ctrdir.rotate_deg( posang + dphi , offset ) ;
        GSkyRegionCircle thisoff( ctrdir , srcrad ) ;

        //  Append offregion to container
        off.append( thisoff ) ;

        dphi += a ;
    }

    //  Save the SkyRegions container
    if ( save and ( fname.size() > 0 ) ) {
        off.save( fname ) ;
    } else {
        std::cout << "Sky (off) region not saved" << std::endl ;       
    }

    //  Return the SkyRegions container
    return off ;

}

int main( int argc , char* argv[] )
{

    // Check the number of arguments in the command line
    if ( ( argc != 17 ) ) {
        
        // If number of arguments is different from 17
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
    double srcrad       = strTodouble( argv[ 13 ] ) ;
    std::string outpath = argv[ 14 ] ;
    std::string obsname = argv[ 15 ] ;
    int threads         = strToint( argv[ 16 ] ) ;

    std::time_t  thistime     = std::time( 0 ) ;
    std:tm *gmt_time          = std::gmtime( &thistime ) ;
    unsigned long int thismjd = get_mjd( gmt_time ) ;
    unsigned long int tmax    = thismjd + interval ;

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

    if ( offset < 1.05 * srcrad ) {
        std::cout << "Sorry, the distance between the centers"
                  << " of the source and the pointing are so"
                  << " close to get background regions"
                  << std::endl ;
        exit( EXIT_FAILURE ) ;
    } else if ( offset > 5.0 ) {
        std::cout << "Sorry, the distance between the centers"
                  << " of the source and the pointing are so"
                  << " far to compute background regions"
                  << std::endl ;
        exit( EXIT_FAILURE ) ;
    }

    // Print number of observation
    print_obs( 10 ) ;

    // Declare GObservations container
    GObservations obslist ;

    // CTA observation container
    std::string id = "0001" ;
    GCTAObservation myobs = single_obs( pntdir , irf , caldb , 0.0 , interval ,
                                        deadc , emin , emax , 
                                        enbins , rad , id  ) ;

    myobs.name( source ) ;
    // Append GCTAObservation myobs to observation_list
    obslist.append( myobs ) ;
    obslist.models( models ) ;
    // Save XML-definition file
    obsname = outpath + "/" + obsname ;
    obslist.save( obsname ) ;

    //  Run simulation
    
    int seed            = ( int ) std::time( 0 ) ;
    std::string logfile = outpath + "/" + source + "ctobssim.log" ;
    std::string outname = outpath + "/" + source + "Events.fits" ;
    ctobssim sim = simobs( obsname , xmlfile , caldb , irf , outname ,
                           seed , pntRA , pntDeC , rad , thismjd , tmax , 
                           emin , emax , deadc , threads , logfile ) ;

    //  CTA-Observation
    GCTAObservation* ctaobs = ( GCTAObservation* ) sim.obs()[ 0 ] ;

    //  Get Energy bounds from observation in the simulation
    GEbounds ebounds  = ctaobs -> ebounds() ;
    GEnergy  thisemin = ebounds.emin() ;
    GEnergy  thisemax = ebounds.emax() ;
    
    //  Create new log Energy bounds
    GEbounds new_ebounds( enbins , thisemin , thisemax ) ;

    //  Get true ebounds
    GEnergy  true_emin( emin , "TeV" ) ;
    GEnergy  true_emax( emax , "TeV" ) ;
    GEbounds true_ebounds( enbins , true_emin , true_emax ) ;

    //  The next section is based in csphagen script
    //  But, I need to test some "new classes" and
    //  I want to avoid compile all the source code
    //  to see that I forget that I was using a pointer
    //  :: :)
    std::cout << "\n\nCalculating background regions" << std::endl ;

    //  Containers for ON & OFF-regions
    GSkyRegions onregions ;
    GSkyRegions offregions ;

    //setting on region
    GSkyRegionCircle srcreg( srcdir , srcrad ) ;
    onregions.append( srcreg ) ;

    //  I will assume that there is only one background region
    //  near to the source center that must be skipped.
    //  Then, the number of background regions to start is 1 + Nskip
    int Nskip  = 1 ;
    int Nmin   = 2 ;
    int Nstart = 1 + Nskip ;
    int Nlim   = 1 + 2 * Nskip ;

    double posang = pntdir.posang_deg( srcdir ) ;

    //  Angular separation between reflected regions
    double alpha = 1.05 * 2.0 * srcrad / offset ;
    
    //Calculation of the number of background regions
    int N = ( int ) ( 2.0 * gammalib::pi / alpha ) ;

    //  If number of bkg regions is less than Nmin, then exit
    if ( N < ( Nmin + Nlim ) ) {
        std::cout << "\tThe number of reflected regions "
                  << "for background estimation is "
                  << "smaller than the minimum required."
                  << "\n\t\t\tAborting..."
                  << std::endl ;
         exit( EXIT_FAILURE ) ;

    }

    std::cout << "Number of background regions "
              << "computed: " << N - Nlim
              << std::endl ;

    //  Angular separation between regions
    double a = 360.0 / N ;

    //  Creating CIRCLE--off-regions
    for ( int s = Nstart  ; s < ( N - Nskip ) ; ++s ) {
        double dphi = s * a ;
        GSkyDir ctrdir( pntdir ) ;
        ctrdir.rotate_deg( posang + dphi , offset ) ;
        GSkyRegionCircle thisoff( ctrdir , srcrad ) ;

        offregions.append( thisoff ) ;
    }

    //  ONOFFModel
    GModels onoffmodels ;

    for( int index=0 ; index < models.size() ; ++index  ) {

        //  Get model from initial models container
        GModel* model = models.at( index ) ;

        //  Get classname, check if is a bkg model,
        //  and append OnOff to bkg model classname
        std::string iname = model -> classname() ;

        if ( iname.find( "CTA" ) != std::string::npos ) {
            
            std::cout << "\t\tCTA Background Model found"
                      << "\n\t\tProcessing..."
                      << std::endl ;

            std::string thisinst = model -> instruments() ;
            thisinst += "OnOff" ;
            model -> instruments( thisinst ) ;
        }

        //  Append model to use for onofmodels container
        model -> tscalc( true ) ;
        onoffmodels.append( *model ) ;
    }

    //  Creating GCTAOnOffObservation :)
    GCTAOnOffObservation onoffobs( *ctaobs , onoffmodels , source , 
                                   true_ebounds , ebounds ,  
                                   onregions , offregions , true ) ;
    
    onoffobs.statistic( "cstat" ) ;
    onoffobs.name( "ONOFF" + source ) ;

    //  This part is to save the relevant fits files
    //  Also, to create the Obs-List XML file
    std::string onname  = outpath + "/" + source + "_pha_on.fits" ;
    std::string offname = outpath + "/" + source + "_pha_off.fits" ;
    std::string arfname = outpath + "/" + source + "_arf.fits" ;
    std::string rmfname = outpath + "/" + source + "_rmf.fits" ;

    GObservations onoffobslist ;
    onoffobslist.append( onoffobs ) ;
    onoffobslist.models( onoffmodels ) ;

    //  Set background and response file names in On spectrum

    onoffobs.on_spec().save( onname , true ) ;
    onoffobs.off_spec().save( offname , true ) ;
    onoffobs.arf().save( arfname , true ) ;
    onoffobs.rmf().save( rmfname , true ) ;
   
    std::string onoffoutmodel = outpath + "/" + source + "onoffmodel.xml" ;
    std::string onoffxmlobs   = outpath + "/" + source + "onoffobs.xml" ;
    onoffobslist.models().save( onoffoutmodel ) ;

    //  Save the xml file with observation-list definitions
    //  I didn't find and obvious way to save directly from
    //  the save method from GObservations.
    //  The last is because, when I tried to access the GPha
    //  object from GCTAOnOffObservation class, the GPha object
    //  is a const-type object, then I can't call any method
    //  trying to modify the object itself :)
    GXml onoffxml ;
    GXmlElement obs_list( "observation_list title=\"observation library\"" ) ;
    GXmlElement onoff_obs( "observation name=\"OnOff" + source + "\" id=\"" 
                     + id + "\" instrument=\"CTAOnOff\" statistic=\"cstat\"" ) ;
    onoff_obs.append( GXmlElement( "parameter name=\"Pha_on\" file=\"" 
                      + onname + "\"") ) ;
    onoff_obs.append( GXmlElement( "parameter name=\"Pha_off\" file=\"" 
                      + offname + "\"") ) ;
    onoff_obs.append( GXmlElement( "parameter name=\"Arf\" file=\"" 
                      + arfname + "\"") ) ;
    onoff_obs.append( GXmlElement( "parameter name=\"Rmf\" file=\"" 
                      + rmfname + "\"") ) ;
    obs_list.append( onoff_obs ) ;
    onoffxml.append( obs_list ) ;
    onoffxml.save( onoffxmlobs ) ;

    //  Save On and Off region files
    std::string onregname  = outpath + "/" + source + "on.reg" ;
    std::string offregname = outpath + "/" + source + "off.reg" ;
    onregions.save( onregname ) ;
    offregions.save( offregname ) ;

    //  So, now it's time to compute the likelihood
    std::string l_outmodel = outpath + "/" + source + "OnOffLike.xml" ;
    std::string l_logfile  = outpath + "/" + source + "ctlike.log" ;
    bool fixspat           = true ;

    ctlike like = obslike( onoffxmlobs , onoffobslist , onoffoutmodel,
                           caldb , irf , l_outmodel , 1.e-4 , 100,
                           fixspat , threads , l_logfile ) ;

    GModels tmodels = like.obs().models() ;
    GModelSky* srcmodel = ( GModelSky* ) tmodels[ source ] ;

    for ( int i = 0; i < srcmodel -> spectral() -> size(); ++i ) {
        GModelPar locpar = srcmodel -> spectral() -> at( i ) ;
        std::string parname = locpar.name() ;
        if ( locpar.is_free() ){
            double value     = locpar.value() ;
            double scale     = locpar.scale() ;
            std::string unit = locpar.unit() ;
            std::cout << parname << ":\t" 
                      << value << " " 
                      << unit << std::endl ;
        }
    }

    GModelPar locpiv = srcmodel -> spectral() -> operator[]( "PivotEnergy" ) ;
    GEnergy   pivotE( locpiv.value() / locpiv.scale() , "TeV" );
    double fitted_flux = srcmodel -> spectral() -> eval( pivotE ) ;
    std::cout << "Fitted flux is: " << fitted_flux << std::endl ;
    double ts = srcmodel -> ts() ;
    std::cout << "TS: " << ts << std::endl ;

    std::cout << like.obs() << std::endl ;

    // Return
    return 0 ;

}


