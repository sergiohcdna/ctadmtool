# ctaAnalysis
This project is dedicated to perform analysis cith ctools and gammalin to compute upper limits for dark matter models in cluster of galaxies.
The code presented here is functional under ctools and gammalib v1.6.3. Also, the code has been tested with development version 1.7 for both, ctools and gammalib.

To install ctools and gammalib you can use conda. Also, you can build from source. The last is particularly useful if you need to test or create new spectral and model classes. Please check [gammalib](http://cta.irap.omp.eu/gammalib/admin/index.html) and [ctools](http://cta.irap.omp.eu/ctools/admin/index.html) pages.

The contents of the project are:

1. **XMLTemplates/**
   * Model templates in XML format (default in gammalib). Currently: IC 310, NGC 1275 and Perseus(Point-source toy model)
2. **additionalData/**
   * ***Sensitivity/***
     - CSV files obtained after running cssens tools to obtain sensitivity for a specified source
   * Files (two-column format) with energy and flux * energy^2  data (in TeV units)
3. **auxfunctions/**
   * ***__init__.py***
     - python package declaration
   * ***auxMan.py***
     - Functions to manage of files, paths and erase lines in screen
   * ***myCTAfuncs.py***
     - ctools and gammalib functions to manage observation container, simulation, and some plotting 
   * ***mycsspec.py***
     - Custom version of csspec tool. This code includes the calculation of TS difference with respect the best parameters computed in likelihood
4. **DMlimits.py**
   * Script for simulation of observation, likelihood fit and calculation of upper limits for dark matter models (Up-to-date)
5. **ONOFFAnalysisGralMean.py**
   * Implementation for simulation of observation, likelihood fit and calculation of upper limits, residuals calculation, spectrum calculation and plotting, for a specified xml-model.
6. **dmexample.py**
   * Example using csdmatter application. See below

- [X] Future work: More comments in code :)
- [X] Future work: Add c++ classes for DM analysis

Feel free to make changes in this file :+1:

## ctools and gammalib setup

You can set the ```GAMMALIB``` and ```CTOOLS``` environment variables as usual. See the documentation about gammalib and ctools. Optionally, you must set the next environment variable

```
	export MYXMLS=/path/to/ctaAnalysis/XMLTemplates
```
This variable is used in some of the XML Model-templates to indicate where a spectrum-file is, for example. You can also use when running DMLimits script to indicate where the input model is

## Example running the csdmatter prototype
csdmatter script is based on cscripts apps. The csdmatter app is a prototype for the dmtool in ctools, so there are several improvements that would be added. When you clone the project, please be sure that you copy or move the files csdmatter.par and csdmatter.txt to $CTOOLS/syspfiles/ and $CTOOLS/share/help/ folders, respectively.

dmexample.py is an example based on other HOW to's from the official documentation. To run dmexample, you only need to type in your term window:

```
$python dmexample.py
```

Please check the file XMLTemplates/Perseus_DMPS_anna_1000GeV_b.xml to change the path to the spectrum file. You can change this file, if you want to try with a different spectrum file

If you find any problem or you have any suggestion, please open an issue :)

## DM Limits Calculation

Implementation of the calculation of exclusion limits on annihilation or decay of DM particles can be found in **DMLimits** script. This is prototype of the analysis chain to be translated into a cscripts class. At this moment, **DMLimits** perform the simulation of an observation based in the input DM-model. If the ```--is_onoff``` argument is present in the command line, then **DMLimits** creates a GCTAOnOff Observation using ```csphagen```. Then, ctlike is used to compute the best-fit parameters. If the ```tscalc``` flag is activated in the XML Model definition, TS is calculated automatically by ctlike. **DMLimits** checks if this flag is activated, if not then compute the TS by copying the observation (obsII) in the ctlike container, remove the source model from obsII, and create a temporary ctlike object to compute the likelihood for the null hypothesis. If TS < 25, then **DMLimits** compute UL to the flux using ```ctulimit```. This process is repeated for some number of realizations indicated with the argument ```--nsims```. All the relevant results are saved to a fits file. The results are:

1. Observation RUNID
2. Number of events in the observation container
3. TS computed
4. Scale factor (ratio of UL differential flux to theoretical flux, both at same reference energy)
5. Logarithm of parameter of interest (cross-section or lifetime)
6. Value and error of free parameters

You can request the help message from **DMLimits** by typping ```python DMLimits.py --help```. The complete message is:

```python3
usage: DMLimits.py [-h] --inmodel PerseusModel.xml --gname Perseus
                   [--coordsys CEL] --ROIra 47.0 --ROIdec 39.0 --ROIradius
                   10.0 --pntra 47.0 --pntdec 39.0 --pntradius 10.0
                   [--caldb prod3b-v2] [--irf North_z20_50h] [--hours 50.0]
                   [--enumbins 20] [--emin 0.01] [--emax 100.0] [--id 0]
                   [--pts 25.0] --nsims 10 [--nthreads 10]
                   [--outpath path/to/file/] [--is_onoff]
                   {anna,decs} ...

This script compute the DM limits using Cirelli et al. Spectrum

positional arguments:
  {anna,decs}           Selecting DM processess to produce gamma-rays
    anna                DM annihilation Process
    decs                DM decay Process

optional arguments:
  -h, --help            show this help message and exit

Source:
  All the relevant information about the source is passed in the inmodel
  argument, where a XML file is passed describing parameters for the
  observation, as position of the source, name, and spatial and spectral
  models for the source. It is assumed that the file provides information
  about the background model used in the observation. For more details, the
  user is referred to the gammalib and ctools documentations :)

  --inmodel PerseusModel.xml
                        File with xml model
  --gname Perseus       Name of the simulation
  --coordsys CEL        Coordinate system used to perform the simulation.
                        Options are: [ GAL , CEL ] ( Galactic , Celestial )
  --ROIra 47.0          First coordinate for center of Region of Interest
                        (according to Coord. Sys.)
  --ROIdec 39.0         Second coordinate for center of Region of Interestn
                        (according to Coord. Sys.)
  --ROIradius 10.0      Radius of the Region of Interest (degrees)

CTA-Intrument options:
  Information about CTA IRF, observation time, etc.

  --pntra 47.0          First coordinate of the pointing direction (according
                        to Coord. Sys.)
  --pntdec 39.0         Second coordinate of the pointing direction (according
                        to Coord. Sys.)
  --pntradius 10.0      Radius of the observation (degrees)
  --caldb prod3b-v2     Database production for the IRF file. Options are: [
                        prod2 , prod3b-v1 , prod3b-v2 ]
  --irf North_z20_50h   Instrument Response Function for CTA
  --hours 50.0          Time for simulation of observation (in hours)
  --enumbins 20         Number of energy bins
  --emin 0.01           Minimum energy for events (TeV)
  --emax 100.0          Maximum energy for events (TeV)
  --id 0                A number identifier for events in simulation.
  --pts 25.0            Value for TS to check if plot point-flux (ts > pts) or
                        upper-limit (ts < pts). Default value: 25.0
  --nsims 10            Number of simulations
  --nthreads 10         Number of threads to speed the calculations
  --outpath path/to/file/
                        Path to save files
  --is_onoff            Boolean to indicate ON/OFF observation type
```

As you can see, it is possible to select between annihilation or decay of DM. The help on this can be displayed using ```python DMLimits.py anna(decs) --help```. The help message is:

```python3
$ python DMLimits.py anna --help
usage: DMLimits.py anna [-h] [--sigmav 1.e-28] [--jfactor 1.e+20]
                        [--channel 11] [--mass 100.0]

optional arguments:
  -h, --help        show this help message and exit
  --sigmav 1.e-28   Annihilation cross-section (cm**3/s)
  --jfactor 1.e+20  Astrophysical J factor (in GeV**2/cm**5)
  --channel 11      Annihilation channel (According to PPPC4DM Tables)
  --mass 100.0      DM mass (in TeV)
```

An example to run the script is:

```python3
$python DMLimits.py --gname Perseus --ROIra 49.946 --ROIdec 41.513 --ROIradius 1.0 --irf North_z20_50h --coordsys CEL --pntra 48.6 --pntdec 40.0 --pntrad 3 --hours 1.0 --emin 0.1 --emax 10.0 --id DMSim --nsims 10 --inmodel $MYXMLS/Perseus_DMPointSource_anna_1000GeV_tau.xml --outpath Perseus_PSDM_annaII --caldb prod3b-v2 --is_onoff anna --mass 10.0 --jfactor 1.2e+19 --channel 11 --sigmav 1.e-28
```

## ON/OFF Analysis

For the (*Wobble*) ON/OFF observation, I implemented two scripts: the first one in python, and the second in C++. An observation in (*Wobble*) ON/OFF mode is where the CTA telescopes are pointing to an alternative sky position lying x degrees from the center of the region of interes (ROI). The python script was written using the csphagen python module, default when ctools is built. In the other hand, the C++ script is currently under construction, because there is not a tool to create ON/OFF observations. The creation of the ON/OFF observation is based in the csphagen python tool, but with several simplifications. The goal of the C++ script is to test new gammalib classes to describe DM scpectral models according to the purposes of this project, and/or calculation of full-likeihood (including nuisance parameters).

### Python Script

The script **ONOFFAnalysisGralMean.py** presents the case for simulation and analysis of an observation in (*Wobble*) ON/OFF mode. You need to specify the path to XML-model file, the name of the source of interest (you can have more than one source described in the template), pointing position, some instrument-related options, and analysis-related options. The script performs n observations of h hours, using IRF and caldb indicated via command-line. After every realization, the code performs a MLE-fit to estimate best parameters and extract the number of events for each bin of energy. Then, an ON/OFF container is created with the average number of events (computed from the n realizations) and a *average* spectrum is computed using ctlike and ctulimit. The script returns a file with average counts before/after MLE fit, the energy bins boundaries, XML model after MLE fit, histograms with parameter distribution for n realizations, file with parameters obtained after MLE fit, a png image with the map centered in the position of source of interest, a png image with the *average* spectrum plot; arf, pha and rmf files obtained after ON/OFF conversion.

To run the script, please be sure to have ctools and gammalib loaded. You can see the help mesagge by typping:

```python3
python3 ONOFFAnalysisGralMean.py --help
```

The complete message is:

```
usage: ONOFFAnalysisGralMean.py [-h] --inmodel PerseusModel.xml --gname
                                Perseus [--coordsys CEL] --ROIra 47.0 --ROIdec
                                39.0 --ROIradius 10.0 --pntra 47.0 --pntdec
                                39.0 --pntradius 10.0 [--caldb prod3b-v1]
                                [--irf North_z20_average_50h] [--hours 50.0]
                                [--enumbins 20] [--emin 0.01] [--emax 100.0]
                                [--id 0] [--pts 25.0] --nsims 10
                                [--outpath path/to/file/] [--is_onoff]
                                [--algorithm SIGNIFICANCE]
                                [--with_sensitivity] [--sensFile CrabSens.txt]
                                [--srcSens Crab] [--addData data.txt]

This script compute the DM limits using Cirelli et al. Spectrum

optional arguments:
  -h, --help            show this help message and exit

Source:
  All the relevant information about the source is passed in the inmodel
  argument, where a XML file is passed describing parameters for the
  observation, as position of the source, name, and spatial and spectral
  models for the source. It is assumed that the file provides information
  about the background model used in the observation. For more details, the
  user is referred to the gammalib and ctools documentations :)

  --inmodel PerseusModel.xml
                        File with xml model
  --gname Perseus       Name of the simulation
  --coordsys CEL        Coordinate system used to perform the simulation.
                        Options are: [ GAL , CEL ] ( Galactic , Celestial )
  --ROIra 47.0          RA for center of Region of Interest (according to
                        Coord. Sys.)
  --ROIdec 39.0         DEC forcenter of Region of Interestn (according to
                        Coord. Sys.)
  --ROIradius 10.0      Radius of the Region of Interest (degrees)

CTA-Intrument options:
  Information about CTA IRF, observation time, etc.

  --pntra 47.0          RA for pointing direction (according to Coord. Sys.)
  --pntdec 39.0         DEC for pointing direction (according to Coord. Sys.)
  --pntradius 10.0      Radius of the observation (degrees)
  --caldb prod3b-v1     Database production for the IRF file. Options are: [
                        prod2 , prod3b-v1 , prod3b-v2 ]
  --irf North_z20_average_50h
                        Instrument Response Function for CTA
  --hours 50.0          Time for simulation of observation (in hours)
  --enumbins 20         Number of energy bins
  --emin 0.01           Minimum energy for events (TeV)
  --emax 100.0          Maximum energy for events (TeV)
  --id 0                A number identifier for events in simulation.
  --pts 25.0            Value for TS to check if plot point-flux (ts > pts) or
                        upper-limit (ts < pts). Default value: 25.0
  --nsims 10            Number of simulations
  --outpath path/to/file/
                        Path to save files
  --is_onoff            Boolean to indicate ON/OFF observation type
  --algorithm SIGNIFICANCE
                        Method used to plot error for residual plot
  --with_sensitivity    Boolean to indicate if spectrum plots need to show
                        sensitivity curves
  --sensFile CrabSens.txt
                        File with data for sensitivity
  --srcSens Crab        Name of Source used to compute the sensitivity
  --addData data.txt    Additional data to plot in spectrum plot
```

An example to run this script looks like:

```python3
python3 ONOFFAnalysisGralMean.py --inmodel XMLTemplates/NGC1275ModelQuiescentState.xml \
--gname NGC1275 \
--ROIra 49.951 \
--ROIdec 41.512 \
--ROIradius 1.0 \
--pntra 48.0 \
--pntdec 40.0 \
--pntradius 3.0 \
--id 0 \
--nsims 100 \
--outpath $HOME/ResultsPerseus/NGC1275QuiescentState100Sim2606202TestTS \
--irf North_z20_50h \
--hours 300.0 \
--caldb prod3b-v2 \
--emin 0.1 \
--with_sensitivity \
--sensFile additionalData/Sensitivity/CrabSens300.0h.txt \
--srcSens Crab \
--addData ''
```

Please note that for the option ```--addData```, you can pass an empty string, in this case no additional data will be plotted in spectrum plot. 

### C++ script

The C++ script is under **cppScripts/** folder. The name of the script is *ONOFFAnalysis.cpp*. I am currently coding the part about ON/OFF observation-container. You can the model via an XML template, and the direction of the pointing in celestial coordinates. If you have more than one source in the XML template you must pass the name of the source you are interested to simulate. As in the *python script* you pass (via command line) information about the energy bounds and bins, calibration database and IRF, radius of ROI (ON/OFF), name of output folder, name to save count-cubes.

To run the C++ script yu must compile *ONOFFAnalysis.cpp*. To do that please follow the next instructions:

1.  After loading your *ctools* and *gammalib* environment, type ```csinfo info```. This would show information about the current installation. 
2.  Put attention to the last four rows. Get the information about GAMMALIB and CTOOLS CFLAGS and LIBS. In my case, I have something like:

```
GAMMALIB  CFLAGS ................. -fopenmp -I/home/shkdna/CTASoft09062020/gamma/include/gammalib 
                                   -I/home/shkdna/CTASoft09062020/gamma/include
```

3.  Change to $HOME/ctaAnalysis/cppScripts
4.  Type:
```cpp
  g++ ONOFFAnalysis.cpp \
  -fopenmp \
  -I/home/shkdna/CTASoft09062020/gamma/include/gammalib \
  -I/home/shkdna/CTASoft09062020/gamma/include \
  -I/home/shkdna/CTASoft09062020/gamma/include/ctools \
  -L/home/shkdna/CTASoft09062020/gamma/lib \
  -lgamma -lctools \
  -o ONOFFAnalysis
```
Change properly the information of the path to folder where you installed gammalib and ctools libs and include.

At this point, everything must be go smoothly. So, you can check all the current options by typing:

```cpp
./ONOFFAnalysis
```

The complete message is:

```cpp
You must check the number of arguments
Usage: ./ONOFFAnalysis <args> CTA-relatedArguments:
	XMLfile  : (string) XML file with model
	pntRA    : (double) R.A. of pointing, in degrees
	pntDeC   : (double) Dec. of pointing, in degrees
	source   : (string) Name of source of interest
	interval : (double) Duration of observation, in hours
	irf      : (string) Intrument response function
	caldb    : (string) Calibration database
	deadc    : (double) Dead time, [0,1]
	emin     : (double) Minimum energy, in TeV
	emax     : (double) Maximum energy, in TeV
	enbins   : (int)    Number of energy bins
	rad      : (double) Radius of region of interest
	srcrad   : (double) Radius of SRC-region of interest
	outpath  : (string) Path to output directory
	obsname  : (string) Name of the XML-file with all the
	                    details about CTAObservation
	threads  : (int)    Number of threads used to
	                    parallelization when possible

```

By the moment, you must pass all the options (sorry, this must be change). Then, to run the script please (I will update this when new changes are added):

```cpp
./ONOFFAnalysis \
/home/shkdna/ctaAnalysis/XMLTemplates/NGC1275ModelQuiescentState.xml \
45.0 \
50.0 \
NGC1275 \
1 \
North_z20_5h \
prod3b-v2 \
0.95 \
0.1 \
10.0 \
40 \
3.0 \
1.0 \
NGC1275Simulation \
NGC1275ObsList.xml \
8
```

The script perform several safety checks at different moments. The most important is to check if the angular distance between the pointing and the source centers are neither too far nor too close to compute the off regions for background estimation. In the python script, this is performed by the csphagen script. In the code:

```cpp
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

```

where ```pntdir``` and ```srcdir``` are the GSkyDir objects for the pointing of the ROI and the coordinates of the source.
