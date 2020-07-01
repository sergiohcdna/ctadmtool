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
4. **ClusterDMlimit.py**
   * Out-of-date script for simulation of observation, likelihood fit and calculation of upper limits for dark matter models in galaxy clusters
5. **ONOFFAnalysisGralMean.py**
   * Implementation for simulation of observation, likelihood fit and calculation of upper limits, residuals calculation, spectrum calculation and plotting, for a specified xml-model.

- [ ] Future work: More comments in code :)
- [ ] Future work: Add c++ classes for DM analysis

Feel free to make changes in this file :+1:

## ctools and gammalib setup

I (Sergio) have two functional versions of ctools and gammalib running in a cluster (@UNAM). I compiled the code from source using: python 3.6.8 (installed with pyenv), gnu 7.2, cmake 3.9, openssl 1.1.1 openmpi  4.0.3 and automake 1.15 (needed if you want to compile the code clonning the project from github, though). The installation is under centos 6. I also have a local (conda) installation in MAC Os-x Mojave. In particular, all the scripts here use python3 syntax :sunglasses:.

## ON/OFF Analysis

The script **ONOFFAnalysisGralMean.py** presents the case for simulation and analysis of an observation in (*Wobble*) ON/OFF mode, where the CTA telescopes are pointing to an alternative sky position (via command-line argumment) lying x degrees from the center of the region of interes (ROI). You need to specify the path to XML-model file, the name of the source of interest (you can have more than one source described in the template), pointing position, some instrument-related options, and analysis-related options. The script performs n observations of h hours, using IRF and caldb indicated via command-line. After every realization, the code performs a MLE-fit to estimate best parameters and extract the number of events for each bin of energy. Then, an ON/OFF container is created with the average number of events (computed from the n realizations) and a *average* spectrum is computed using ctlike and ctulimit. The script returns a file with average counts before/after MLE fit, the energy bins boundaries, XML model after MLE fit, histograms with parameter distribution for n realizations, file with parameters obtained after MLE fit, a png image with the map centered in the position of source of interest, a png image with the *average* spectrum plot; arf, pha and rmf files obtained after ON/OFF conversion.

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
