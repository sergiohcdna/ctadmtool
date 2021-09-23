# ctaAnalysis
This project is dedicated to perform analysis with ctools and gammalib to compute upper limits of annihilation cross-section and decay lifetime of dark matter particles using CTA observations.

The code presented here is functional at least for versions of ctools and gammalib greater than 1.6.3.

To install ctools and gammalib you can use conda. Also, you can build from source. Please check [gammalib](http://cta.irap.omp.eu/gammalib/admin/index.html) and [ctools](http://cta.irap.omp.eu/ctools/admin/index.html) pages.

The core of the analysis tool is **csdmatter**. This class is based in **csscripts** installed with ctools. One of the key features is the management of dark-matter spectrum via the *gammalib* **GModelSpectralTable** class. This has the advantage to only use one fits file to describe the spectrum for candidates with masses in a determined range, dedicated annihilation channels (following the convention from the PPPC4DMID project) and particular energy ranges. You can check the script *create_dmtable* under the **examples/** folder to check how to create a spectrum fits file. The *create_dmtable* also includes the estimation of EBL atenuattion and you can indicate the annihilation cross-section (*logsigmav* parameter) and astropysical jfactor (*logastfactor* parameter) to compute the overall normalization of the spectrum. You can change the value of cross-section and jfactor according to the requirements of the source (as mentioned before, the cross-section and jfactor only are used to compute the normalization of the spectra). Then, you can use only one spectrum file for a family of targets at the same redshift. The *create_dmtable* script uses two internal classes to interpolate the spectrum at energy and masses indicated by the user via command-line.

The parameters in the spectrum fits file are the following:

1. Mass. The mass of the dark matter candidate (in GeV). [**Default: fixed**]. In the case of a hypothetical dark-matter signal, you should free this parameter.
2. Channel. Annihilation channel. I took the same number convetion as in the PPPC4DMID project. [**Default: fixed**]. This variable should not set free during the analysis.
3. Normalization. Overall normalization of the dark-matter spectra. [**Default: free**]

To use the *create_dmtable* script you can use:

```bash
$ python create_dmtable.py --srcname AquariusII --jfactor 1.8621e+18 --sigmav 3.0e-26 --z 0 --mpoints 100 --emin 30.0 --nebins 450
```

## Current known issues:

1. **csdmatter** fails to correctly process *GCTAOnOffObservations*

If you found any issue during the test of use of this project please open an Issue.

## Description

The contents of the project are:

1.  **old**
2.  **ctaAnalysis**
  * *data*
  * *dmspectrum*
  * *pfiles*
  * *tools*
  * **csdmatter.py**
3.  **examples**
4.  **LICENSE**
5.  **MANIFEST.in**
6.  **README.md**
7.  **setup.cfg**
8.  **setup.py**

## ctools and gammalib setup

You can set the ```GAMMALIB``` and ```CTOOLS``` environment variables as usual. See the documentation about gammalib and ctools.

##  The ctaAnalysis python package

The ctaAnalysis is a python package to make life a little easier. There are two subpackages inside the package:

1. **dmspectrum**
2. **tools**

And the relevant classes are:

1.  *dmspectra*
  - To compute the number of photons produced during annihilation of dark matter particles using PPPC4DMID tables. Here, it is also included EBL attenuation using ebl-table project
2.  *dmflux*
  - To compute gamma-ray flux using ref. cross-section times total astrophysical factor times spectra of photons produced during annihilation
3.  *createmodels*
  - To generate gammalib.GModels during execution time without generating XML files before. (**Deprecated**. As noted above, the management of dark-matter spectrum is via the *gammalib* class **GModelSpectralTable**. This will be removed in new versions of the package)

There are also some files in *data* and *pfiles* folders:

1.  *data*
  - Tables from PPPC4DMID project
  - Events cube from a CTA simulation of Perseus region for 10h
2.  *pfiles*
  - Parameter file of **csdmatter** app
  - Help of **csdmatter** app

Finally, **ctaAnalysis** contains also the **csdmatter** app, based on cscripts.

##  Installation

In order to hava **ctaAnalysis** package availabe in your system you must to be sure that *ctools* and *gammalib* are loaded. Then to install **ctaAnalysis** you have two options:

1. Cloning:
  - `$ git clone git@github.com:sergiohcdna/ctaAnalysis.git`
  - `$ cd ctaAnalysis`
  - `$ python -m pip install .`

2. Using pip directly:
  - `$ python -m pip install git@github.com:sergiohcdna/ctaAnalysis.git`

Please note, that, if you want to contribute to the development of **csdmatter** and related classes, you must use the first option. Additionally, you can create a branch.

##  The csdmatter app

The csdmatter app is based on how the cscripts are implemented within ctools. The csdmatter computes the upper-limits (at this moment, just) for annihilation cross-section for a famlily of mass points of dark matter particles. You can refer to *pfiles/csdmatter.par* and *pfiles/csdmatter.txt* to check the full list of input parameters, and the help of the app.

### Run the csdmatter app

Because the script is not part of the default cscripts, you must execute it using inside a python script.

### DM Limits Calculation

To compute the ULs, **csdmatter** compute the ratio between the expected flux and the upper-limit flux obtained during the fit. The ratio is used to compute the exclusion limit using the reference annihilation cross-section used to compute the normalization of the dark matter spectrum. The ULs are computed for the number of masses points (*mnumpoints* parameter) in the range of masses indicated by the user (via *mmin* and *mmax* parameters)


### Results

As already it was mentioned, the **csdmatter** app computes the upper-limits for a family of mass points (you can specify as many points as you need using the input parameter *mnumpoints*). These masses corresponds to differentes dark matter particles. They are separated logarithmically. For every mass value, the **csdmatter** app generate the corresponding GModel, and compute the upper-limit. Then, for every mass point, the following results are saved:

1. MinEnergy
2. MaxEnergy
3. Mass
4. Flux
5. EFlux
6. LogL
7. TS
8. UpperLimit
9. ULCrossSection
10. RefCrossSection
11. ScaleFactor
12. 
At the end, the results for all mass points are saved into a fits file.

If you use the *run* method, the results table can be accessed via the *dmatter_fits* method. For example:

```python
from ctaAnalysis.csdmatter import csdmatter

thistool = csdmatter()
# After all the parameter initialization

thistools.run()
results = thistool.dmatter_fits()
print(results)
```

##  Examples folder

Related to the **ctaAnalysis** package, you can find examples on how to use:

1.  dmspectra class
2.  dmflux class
3.  csdmatter app
4.  *create_dmtable*. How to create dark matter spectra fits file using **GModelSpectralTable**

### Old

There is afolder named **old**. The contents will be deprecated for posteriors versions of the **csdmatter** app. If we need some of the scripts within **old**, then would be moved to the current structure of the **ctaAnalysis package**

##  Final comments

Please feel free to check the code and test. If you find any issue, bug, problem or you have a suggestion, open an issue.
