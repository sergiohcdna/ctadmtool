# ctaAnalysis
This project is dedicated to perform analysis with ctools and gammalib to compute upper limits of annihilation cross-section and decay lifetime of dark matter particles using CTA observations.

The code presented here is functional at least for versions of ctools and gammalib greater than 1.6.3.

To install ctools and gammalib you can use conda. Also, you can build from source. Please check [gammalib](http://cta.irap.omp.eu/gammalib/admin/index.html) and [ctools](http://cta.irap.omp.eu/ctools/admin/index.html) pages.

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
  - To generate gammalib.GModels during execution time without generating XML files before

There are also some files in *data* and *pfiles* folders:

1.  *data*
  - Tables from PPPC4DMID project
  - Events cube from a CTA simulation of Perseus region for 10h
2.  *pfile*
  - Parameter file of **csdmatter** app
  - Help of **csdmatter** app

Finally, **ctaAnalysis** contains also the **csdmatter** app, based on cscripts

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

By default, the models to describe the gamma-ray flux from dark matter annihilaions are parametrized using FileFunction model type. For this models, the only parameter is the *Prefractor* that accounts for the normalization of the gamma-ray differential-flux. Then, in order to compute the upper-limit of the cross-section, we compute a scale factor that is the ratio of the theoretical gamma-ray flux at some reference energy and the gamma-ray flux computed during the ctulimit instance. The **csdmatter** app ask for a reference cross-section to compute the gamma-ray flux so, using this reference value and the scale factor described previously.

### Results

As already it was mentioned, the **csdmatter** app computes the upper-limits for a family of mass points (you can specify as many points as you need using the input parameter *mnumpoints*). These masses corresponds to differentes dark matter particles. They are separated logarithmically. For every mass value, the **csdmatter** app generate the corresponding GModel, and compute the upper-limit. Then, for every mass point, the following results are saved:

1.  Reference Energy
2.  Flux computed during Likelihood fit
3.  Error on the flux
4.  TS
5.  UL on flux
6.  Reference cross-section
7.  Scale factor
8.  UL of the cross-section

At the end, the results for all mass points are saved into a fits file

##  Examples folder

Related to the **ctaAnalysis** package, you can find examples on how to use:

1.  dmspectra class
2.  dmflux class
3.  csdmatter app

### Old

There is afolder named **old**. The contents will be deprecated for posteriors versions of the **csdmatter** app. If we need some of the scripts within **old**, then would be moved to the current structure of the **ctaAnalysis package**

##  Final comments

Please feel free to check the code and test. If you find any issue, bug, problem or you have a suggestion, open an issue.