# ctadmtool Package

This project is dedicated to perform analysis with `ctools` and `gammalib` to compute upper limits of annihilation cross-section and decay lifetime of dark matter particles using CTA observations.

The code presented here is functional at least for versions of `ctools` and `gammalib` greater than 1.6.3.

To install `ctools` and `gammalib` you can use conda. Also, you can build from source. Please check [gammalib](http://cta.irap.omp.eu/gammalib/admin/index.html) and [ctools](http://cta.irap.omp.eu/ctools/admin/index.html) pages.

The core of the analysis tool is `csdmatter`. This class is based in `csscripts` installed with ctools. One of the key features is the management of dark-matter spectrum via the `gammalib` `GModelSpectralTable` class. This has the advantage to only use one fits file to describe the spectrum for candidates with masses in a determined range, dedicated annihilation channels (following the convention from the [PPPC4DMID](http://www.marcocirelli.net/PPPC4DMID.html) project) and particular energy ranges. You can check the jupyter notebook `tablemodel` under the **notebooks/** folder to check how to create a spectrum fits file. The `tablemodel` notebook uses `dmtable` class to interpolate the spectrum for different masses, energy and channels using `dmspectrum` class. The spectrum is save into a fits table that is esaily ingested by `ctools` and `gammalib`.

The fits table has three parameters to describe the spectrum:

1. Mass. The mass of the dark matter candidate (in GeV). [**Default: fixed**]. In the case of a hypothetical dark-matter signal, you should free this parameter.
2. Channel. Annihilation channel. I took the same number convetion as in the PPPC4DMID project. [**Default: fixed**]. This variable should not set free during the analysis.
3. Normalization. Overall normalization of the dark-matter spectra. [**Default: free**]

Aditionally, you can select whether or not to include electroweak (EW) corrections in the spectrum. By default, this parameter is set to `True`.

## Current known issues:

1. Compute gamma-ray spectrum of decay to W's of WIMPs with masses below 160 GeV

If you found any issue during test or use of this project, please open an Issue.

## Description

The contents of the project are:

2.  **ctadmtool**
  * *data*
  * *dmspectrum*
  * *pfiles*
  * *tools*
  * **csdmatter.py**
4.  **LICENSE**
5.  **MANIFEST.in**
6.  **README.md**
7.  **setup.cfg**
8.  **setup.py**

## ctools and gammalib setup

You can set the ```GAMMALIB``` and ```CTOOLS``` environment variables as usual. See the documentation about gammalib and ctools.

##  The ctadmtool python package

`ctadmtool` is a python package to compute exclusion limits for model-independent dark matter searches with CTA. The package is an effort to have a common set of tools and use as basic example for analysis. There are two subpackages:

1. **dmspectrum**
2. **tools**

The relevan classes (about physics) are:

1.  *dmspectra*
  - To compute the number of photons produced during annihilation or decay of dark matter particles using PPPC4DMID tables. Here, it is also included EBL attenuation using [ebl-table project](https://github.com/me-manu/ebltable)
2.  *dmflux_table*
  - To generate `GModelSpectralTable` models for annihilation or decay of dark matter particles.

There are also some files in *data* and *pfiles* folders:

1.  *data*
  - Tables from PPPC4DMID project
2.  *pfiles*
  - Parameter file of **csdmatter** app
  - Help of **csdmatter** app

Finally, **ctadmtool** contains also the **csdmatter** app, based on cscripts.

##  The csdmatter app

The csdmatter app is based on how the cscripts are implemented within ctools. The csdmatter computes the upper-limits (at this moment, just) for annihilation cross-section for a famlily of mass points of dark matter particles. You can refer to [*pfiles/csdmatter.par*](ctadmtool/pfiles/csdmatter.par) and [*pfiles/csdmatter.txt*](ctadmtool/pfiles/csdmatter.txt) to check the full list of input parameters, and the help of the app.

##  Installation

To have **ctadmtool** package availabe in your system you must to be sure that *ctools* and *gammalib* are loaded. Then to install **ctadmtool** you have two options:

1. Cloning:
  - `$ git clone git@github.com:sergiohcdna/ctadmtool.git`
  - `$ cd ctadmtool`
  - `$ python -m pip install .`

2. Using pip directly:
  - `$ python -m pip install git@github.com:sergiohcdna/ctadmtool.git`

Please note, that, if you want to contribute to the development of **csdmatter** and related classes, you must use the first option. Additionally, you can create a branch.

**Note:** There is `Deprecation` message when usin `python -m pip install .` referring that the local packages will be building in-place. You can take a look at the [discussion](https://github.com/pypa/pip/issues/7555). To avoid the deprecation message and test this new feature you must use:

```
$ python -m pip install . --use-feature=in-tree-build
```

**Note:** If you are updating a previous installation of `ctadmtool`, please be sure that you don't have any old `csdmatter.par` files in other locations. This will create a problem if new parameters were added to the version you are trying to install. A manual solution is to erase the `csdmatter.par` files from the `$CTOOLS/syspfiles` and your `$HOME`:

```bash
$ rm $CTOOLS/syspfiles/csdmatter.par
$ rm $HOME/pfiles/csdmatter.par
```

After this, you can install `ctadmtool`.

### Running the csdmatter app

Because the script is not part of the default cscripts, you must execute it using inside a python script.

### DM Limits Calculation

To compute the ULs, **csdmatter** compute the ratio between the expected integrated flux and the upper-limit integrated flux obtained during the fit. Both fluxes are computed between minimum energy `emin` and 95% of the available energy in the process. The ratio is used to compute the exclusion limit using the reference annihilation cross-section or the decay lifetime used to compute the normalization of the dark matter spectrum. The ULs are computed for masses points (*mnumpoints* parameter) in the range of masses indicated by the user (via *mmin* and *mmax* parameters)

### Results

As already it was mentioned, the **csdmatter** app computes the upper-limits for a family of mass points (you can specify as many points as you need using the input parameter *mnumpoints*). These masses corresponds to differentes dark matter particles. They are separated logarithmically. For every mass value, the **csdmatter** app generate the corresponding GModel, and compute the upper-limit. Then, for every mass point, the following results are saved:

1. MinEnergy &#10140; Minimum Energy
2. MaxEnergy &#10140; Maximum Energy
3. Mass &#10140; Mass of the dark matter candidate
4. Flux &#10140; Flux at reference energy
5. ErrFlux &#10140; Error associated to flux
6. E2Flux &#10140; Energy squared times flux
7. E2ErrFlux &#10140; Energy squared times flux error
6. LogL &#10140; Log-Likelihood obtained during the fit
7. TS &#10140; Test Statistic
8. UpperLimit &#10140; Upper limit to the flux
11. ScaleFactor &#10140; Ratio between theoretical and upper-limit flux
12. ULLifetime or ULCrossSection &#10140; Exclusion Limit to Lifetime or Cross-section
13. RefLifetime or RefCrossSection &#10140; Reference values used to compute the dark matter flux

At the end, the results for all mass points are saved into a fits file.

If you use the `run` method, the results table can be accessed via the `dmatter_fits` method. For example:

```python
from ctadmtool.csdmatter import csdmatter

thistool = csdmatter()
# After all the parameter initialization

thistools.run()
results = thistool.dmatter_fits()
print(results)
```

You can take a look at the jupyter notebooks `get_decayllimits` and `get_dmulimits` to learn how to run the `csdmatter` tool to get exclusion limits for a *Toy dark-matter halo*.

## Notebooks

There are several jupyter notebooks to show you how to use the package. The notebooks are:

- `dmspectrum` &#10140; Shows how to use the`dmspectrum` class to get dark matter spectrum for annihilation or decay
- `tablemodel` &#10140; Shows how to use `dmtable` class to create `GModelSpectralTable` models for dark matter
- `plottingdmSpectra`&#10140; Shows how to plot spectra from table models
- `dmflux` &#10140; Shows how to compute gamma-ray fluxes using `dmspectrum` class. comparison with *Clumpy* gamma-ray flux for decay.
- `dmsimulation` &#10140; Simulate the a dark-matter signal using `GModelSpectralTable` models
- `get_dmulimits` &#10140; Shows how to use `csdmatter` tool to compute exclusion limits for annihilation of dark matter particles
- `getdecayllimits` &#10140; Shows how to use `csdmatter` tool to compute exclusion limits for decay of dark matter particles

##  Final comments

Please feel free to check the code and test. If you find any issue, bug, problem or you have a suggestion, open an issue.
