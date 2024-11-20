# galaxyEmulator

Galaxy emulation from [IllustrisTNG](https://www.tng-project.org/) using [SKIRT](https://skirt.ugent.be/root/_home.html) for various observations.

## Dependence
Python verison:  
`python==3.11.10`  

Python packages:  
`numpy==1.26.4`  
`scipy==1.13.1`  
`matplotlib==3.9.2`  
`astropy==5.3`  
`skimage==0.24.0`  
`joblib==1.4.2`  
`matplotlib_scalebar==0.8.1`  
`h5py==3.9.0`  
`termcolor==2.4.0`  
`photutils==2.0.2`

Other packages:  
[`illustris_python`](https://github.com/illustristng/illustris_python)  
[`SKIRT`](https://skirt.ugent.be/root/_installation_guide.html)  

## Recommended folder tree
```bash
├── galaxyEmulator/  
│   ├── Data/  
│   ├── galaxyEmulator/  
│   ├── Notebooks/  
│   └── (workspace)/  
├── TNG-100/  
│   ├── groups_N/  
│   │   ├── fof_subhalo_tab_N.n.hdf5  
│   │   └── ...  
│   └── snapdir_N/  
│       ├── snap_N.n.hdf5  
│       └── ...  
└── postprocessing/  
    └── offsets/  
        ├── offsets_N.hdf5  
        └── ...  
```

## Filters and PSFs
Filters and PSFs of specific survey are saved in `Data/filters/(survey)` and `Data/PSFs/(survey)` directories.  
The format of filters are recommended to be plain text file with extension `.fil`, including two columns: wavelength (in angstrom) and throughputs.  
PSFs are recommended to be `numpy.array`, and can be opened by `numpy.loadtxt` or `numpy.load`.  
Please make sure that the filters and PSFs exist in correct directories and format before running galaxyEmulator.

## Usage
### Initialization
```Python
python init.py --workspace=workspace --surveys="CSST,JWST"
```
`config.ini`, `config_CSST.ini` and `config_JWST.ini` will be created in `workspace` directory. Please freely edit them as wish.  

if `surveys` are not specified, `postPostprocessing` will be set as False. However, for consistency, `PostProcess` class must be initialized, and only data cube files will be saved.  

Note:  
Currently, we only upload throughputs and PSFs for CSST, and other surveys will be added later. (Nov. 14, 2024)  

### Run
Enter workspace, and create a python file named emulator.py
```Python
# emulator.py

import sys
sys.path.append('..')

from galaxyEmulator.config import Configuration
from galaxyEmulator.preprocess import PreProcess
from galaxyEmulator.postprocess import PostProcess

config = Configuration() # initialize Configuration class
conf = config.get_config() # read config from current directory

preprocess = PreProcess(conf) # initialize PreProcess class
preprocess.get_subhalos() # get subhalos 
subhaloIDs = preprocess.get_subhaloIDs() # get subhaloIDs

for ID in subhaloIDs:
    preprocess.subhalo(subhaloID=ID) # get properties of subhalo
    preprocess.prepare() # prepare for simulation

    preprocess.runSKIRT() # run SKIRT

    postprocess = PostProcess(properties=preprocess.properties, config=conf) # initialize PostProcess class
    postprocess.runPostprocess() # run postprocessing
```  
then, `python emulator.py`.  

Or you can interactively run by specifying the subhaloIDs in jupyter as illustrated in [Notebooks/tutorial.ipynb](https://github.com/xczhou-astro/galaxyEmulator/blob/main/Notebooks/tutorial.ipynb).  
## Classes
### Configuration
```Python
config = Configuration(surveys=None)
```
surveys: `str (N,)` or `list`, considered surveys. If None, configurations will be read from current directory, or configuration files will be created from templates.  
```Python
conf = config.get_config()
```  
return configurations.  
```
config.add_survey(surveys)
```
surveys: `str (N,)` or `list`, new surveys to be added. Call `get_config()` to update configurations.  
```
config.remove_survey(surveys)
```
surveys: `str (N,)` or `list`, surveys will be removed. Call `get_config()` to update configurations.  
```
config.save_config()
```
manually save configurations to `config.ini`.  
```
config.check_config()
```
manually check configurations.  

### PreProcess
```Python
preprocess = PreProcess(config)
```  
```Python
preprocess.get_subhalos()
```  
Subhalos in stellar mass range in `config.ini` will be read.  
```Python
subhaloIDs = preprocess.get_subhaloIDs()
```  
return subhaloIDs for subhalos obtained in `get_subhalos()`.  
```Python
stellarMasses =preprocess.get_stellarMasses()
```
return stellar masses for subhalos obtained in `get_subhalos()`.
```Python
preprocess.subhalo(subhaloID)
```
`subhaloID`: `int`: Intialize simulation for subhalo with subhaloID.  
```Python
preprocess.prepare()
```
Preparation for simulation, including retrieving particles and creating .ski file.  
```Python
properties = preprocess.get_properties()
```
Return properties used for postprocessing.  
```Python
preprocess.runSKIRT()
```
Run SKIRT.  

### PostProcess
```Python
postprocess = PostProcess(properties, config)
```
`properties`: `dict`: properties from `PreProcess`, can be obtained by `PreProcess.get_properties()` or `PreProcess.properties`.  
`config`: `dict`: configurations.  
```Python
postprocess.runPostProcess(showImages=False)
```  
`showImages`: `bool`: if show galaxy images and SEDs. Run postprocessing for galaxy simulation.  

## Config.ini
`dataDir`:  
`str`, Data directory of galaxyEmulator.  

`filePath`:  
`str`, Directory for TNG simulation.  

`workingDir`:  
`str`, Directory for execution of SKIRT.  

`simulationMode`:  
`str`, `ExtinctionOnly` or `DustEmission`, Simulation Mode, `DustEmission` should be used when near- or mid-infrared bands are considered.  

`includeDust`:  
`bool`, If include dusts generated from Gas particles; must be True if `simulationMode=DustEmission`  

`dustEmissionType`:  
`str`, `Equilibrium` or `Stochastic`, Dust emission type.  

`dustModel`:  
`str`, `ZubkoDustMix`, `DraineLiDustMix` or `ThemisDustMix`, Dust model.  

`minWavelength`, `maxWavelength`:  
`float`, Considered wavelength range, in micron; `maxWavelength` should be higher than maximum wavelength of filters considered.  

`boxLengthScale`:  
`float`, Determine the boxsize to retrieve particles;  

`maxBoxLength`:  
`float`, Maximum boxsize, in kpc; `boxisze = min(halfStellarMassRadius * boxLengthScale, maxBoxLength)`.  

`wavelengthGrid`:  
`str`, `Linear` or `log`, Wavelength grid type for SKIRT.  

`numWavelengths`:  
`int`, Wavelength bins for SKIRT.  

`minLevel`, `maxLevel`:  
`int`: Octree min/max level refinement for dust calculation.  

`numPackets`:  
`float`: number of photo packets launched during simulation; determine the level of signal to noise in the results.  

`SEDFamily`:  
`str`, `BC03` or `FSPS`, SED family for quenched star particles.  

`initialMassFunction`:  
`str`, `Chabrier` or `Salpeter` for `BC03`; `Chabrier`, `Salpeter` or `Kroupa` for `FSPS`, Intial mass function for SED family.  

`minStellarMass`, `maxStellarMass`:  
`float` or `inf`, Stellar mass range for subhalos, in 1e10 Msun; `inf` for infinite.  

`numViews`:  
`int`, Number of instrument views for datacube generation by SKIRT.  

`randomViews`:  
`bool`, If specify views by randoms.  

`inclinations`, `azimuths`:  
`float (numViews,)`, Views; must be provided if `randomViews=False`; inclinations: 0 ~ 180, azimuths: -360 ~ 360.  

`FoVboxLengthRatio`:  
`float`, Ratio for field of view, `Fov = Boxsize * FoVboxLengthRatio`.  

`postProcessing`:  
`bool`, If run postprocessing.  

`saveDataCube`:  
`bool`, If save data cubes; must be True if `postProcessing=False`.  

`spatialResol`:  
`float`, Spatial resolution for data cube, in pc; must be provided if `postProcessing=False`.  

`surveys`:  
`str (N,)`, Considered surveys; must be provided if `postProcessing=True`.  

`displaySED`:  
`bool`, If display SED.  

`displaySEDxlogscale`:  
`bool`, If display SED in logscale for wavelength.  

`snapNum`:  
`int`,  Snapshot ID of TNG simulation.  

`fixedRedshift`:  
`float`, Redshift of snapshot ID; should be larger than 0 when `snapNum=99` to avoid error.  

`numThreads`:  
`int`, Number of Threads to run SKIRT. The perforamnce will not be improved if `numThreads > 24`.  

`recordComponents`:  
`bool`, If including transparent, primary direct, primary scattered, secondarytransparent, secondarydirect, secondaryscattered components for data cube and SED; Consumption of memory increases to 7 times if True.  

`ageThreshold`:  
`float`, Age threshold for discriminating star-forming and quenched star particles, in Myr.  

`logCompactnessMean`, `logCompactnessStd`:  
`float`,  logCompactness for star-forming particles, sampled from normal distribution.  

`logPressure`:  
`float`, `log10[(Pressure/k_B)/cm^-3 K] = logPressure` for star-forming particles.  

`PDRClearingTimescale`:  
`float`, covering factor is `f = e^(-t / PDRClearingTimescale)`, where t is the age.  

`temperatureThreshold`:  
`float`, Gas particles lower than `temperatureThreshold` will be considered as dusts.  

`massFraction`:  
`float`, Fraction of the metallic gas locked up in dust.  

`DISMModel`:  
`str`, `Camps_2016` or `Torrey_2012`, dust-containing ISM (DISM) model.  

`numSilicateSizes`, `numGraphiteSizes`, `numPAHSizes`, `numHydrocarbonSizes`:  
`int`, Number of bins for dust grains.  

## Config_\[survey\].ini
Config_\[survey\].ini is generated if `postProcessing=True` and `surveys` are provided.  

`filters`:  
`str (N,)`, Considered filters for survey.  

`resolFromPix`:  
`bool`, If use resolution derived from pixel scales; instrumental effects can only be added if `resolFromPix=True`.  

`resolution`:  
`float`, Spatial resolution, in pc; override `spatialResol` in config.ini.  

`pixelScales`:  
`float (N,)` or `float (1,)`, Pixel scales for considered filters, in arcsec; provide one value if homegenous, otherwise must be consistent with filters.  

`numExposure`:  
`int (N,)` or `int (1,)` Number exposures or number of filters.  

`exposureTime`:  
`float (N,)` or `int (1,)`, Exposure time, in second.  

`aperture`:  
`float`, Aperture size for instrument, in meter.  

`includePSF`:  
`bool`, If include PSF effects; can be True if `resolFromPix=True`, otherwise must be False.  

`PSFFromFile`:  
`bool`, If PSFs are from files; pixel scales of PSF array must be consistent with `pixelScales`.  

`PSFFWHM`:  
`float (N,)` or `float (1,)`, FWHM to create Gaussian PSF kernels, in arcsec.  

`includeBkg`:  
`bool`, If include background; can be True if `resolFromPix=True`, otherwise must be False.  

`bkgNoise`:  
`float (N,)`, Sigmas for Gaussian background for each filter.  

`imgDisplay`:  
`bool`, If display galaxy images.  

`RGBImg`:  
`bool`, If display rgb image created by 3 bands.  
RGB composition by [HumVI](https://github.com/drphilmarshall/HumVI), please use with caution.  
For better illustration, several parameters need to be adjusted in `galaxyEmulator.utils.convert_to_rgb`.  

`RGBFilters`:  
`str (3,)`, 3 bands to create rgb image.  

`displayFilter`:  
`str`, Band image to display if `RGBImg=False`.  
