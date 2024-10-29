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

## Usage
### Initialization
```Python
python init.py -w=workspace
```
config.ini will be created in workspace, edit this file then:
```
python init.py -w=workspace
```
config_\[survey\].ini will be created if `postprocessing` is True and `surveys` are given in config.ini

### Run
```Python
import sys
sys.path.append('.')

from galaxyEmulator.config import Configuration
from galaxyEmulator.preprocess import PreProcess
from galaxyEmulator.postprocess import PostProcess

config = Configuration('config.ini')
conf = config.get_config()

preprocess = PreProcess(conf)
preprocess.get_subhalos()
preprocess.get_subhaloIDs()

preprocess.subhalo(subhaloID=2)
preprocess.prepare()

preprocess.runSKIRT()

postprocess = PostProcessing(properties=preprocess.properties, conf)
postprocess.runPostProcess()
```

## Classes
### Configuration
```Python
config = Configuration(main_config_file)
```
`main_config_file`: `str`, filename for config.ini  
```Python
conf = config.get_config()
```  
`conf`: `dict`, configurations including config.ini and config_\[survey\].ini

### PreProcess
```Python
preprocess = PreProcess(config)
```  
```Python
preprocess.get_subhalos()
```  
Subhalos in stellar mass range in config.ini will be read.  
```Python
preprocess.get_subhaloIDs()
```  
return subhaloIDs for subhalos obtained in `get_subhalos`.  
```Python
preprocess.get_stellarMasses()
```
return stellar masses for subhalos obtained in `get_subhalos`.
```Python
preprocess.subhalo(subhaloID)
```
`subhaloID`: `int`: ID of subhalo used for galaxy simulation.  
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

### PostProcessing
```Python
postprocess = PostProcessing(properties, conf)
```
`properties`: `dict`: properties from `PreProcessing`.  
```Python
postprocess.runPostProcess()
```
Run postprocessing for subhalo obtained before.  