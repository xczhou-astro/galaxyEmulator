# README #

galaxyEmulator for space surveys using SKIRT project

## Dependence
numpy

scipy

matplotlib

pandas

astropy

termcolor

[illustris_python](https://github.com/illustristng/illustris_python)

[SKIRT](https://skirt.ugent.be/root/_home.html)

## Usage

### Add filters
`python add_filters.py CSST.NUV,CSST.u,CSST.g --directory='CSST_filters`

throughput files should be readable by `np.loadtxt()`, and a header is required indicating the unit of wavelength

or 

You can add filters directly to `Data/filters/`, with format like `CSST/filter.fil`

### Initialization

`python init.py `


### Edit config.ini in workspace

`cd workspace`

edit `config.ini`

### Run emulator

`python run_emulator.py`

### Illustris-TNG data structure
```
-- TNG-100

    -- groups_{snapnum}
    
        -- fof_subhalo_tab_{snapnum}.{n}.hdf5

    -- snapdir_{snapnum}
    
        -- snap_{snapnum}.{n}.hdf5


-- postprocessing

    -- offsets
    
        -- offsets_{snapnum}.hdf5
```
