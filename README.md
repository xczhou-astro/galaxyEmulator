# README #

galaxyEmulator for space surveys using SKIRT project

## Dependence

numpy

scipy

matplotlib

pandas

astropy

termcolor

and 

[illustris_python](https://github.com/illustristng/illustris_python)

and 

[SKIRT](https://skirt.ugent.be/root/_home.html)

## Usage

### Add filters
`python add_filters.py CSST.NUV,CSST.u,CSST.g --directory='CSST_filters'`

throughput files should be readable by np.loadtxt(), and a header is required indicating the unit of wavelength

or 

You can add filters directly to Data/filters/, with format like CSST/filter.fil

### Initialization

`python init.py [--workspace=workspace]`


### Edit config.ini in workspace

`cd workspace`

edit `config.ini`

and 

`python run_emulator.py`
