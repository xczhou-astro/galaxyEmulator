import numpy as np
import os
import argparse
from termcolor import colored
import sys
import re
from galaxyEmulator.utils import *

parser = argparse.ArgumentParser()
parser.add_argument('--directory', type=str, default='.',
                   help='where to find these files, default in current directory')
parser.add_argument('filters', type=str, 
                    help='throughput files in {survey}.{filter) format')
args = parser.parse_args()

directory = args.directory
filters = args.filters

filters = split(filters)

surveys = [fil.split('.')[0] for fil in filters]
bands = [fil.split('.')[1] for fil in filters]

dirls = os.listdir(directory)
dirls = [name for name in dirls if not name.startswith('.')]

filenames = []
not_found = []
for survey, bd in zip(surveys, bands):
    pattern = rf'\b{bd}(?![\*+])'
    thr_name = [ls for ls in dirls if (bool(re.search(pattern, ls)))]
    if len(thr_name):
        filenames += thr_name
    else:
        not_found.append([survey, bd])


filter_files = [os.path.join(directory, filename) for filename in filenames]


if len(not_found):
    for info in not_found:
        print(colored(f'filter {info[0]}.{info[1]} not found', 'red'))
    sys.exit(0)


def save_filter(filter_file, filter_name):
         
    wavelength_scale = get_wavelength_scale(filter_file)
    
    transmission = np.loadtxt(filter_file)
        
    transmission[:, 0] = transmission[:, 0] * wavelength_scale
    # add a new feature calculating the pivot wavelength and save
    # used to compare with the maxWavelength for skirt
    
    strings = filter_name.split('.')
    survey = strings[0]
    band = strings[1]
    
    os.makedirs(f'Data/filters/{survey}', exist_ok=True)
    np.savetxt(f'Data/filters/{survey}/{band}.fil', transmission, 
               header='angstrom throughputs')
    
    print(f'{survey}.{band} successfully saved!')
        
for filter_file, filter_name in zip(filter_files, filters):
    save_filter(filter_file, filter_name)

