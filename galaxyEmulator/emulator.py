from .config import *
from .utils import *
from .prepare import *
from .postprocess import *
import sys
import os
import numpy as np
from astropy.cosmology import Planck15
import h5py
from termcolor import colored

def emulator(config_file):
    
    config = get_config(config_file)
    issue_flag = check_config(config)
    if issue_flag > 0:
        print('Please edit config.ini')
        sys.exit(0)
    
    snapnum = np.int32(config['snapNum'])
    
    snap = h5py.File(os.path.join(config['filePath'], 
                     f'snapdir_{snapnum:03d}/snap_{snapnum:03d}.0.hdf5'), 'r')
    snapz = dict(snap['Header'].attrs.items())['Redshift']
    
    subhalos = get_subhalos(config['filePath'], snapnum)
    minStellarMass = np.float32(config['minStellarMass'])
    if config['maxStellarMass'] == 'inf':
        maxStellarMass = np.inf
    else:
        maxStellarMass = np.float32(config['maxStellarMass'])
    
    stellarMass = subhalos['SubhaloMassType'][:, 4] / Planck15.h
    
    subhalo_indices = np.where((stellarMass > minStellarMass) \
                            & (stellarMass < maxStellarMass) \
                            & (subhalos['SubhaloFlag'] == 1)\
                            & (subhalos['SubhaloParent'] == 0))[0]
    
    subhaloNums = subhalo_indices.shape[0]
    print(f'Number of subhalos is {subhaloNums} in snapshot {snapnum} in stellar mass range {minStellarMass} to {maxStellarMass} [10^10 M_sun]')
    if config['numGeneration'] == 'All':
        print('Generate all galaxies')
        indices = subhalo_indices
    else:
        numGen = np.int32(config['numGeneration'])
        print(f'Generate {numGen} galaxies')
        rand = np.random.choice(subhaloNums, numGen, replace=False)
        indices = subhalo_indices[rand]
        
    for subhaloID in indices:
        print(f'Processing subhalo {subhaloID}')
        print(f'Stellar Mass: {stellarMass[subhaloID]} [10^10 M_sun]')
        
        boxLength = subhalos['SubhaloHalfmassRadType'][:, 4][subhaloID] \
            * np.float32(config['boxLengthScale']) / Planck15.h
            
        boxLength = np.min([boxLength, np.float32(config['maxBoxLength'])])
        print('boxLength', np.around(boxLength, 2))
        
        workingDir = config['workingDir']
        os.makedirs(workingDir, exist_ok=True)
        print('Getting stellar or gas particles')
        particles_from_tng(subhaloID, snapnum, snapz, subhalos, boxLength, config)
        fixedRedshift = np.float32(config['fixedRedshift'])
        print('Creating .ski file')
        property = modify_ski_file(fixedRedshift, boxLength, config)
        print('Beginning skirt program')
        flag = run_skirt(config)
        if flag:
            print(colored('run SKIRT error, exit', 'red'))
            sys.exit(0)
        print(f'Postprocessing subhalo {subhaloID}')
        postprocess(subhaloID, property, config)
        print(f'Finshing processing subhalo {subhaloID}')