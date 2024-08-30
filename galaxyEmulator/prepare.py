import os
import numpy as np
from astropy.cosmology import Planck15
import astropy.units as u
import astropy.constants as const
import subprocess
import illustris_python as ill
from .utils import *

class Galaxy:
    '''
    Convert coordinates of subhalo particles in 
    the central frame of the subhalo and store
    coordinates and mass in proper units
    '''
    def __init__(self, position, mass, subhaloPos, h=0.6774, a=1):
        
        self.pos = (position - subhaloPos) * a / h # in kpc
        self.mass = (mass * 10**10) / h # in M_sun
        
def get_z_from_snapnum(info_TNG, snapnum):
    '''
    info_TNG: the information of TNG-100 containing redshift and snapshot number
    snapnum: the snapshot ID
    
    return: 
    snap_z: the redshift of snapshot
    '''
    
    snap_z = info_TNG[info_TNG['Snapshot'] == snapnum]['Redshift']
    snap_z = snap_z.values[0]
    return snap_z

def get_snapnum(info_TNG, z):
    '''
    info_TNG: the information of TNG-100 containing redshift and snapshot number
    z: the redshift of mock galaxy

    return: 
    snapnum: the snapshot ID
    snap_z: the snapshot redshift
    '''
    snap_zs = info_TNG['Redshift']
    snap_nums = info_TNG['Snapshot']

    idx = np.where(info_TNG['Redshift'] < z)[0][0]
    snapnum = snap_nums.iloc[idx]
    snap_z = snap_zs.iloc[idx]

    return snapnum, snap_z

def get_subhalos(filePath, snapnum):
    
    '''
    filePath: filePath of TNG simulation files
    snapnum: the snapshot ID
    
    return:
    subhalos: information of subhalos in dict
    '''

    subhalos = ill.groupcat.loadSubhalos(filePath, snapnum)
    return subhalos

def particles_from_tng(idx, snapnum, 
                       snap_z, subhalos, boxLength,
                       config):
    '''
    Generate star and dust particles from snapshot
    
    idx: the index of subhalo
    snapnum: snapshot number
    snap_z: redshift of the snapshot
    subhalos: subhalos read by func: get_subhalos
    boxLength: box length for selecting particles
    config: configuration in config.ini
    '''
    region = boxLength / 2. # in kpc
    
    a = 1 / (1 + snap_z)
    h = Planck15.h
    filePath = config['filePath']
    workingDir = config['workingDir']
    
    fields = ['GFM_InitialMass', 'GFM_Metallicity', 'GFM_StellarFormationTime',
              'Coordinates']
    starPart = ill.snapshot.loadSubhalo(filePath, snapnum, idx, 'star', fields)
    mask = np.where(starPart['GFM_StellarFormationTime'] > 0)[0]
    for key in starPart.keys():
        if key == 'count':
            pass
        else:
            starPart[key] = starPart[key][mask]

    starPart['age'] = fage(snap_z) - fage(1/starPart['GFM_StellarFormationTime'] - 1) # in Myr
    g = Galaxy(starPart['Coordinates'], starPart['GFM_InitialMass'], 
               subhalos['SubhaloPos'][idx], h, a)
    
    ageThreshold = np.float32(config['ageThreshold'])
    
    # extract star-forming star particles
    res = {}
    mask = np.where((np.abs(g.pos[:, 0]) < region)\
                    & (np.abs(g.pos[:, 1]) < region)\
                    & (np.abs(g.pos[:, 2]) < region)\
                    & (starPart['age'] < ageThreshold))[0]
    size = mask.shape[0]
    print('Star-forming star particles:', size)
    res['x'] = g.pos[:, 0][mask] # in kpc
    res['y'] = g.pos[:, 1][mask]
    res['z'] = g.pos[:, 2][mask]
    res['smoothLength'] = np.array([0.74] * size) # softening length from TNG website
    res['sfr'] = g.mass[mask] / np.float32(eval(config['ratioSFR']))
    res['Z'] = starPart['GFM_Metallicity'][mask]
    res['compactness'] = np.array([10**np.float32(config['logCompactness'])] * size)
    pressure = (10**np.float32(config['logPressure']) * const.k_B) * (u.J * u.cm**-3).to(u.J * u.m**-3)
    pressure = pressure.value
    res['pressure'] = np.array([pressure] * size) 
    res['covering'] = np.array([np.float32(config['coveringFactor'])] * size)
    
    header = 'star forming particles\n' \
                + '\n' \
                + 'Column 1: x-coordinate (kpc)\n' \
                + 'Column 2: y-coordinate (kpc)\n' \
                + 'Column 3: z-coordinate (kpc)\n' \
                + 'Column 4: smoothing length (kpc)\n' \
                + 'Column 5: star formation rate (Msun/yr)\n' \
                + 'Column 6: metallicity (1)\n' \
                + 'Column 7: compactness (1)\n' \
                + 'Column 8: pressure (J/m3)\n' \
                + 'Column 9: coveringFactor (1)\n'

    info = np.column_stack((res['x'], res['y'], res['z'],
                            res['smoothLength'], res['sfr'],
                            res['Z'], res['compactness'],
                            res['pressure'], res['covering'])) # 9 params

    np.savetxt(workingDir + '/starforming_stars.txt', info, header=header)
    
    # extract quenched star particles
    res = {}
    mask = np.where((np.abs(g.pos[:, 0]) < region)\
                    & (np.abs(g.pos[:, 1]) < region)\
                    & (np.abs(g.pos[:, 2]) < region)\
                    & (starPart['age'] > ageThreshold))[0]
    size = mask.shape[0]
    print('Quenched star particles:', size)
    res['x'] = g.pos[:, 0][mask] # in kpc
    res['y'] = g.pos[:, 1][mask] # in kpc
    res['z'] = g.pos[:, 2][mask] # in kpc
    res['smoothLength'] = np.array([0.74] * size) # in kpc (need investigation)
    res['initialMass'] = g.mass[mask]
    res['Z'] = starPart['GFM_Metallicity'][mask]
    res['age'] = starPart['age'][mask]

    header = 'quenched particles\n' \
                + '\n' \
                + 'Column 1: x-coordinate (kpc)\n'\
                + 'Column 2: y-coordinate (kpc)\n'\
                + 'Column 3: z-coordinate (kpc)\n'\
                + 'Column 4: smoothing length (kpc)\n'\
                + 'Column 5: initial mass (Msun)\n'\
                + 'Column 6: metallicity (1)\n'\
                + 'Column 7: age (Myr)\n'
    
    info = np.column_stack((res['x'], res['y'], res['z'],
                            res['smoothLength'], res['initialMass'],
                            res['Z'], res['age']))
    np.savetxt(workingDir + '/quenched_stars.txt', info, header=header)
    
    res = {}
    if config['includeDust'] == True:
        try:
            fields = ['GFM_Metallicity', 'Coordinates', 'Masses',
                      'InternalEnergy', 'StarFormationRate', 'ElectronAbundance']
            gasPart = ill.snapshot.loadSubhalo(filePath, snapnum, idx, 'gas', fields)
            gasPart['Temperature'] = u2temp(gasPart['InternalEnergy'], 
                                            gasPart['ElectronAbundance'])
        
            gas = Galaxy(gasPart['Coordinates'], gasPart['Masses'], 
                         subhalos['SubhaloPos'][idx], h, a)
            mask = np.where((np.abs(gas.pos[:, 0]) < region)\
                            & (np.abs(gas.pos[:, 1]) < region)\
                            & (np.abs(gas.pos[:, 2]) < region)\
                            & (gasPart['StarFormationRate'] > 0)\
                            & (gasPart['Temperature'] < np.float32(config['temperatureThreshold'])))[0]
            size = mask.shape[0]
            print('DUST particles:', size)
            res['x'] = gas.pos[:, 0][mask] # in kpc
            res['y'] = gas.pos[:, 1][mask]
            res['z'] = gas.pos[:, 2][mask]
            res['smoothLength'] = np.array([0.18] * size)
            res['mas'] = gas.mass[mask]
            res['Z'] = gasPart['GFM_Metallicity'][mask]
        except:
            res['x'] = np.array([])
            res['y'] = np.array([])
            res['z'] = np.array([])
            res['smoothLength'] = np.array([])
            res['mas'] = np.array([])
            res['Z'] = np.array([])

    else:
        res['x'] = np.array([])
        res['y'] = np.array([])
        res['z'] = np.array([])
        res['smoothLength'] = np.array([])
        res['mas'] = np.array([])
        res['Z'] = np.array([])

    header = 'Gas particles\n' \
                + '\n' \
                + 'Column 1: x-coordinate (kpc)\n'\
                + 'Column 2: y-coordinate (kpc)\n'\
                + 'Column 3: z-coordinate (kpc)\n'\
                + 'Column 4: smoothing length (kpc)\n'\
                + 'Column 5: gas mass (Msun)\n'\
                + 'Column 6: metallicity (1)\n'

    info = np.column_stack((res['x'], res['y'], res['z'],
                            res['smoothLength'], res['mas'],
                            res['Z']))
    np.savetxt(workingDir + '/dusts.txt', info, header=header)
    
def modify_ski_file(z, boxLength, config):
    '''
    z: redshift of source
    boxLength: boxLength for calculating radiative transfer, in kpc
    workingDir: working directory for skirt
    config: configurations in config.ini
    
    return:
    property: property useful for next step
    '''
    property = {}
    mode = config['simulationMode']
    ski_file = f'../Data/ski_templates/{mode}.ski.template'
    
    with open(ski_file, 'r') as file:
        data = file.read()
        
    # if config['generationMode'] == 'Local':
    #     universe = '<FlatUniverseCosmology redshift="0.008" reducedHubbleConstant="0.6774" matterDensityFraction="0.3075"/>'
    #     data = data.replace(universe, '<LocalUniverseCosmology/>')
    #     distance = np.float32(config['localDistance'])
    #     data = data.replace('distance="0 Mpc"', f'distance="{distance} Mpc"')
    #     property['redshift'] = 0
    #     property['lumiDis'] = distance

    # else:
    #     data = data.replace('redshift="0.008"', f'redshift="{z}"')
    #     property['redshift'] = z
    #     distance = Planck15.luminosity_distance(z).value
    #     property['lumiDis'] = distance
        
    data = data.replace('redshift="0.008"', f'redshift="{z}"')
    property['redshift'] = z
    distance = Planck15.luminosity_distance(z).value
    property['lumiDis'] = distance
    
    SEDFamily = {'BC03': 'BruzualCharlotSEDFamily',
                 'FSPS': 'FSPSSEDFamily'}
    SEDFamily = SEDFamily[config['SEDFamily']]

    if SEDFamily == 'FSPSSEDFamily':
        data = data.replace('resolution="High"', '')
    
    initialMassFunction = config['initialMassFunction']
    
    data = data.replace('BruzualCharlotSEDFamily', SEDFamily)
    data = data.replace('Chabrier', initialMassFunction)
    
    numPackets = np.float32(config['numPackets'])
    data = data.replace('numPackets="1e7"', f'numPackets="{numPackets}"')
    
    minWavelength = np.float32(config['minWavelength'])
    data = data.replace('minWavelength="0.01 micron"', f'minWavelength="{minWavelength} micron"')
    
    maxWavelength = np.float32(config['maxWavelength'])
    data = data.replace('maxWavelength="1.2 micron"', f'maxWavelength="{maxWavelength} micron"')

    minLevel = np.int32(config['minLevel'])
    maxLevel = np.int32(config['maxLevel'])
    data = data.replace('minLevel="8"', f'minLevel="{minLevel}"')
    data = data.replace('maxLevel="12"', f'maxLevel="{maxLevel}"')
    
    spatialRange = boxLength / 2. * 10**3 # in pc
    data = data.replace('minX="-5e4 pc"', f'minX="{-spatialRange} pc"')
    data = data.replace('maxX="5e4 pc"', f'maxX="{spatialRange} pc"')
    data = data.replace('minY="-5e4 pc"', f'minY="{-spatialRange} pc"')
    data = data.replace('maxY="5e4 pc"', f'maxY="{spatialRange} pc"')
    data = data.replace('minZ="-5e4 pc"', f'minZ="{-spatialRange} pc"')
    data = data.replace('maxZ="5e4 pc"', f'maxZ="{spatialRange} pc"')
    
    wavelengthGrid = config['wavelengthGrid']
    grid_type = {'Linear': 'LinWavelengthGrid',
            'Log': 'LogWavelengthGrid'}
    grid_type = grid_type[wavelengthGrid]
    data = data.replace('LinWavelengthGrid', grid_type)
    
    numWavelengths = config['numWavelengths']
    data = data.replace('numWavelengths="1000"', f'numWavelengths="{numWavelengths}"')
    
    numViews = np.int32(config['numViews'])
    randomViews = config['randomViews']
    if randomViews == True:
        inclinations = np.random.uniform(0, 180, numViews)
        azimuths = np.random.uniform(-360, 360, numViews)
    else:
        
        inclinations = split(config['inclinations'], float)
        azimuths = split(config['azimuths'], float)

    # property['inclinations'] = inclinations.tolist()
    # property['azimuths'] = azimuths.tolist()
    property['inclinations'] = inclinations
    property['azimuths'] = azimuths
    
    begin_str = '<FullInstrument'
    end_str = '</FullInstrument>'
    offset = len(end_str)
    
    idx_begin = data.index(begin_str)
    idx_end = data.index(end_str) + offset
    
    instrumentInfo = data[idx_begin: idx_end]
    remainingInfo = data[idx_end:]
    
    data = data.replace(instrumentInfo, '')
    
    instrumentInfo = '\n' + instrumentInfo + '\n'
    fovscale = np.float32(config['FoVboxLengthRatio'])
    fov = fovscale * boxLength * 10**3 # in pc
    property['FoV'] = fov
    
    numPixels = int(fov / np.float32(config['resolution']))
    
    insert_begin_idx = idx_begin
    for i, (inclination, azimuth) in enumerate(zip(inclinations, azimuths)):
        info = instrumentInfo.replace('view', f'view_{i:02d}')
        info = info.replace('inclination="0 deg"', f'inclination="{inclination} deg"')
        info = info.replace('azimuth="0 deg"', f'azimuth="{azimuth} deg"')
        info = info.replace('fieldOfViewX="1e5 pc"', f'fieldOfViewX="{fov} pc"')
        info = info.replace('fieldOfViewY="1e5 pc"', f'fieldOfViewY="{fov} pc"')
        info = info.replace('numPixelsX="1000"', f'numPixelsX="{numPixels}"')
        info = info.replace('numPixelsY="1000"', f'numPixelsY="{numPixels}"')
        
        data = data[:insert_begin_idx] + info
        insert_begin_idx = insert_begin_idx + len(info)

    data = data + remainingInfo
    
    workingDir = config['workingDir']
    with open(workingDir + '/skirt.ski', 'w') as file:
        file.write(data)

    return property

def run_skirt(config):
    
    base = os.getcwd()
    os.chdir(config['workingDir'])
    numThreads = int(config['numThreads'])
    if numThreads > 24:
        numThreads = 24
    command = f'skirt -t {numThreads} skirt.ski'
    process = subprocess.Popen(command, shell=True)
    flag = process.wait()
    
    os.chdir(base)

    return flag
    