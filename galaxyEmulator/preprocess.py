import numpy as np
import h5py
import os
import illustris_python as ill
from astropy.cosmology import Planck15
import astropy.units as u
import subprocess
from .utils import *

class Galaxy:

    def __init__(self, position, mass, subhaloPos, a, h):
        self.pos = position * a / h - subhaloPos # kpc
        self.mass = mass * 10**10 # Msun
        

class PreProcess:
    
    def __init__(self, config):
        self.config = config
        self.snapnum = np.int32(self.config['snapNum'])
        self.cosmology = Planck15
        self.snapz = self.__get_snapz()
        self.h = self.cosmology.h
        self.a = 1 / (1 + self.snapz)
        self.workingDir = self.config['workingDir']
        self.dataDir = self.config['dataDir']

    # def __get_snapz(self):
    #     tng = self.config['filePath'].split('/')[-1]
    #     aux_filename = self.config['dataDir'] + '/auxiliary/' + tng + '.txt'
    #     aux = np.loadtxt(aux_filename, skiprows=12, usecols=(0, 1, 2))
    #     snapz = aux[:, 1][list(np.int32(aux[:, 0])).index(self.snapnum)]
    #     return snapz

    def __get_snapz(self):
        snap = h5py.File(os.path.join(self.config['filePath'],
                                      f'snapdir_{self.snapnum:03d}/snap_{self.snapnum:03d}.0.hdf5'), 'r')
        snapz = dict(snap['Header'].attrs.items())['Redshift']
        return snapz
        
    def __read_subhalos(self):
        snap_subhalos = ill.groupcat.loadSubhalos(self.config['filePath'], self.snapnum)
        return snap_subhalos

    def get_subhalos(self):

        minStellarMass = np.float32(self.config['minStellarMass'])
        if self.config['maxStellarMass'] == 'inf':
            maxStellarMass = np.inf
        else:
            maxStellarMass = np.float32(self.config['maxStellarMass'])
        
        snap_subhalos = self.__read_subhalos()
        stellarMass = snap_subhalos['SubhaloMassType'][:, 4] / self.h # 1e10 Msun

        subhalo_indices = np.where((stellarMass > minStellarMass) \
                                    & (stellarMass < maxStellarMass) \
                                    & (snap_subhalos['SubhaloFlag'] == 1)\
                                    & (snap_subhalos['SubhaloParent'] == 0))[0]
        self.subhaloIDs = subhalo_indices
        self.stellarMasses = stellarMass[self.subhaloIDs]
        halfMassRad = snap_subhalos['SubhaloHalfmassRadType'][:, 4] * self.a / self.h # kpc
        self.halfMassRadii = halfMassRad[self.subhaloIDs]
        self.subhaloNum = self.subhaloIDs.shape[0]
        subhaloPos = snap_subhalos['SubhaloPos'] * self.a / self.h
        self.subhaloPos = subhaloPos[self.subhaloIDs]

        print(f'{self.subhaloNum} subhalos in snapshot {self.snapnum} in stellar mass from {minStellarMass} to {maxStellarMass} [10^10 M_sun]')


    def get_subhaloIDs(self):
        return self.subhaloIDs

    def get_stellarMasses(self):
        return self.stellarMasses * 1e10 * u.Msun

    def subhalo(self, subhaloID):

        self.idx = list(self.subhaloIDs).index(subhaloID)
        self.id = self.subhaloIDs[self.idx]
        self.mass = self.stellarMasses[self.idx]
        self.radius = self.halfMassRadii[self.idx]
        self.pos = self.subhaloPos[self.idx]

        print(f'Stellar Mass of Subhalo {self.id} is {self.mass} [10^10 M_sun]')

    def __get_particles(self):
        
        print('Retrieving Stellar and Gas particles.')

        self.boxLength = self.radius * np.float32(self.config['boxLengthScale'])
        self.boxLength = np.min([self.boxLength, np.float32(self.config['maxBoxLength'])])

        region = self.boxLength / 2.

        # a = 1 / (1 + self.snapz)

        fields = ['GFM_InitialMass', 'GFM_Metallicity', 'GFM_StellarFormationTime',
                'Coordinates']
        
        starPart = ill.snapshot.loadSubhalo(self.config['filePath'], self.snapnum,
                                            self.id, 'star', fields)
        mask = np.where(starPart['GFM_StellarFormationTime'] > 0)[0]
        for key in starPart.keys():
            if key == 'count':
                pass
            else:
                starPart[key] = starPart[key][mask]
        
        starPart['age'] = fage(self.snapz) - fage(1/starPart['GFM_StellarFormationTime'] - 1) # in Myr
        g = Galaxy(starPart['Coordinates'], starPart['GFM_InitialMass'], self.pos, a=self.a, h=self.h)


        ageThreshold = np.float32(self.config['ageThreshold'])

        part = {}
        mask = np.where((np.abs(g.pos[:, 0]) < region)\
                    & (np.abs(g.pos[:, 1]) < region)\
                    & (np.abs(g.pos[:, 2]) < region)\
                    & (starPart['age'] < ageThreshold))[0]
        size = mask.shape[0]
        print('Star-forming star particles:', size)

        part['x'] = g.pos[:, 0][mask] # in kpc
        part['y'] = g.pos[:, 1][mask]
        part['z'] = g.pos[:, 2][mask]
        part['smoothLength'] = np.array([0.74] * size) # softening length from TNG website
        part['sfr'] = g.mass[mask] / np.float32(eval(self.config['ratioSFR']))
        part['Z'] = starPart['GFM_Metallicity'][mask]
        # res['compactness'] = np.array([10**np.float32(config['logCompactness'])] * size) 
        part['compactness'] = np.array([np.float32(self.config['logCompactness'])] * size)
        pressure = (10**np.float32(self.config['logPressure']) * const.k_B) * (u.J * u.cm**-3).to(u.J * u.m**-3) # J * m**3 == Pa
        pressure = pressure.value
        part['pressure'] = np.array([pressure] * size) 
        part['covering'] = np.array([np.float32(self.config['coveringFactor'])] * size)

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

        info = np.column_stack((part['x'], part['y'], part['z'],
                            part['smoothLength'], part['sfr'],
                            part['Z'], part['compactness'],
                            part['pressure'], part['covering'])) # 9 params

        np.savetxt(self.workingDir + '/starforming_stars.txt', info, header=header)

        part = {}
        mask = np.where((np.abs(g.pos[:, 0]) < region)\
                        & (np.abs(g.pos[:, 1]) < region)\
                        & (np.abs(g.pos[:, 2]) < region)\
                        & (starPart['age'] > ageThreshold))[0]
        size = mask.shape[0]
        print('Quenched star particles:', size)
        part['x'] = g.pos[:, 0][mask] # in kpc
        part['y'] = g.pos[:, 1][mask] # in kpc
        part['z'] = g.pos[:, 2][mask] # in kpc
        part['smoothLength'] = np.array([0.74] * size) # in kpc (need investigation)
        part['initialMass'] = g.mass[mask]
        part['Z'] = starPart['GFM_Metallicity'][mask]
        part['age'] = starPart['age'][mask]

        header = 'quenched particles\n' \
                    + '\n' \
                    + 'Column 1: x-coordinate (kpc)\n'\
                    + 'Column 2: y-coordinate (kpc)\n'\
                    + 'Column 3: z-coordinate (kpc)\n'\
                    + 'Column 4: smoothing length (kpc)\n'\
                    + 'Column 5: initial mass (Msun)\n'\
                    + 'Column 6: metallicity (1)\n'\
                    + 'Column 7: age (Myr)\n'
        
        info = np.column_stack((part['x'], part['y'], part['z'],
                                part['smoothLength'], part['initialMass'],
                                part['Z'], part['age']))
        np.savetxt(self.workingDir + '/quenched_stars.txt', info, header=header)

        part = {}
        if self.config['includeDust']:
            try:
                fields = ['GFM_Metallicity', 'Coordinates', 'Masses',
                        'InternalEnergy', 'StarFormationRate', 'ElectronAbundance']
                gasPart = ill.snapshot.loadSubhalo(self.config['filePath'], 
                        self.snapnum, self.id, 'gas', fields)
                gasPart['Temperature'] = u2temp(gasPart['InternalEnergy'], 
                                                gasPart['ElectronAbundance'])
            
                gas = Galaxy(gasPart['Coordinates'], gasPart['Masses'], 
                            self.pos, a=self.a, h=self.h)
                mask = np.where((np.abs(gas.pos[:, 0]) < region)\
                                & (np.abs(gas.pos[:, 1]) < region)\
                                & (np.abs(gas.pos[:, 2]) < region)\
                                & (gasPart['StarFormationRate'] > 0)\
                                & (gasPart['Temperature'] < np.float32(self.config['temperatureThreshold'])))[0]
                size = mask.shape[0]
                print('Dust particles:', size)
                part['x'] = gas.pos[:, 0][mask] # in kpc
                part['y'] = gas.pos[:, 1][mask]
                part['z'] = gas.pos[:, 2][mask]
                part['smoothLength'] = np.array([0.18] * size)
                part['mass'] = gas.mass[mask]
                part['Z'] = gasPart['GFM_Metallicity'][mask]
            except:
                part['x'] = np.array([])
                part['y'] = np.array([])
                part['z'] = np.array([])
                part['smoothLength'] = np.array([])
                part['mass'] = np.array([])
                part['Z'] = np.array([])
        else:
            part['x'] = np.array([])
            part['y'] = np.array([])
            part['z'] = np.array([])
            part['smoothLength'] = np.array([])
            part['mass'] = np.array([])
            part['Z'] = np.array([])

        info = np.column_stack((part['x'], part['y'], part['z'],
                                part['smoothLength'], part['mass'],
                                part['Z']))

        np.savetxt(self.workingDir + '/dusts.txt', info, header=header)

    def __get_properties(self, z):
        properties = {}
        properties['subhaloID'] = self.id
        properties['redshift'] = z
        distance = self.cosmology.luminosity_distance(z).value
        properties['lumiDis'] = distance
        if self.config['postProcessing']:
            surveys = split(self.config['surveys'])
            spatialResols = []
            for survey in surveys:
                filters = split(self.config[f'filters_{survey}'])
                numfilters = len(filters)
                if self.config[f'resolFromPix_{survey}']:
                    pixelScales = split(self.config[f'pixelScales_{survey}'], float)
                    pixelScales = extend(pixelScales, numfilters)
                    resol = [distance * 10**6 * np.deg2rad(ang / 3600) for ang in pixelScales]
                    ps_in_sr = [(ps**2 * u.arcsec**2).to(u.sr).value for ps in pixelScales]
                    
                else:
                    resol = split(self.config[f'resolution_{survey}'], float)
                    resol = extend(resol, numfilters)
                    pixelScales = [np.rad2deg(res / (distance * 10**6)) * 3600 for res in resol]
                    ps_in_sr = [(ps**2 * u.arcsec**2).to(u.sr).value for ps in pixelScales]
                
                properties[f'resolution_{survey}'] = resol
                properties[f'angleRes_{survey}'] = pixelScales
                properties[f'ps_in_sr_{survey}'] = ps_in_sr

                spatialResols += resol

                numExposure = split(self.config[f'numExposure_{survey}'], int)
                numExposure = extend(numExposure, numfilters)
                properties[f'numExposure_{survey}'] = numExposure

                exposureTime = split(self.config[f'exposureTime_{survey}'], float)
                exposureTime = extend(exposureTime, numfilters)
                properties[f'exposureTime_{survey}'] = exposureTime
                
                aperture = np.float32(self.config[f'aperture_{survey}'])
                properties[f'aperture_{survey}'] = aperture
                
                properties[f'numfilters_{survey}'] = numfilters
                properties[f'filters_{survey}'] = filters
            
            spatialResol = np.min(spatialResols)

            for survey in surveys:
                ratios = [spatialResol / res for res in properties[f'resolution_{survey}']]
                properties[f'ratios_{survey}'] = ratios

        else:
            spatialResol = np.float32(self.config['spatialResol'])

        properties['baseRes'] = spatialResol

        numViews = np.int32(self.config['numViews'])
        properties['numViews'] = numViews
        randomViews = self.config['randomViews']
        if randomViews:
            inclinations = np.random.uniform(0, 180, numViews).tolist()
            azimuths = np.random.uniform(-360, 360, numViews).tolist()
        else:
            inclinations = split(self.config['inclinations'], float)
            azimuths = split(self.config['azimuths'], float)
        
        properties['inclinations'] = inclinations
        properties['azimuths'] = azimuths
        fovscale = np.float32(self.config['FoVboxLengthRatio'])
        fov = fovscale * self.boxLength * 10**3
        properties['FoV'] = fov
        
        return properties


    def __create_ski(self):

        print('Creating .ski file.')

        mode = self.config['simulationMode']

        ski_file = os.path.join(self.dataDir,  f'ski_templates/{mode}_template.ski')

        with open(ski_file, 'r') as file:
            data = file.read()

        z = np.float32(self.config['fixedRedshift'])

        data = data.replace('redshift="0.008"', f'redshift="{z}"')

        self.properties = self.__get_properties(z=z)

        SEDFamily = {'BC03': 'BruzualCharlotSEDFamily',
                 'FSPS': 'FSPSSEDFamily'}
        SEDFamily = SEDFamily[self.config['SEDFamily']]

        if SEDFamily == 'FSPSSEDFamily':
            data = data.replace('resolution="High"', '')

        if mode == 'DustEmission':
            dustEmissionType = self.config['dustEmissionType']
            data = data.replace('dustEmissionType="Equilibrium"',
                                f'dustEmissionType="{dustEmissionType}"')
        
        initialMassFunction = self.config['initialMassFunction']
        
        data = data.replace('BruzualCharlotSEDFamily', SEDFamily)
        data = data.replace('Chabrier', initialMassFunction)
        
        numPackets = np.float32(self.config['numPackets'])
        data = data.replace('numPackets="1e7"', f'numPackets="{numPackets}"')
        
        minWavelength = np.float32(self.config['minWavelength'])
        data = data.replace('minWavelength="0.01 micron"', f'minWavelength="{minWavelength} micron"')
        
        maxWavelength = np.float32(self.config['maxWavelength'])
        data = data.replace('maxWavelength="1.2 micron"', f'maxWavelength="{maxWavelength} micron"')

        dustConfig = '<ZubkoDustMix numSilicateSizes="15" numGraphiteSizes="15" numPAHSizes="15"/>'
        dustModel = self.config['dustModel']
        numSilicateSizes = np.int32(self.config['numSilicateSizes'])
        numGraphiteSizes = np.int32(self.config['numGraphiteSizes'])
        numPAHSizes = np.int32(self.config['numPAHSizes'])
        numHydrocarbonSizes = np.int32(self.config['numHydrocarbonSizes'])
        if dustModel == 'ThemisDustMix':
            dustConfigNew = f'<{dustModel} numHydrocarbonSizes="{numHydrocarbonSizes}" numSilicateSizes="{numSilicateSizes}"/>'
        else:
            dustConfigNew = f'<{dustModel} numSilicateSizes="{numSilicateSizes}" numGraphiteSizes="{numGraphiteSizes}" numPAHSizes="{numPAHSizes}"/>'
        data = data.replace(dustConfig, dustConfigNew)

        minLevel = np.int32(self.config['minLevel'])
        maxLevel = np.int32(self.config['maxLevel'])
        data = data.replace('minLevel="8"', f'minLevel="{minLevel}"')
        data = data.replace('maxLevel="12"', f'maxLevel="{maxLevel}"')
        
        spatialRange = self.boxLength / 2. * 10**3 # in pc
        data = data.replace('minX="-5e4 pc"', f'minX="{-spatialRange} pc"')
        data = data.replace('maxX="5e4 pc"', f'maxX="{spatialRange} pc"')
        data = data.replace('minY="-5e4 pc"', f'minY="{-spatialRange} pc"')
        data = data.replace('maxY="5e4 pc"', f'maxY="{spatialRange} pc"')
        data = data.replace('minZ="-5e4 pc"', f'minZ="{-spatialRange} pc"')
        data = data.replace('maxZ="5e4 pc"', f'maxZ="{spatialRange} pc"')
        
        wavelengthGrid = self.config['wavelengthGrid']
        grid_type = {'Linear': 'LinWavelengthGrid',
                'Log': 'LogWavelengthGrid'}
        grid_type = grid_type[wavelengthGrid]
        data = data.replace('LinWavelengthGrid', grid_type)
        
        numWavelengths = np.int32(self.config['numWavelengths'])
        data = data.replace('numWavelengths="1000"', f'numWavelengths="{numWavelengths}"')

        begin_str = '<FullInstrument'
        end_str = '</FullInstrument>'
        offset = len(end_str)
        
        idx_begin = data.index(begin_str)
        idx_end = data.index(end_str) + offset
        
        instrumentInfo = data[idx_begin: idx_end]
        remainingInfo = data[idx_end:]
        
        data = data.replace(instrumentInfo, '')
        
        instrumentInfo = '\n' + instrumentInfo + '\n'
        fovscale = np.float32(self.config['FoVboxLengthRatio'])
        fov = fovscale * self.boxLength * 10**3 # in pc

        spatialResol = self.properties['baseRes']
        inclinations = self.properties['inclinations']
        azimuths = self.properties['azimuths']
        numViews = self.properties['numViews']

        numPixels = int(fov / spatialResol)

        recordComponents = 'true' if self.config['recordComponents'] else 'false'
        
        insert_begin_idx = idx_begin
        for i, (inclination, azimuth) in enumerate(zip(inclinations, azimuths)):
            info = instrumentInfo.replace('view', f'view_{i:02d}')
            info = info.replace('inclination="0 deg"', f'inclination="{inclination} deg"')
            info = info.replace('azimuth="0 deg"', f'azimuth="{azimuth} deg"')
            info = info.replace('fieldOfViewX="1e5 pc"', f'fieldOfViewX="{fov} pc"')
            info = info.replace('fieldOfViewY="1e5 pc"', f'fieldOfViewY="{fov} pc"')
            info = info.replace('numPixelsX="1000"', f'numPixelsX="{numPixels}"')
            info = info.replace('numPixelsY="1000"', f'numPixelsY="{numPixels}"')
            info = info.replace('recordComponents="false"', f'recordComponents="{recordComponents}"')
            
            
            data = data[:insert_begin_idx] + info
            insert_begin_idx = insert_begin_idx + len(info)

        data = data + remainingInfo
        with open(self.workingDir + '/skirt.ski', 'w') as file:
            file.write(data)

        print('------estimate memory usage------')
        print(f'numViews: {numViews}')
        print(f'numSpatialPixels: {numPixels}')
        print(f'numWavelengthPixels: {numWavelengths}')
        factor = 7 if self.config['recordComponents'] else 1
        numPixels = np.float64(numPixels) # avoid overflow
        dataCubeSize = np.int64(numPixels ** 2 * numWavelengths * numViews)
        dataSize_in_GB = np.around(dataCubeSize * 8 * factor * 1e-9, 3)
        print(f'Estimated memory usage: {dataSize_in_GB} GB')

    def prepare(self):
        os.makedirs(self.workingDir, exist_ok=True)
        self.__get_particles()
        self.__create_ski()

    def runSKIRT(self):
        self.__run_skirt()

    def __run_skirt(self):

        print('Running SKIRT')

        base = os.getcwd()
        os.chdir(self.workingDir)
        numThreads = int(self.config['numThreads'])
        if numThreads > 24:
            numThreads = 24
        command = f'skirt -t {numThreads} skirt.ski'
        process = subprocess.Popen(command, shell=True)
        flag = process.wait()

        os.chdir(base)
        return flag
    
    def get_properties(self):
        return self.properties

    # def runSKIRT(self):
    #     self.__get_particles()
    #     self.__create_ski()
    #     self.__run_skirt()