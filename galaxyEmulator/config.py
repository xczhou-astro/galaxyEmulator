import os
from termcolor import colored
import numpy as np
from .utils import *

class Configuration:

    def __init__(self, main_config_file='../Data/config/config.ini'):
        self.main_config = self.__read_config(main_config_file)
        
        survey_config_template = os.path.join(self.main_config['dataDir'], 'config/config_survey.ini')
        self.survey_config_template = self.__read_config(survey_config_template)
        
        self.surveys = self.__retrieve_surveys()
        
        self.flag_count = 0

    def __read_config(self, file, instrument=None):
        if instrument is not None:
            suffix = f'_{instrument}'
        else:
            suffix = ''

        config = {}
        with open(file) as fp:
            for line in fp:
                line = line.strip()
                if line.startswith('#') or len(line) == 0:
                    continue
                key, val = "".join(line.split()).split('=')

                if len(val) == 0:
                    continue
                key = key + suffix
                config[key] = val

                if val == 'True':
                    config[key] = True
                elif val == 'False':
                    config[key] = False
                else:
                    pass
                
        return config
    
    def __retrieve_surveys(self):
        if 'postProcessing' in self.main_config and self.main_config['postProcessing']:
            surveys = split(self.main_config['surveys'])
            return surveys
        else:
            return None

    def __create_survey_config(self):
        survey_config = {}
        if self.surveys is not None:
            for survey in self.surveys:
                for key, value in self.survey_config_template.items():
                    survey_config[key + f'_{survey}'] = value
            
        else:
            pass
        
        return survey_config

    
    def get_config(self):
        if self.surveys is not None:
            for survey in self.surveys:
                self.config = self.main_config
                if os.path.exists(f'config_{survey}.ini'):
                    survey_config = self.__read_config(f'config_{survey}.ini', instrument=survey)
                    self.config = self.config | survey_config
                else:
                    survey_config = self.__create_survey_config()
                    self.config = self.config | survey_config
                    self.__save_config()
        else:
            self.config = self.main_config
            
        self.check_config()
        return self.config
    
    
    def __save_config(self):

        keys = list(self.config.keys())
        if len(self.surveys):
            for survey in self.surveys:
                keys_survey = [key for key in keys if survey in key]
                filename = f'config_{survey}.ini'
                with open(filename, 'w') as f:
                    for key in keys_survey:
                        key_strip = key.strip(f'_{survey}')
                        f.write(key_strip + ' = ' + str(self.config[key]))
                        f.write('\n')
                        keys.remove(key)

        main_config_filename = 'config.ini'
        with open(main_config_filename, 'w') as f:
            for key in keys:
                f.write(key + ' = ' + str(self.config[key]))
                f.write('\n')

    def issue(self, message):
        print(colored(message, 'red'))

    def exist(self, key, threshold=None):
        if key in self.config:
            if threshold is not None:
                values = split(self.config[f'{key}'])
                values = np.array([np.float32(val) for val in values])
                check = (values > threshold).all()
                if check == False:
                    self.issue(f'{key} value must > {threshold}.')
                    self.flag_count += 1
                    return False
                else:
                    return True
            else:
                return True
        
        else:
            self.issue(f'{key} is not provided.')
            self.flag_count += 1
            return False
        
    def exist_return(self, key, threshold=None):
        value = self.exist(key, threshold=threshold)
        if value:
            return self.config[key]
        else:
            return False
        
    def match(self, base_key, key, homo=False):
        
        if self.exist(key):
            
            base = split(self.config[base_key])
            values = split(self.config[key])
            if len(values) == len(base):
                return True
            elif (len(values) == 1) & (homo):
                return True
            else:
                self.issue(f'number of {key} inconsistent with number of {base_key}.')
                self.flag_count += 1
                return False
    
    def logic(self, base_key, key, relation='match'):
        if relation == 'match':
            if self.config[base_key] != self.config[key]:
                self.issue(f'{key} must be {self.config[base_key]} if {base_key} is {self.config[base_key]}.')
                self.flag_count += 1
                return False
        elif relation == 'diff':
            if self.config[base_key] == self.config[key]:
                self.issue(f'{key} must be {not self.config[base_key]} if {base_key} is {self.config[base_key]}.')
                self.flag_count += 1
                return False
        else:
            return True

    def check_config(self):

        # print(colored('Conflicts on config are indicated in', 'green'), colored('RED', 'red'))

        if self.exist('filePath'):
            if not os.path.exists(self.config['filePath']):
                self.issue('filePath not found.')
                self.flag_count += 1
        
            posprocessing_path = '/'.join(self.config['filePath'].split('/')[:-1]) + '/postprocessing'
            if not os.path.exists(posprocessing_path):
                self.issue('postprocessing path not found.')
                self.flag_count += 1
                
        self.exist('workingDir')
        
        simulationMode = self.exist_return('simulationMode')
        if simulationMode:
            if (simulationMode != 'ExtinctionOnly') & (simulationMode != 'DustEmission'):
                self.issue('simulationMode unrecognized.')
                self.flag_count += 1

            elif simulationMode == 'DustEmission':
                includeDust = self.exist_return('includeDust')
                if not includeDust:
                    self.issue('includeDust must be True if simulationMode is DustEmission.')
                    self.flag_count += 1
                else:
                    dustEmissionType = self.exist_return('dustEmissionType')
                    if dustEmissionType:
                        if (dustEmissionType != 'Equilibrium') & (dustEmissionType != 'Stochastic'):
                            self.issue('dustEmissionType unrecognized.')
                            self.flag_count += 1
            else:
                pass
            
        
        dustModel = self.exist_return('dustModel')
        if dustModel:
            if (dustModel != 'ZubkoDustMix') & (dustModel != 'DraineLiDustMix') & (dustModel != 'ThemisDustMix'):
                self.issue('dustModel unrecognized.')
                self.flag_count += 1
                
        self.exist('minWavelength')
        # exist('maxWavelength')
        self.exist('boxLengthScale')
        self.exist('maxBoxLength')
                
        wavelengthGrid = self.exist_return('wavelengthGrid')
        if wavelengthGrid:
            if (wavelengthGrid != 'Linear') & (wavelengthGrid != 'Log'):
                self.issue('wavelengthGrid unrecognized.')
                self.flag_count += 1
                
        self.exist('numWavelengths')
        self.exist('minLevel')
        self.exist('maxLevel')
        self.exist('numPackets')
        
        if self.exist('numViews'):
            if not self.exist_return('randomViews'):
                if self.exist('inclinations') & self.exist('azimuths'):
                    numViews = np.int32(self.config['numViews'])
                    inclinations = split(self.config['inclinations'], float)
                    azimuths = split(self.config['azimuths'], float)
                    
                    if (len(inclinations) != numViews) | (len(azimuths) != numViews):
                        self.issue('number of inclinations or azimuth inconsistent with numViews.')
                        self.flag_count += 1
                    
                    inc_in_cond = [True if 0 <= inc <= 180 else False for inc in inclinations]
                    azi_in_cond = [True if -360 <= azi <= 360 else False for azi in azimuths]

                    if not all(inc_in_cond):
                        self.issue('inclinations must be in 0 to 180 degree.')
                        self.flag_count += 1

                    if not all(azi_in_cond):
                        self.issue('azimuths must be in -360 to 360 degree.')
                        self.flag_count += 1
                    
        SEDFamily = self.exist_return('SEDFamily')
        if SEDFamily:
            if (SEDFamily != 'BC03') & (SEDFamily != 'FSPS'):
                self.issue('SEDFamily unrecognized.')
                self.flag_count += 1
                
            elif SEDFamily == 'BC03':
                initialMassFunction = self.exist_return('initialMassFunction')
                if initialMassFunction:
                    if (initialMassFunction != 'Chabrier') & (initialMassFunction != 'Salpeter'):
                        self.issue('initialMassFunction unrecognized for BC03 SEDFamily.')
                        self.flag_count += 1
                    
            elif SEDFamily == 'FSPS':
                initialMassFunction = self.exist_return('initialMassFunction')
                if initialMassFunction:
                    if  (initialMassFunction != 'Chabrier') & (initialMassFunction != 'Salpeter') & (initialMassFunction != 'Kroupa'):
                        self.issue('initialMassFunction unrecognized for FSPS SEDFamily.')
                        self.flag_count += 1
                        
        self.exist('minStellarMass')
        self.exist('maxStellarMass')
        self.exist('FoVboxLengthRatio')
                        
        if self.exist_return('displaySED'):
            self.exist('displaySEDxlogscale')
        
        if not self.exist_return('postProcessing'):
            self.exist('spatialResol')
            self.logic('postProcessing', 'saveDataCube', 'diff')
            
        else:
            self.exist('saveDataCube')
            surveys = self.exist_return('surveys')
            pivots = []
            if surveys:
                surveys = split(surveys)
                for survey in surveys:
                    filters = self.exist_return(f'filters_{survey}')
                    if filters:
                        filters = split(filters)
                        for filter in filters:
                            if not os.path.exists(f'../Data/filters/{survey}/{filter}.fil'):
                                self.issue(f'Throughput file {survey}/{filter} not found')
                                # self.issue('Please use add_filters.py or directly add them in Data/filters.')
                                self.flag_count += 1
                            else:
                                pivots.append(calc_pivot(survey, filter))
                                
                        if self.exist(f'numExposure_{survey}'):
                            self.match(f'filters_{survey}', f'numExposure_{survey}', True)
                        
                        if self.exist(f'exposureTime_{survey}'):
                            self.match(f'filters_{survey}', f'exposureTime_{survey}', True)
                        
                        self.exist(f'aperture_{survey}')
                        
                        if self.exist_return(f'resolFromPix_{survey}'):
                            pixelScales = self.exist_return(f'pixelScales_{survey}', 0)
                            if pixelScales:
                                self.match(f'filters_{survey}', f'pixelScales_{survey}', True)
                            
                            if self.exist_return(f'includePSF_{survey}'):
                                PSFFromFile = self.exist_return(f'PSFFromFile_{survey}')
                                if not PSFFromFile:
                                    PSFFWHM = self.exist_return(f'PSFFWHM_{survey}')
                                    if PSFFWHM:
                                        self.match(f'filters_{survey}', f'PSFFWHM_{survey}', True)
                                        PSFModel = self.exist_return(f'PSFModel_{survey}')
                                        if (PSFModel != 'Moffat') & (PSFModel != 'Gaussian'):
                                            self.issue('PSFModel unrecognized.')
                                            self.flag_count += 1


                                else:
                                    if os.path.exists(f'../Data/PSFs/{survey}'):
                                        psffiles = os.listdir(f'../Data/PSFs/{survey}')
                                        psffilters = [name.split('.')[0] for name in psffiles]
                                        for filter in filters:
                                            if not filter in psffilters:
                                                self.issue(f'PSF file {survey}/{filter} not found.')
                                                self.issue(f'Please add PSF file in Data/PSFs.')
                                                self.flag_count += 1
                                    else:
                                        self.issue(f'PSF folder of {survey} not found.')
                                        self.flag_count += 1
                            
                            if self.exist_return(f'includeBkg_{survey}'):
                                bkgNoise = self.exist_return(f'bkgNoise_{survey}')
                                if bkgNoise:
                                    self.match(f'filters_{survey}', f'bkgNoise_{survey}')
                            
                            if self.exist_return(f'imgDisplay_{survey}'):
                                if self.exist_return(f'RGBImg_{survey}'):
                                    RGBFilters = self.exist_return(f'RGBFilters_{survey}')
                                    RGBFilters = split(RGBFilters)
                                    if not all([True if fil in filters else False for fil in RGBFilters]):
                                        self.issue(f'RBGFilters not in {survey} filter set.')
                                        self.flag_count += 1
                                else:
                                    displayFilter = self.exist_return(f'displayFilter_{survey}')
                                    if displayFilter:
                                        if not displayFilter in filters:
                                            self.issue(f'displayFilter not in {survey} filter set.')
                                            self.flag_count += 1
                        else:
                            resolution = self.exist_return(f'resolution_{survey}', 0)
                            if resolution:
                                self.match(f'filters_{survey}', f'resolution_{survey}', True)
                                
                            if self.exist_return(f'imgDisplay_{survey}'):
                                if self.exist_return(f'RGBImg_{survey}'):
                                    RGBFilters = self.exist_return(f'RGBFilters_{survey}')
                                    RGBFilters = split(RGBFilters)
                                    if not all([True if fil in filters else False for fil in RGBFilters]):
                                        self.issue(f'RBGFilters not in {survey} filter set.')
                                        self.flag_count += 1
                                else:
                                    displayFilter = self.exist_return(f'displayFilter_{survey}')
                                    if displayFilter:
                                        if not displayFilter in filters:
                                            self.issue(f'displayFilter not in {survey} filter set.')
                                            self.flag_count += 1
                                            
            if len(pivots) & (self.exist('maxWavelength')):
                max_pivot = np.max(pivots)
                maxWavelength = np.float32(self.config['maxWavelength']) * 10**4
                
                if max_pivot > maxWavelength:
                    self.issue('maxWavelength smaller than max wavelength for filters.')
                    self.flag_count += 1
                
                if (maxWavelength > 2 * 10**4) & (self.config['simulationMode'] != 'DustEmission'):
                    self.issue('filters reaching infrared, simulationMode should be DustEmission.')
                    self.flag_count += 1
            
                                            
        self.exist('snapNum')
        self.exist('fixedRedshift', 0)
        
        # numGeneration = self.exist_return(f'numGeneration')
        # if numGeneration:
        #     if numGeneration == 'all':
        #         pass
        #     elif numGeneration.isdigit():
        #         number = np.int32(numGeneration)
        #         if 'subhaloIDs' in self.config:
        #             subhaloIDs = split(self.config['subhaloIDs'])
        #             if number != len(subhaloIDs):
        #                 self.issue('number of subhaloIDs inconsistent with numGeneration.')
        #                 self.flag_count += 1
        #     else:
        #         self.issue('numGeneration unrecognized.')
        #         self.flag_count += 1
        
        if self.exist('numThreads'):
            if np.int32(self.config['numThreads']) > 24:
                self.issue('No speed improvement with numThreads larger than, 24, falling back to 24.')

        self.exist('recordComponents')
        self.exist('ageThreshold')
        self.exist('ratioSFR')
        self.exist('logCompactness')
        self.exist('logPressure')
        self.exist('coveringFactor')
        self.exist('temperatureThreshold')
        self.exist('numSilicateSizes')
        self.exist('numGraphiteSizes')
        self.exist('numPAHSizes')
        self.exist('numHydrocarbonSizes')

        if self.flag > 0:
            print(colored('Conflicts on config are indicated in', 'green'), colored('RED', 'red'))
            
        return self.flag_count