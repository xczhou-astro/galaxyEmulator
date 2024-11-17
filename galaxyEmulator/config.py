import os
from termcolor import colored
import numpy as np
from .utils import *
from pathlib import Path

class Configuration:

    def __init__(self, surveys=None):
        if surveys is not None:
            if isinstance(surveys, str):
                self.surveys = split(surveys)
            else:
                self.surveys = surveys
        else:
            self.surveys = []
        
        dataDir = '../Data'
        self.main_config_template = self.__read_config(
            os.path.join(dataDir, 'config/config_main.ini')
        )
        self.survey_config_template = self.__read_config(
            os.path.join(dataDir,'config/config_survey.ini')
        )
        
        # self.flag_count = 0

    def __load_config(self):
        self.config = self.__read_config('config.ini')
        config_survey_path = Path('.')
        config_survey_files = list(config_survey_path.glob('config_*.ini'))
        detected_surveys = [
            file.name.split('_')[1].split('.')[0] for file in config_survey_files
        ]

        for survey in detected_surveys:
            if survey not in self.surveys:
                self.add_survey(survey)
            config_survey = self.__read_config(f'config_{survey}.ini', survey)

            self.config = self.config | config_survey
        
        self.__modify_main_config()

    def add_survey(self, surveys):

        if surveys is not None:
            if isinstance(surveys, str):
                surveys_add = split(surveys)
            else:
                surveys_add = surveys
        else:
            surveys_add = []

        added_surveys = []
        for survey in surveys_add:
            if survey in self.surveys:
                print(f'{survey} already exists.')
            else:
                added_surveys.append(survey)

        self.surveys += added_surveys

    def  __read_config(self, file, insrument=None):
        if insrument is not None:
            suffix = f'_{insrument}'
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
    
    def __modify_main_config(self):
        if len(self.surveys) == 0:
            self.config['postProcessing'] = 'False'
            self.config['surveys'] = ''
        else:
            self.config['postProcessing'] = 'True'
            self.config['surveys'] = ','.join(self.surveys)

    def __save_main_config(self, main_config):
        keys = list(main_config.keys())
        with open('config.ini', 'w') as f:
            for key in keys:
                f.write(f'{key} = {main_config[key]}')
                f.write('\n')

    def __create_survey_config(self, survey):
        survey_config = {}
        for key, val in self.survey_config_template.items():
            survey_config[key + f'_{survey}'] = val
        return survey_config
    
    def __save_survey_config(self, survey_config, survey):
        keys = list(survey_config.keys())
        with open(f'config_{survey}.ini', 'w') as f:
            for key in keys:
                key_strip = key.split(f'_{survey}')[0]
                f.write(f'{key_strip} = {survey_config[key]}')
                f.write('\n')
 
    def __update_config(self):
        if len(self.surveys) == 0:
            if os.path.exists('config.ini'):
                self.__load_config()
                # self.__modify_main_config()
            else:
                self.config = self.main_config_template
                self.__modify_main_config()
                self.__save_main_config(self.config)
        else:
            if os.path.exists('config.ini'):
                self.__load_config() # self.surveys will be updated
                for survey in self.surveys:
                    if not os.path.exists(f'config_{survey}.ini'):
                        survey_config = self.__create_survey_config(survey)
                        self.__save_survey_config(survey_config, survey)
                        self.config = self.config | survey_config
            
            else:
                self.config = self.main_config_template
                self.__modify_main_config()
                self.__save_main_config(self.config)
                for survey in self.surveys:
                    survey_config = self.__create_survey_config(survey)
                    self.__save_survey_config(survey_config, survey)
                    self.config = self.config | survey_config
    
    def get_config(self):
        self.__update_config()
        
        flags = self.check_config()
        return self.config
    
    def init(self):
        self.__update_config()
        flags = self.check_config()
        self.save_config()
    
    def save_config(self, directory='.'):
        keys = list(self.config.keys())
        if self.surveys is not None:
            for survey in self.surveys:
                keys_survey = [key for key in keys if survey in key]
                filename = os.path.join(directory, f'config_{survey}.ini')
                with open(filename, 'w') as f:
                    for key in keys_survey:
                        key_strip = key.strip(f'_{survey}')
                        f.write(key_strip + ' = ' + str(self.config[key]))
                        f.write('\n')
                        keys.remove(key)
        
        main_config_filename = os.path.join(directory, 'config.ini')
        with open(main_config_filename, 'w') as f:
            for key in keys:
                f.write(key + ' = ' + str(self.config[key]))
                f.write('\n')

    def __issue(self, message):
        print(colored(message, 'red'))

    def __exist(self, key, threshold=None):
        if key in self.config:
            if threshold is not None:
                values = split(self.config[f'{key}'])
                values = np.array([np.float32(val) for val in values])
                check = (values > threshold).all()
                if check == False:
                    self.__issue(f'{key} value must > {threshold}.')
                    self.flag_count += 1
                    return False
                else:
                    return True
            else:
                return True
        
        else:
            self.__issue(f'{key} is not provided.')
            self.flag_count += 1
            return False
        
    def __exist_return(self, key, threshold=None):
        value = self.__exist(key, threshold=threshold)
        if value:
            return self.config[key]
        else:
            return False
        
    def __match(self, base_key, key, homo=False):
        
        if self.__exist(key):
            
            base = split(self.config[base_key])
            values = split(self.config[key])
            if len(values) == len(base):
                return True
            elif (len(values) == 1) & (homo):
                return True
            else:
                self.__issue(f'number of {key} inconsistent with number of {base_key}.')
                self.flag_count += 1
                return False
    
    def __logic(self, base_key, key, relation='match'):
        if relation == 'match':
            if self.config[base_key] != self.config[key]:
                self.__issue(f'{key} must be {self.config[base_key]} if {base_key} is {self.config[base_key]}.')
                self.flag_count += 1
                return False
        elif relation == 'diff':
            if self.config[base_key] == self.config[key]:
                self.__issue(f'{key} must be {not self.config[base_key]} if {base_key} is {self.config[base_key]}.')
                self.flag_count += 1
                return False
        else:
            return True

    def check_config(self):

        self.flag_count = 0

        print(colored('Conflicts in config are indicated in', 'green'), colored('RED.', 'red'))


        if self.__exist('filePath'):
            if not os.path.exists(self.config['filePath']):
                self.__issue('filePath not found.')
                self.flag_count += 1
        
            posprocessing_path = '/'.join(self.config['filePath'].split('/')[:-1]) + '/postprocessing'
            if not os.path.exists(posprocessing_path):
                self.__issue('postprocessing path not found.')
                self.flag_count += 1
                
        self.__exist('workingDir')
        
        simulationMode = self.__exist_return('simulationMode')
        if simulationMode:
            if (simulationMode != 'ExtinctionOnly') & (simulationMode != 'DustEmission'):
                self.__issue('simulationMode unrecognized.')
                self.flag_count += 1

            elif simulationMode == 'DustEmission':
                includeDust = self.__exist_return('includeDust')
                if not includeDust:
                    self.__issue('includeDust must be True if simulationMode is DustEmission.')
                    self.flag_count += 1
                else:
                    dustEmissionType = self.__exist_return('dustEmissionType')
                    if dustEmissionType:
                        if (dustEmissionType != 'Equilibrium') & (dustEmissionType != 'Stochastic'):
                            self.__issue('dustEmissionType unrecognized.')
                            self.flag_count += 1
            else:
                pass
            
        
        dustModel = self.__exist_return('dustModel')
        if dustModel:
            if (dustModel != 'ZubkoDustMix') & (dustModel != 'DraineLiDustMix') & (dustModel != 'ThemisDustMix'):
                self.__issue('dustModel unrecognized.')
                self.flag_count += 1
                
        self.__exist('minWavelength')
        # exist('maxWavelength')
        self.__exist('boxLengthScale')
        self.__exist('maxBoxLength')
                
        wavelengthGrid = self.__exist_return('wavelengthGrid')
        if wavelengthGrid:
            if (wavelengthGrid != 'Linear') & (wavelengthGrid != 'Log'):
                self.__issue('wavelengthGrid unrecognized.')
                self.flag_count += 1
                
        self.__exist('numWavelengths')
        self.__exist('minLevel')
        self.__exist('maxLevel')
        self.__exist('numPackets')
        
        if self.__exist('numViews'):
            if not self.__exist_return('randomViews'):
                if self.__exist('inclinations') & self.__exist('azimuths'):
                    numViews = np.int32(self.config['numViews'])
                    inclinations = split(self.config['inclinations'], float)
                    azimuths = split(self.config['azimuths'], float)
                    
                    if (len(inclinations) != numViews) | (len(azimuths) != numViews):
                        self.__issue('number of inclinations or azimuth inconsistent with numViews.')
                        self.flag_count += 1
                    
                    inc_in_cond = [True if 0 <= inc <= 180 else False for inc in inclinations]
                    azi_in_cond = [True if -360 <= azi <= 360 else False for azi in azimuths]

                    if not all(inc_in_cond):
                        self.__issue('inclinations must be in 0 to 180 degree.')
                        self.flag_count += 1

                    if not all(azi_in_cond):
                        self.__issue('azimuths must be in -360 to 360 degree.')
                        self.flag_count += 1
                    
        SEDFamily = self.__exist_return('SEDFamily')
        if SEDFamily:
            if (SEDFamily != 'BC03') & (SEDFamily != 'FSPS'):
                self.__issue('SEDFamily unrecognized.')
                self.flag_count += 1
                
            elif SEDFamily == 'BC03':
                initialMassFunction = self.__exist_return('initialMassFunction')
                if initialMassFunction:
                    if (initialMassFunction != 'Chabrier') & (initialMassFunction != 'Salpeter'):
                        self.__issue('initialMassFunction unrecognized for BC03 SEDFamily.')
                        self.flag_count += 1
                    
            elif SEDFamily == 'FSPS':
                initialMassFunction = self.__exist_return('initialMassFunction')
                if initialMassFunction:
                    if  (initialMassFunction != 'Chabrier') & (initialMassFunction != 'Salpeter') & (initialMassFunction != 'Kroupa'):
                        self.__issue('initialMassFunction unrecognized for FSPS SEDFamily.')
                        self.flag_count += 1
                        
        self.__exist('minStellarMass')
        self.__exist('maxStellarMass')
        self.__exist('FoVboxLengthRatio')
                        
        if self.__exist_return('displaySED'):
            self.__exist('displaySEDxlogscale')
        
        if not self.__exist_return('postProcessing'):
            self.__exist('spatialResol')
            self.__logic('postProcessing', 'saveDataCube', 'diff')
            
        else:
            self.__exist('saveDataCube')
            surveys = self.__exist_return('surveys')
            pivots = []
            if surveys:
                surveys = split(surveys)
                for survey in surveys:
                    filters = self.__exist_return(f'filters_{survey}')
                    if filters:
                        filters = split(filters)
                        for filter in filters:
                            if not os.path.exists(f'../Data/filters/{survey}/{filter}.fil'):
                                self.__issue(f'Throughput file {survey}/{filter} not found! Please specify correct filters!')
                                # self.__issue('Please use add_filters.py or directly add them in Data/filters.')
                                self.flag_count += 1
                            else:
                                pivots.append(calc_pivot(survey, filter))
                                
                        if self.__exist(f'numExposure_{survey}'):
                            self.__match(f'filters_{survey}', f'numExposure_{survey}', True)
                        
                        if self.__exist(f'exposureTime_{survey}'):
                            self.__match(f'filters_{survey}', f'exposureTime_{survey}', True)
                        
                        self.__exist(f'aperture_{survey}')
                        
                        if self.__exist_return(f'resolFromPix_{survey}'):
                            pixelScales = self.__exist_return(f'pixelScales_{survey}', 0)
                            if pixelScales:
                                self.__match(f'filters_{survey}', f'pixelScales_{survey}', True)
                            
                            if self.__exist_return(f'includePSF_{survey}'):
                                PSFFromFile = self.__exist_return(f'PSFFromFile_{survey}')
                                if not PSFFromFile:
                                    PSFFWHM = self.__exist_return(f'PSFFWHM_{survey}')
                                    if PSFFWHM:
                                        self.__match(f'filters_{survey}', f'PSFFWHM_{survey}', True)
                                        # PSFModel = self.__exist_return(f'PSFModel_{survey}')
                                        # if (PSFModel != 'Moffat') & (PSFModel != 'Gaussian'):
                                        #     self.__issue('PSFModel unrecognized.')
                                        #     self.flag_count += 1


                                else:
                                    if os.path.exists(f'../Data/PSFs/{survey}'):
                                        psffiles = os.listdir(f'../Data/PSFs/{survey}')
                                        psffilters = [name.split('.')[0] for name in psffiles]
                                        for filter in filters:
                                            if not filter in psffilters:
                                                self.__issue(f'PSF file {survey}/{filter} not found.')
                                                self.__issue(f'Please add PSF file in Data/PSFs.')
                                                self.flag_count += 1
                                    else:
                                        self.__issue(f'PSF folder of {survey} not found.')
                                        self.flag_count += 1
                            
                            if self.__exist_return(f'includeBkg_{survey}'):
                                bkgNoise = self.__exist_return(f'bkgNoise_{survey}')
                                if bkgNoise:
                                    self.__match(f'filters_{survey}', f'bkgNoise_{survey}')
                            
                            if self.__exist_return(f'imgDisplay_{survey}'):
                                if self.__exist_return(f'RGBImg_{survey}'):
                                    RGBFilters = self.__exist_return(f'RGBFilters_{survey}')
                                    RGBFilters = split(RGBFilters)
                                    if not all([True if fil in filters else False for fil in RGBFilters]):
                                        self.__issue(f'RBGFilters not in {survey} filter set.')
                                        self.flag_count += 1
                                else:
                                    displayFilter = self.__exist_return(f'displayFilter_{survey}')
                                    if displayFilter:
                                        if not displayFilter in filters:
                                            self.__issue(f'displayFilter not in {survey} filter set.')
                                            self.flag_count += 1
                        else:
                            resolution = self.__exist_return(f'resolution_{survey}', 0)
                            if resolution:
                                self.__match(f'filters_{survey}', f'resolution_{survey}', True)
                                
                            if self.__exist_return(f'imgDisplay_{survey}'):
                                if self.__exist_return(f'RGBImg_{survey}'):
                                    RGBFilters = self.__exist_return(f'RGBFilters_{survey}')
                                    RGBFilters = split(RGBFilters)
                                    if not all([True if fil in filters else False for fil in RGBFilters]):
                                        self.__issue(f'RBGFilters not in {survey} filter set.')
                                        self.flag_count += 1
                                else:
                                    displayFilter = self.__exist_return(f'displayFilter_{survey}')
                                    if displayFilter:
                                        if not displayFilter in filters:
                                            self.__issue(f'displayFilter not in {survey} filter set.')
                                            self.flag_count += 1
                                            
            if len(pivots) & (self.__exist('maxWavelength')):
                max_pivot = np.max(pivots)
                maxWavelength = np.float32(self.config['maxWavelength']) * 10**4
                
                if max_pivot > maxWavelength:
                    self.__issue('maxWavelength smaller than max wavelength for filters.')
                    self.flag_count += 1
                
                if (maxWavelength > 2 * 10**4) & (self.config['simulationMode'] != 'DustEmission'):
                    self.__issue('filters reaching infrared, simulationMode should be DustEmission.')
                    self.flag_count += 1
            
                                            
        self.__exist('snapNum')
        self.__exist('fixedRedshift', 0)
        
        # numGeneration = self.__exist_return(f'numGeneration')
        # if numGeneration:
        #     if numGeneration == 'all':
        #         pass
        #     elif numGeneration.isdigit():
        #         number = np.int32(numGeneration)
        #         if 'subhaloIDs' in self.config:
        #             subhaloIDs = split(self.config['subhaloIDs'])
        #             if number != len(subhaloIDs):
        #                 self.__issue('number of subhaloIDs inconsistent with numGeneration.')
        #                 self.flag_count += 1
        #     else:
        #         self.__issue('numGeneration unrecognized.')
        #         self.flag_count += 1
        
        if self.__exist('numThreads'):
            if np.int32(self.config['numThreads']) > 24:
                self.__issue('No speed improvement with numThreads larger than, 24, falling back to 24.')

        self.__exist('recordComponents')
        self.__exist('ageThreshold')
        self.__exist('logCompactnessMean')
        self.__exist('logCompactnessStd')
        self.__exist('massFraction')
        self.__exist('PDRClearingTimescale')
        self.__exist('logPressure')
        self.__exist('temperatureThreshold')
        DISMModel = self.__exist_return('DISMModel')
        if DISMModel:
            if (DISMModel != 'Camps_2016') & (DISMModel != 'Torrey_2012'):
                self.__issue('DISMModel unrecognized.')
                self.flag_count += 1
        
        self.__exist('numSilicateSizes')
        self.__exist('numGraphiteSizes')
        self.__exist('numPAHSizes')
        self.__exist('numHydrocarbonSizes')

        if self.flag_count == 0:
            print('No conflicts in config')
        
        self.config['flags'] = self.flag_count

        return self.flag_count