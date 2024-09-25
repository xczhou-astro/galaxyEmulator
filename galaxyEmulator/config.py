import os
import numpy as np
import re
from termcolor import colored
from .utils import *

def get_config(default_config_file='../Data/config_template.ini'):
    config = {}
    with open(default_config_file) as fp:
        for line in fp:
            if line.startswith('#') or line == '\n':
                continue

            key, val = "".join(line.strip().split()).split('=')
            # key, val = line.strip().split(' = ', 1)
            if len(val) == 0:
                continue
                
            config[key] = val
            
            if val == 'True':
                config[key] = True
            elif val == 'False':
                config[key] = False
            else:
                pass
            
    return config

def save_config(config):
    
    keys = list(config.keys())
    
    with open('config.ini', 'w') as f:
        for key in keys:
            f.write(key + ' = ' + str(config[key]))
            
    print(f'config.ini saved')
    return config
    
# def issue(flag, message=None):
#     '''
#     flag must be passed as string
#     '''
#     if flag in locals() or flag in globals():
#         print(colored(message, 'red'))
#     else:
#         pass
    

def check_config(config):
    
    def issue(message):
        print(colored(message, 'red'))
        
    print(colored('Configuration conflicts are indicated in', 'green'), colored('RED', 'red'))
    # print(colored('Configuration conflicts are indicated in RED', 'green'))
    
    global flag_count
    flag_count = 0
    
    if not os.path.exists(config['filePath']):
        issue('filePath not found.')
        flag_count += 1
    
    posprocessing_path = '/'.join(config['filePath'].split('/')[:-1]) + '/postprocessing'
    if not os.path.exists(posprocessing_path):
        issue('postprocessing path not found.')
        flag_count += 1
        
    if (config['simulationMode'] != 'ExtinctionOnly') \
            & (config['simulationMode'] != 'DustEmission'):
        issue('simulationMode unrecognized.')
        flag_count += 1
        
    if (config['simulationMode'] == 'DustEmission') \
            & (config['includeDust'] == False):
        issue('includeDust must be True if simulationMode is DustEmission.')
        flag_count += 1

    if (config['simulationMode'] == 'DustEmission'):
        if (config['dustEmissionType'] != 'Equilibrium')\
            & (config['dustEmissionType'] != 'Stochastic'):
            issue('dustEmissionType unrecognized.')
            flag_count += 1

    if (config['dustModel'] != 'ZubkoDustMix') \
            & (config['dustModel'] != 'DraineLiDustMix')\
                & (config['dustModel'] != 'ThemisDustMix'):
        issue('dustModel unrecognized.')
        flag_count += 1
        
    if (config['wavelengthGrid'] != 'Linear') \
            & (config['wavelengthGrid'] != 'Log'):
        issue('wavelengthGrid unrecognized.')
        flag_count += 1
    
    if not config['randomViews']:
        numViews = np.int32(config['numViews'])
        inclinations = split(config['inclinations'], float)
        azimuths = split(config['azimuths'], float)
        
        if (len(inclinations) != numViews) | (len(azimuths) != numViews):
            issue('number of inclinations or azimuth inconsistent with numViews.')
            flag_count += 1

        inc_in_cond = [True if 0 <= inc <= 180 else False for inc in inclinations]
        azi_in_cond = [True if -360 <= azi <= 360 else False for azi in azimuths]

        if not all(inc_in_cond):
            issue('inclinations must be in 0 to 180 degree.')
            flag_count += 1

        if not all(azi_in_cond):
            issue('azimuths must be in -360 to 360 degree.')
            flag_count += 1
            
    if config['SEDFamily'] == 'BC03':
        if (config['initialMassFunction'] != 'Chabrier') \
                & (config['initialMassFunction'] != 'Salpeter'):
            issue('initialMassFunction unrecognized for BC03 SEDFamily')
            flag_count += 1

    elif config['SEDFamily'] == 'FSPS':
        if (config['initialMassFunction'] != 'Chabrier') \
                & (config['initialMassFunction'] != 'Salpeter') \
                    & (config['initialMassFunction'] != 'Kroupa'):
            issue('initialMassFunction unrecognized for FSPS SEDFamily')
            flag_count += 1
    else:
        issue('SEDFamily unrecognized.')
        flag_count += 1

    def exist(key):
        if key in config:
            return True
        else:
            global flag_count
            issue(f'{key} is not provided.')
            flag_count += 1

    def match(base_key, key):
        if exist(base_key):
            base = split(config[base_key])
            values = split(config[key])
            if len(values) == len(base):
                return True
            else:
                global flag_count
                issue(f'number of {key} inconsistent with number of {base_key}.')
                flag_count += 1

    def logic(base_key, key, relation='match'):
        global flag_count
        if relation == 'match':
            if config[base_key] != config[key]:
                issue(f'{key} must be {config[base_key]} if {base_key} is {config[base_key]}.')
                flag_count += 1
        elif relation == 'diff':
            if config[base_key] == config[key]:
                issue(f'{key} must be {not config[base_key]} if {base_key} is {config[base_key]}.')
                flag_count += 1
        else:
            pass
            

    if config['postProcessing']:
        exist('saveDatacube')
        if exist('filters'):
            filters = split(config['filters'])
            surveys = [name.split('.')[0] for name in filters]
            filterNames = [name.split('.')[1] for name in filters]

            pivots = []
            for survey, filter in zip(surveys, filterNames):
                if os.path.exists(f'../Data/filters/{survey}'):
                    filfiles = os.listdir(f'../Data/filters/{survey}')
                    filfilters = [name.split('.')[0] for name in filfiles]
                    if not filter in filfilters:
                        issue(f'Throughput {survey}.{filter} not found.')
                        issue('Please use add_filters.py or directly add them in Data/filters.')
                        flag_count += 1
                    else:
                        pivots.append(calc_pivot(survey, filter))
                        
                else:
                    issue(f'Throughput {survey}.{filter} not found.')
                    issue('Please use add_filters.py or directly add them in Data/filters')
                    flag_count += 1

            if len(pivots):
                max_pivot = np.max(pivots)
                maxWavelength = np.float32(config['maxWavelength']) * 10**4

                if max_pivot > maxWavelength:
                    issue('maxWavelength smaller than max wavelength for filters.')
                    flag_count += 1
                
                if (maxWavelength > 2 * 10**4) & (config['simulationMode'] != 'DustEmission'):
                    issue('filters reaching infrared, simulationMode should be DustEmission.')
                    flag_count += 1

            for survey, filter in zip(surveys, filterNames):
                idx_survey = [i for i, sur in enumerate(surveys) if sur == survey]
                survey_filters = [filterNames[i] for i in idx_survey]

                exist(f'exposureTime_{survey}')
                exist(f'mirrorSize_{survey}')

                if config['resolFromPix']:
                    if config['includePSF']:
                        if config['PSFFromFile']:
                            if os.path.exists(f'../Data/PSFs/{survey}'):
                                psffiles = os.listdir(f'../Data/PSFs/{survey}/')
                                psffilters = [name.split('.')[0] for name in psffiles]
                                if not filter in psffilters:
                                    issue(f'PSF {survey}.{filter} not found.')
                                    issue('Please add PSF file in Data/PSFs.')
                                    flag_count += 1
                            else:
                                issue(f'PSF {survey}.{filter} not found.')
                                issue('Please add PSF file in Data/PSFs')
                                flag_count += 1
                        else:
                            if exist('PSFFWHM'):
                                match('filters', 'PSFFWHM')
                    
                    if config['includeBackground']:
                        if exist('backgroundSigma'):
                            match('filters', 'backgroundSigma')
                    
                    if config[f'imgDisplay_{survey}']:
                        if config[f'RGBImg_{survey}']:
                            if exist(f'RGBFilters_{survey}'):
                                RGBFilters = split(config[f'RGBFilters_{survey}'])
                                if not all([True if fil in survey_filters else False for fil in RGBFilters]):
                                    issue(f'RGBFilters not in filter set for {survey}.')
                                    flag_count += 1

                        else:
                            if exist(f'displayFilter_{survey}'):
                                if not config[f'displayFilter_{survey}'] in survey_filters:
                                    issue(f'displayFilter not in filter set for {survey}.')
                                    flag_count += 1
                        
                        exist(f'displaySEDxlogscale_{survey}')

                else:
                    logic('resolFromPix', 'includePSF', 'match')
                    logic('resolFromPix', 'includeBackground', 'match')
                    exist('resolution')
                    if config[f'imgDisplay_{survey}']:
                        if config[f'RGBImg_{survey}']:
                            if exist(f'RGBFilters_{survey}'):
                                RGBFilters = split(config[f'RGBFilters_{survey}'])
                                if not all([True if fil in survey_filters else False for fil in RGBFilters]):
                                    issue(f'RGBFilters not in filter set for {survey}.')
                                    flag_count += 1

                        else:
                            if exist(f'displayFilter_{survey}'):
                                if not config[f'displayFilter_{survey}'] in survey_filters:
                                    issue(f'displayFilter not in filter set for {survey}.')
                                    flag_count += 1
                        
                        exist(f'displaySEDxlogscale_{survey}')

        if exist('pixelScales'):
            match('filters', 'pixelScales')

        if exist('numExposure'):
            match('filters', 'numExposure')
    
    else:
        logic('postProcessing', 'saveDatacube', 'diff')
                    
    fixedRedshift = np.float32(config['fixedRedshift'])
    if (config['snapNum'] == '99') & (fixedRedshift == 0):
        issue('fixedRedshift should be larger than 0 if snapNum == 99.')
        flag_count += 1
        
    if np.int32(config['numThreads']) > 24:
        issue('No speed improvement with numThreads larger than, 24, falling back to 24.')
    
    return flag_count