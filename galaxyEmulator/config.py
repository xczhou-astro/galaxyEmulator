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
    flag_count = 0
    issue_flags = {}
    
    if not os.path.exists(config['filePath']):
        issue_flags['TNG_path_flag'] = 1
        TNG_path_flag = 1
        flag_count += 1
    
    posprocessing_path = '/'.join(config['filePath'].split('/')[:-1]) + '/postprocessing'
    if not os.path.exists(posprocessing_path):
        postprocessing_path_flag = 1
        issue_flags['postprocessing_path_flag'] = 1
        flag_count += 1
        
    if (config['simulationMode'] != 'ExtinctionOnly') \
            & (config['simulationMode'] != 'DustEmission'):
        issue_flags['simulationMode_flag'] = 1
        simulationMode_flag = 1
        flag_count += 1
        
    if (config['simulationMode'] == 'DustEmission') \
            & (config['includeDust'] == False):
        issue_flags['includeDust'] = 1
        includeDust_flag = 1
        flag_count += 1
        
    if (config['wavelengthGrid'] != 'Linear') \
            & (config['wavelengthGrid'] != 'Log'):
        issue_flags['wavelengthGrid_flag'] = 1
        wavelengthGrid_flag = 1
        flag_count += 1
    
    if not config['randomViews']:
        numViews = np.int32(config['numViews'])
        inclinations = split(config['inclinations'], float)
        azimuths = split(config['azimuths'], float)
        
        if (len(inclinations) != numViews) | (len(azimuths) != numViews):
            issue_flags['angle_flag'] = 1
            angle_flag = 1
            flag_count += 1
            
    if config['SEDFamily'] == 'BC03':
        if (config['initialMassFunction'] != 'Chabrier') \
                & (config['initialMassFunction'] != 'Salpeter'):
            issue_flags['BC03_MF_flag'] = 1
            BC03_MF_flag = 1
            flag_count += 1

    elif config['SEDFamily'] == 'FSPS':
        if (config['initialMassFunction'] != 'Chabrier') \
                & (config['initialMassFunction'] != 'Salpeter') \
                    & (config['initialMassFunction'] != 'Kroupa'):
            issue_flags['FSPS_MF_flag'] = 1
            FSPS_MF_flag = 1
            flag_count += 1
    else:
        issue_flags['SEDFamily_flag'] = 1
        SEDFamily_flag = 1
        flag_count += 1
        
    filters = split(config['filters'])
    surveys = [name.split('.')[0] for name in filters]
    filterNames = [name.split('.')[1] for name in filters]
    
    non_existing_filters = []
    pivots = []
    for survey, filter in zip(surveys, filterNames):
        exist = os.path.exists(f'../Data/filters/{survey}/{filter}.fil')
        if not exist:
            non_existing_filters.append([survey, filter])
            flag_count += 1
        else:
            pivots.append(calc_pivot(survey, filter))
    

    non_existing_psfs = []
    if config['includePSF']:
        if config['PSFFromFile']:
            for survey, filter in zip(surveys, filterNames):
                exist = os.path.exists(f'../Data/PSFs/{survey}/{filter}.npy')
                if not exist:
                    non_existing_psfs.append([survey, filter])
                    flag_count += 1
        else:
            fwhms = split(config['PSFFWHM'], float)
            if len(fwhms) != len(filters):
                issue_flags['numPSF_flag'] = 1
                numPSF_flag = 1
                flag_count += 1
                
    if config['includeBackground']:
        sigmas = split(config['backgroundSigma'], float)
        if len(sigmas) != len(filters):
            issue_flags['numBackground_flag'] = 1
            numBackground_flag = 1
            flag_count += 1
            
            
    max_pivot = np.max(pivots)
    maxWavelength = np.float32(config['maxWavelength']) * 10**3 * 10 # compare in angstrom

    if max_pivot > maxWavelength:
        issue_flags['maxWavelength_flag'] = 1
        maxWavelength_flag = 1
        flag_count += 1
    
    if (maxWavelength > 2 * 10**4) & (config['simulationMode'] != 'DustEmission'):
        issue_flags['reachInfrared_flag'] = 1
        reachInfrared_flag = 1
        flag_count += 1
    
    surveys, numfilters = np.unique(surveys, return_counts=True)
    
    survey_flags = {}
    for survey in surveys:
        survey_filters = [filter for filter in filters if bool(re.match(survey, filter, re.I))]
        filterNames = [name.split('.')[1] for name in survey_filters]

        imgDisplay_key_exist = f'imgDisplay_{survey}' in config

        if imgDisplay_key_exist:
        
            if config[f'imgDisplay_{survey}']:

                rgb_or_one_exist = f'RGBImg_{survey}' in config

                if rgb_or_one_exist:
                
                    if config[f'RGBImg_{survey}']:
                        
                        RGBFilters = split(config[f'RGBFilters_{survey}'])
                        
                        if not len(RGBFilters) == 3:
                            filtersNot3_flag = 1
                            survey_flags[f'filtersNot3_{survey}'] = filtersNot3_flag
                            flag_count += 1
                        else:
                            isinfilters = all([fil in filterNames for fil in RGBFilters])
        
                            if not isinfilters:
                                RGBNotInFilters_flag = 1
                                survey_flags[f'RGBNotInFilters_{survey}'] = RGBNotInFilters_flag
                                flag_count += 1
                                
                    else:
                        displayFilter = config[f'displayFilter_{survey}']
                        if not displayFilter in filterNames:
                            filterNotIn_flag = 1
                            survey_flags[f'filterNotIn_{survey}'] = filterNotIn_flag
                            flag_count += 1
                else:
                    survey_flags[f'rgb_or_one_exist_{survey}'] = 1
                    flag_count += 1

        else:
            survey_flags[f'imgDisplay_exist_{survey}'] = 1
            flag_count += 1
                    
    fixedRedshift = np.float32(config['fixedRedshift'])
    if (config['snapNum'] == '99') & (fixedRedshift == 0):
        issue_flags['fixedRedshift_flag'] = 1
        fixedRedshift_flag = 1
        flag_count += 1
        
    if np.int32(config['numThreads']) > 24:
        issue_flags['numThreads_flag'] = 1
        numThreads_flag = 1
    
    print('config:')
    for key in config.keys():
        print(f'{key} = {config[key]}')

    def issue(flag, message):
        if flag in issue_flags:
            print(colored(message, 'red'))

    if flag_count > 0:
        print(colored('Found several conflictions in config.ini, please revise!!!', 'green'))
    issue('TNG_path_flag', 'filePath not found')
    issue('postprocessing_path_flag', 'postprocessing path not found')
    issue('simulationMode_flag', 'simulationMode unrecognized')
    issue('includeDust_flag', 'includeDust should be True when simulationMode is DustEmission')
    issue('wavelengthGrid_flag', 'wavelengthGrid unrecognized')
    issue('angle_flag', 'number of inclinations or azimuth inconsistent with numViews')
    issue('BC03_MF_flag', 'initialMassFunction unrecognized for BC03 SEDFamily')
    # add checking for fsps family directory
    issue('FSPS_MF_flag', 'initialMassFunction unrecognized for FSPS SEDFamily')
    issue('SEDFamily_flag', 'SEDFamily unrecognized')
    for non_existing_fil in non_existing_filters:
        print(colored(f'Throughput {non_existing_fil[0]}.{non_existing_fil[1]} not found', 'red'))
        print(colored('use add_filters.py or directly add them in Data/filters', 'red'))
    for non_existing_psf in non_existing_psfs:
        print(colored(f'PSF {non_existing_psf[0]}.{non_existing_psf[1]} not found', 'red'))
        print(colored('Please add PSF file in Data/PSFs', 'red'))
    issue('numPSF_flag', 'number of PSFFWHM inconsistent with number of filters')
    issue('numBackground_flag', 'number of backgroundSigma inconsistent with number of filters')
    issue('maxWavelength_flag', 'maxWavelength smaller than max wavelength for filters, should be extented')
    issue('reachInfrared_flag', 'filters reaching infared, simulationMode should be DustEmission')
    for survey in surveys:
        if f'imgDisplay_exist_{survey}' in survey_flags:
            print(colored(f'imgDisplay_{survey} not found', 'red'))
        if f'rgb_or_one_exist_{survey}' in survey_flags:
            print(colored(f'RGBImg_{survey} not found', 'red'))
        if f'filtersNot3_{survey}' in survey_flags:
            print(colored(f'RGB image cannot be created using RGBFilters_{survey} set', 'red'))
        if f'RGBNotInFilters_{survey}' in survey_flags:
            print(colored(f'RGBFilters_{survey} not found in considered filters', 'red'))
        if f'filterNotIn_{survey}' in survey_flags:
            print(colored(f'displayFilter_{survey} not found in considred filters', 'red'))
    issue('fixedRedshift_flag', 'fixedRedshift should larger than 0 when snapNum == 99')
    issue('numThreads_flag', 'no speed improvement with numThreads larger than 24, falling back to 24')
    
    return flag_count
    # need to fill more checks
    