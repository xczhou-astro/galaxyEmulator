import os
import argparse
import shutil

parser = argparse.ArgumentParser()

parser.add_argument('-w', '--workspace', type=str, default='workspace')

args = parser.parse_args()

workspace = args.workspace

os.makedirs(f'{workspace}',exist_ok=True)

src = 'Data/config/config_main.ini'
tar = os.path.join(workspace, 'config.ini')

def read_config(file):
    config = {}
    with open(file) as fp:
        for line in fp:
            line = line.strip()
            if line.startswith('#') or len(line) == 0:
                continue
            key, val = "".join(line.split()).split('=')
            
            if len(val) == 0:
                continue
            
            key = key
            
            config[key] = val
            
            if val == 'True':
                config[key] = True
            elif val == 'False':
                config[key] = False
            else:
                pass
    return config

if not os.path.exists(tar):
    shutil.copy(src, tar)
    print(f'Initialization finished! Please edit config.ini in {workspace}/')

else:
    config = read_config(tar)
    src = 'Data/config/config_survey.ini'
    if config['postProcessing']:
        surveys = [inc for inc in "".join(config['surveys'].split()).split(',')]
        for survey in surveys:
            tar = os.path.join(workspace, f'config_{survey}.ini')
            shutil.copy(src, tar)
            print(f'config_{survey} created! Please edit config_{survey} in {workspace}/')
    
    