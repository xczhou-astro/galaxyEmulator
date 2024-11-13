import os
import argparse
import shutil
from galaxyEmulator.utils import split
from galaxyEmulator.config import Config

parser = argparse.ArgumentParser()

parser.add_argument('-w', '--workspace', type=str, default='workspace')
parser.add_argument('-s', '--surveys', type=str, default='CSST')

args = parser.parse_args()

workspace = args.workspace
surveys = args.surveys

os.makedirs(f'{workspace}',exist_ok=True)

configuration = Config('Data', surveys=surveys)
conf = configuration.get_config()

print('Configuration files are created.')