import os
import argparse
from galaxyEmulator.utils import split
from galaxyEmulator.config import Configuration

parser = argparse.ArgumentParser()

parser.add_argument('-w', '--workspace', type=str, default='workspace')
parser.add_argument('-s', '--surveys', type=str, default='None')

args = parser.parse_args()

workspace = args.workspace
surveys = args.surveys

os.makedirs(f'{workspace}',exist_ok=True)

os.chdir(f'{workspace}')

configuration = Configuration(surveys=surveys)
conf = configuration.init()

print(f'Configuration files are created in {workspace}. Please edit them!')