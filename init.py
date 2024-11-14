import os
import argparse
from galaxyEmulator.utils import split
from galaxyEmulator.config import Config

parser = argparse.ArgumentParser()

parser.add_argument('-w', '--workspace', type=str, default='workspace')
parser.add_argument('-s', '--surveys', type=str, default='None')

args = parser.parse_args()

workspace = args.workspace
surveys = args.surveys

os.makedirs(f'{workspace}',exist_ok=True)

configuration = Config(surveys=surveys, dataDir='Data')
conf = configuration.init(workspace)

print(f'Configuration files are created in {workspace}. Please edit them!')