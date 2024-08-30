import os
import argparse
import shutil

parser = argparse.ArgumentParser()

parser.add_argument('--workspace', type=str, default='workspace')

args = parser.parse_args()

workspace = args.workspace

os.makedirs(workspace, exist_ok=True)
src = 'Data/config_templates.ini'
dst = os.path.join(workspace, 'config.ini')

shutil.copy(src, dst)

src = 'Data/run_emulator.py'
dst = os.path.join(workspace, 'run_emulator.py')

shutil.copy(src, dst)

print(f'Initialization complete! Please edit the config.ini in {workspace}!')
