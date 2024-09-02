import sys
import os

galaxyEmulator_path = '..'
sys.path.append(galaxyEmulator_path)
if not os.path.exists(os.path.join(galaxyEmulator_path, 'galaxyEmulator')):
    print('Please edit the system path')
    sys.exit(0)

from galaxyEmulator.emulator import emulator

if __name__ == '__main__':
    
    emulator('config.ini')