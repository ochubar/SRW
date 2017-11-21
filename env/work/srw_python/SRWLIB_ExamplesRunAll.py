#############################################################################
# SRWLIB Examples: Runs all SRW Examples
# v 0.01
#############################################################################

from __future__ import print_function #Python 2.7 compatibility
import os

print('SRWLIB Python Examples: run all')
print('***Running all SRW Python examples one by one***')

exFileNames = [
    'SRWLIB_Example01.py',
    'SRWLIB_Example01_kick_matr.py',
    'SRWLIB_Example02.py',
    'SRWLIB_Example03.py',
    'SRWLIB_Example04.py',
    'SRWLIB_Example05.py',
    'SRWLIB_Example06.py',
    'SRWLIB_Example06_PETRA.py',
    'SRWLIB_Example07.py',
    'SRWLIB_Example08.py',
    'SRWLIB_Example09.py',
    'SRWLIB_Example10.py',
    'SRWLIB_Example11.py',
    'SRWLIB_Example12.py',
    'SRWLIB_Example13.py',
    'SRWLIB_Example14.py',
    'SRWLIB_Example15.py',
    'SRWLIB_Example16.py',
    ]

for i in range(len(exFileNames)):
    scriptName = exFileNames[i]
    print('Starting', scriptName + '.', 'When done testing, close its windows to continue.')
    os.system('python ' + scriptName)
