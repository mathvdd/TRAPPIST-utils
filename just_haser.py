#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 16 19:34:46 2022

@author: Mathieu
"""
import trap_reduction
import trap_plot
from trapconfig import param
import os
import shutil

comet = 'CK19L030'
night = '2022-01-08'
obs = 'TN'

reduced_dir = os.path.join(param['reduced'], comet, night.replace('-','') + obs)
profiles_dir = os.path.join(reduced_dir, "profiles")
haser_dir = os.path.join(reduced_dir, "haser")
garbage_dir = os.path.join(reduced_dir, "probably_garbage")

if os.path.exists(param['tmpout']):
    shutil.rmtree(param['tmpout'])
os.mkdir(param['tmpout'])

print('--- copying files ---')
for path, subdirs, files in os.walk(profiles_dir):
    for file in files:
        shutil.copy(os.path.join(path, file), os.path.join(param['tmpout'], file))
shutil.copy(os.path.join(haser_dir, 'inputhaser-BC'), os.path.join(param['tmpout'], 'inputhaser-BC'))
shutil.copy(os.path.join(garbage_dir, 'ephem.brol'), os.path.join(param['tmpout'], 'ephem.brol'))

conda = True if param['conda'] == 'True' else False #wether to use 'source activate to launch cl or not'
print('--- launching hasercalctest ---')
while True:
    trap_reduction.clhasercalctest(param['iraf'], 'yes', conda=conda)
    trap_plot.plot_haserprofile(param['tmpout'],comet_name=comet)
    while True:
        inp = input('relaunch hasercalctext (r) or overwrite results in reduced dir (c)? [r/c]')
        if inp == 'r' or inp == 'R' or inp == 'c' or inp == 'C':
            break
    if inp == 'r' or inp == 'R':
        print('relaunching hasercalctest')
        continue
    elif inp == 'c' or inp == 'C':
        print('overwriting results')
        for path2, subdirs2, files2 in os.walk(param['tmpout']):
            for file in files2:
                if 'haser' in file and 'tmp' not in file:
                    shutil.copy(os.path.join(path2, file), os.path.join(haser_dir, file))
                    print('copied', file, "in reduced dir")
        break