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
import datetime
import pandas as pd

comet = 'CK21A010'
Qfitlim = (3.5, 4.1)
fc = {'OH':5,
      'NH':20,
      'CN':25,
      'C3':190,
      'C2':170}
# fc_default = {'OH':5,
#       'NH':20,
#       'CN':25,
#       'C3':190,
#       'C2':170}

def rewrite_fc_in_haserinput(fc):
    inputhaser_path = os.path.join(param["tmpout"], 'inputhaser-BC')
    inputhaser = pd.read_csv(inputhaser_path, header=None, sep='\s+')
    for index, row in inputhaser.iterrows():
        rad = pd.read_csv(os.path.join(param["tmpout"], 'rad_' + row[0]), header=None, sep='\s+')
        filt = rad.iloc[0,14]
        inputhaser.iloc[index,2] = fc.get(filt)
    inputhaser.to_csv(os.path.join(param["tmpout"], 'inputhaser-BC'), index=False, header=False, sep = ' ')

def haser_reduce_1night(comet, night, obs, Qfitlim, check=True):
    reduced_dir = os.path.join(param['reduced'], comet, night.replace('-','') + obs)
    profiles_dir = os.path.join(reduced_dir, "profiles")
    haser_dir = os.path.join(reduced_dir, "haser")
    garbage_dir = os.path.join(reduced_dir, "probably_garbage")
    print(f'Initiating haser for {comet}, {night}, {obs} with range 10E{Qfitlim} km')
    if os.path.exists(param['tmpout']):
        shutil.rmtree(param['tmpout'])
    os.mkdir(param['tmpout'])
    # print('--- copying files ---')
    for path, subdirs, files in os.walk(profiles_dir):
        for file in files:
            shutil.copy(os.path.join(path, file), os.path.join(param['tmpout'], file))
    shutil.copy(os.path.join(haser_dir, 'inputhaser-BC'), os.path.join(param['tmpout'], 'inputhaser-BC'))
    shutil.copy(os.path.join(garbage_dir, 'ephem.brol'), os.path.join(param['tmpout'], 'ephem.brol'))
    rewrite_fc_in_haserinput(fc)
    
    conda = True if param['conda'] == 'True' else False #wether to use 'source activate to launch cl or not'
    print('--- launching hasercalctest ---')
    if check == True:
        while True:
            trap_reduction.clhasercalctest(param['iraf'], arg='yes', Qproflow=Qfitlim[0], Qprofhigh=Qfitlim[1], conda=conda)
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
    else:
        trap_reduction.clhasercalctest(param['iraf'], arg='yes', Qproflow=Qfitlim[0], Qprofhigh=Qfitlim[1], conda=conda)
        trap_plot.plot_haserprofile(param['tmpout'],comet_name=comet)
        for path2, subdirs2, files2 in os.walk(param['tmpout']):
            for file in files2:
                if 'haser' in file and 'tmp' not in file:
                    shutil.copy(os.path.join(path2, file), os.path.join(haser_dir, file))
                    print('copied', file, "in reduced dir")

comet = 'CK21A010'
night = '2021-12-31'
obs = 'TS'     
haser_reduce_1night(comet, night, obs, Qfitlim)

# if __name__ == "__main__":
#     dt = datetime.datetime.now()
#     comet_dir = os.path.join(param['reduced'], comet)
#     for path, subdirs, files in os.walk(comet_dir):
#         if (path.split('/')[-1] == 'haser') and os.path.isfile(
#                 os.path.join(path, 'inputhaser-BC')):
#             night = (path.split('/')[-2][:4] +'-'+ path.split('/')[-2][4:6] +'-'+ path.split('/')[-2][6:8])
#             obs = path.split('/')[-2][8:10]
#             haser_reduce_1night(comet, night, obs, Qfitlim, check=False)
#     print('Executed in ', datetime.datetime.now() - dt)
