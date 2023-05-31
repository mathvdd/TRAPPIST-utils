#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 16 12:22:20 2023

@author: Mathieu Vander Donckt
mathieu.vanderdonckt@uliege.be
"""

import trap_reduction
import trap_plot
from trapconfig import param
import os
import shutil
import datetime
import pandas as pd
from just_haser import rewrite_fc_in_haserinput

working_dir = '/home/Mathieu/Documents/TRAPPIST/reduced_data/0073P'
Qfitlim = (3.6, 4.1) # limit in log10 km for the range over which Q is fitted
# Qfitlim = (None,None)
# fc = None
fc = {'OH':19,
      'NH':24,
      'CN':30,
      'C3':248,
      'C2':170}
kitty=None

if __name__ == "__main__":
    
    dirlist = []
    for path, subdirs, files in sorted(os.walk(working_dir)):
        if (path.split('/')[-1] == 'centering') and os.path.isfile(
                os.path.join(path, 'centerlist')):
            dirlist.append(path)
    
    print(f'Found {len(dirlist)} directories')
    print(f'Qfitlim = {Qfitlim}')
    print(f'fc = {fc}')
    input(f'!!! this will overwrite any previous reduction in {working_dir}\nPress any key to continue')
    
    dt = datetime.datetime.now()
    conda = True if param['conda'] == 'True' else False
    
    dirlist = []
    for path, subdirs, files in sorted(os.walk(working_dir)):
        if (path.split('/')[-1] == 'centering') and os.path.isfile(
                os.path.join(path, 'centerlist')):
            dirlist.append(path)
    
    count = 0
    for path in sorted(dirlist):
        count += 1
        print(f'Current dir is {path}')
        print(f'{count} / {len(dirlist)}')
        root = path[:-10]
        comet = path.split('/')[-3]
        night = (path.split('/')[-2][:4] +'-'+ path.split('/')[-2][4:6] +'-'+ path.split('/')[-2][6:8])
        obs = path.split('/')[-2][8:10]
        if obs == 'TN':
            pixsize = 0.60
        elif obs == 'TS':
            pixsize = 0.65
        print('--- working with', root,'---')
        
        print('--- Initating tmpdata and tmpout dirs ---')
        if os.path.exists(param['tmpdata']):
            shutil.rmtree(param['tmpdata'])
        os.mkdir(param['tmpdata'])
        if os.path.exists(param['tmpout']):
            shutil.rmtree(param['tmpout'])
        os.mkdir(param['tmpout'])
        
        print('--- copying previous reduction to tmpout ---')
        shutil.copy(os.path.join(root, 'centering', 'centerlist'), os.path.join(param['tmpout']))
        for path2, subdirs2, files2 in os.walk(os.path.join(root, 'images')):
            for file in files2:
                if file.endswith('.fits') and file.startswith('TRAP'):
                    shutil.copy(os.path.join(path2, file), os.path.join(param['tmpout']))
        try:
            shutil.copy(os.path.join(root, 'probably_garbage', 'calib.dat'), os.path.join(param['tmpout']))
        except:
            import just_haser
            just_haser.redo_calib_dat(param['tmpout'],comet=comet)
        shutil.copy(os.path.join(root, 'probably_garbage', 'eph.dat'), os.path.join(param['tmpout']))
        shutil.copy(os.path.join(root, 'probably_garbage', 'ephem.brol'), os.path.join(param['tmpout']))
        try:
            shutil.copy(os.path.join(root, 'haser', 'inputhaser-BC'), os.path.join(param['tmpout']))
            shutil.copy(os.path.join(root, 'haser', 'outputhaser-BC'), os.path.join(param['tmpout']))            
            shutil.copy(os.path.join(root, 'probably_garbage', 'tmp.txt'), os.path.join(param['tmpout']))
            haser = True
        except:
            print('no haser files')
            haser = False

        print('--- starting reduction ---')
        centerlist = pd.read_csv(os.path.join(param['tmpout'], 'centerlist'),header=None, sep=' '
                                              ,usecols=[0,2,3,5,10],
                                              names=['file', 'xcent', 'ycent', 'filt', 'ctnmethod'])
        centerlist.drop_duplicates(subset=['file'],keep='last',inplace=True)
        print(centerlist)
        for index, item in centerlist.iterrows():
            trap_reduction.clafrhocalcext(param['iraf'], str(pixsize), item['file'], str(item['xcent']), str(item['ycent']), str(0), conda=conda)
        trap_plot.plot_centering_profile(param['tmpout'], comet_name=comet,kitty=kitty)
        trap_reduction.clean_afrhotot(param['tmpout'])

        # haser is False
        print('--- launching hasercalctest ---')
        if haser== True:
            haserinput = pd.read_csv(os.path.join(param['tmpout'], 'inputhaser-BC'), header=None, sep="\s+")
            if len(haserinput) > 0: 
                if fc != None:
                    rewrite_fc_in_haserinput(fc)
                haseroutput = pd.read_csv(os.path.join(param['tmpout'], 'outputhaser-BC'), header=None, sep="\s+")
                if Qfitlim[0] == None:
                    Qfitlim[0] = haseroutput[16][0]
                if Qfitlim[1] == None:
                    Qfitlim[1] = haseroutput[17][0]
                trap_reduction.clhasercalctest(param['iraf'], arg='yes', Qproflow=Qfitlim[0], Qprofhigh=Qfitlim[1], conda=conda)
                trap_plot.plot_haserprofile(param['tmpout'],comet_name=comet,kitty=kitty)
        else:
            print('No haserinput filters found, skipping Haser') 

        print('--- copying files in reduced directory ---')
        images_dir = os.path.join(root, "images")
        profiles_dir = os.path.join(root, "profiles")
        afrho_dir = os.path.join(root, "afrho")
        haser_dir = os.path.join(root, "haser")
        garbage_dir = os.path.join(root, "probably_garbage")
        centering_dir = os.path.join(root, "centering")
        for path2, subdirs2, files2 in os.walk(param['tmpout']):
            for file in files2:
                if file.endswith('.fits') and 'tmp' not in file and 'coms' not in file:
                    # shutil.copy(os.path.join(path2, file), os.path.join(images_dir, file))
                    print('not copied', file, "in reduced dir")
                elif 'rad' in file and 'tmp' not in file:
                    shutil.copy(os.path.join(path2, file), os.path.join(profiles_dir, file))
                    print('copied', file, "in reduced dir")
                elif 'afrho' in file and 'tmp' not in file:
                    shutil.copy(os.path.join(path2, file), os.path.join(afrho_dir, file))
                    print('copied', file, "in reduced dir")
                elif 'haser' in file and 'tmp' not in file:
                    shutil.copy(os.path.join(path2, file), os.path.join(haser_dir, file))
                    print('copied', file, "in reduced dir")
                elif 'center' in file and 'tmp' not in file:
                    if file.endswith('.png'):
                        shutil.copy(os.path.join(path2, file), os.path.join(centering_dir, file))
                        print('copied', file, "in reduced dir")
                    else:
                        print('not copied', file, "in reduced dir")
                elif 'palist' in file and 'tmp' not in file:
                    # shutil.copy(os.path.join(path2, file), os.path.join(centering_dir, file))
                    print('not copied', file, "in reduced dir")
                else:
                    shutil.copy(os.path.join(path2, file), os.path.join(garbage_dir, file))
                    print('copied', file, "in reduced dir")
    print('Executed in ', datetime.datetime.now() - dt)
