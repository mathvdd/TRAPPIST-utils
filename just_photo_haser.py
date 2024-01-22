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

working_dir = '/home/Mathieu/Documents/TRAPPIST/reduced_data/0012P/20231218TN/centering'
# Qfitlim = (3.6, 4.1) # limit in log10 km for the range over which Q is fitted
Qfitlim = [None,None]
fc = None
# fc = {'OH':19,
#       'NH':24,
#       'CN':30,
#       'C3':248,
#       'C2':170}
kitty=False
check_centering = True


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
        shutil.copy(os.path.join(root, 'centering', 'center_comment'), os.path.join(param['tmpout']))
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
        
        def add_comment(row, comment):
            row_comment = comment.loc[comment['file'] == row['file'], 'comment']
            return row_comment.to_string(index = False)
        if check_centering:
            while True:
                centerlist = pd.read_csv(os.path.join(param['tmpout'], 'centerlist'),header=None, sep=' '
                                         ,usecols=[0,2,3,5,10],
                                         names=['file', 'xcent', 'ycent', 'filt', 'ctnmethod'])
                comment = pd.read_csv(os.path.join(param['tmpout'], 'center_comment'),header=None, sep=','
                                         ,usecols=[0,1,2],
                                         names=['file', 'filt', 'comment'])
                centerlist['comment'] = centerlist.apply(lambda row: add_comment(row, comment), axis=1)
                print(centerlist)
                while True:
                    solocomete = False
                    inp = input('\nCheck individual images for centering.\
                    \n   - Relaunch afrhocalcext for all images (r)\
                    \n   - Relaunch afrhocalcext for an individual file (IMINDEX XCENTER YCENTER BOXSIZE)\
                    \n     optional: YMIN YMAX\
                    \n   - Bypass (b)\
                    \n   - Add a comment (c IMINDEX comment)\
                    \n   :')
                    if inp == 'r' or inp == 'R' or inp == 'b' or inp == 'B' or inp.split(' ')[0] == 'c' or inp.split(' ')[0] == 'C':
                        break
                    elif len(inp.split(' ')) == 4 or len(inp.split(' ')) == 6: #check of good format for a one file reduction
                        solocomete = True
                        try:
                            FILE = centerlist.iloc[[inp.split(' ')[0]]].file.values[0]
                        except:
                            solocomete = False
                            print('IMINDEX wrong format')
                        try:
                            XCENTER = float(inp.split(' ')[1])
                            if XCENTER ==0:
                                XCENTER = centerlist.iloc[[inp.split(' ')[0]]].xcent.values[0]
                        except:
                            solocomete = False
                            print('XCENTER wrong format')
                        try:
                            YCENTER = float(inp.split(' ')[2])
                            if YCENTER ==0:
                                YCENTER = centerlist.iloc[[inp.split(' ')[0]]].ycent.values[0]
                        except:
                            solocomete = False
                            print('YCENTER wrong format')
                        try:
                            BOXSIZE = int(inp.split(' ')[3])
                            if BOXSIZE !=0 and BOXSIZE < 3:
                                print('Use a BOXSIZE value <= 3 or 0 to use the exact values of XCENTER and YCENTER')
                                continue
                        except:
                            solocomete = False
                            print('BOXSIZE wrong format')

                        if len(inp.split(' ')) == 6:
                            try:
                                ZMIN = float(inp.split(' ')[4])
                            except:
                                solocomete = False
                                print('ZMIN wrong format')
                            try:
                                ZMAX = float(inp.split(' ')[5])
                            except:
                                solocomete = False
                                print('ZMAX wrong format')
                        else:
                            ZMIN = None
                            ZMAX = None
                                
                        if solocomete == True: # break if good format
                            break
                    else:
                        print('wrong input')
                if solocomete == True:
                    trap_reduction.clafrhocalcext(param['iraf'], str(pixsize), FILE, str(XCENTER), str(YCENTER), str(BOXSIZE), conda=conda)
                    trap_plot.plot_centering_profile(param['tmpout'], solocomet=True, comet_name=comet,zmin=ZMIN,zmax=ZMAX, kitty=kitty)
                    if kitty == True:
                                try:
                                    os.system(f'kitty +kitten icat "{param["tmpout"] + "/" +  FILE[:-5] + "_centering.png"}"')
                                except:
                                    print('kitty command failed')
                    continue
                if inp == 'r' or inp == 'R':
                    print('relaunching afrhocalcext')
                    trap_reduction.clafrhocalcext(param['iraf'], pixsize[1], str(0), str(0), str(0), str(0), conda=conda)
                    trap_plot.plot_centering_profile(param['tmpout'], comet_name=comet,kitty=kitty)
                    if kitty == True:
                                try:
                                    os.system(f'for f in {param["tmpout"]}/*_centering.png ; do kitty +kitten icat "$f" ; done')
                                except:
                                    print('kitty command failed')
                    continue
                elif inp == 'b' or inp == 'B':
                    break
                elif inp.split(' ')[0] == 'c' or inp.split(' ')[0] == 'C':
                    try:
                        FILE = centerlist.iloc[[inp.split(' ')[1]]].file.values[0]
                    except:
                        print('Comment wrong format')
                        continue
                    str_comment = inp.split(' ', 2)[-1]
                    comment.loc[comment['file'] == FILE, 'comment'] = str_comment
                    comment.to_csv(os.path.join(param['tmpout'], 'center_comment'), index=False, sep=",", header=False)
        
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
