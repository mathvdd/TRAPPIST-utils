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
import get_ephem

# comet = '0073P'
# Qfitlim = (4, 4.6)
# fc = {'OH':5,
#       'NH':20,
#       'CN':25,
#       'C3':190,
#       'C2':170}
# fc_default = {'OH':5,
#       'NH':20,
#       'CN':25,
#       'C3':190,
#       'C2':170}

def rewrite_fc_in_haserinput(fc):
    inputhaser_path = os.path.join(param["tmpout"], 'inputhaser-BC')
    try:
        inputhaser = pd.read_csv(inputhaser_path, header=None, sep='\s+',names=[0,1,2,3])
    except:
        print('WARNING: could not read', inputhaser_path)
        print('trying with error_bad_lines=False')
        inputhaser = pd.read_csv(inputhaser_path, header=None, sep='\s+',names=[0,1,2,3]
                                 ,error_bad_lines=False)
        input(inputhaser)
    for index, row in inputhaser.iterrows():
        print(row)
        rad = pd.read_csv(os.path.join(param["tmpout"], 'rad_' + row[0]), header=None, sep='\s+')
        filt = rad.iloc[0,14]
        inputhaser.iloc[index,2] = fc.get(filt)
        inputhaser.iloc[index,3] = 0
    inputhaser.to_csv(os.path.join(param["tmpout"], 'inputhaser-BC'), index=False, header=False, sep = ' ')

def redo_calib_dat(imdir, comet):
        fitstable = trap_reduction.get_fitstable(imdir)
        # fitstable = fitstable.loc[fitstable['type'].isin(['LIGHT', 'Light Frame'])
        #                           & ~fitstable['file'].isin(['coms.fits', 'comsplus.fits'])]
        # print(fitstable)
        print('--- generating ephemeris ---')
        ephemeris = get_ephem.ephemeris()
        ephemeris.retrieve_param_from_fits(imdir)
        ephemeris.query_input(unique_target=False, target=comet, convert_MPC_Horizon=True)
        ephemeris.generate_ephem_files(param['tmpout'])
        print('--- generating calib file ---')
        ZP_warning = trap_reduction.generate_ZP(param['calib'], ephemeris, fitstable, output_dir=param['tmpout'])
        if ZP_warning == True:
            input("press enter")

def haser_reduce_1night(comet, night, obs, Qfitlim, check=True, redo_ZP=False,kitty=False):
    reduced_dir = os.path.join(param['reduced'], comet, night.replace('-','') + obs)
    profiles_dir = os.path.join(reduced_dir, "profiles")
    images_dir = os.path.join(reduced_dir, "images")
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
    
    if redo_ZP == True:
        redo_calib_dat(imdir= images_dir, comet=comet)
    elif redo_ZP == False:
        shutil.copy(os.path.join(garbage_dir, 'calib.dat'), os.path.join(param['tmpout'], 'calib.dat'))
    # input('end')
    conda = True if param['conda'] == 'True' else False #wether to use 'source activate to launch cl or not'
    print('--- launching hasercalctest ---')
    if check == True:
        while True:
            trap_reduction.clhasercalctest(param['iraf'], arg='yes', Qproflow=Qfitlim[0], Qprofhigh=Qfitlim[1], conda=conda)
            trap_plot.plot_haserprofile(param['tmpout'],comet_name=comet,kitty=kitty)
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
                shutil.copy(os.path.join(param['tmpout'], 'calib.dat'), os.path.join(garbage_dir, 'calib.dat'))
                print('copied calib.dat in garbage_dir')
                break
    else:
        trap_reduction.clhasercalctest(param['iraf'], arg='yes', Qproflow=Qfitlim[0], Qprofhigh=Qfitlim[1], conda=conda)
        trap_plot.plot_haserprofile(param['tmpout'],comet_name=comet,kitty=kitty)
        for path2, subdirs2, files2 in os.walk(param['tmpout']):
            for file in files2:
                if 'haser' in file and 'tmp' not in file:
                    shutil.copy(os.path.join(path2, file), os.path.join(haser_dir, file))
                    print('copied', file, "in reduced dir")
        shutil.copy(os.path.join(param['tmpout'], 'calib.dat'), os.path.join(garbage_dir, 'calib.dat'))
        print('copied calib.dat in garbage_dir')


# night = '2022-01-02'
# obs = 'TS'     
# haser_reduce_1night(comet, night, obs, Qfitlim)

if __name__ == "__main__":
    kitty = None
    redo_ZP=False
    
    Qfitlim = (3.4, 4.3)
    fc = {'OH':16, #2019 L3
          'NH':22,
          'CN':44,
          'C3':235,
          'C2':221,
          'CO+':56,
          'H2O':129}
    
    # Qfitlim = (4.4, 5.0)
    # fc = {'OH':16, #2019 L3
    #       'NH':22,
    #       'CN':31,
    #       'C3':237,
    #       'C2':185,
    #       'CO+':56,
    #       'H2O':129}
    comet = 'CK21A010'
    dt = datetime.datetime.now()
    comet_dir = os.path.join(param['reduced'], comet)
    # comet_dir = '/home/Mathieu/Documents/TRAPPIST/reduced_data/CK19L030/20230220TS'
    for path, subdirs, files in os.walk(comet_dir):
        if (path.split('/')[-1] == 'haser') and os.path.isfile(
                os.path.join(path, 'inputhaser-BC')):
            night = (path.split('/')[-2][:4] +'-'+ path.split('/')[-2][4:6] +'-'+ path.split('/')[-2][6:8])
            obs = path.split('/')[-2][8:10]
            haser_reduce_1night(comet, night, obs, Qfitlim, check=False, redo_ZP=redo_ZP,kitty=kitty)
            if kitty == True:
                try:
                    os.system(f'for f in {param["tmpout"]}/*_haserprofile.png ; do kitty +kitten icat "$f" ; done')
                except:
                    print('kitty command failed')
                        
    print('Executed in ', datetime.datetime.now() - dt)
