#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 28 15:07:56 2021

@author: Mathieu Vander Donckt
mathieu.vanderdonckt@uliege.be

Append login.cl with tasks and TRAPDAT paths
comment stty xgterm in login.cl to be able to launch scripts from terminal
with iraf installed in iraf27 environment
login.cl needs to be in the folder where iraf is launched
change line 14 of renamefits to accept direct paths
trapccd in iraf folder
change line 44 of hasercalctewt to redirect to the script directory
"""
import get_ephem
import directory_structure
import os
import trap_reduction
import shutil
import trap_plot
from datetime import datetime
import pandas as pd

dt = datetime.now()

ds = directory_structure.directory_structure()
ephemeris = None # necessary not to ask for the target name every time

# ds.raw = os.path.join(ds.raw, "2020T2darktest")
# ds.reduced = os.path.join(ds.reduced, "2020T2darktest")
ds.raw = "/home/Mathieu/Documents/TRAPPIST/raw_data/CK20T020/TN/20210418"
ds.reduced= "/home/Mathieu/Documents/TRAPPIST/reduced_data/CK20T020"

for path, subdirs, files in os.walk(ds.raw):
    if 'Calibration' not in path and 'Autoflat' not in path and (any('.fits' in file for file in files) or any('.fts' in file for file in files)):
        raw_dir = path
        print("working with", raw_dir)
        while True: 
            print('--- Renaming .fts to .fits ---')
            trap_reduction.renameftsfits(raw_dir)
            print('--- checking for calib files ---')
            fitstable = trap_reduction.get_fitstable(raw_dir)
            check_calib_warning = trap_reduction.check_calib(fitstable)
            if check_calib_warning == True:
                inp = input("add calib files and press enter (b for bypass)")
                if inp == 'b' or inp =="B":
                    break
            else:
                break
        # continue
        print('--- Initating tmpdata and tmpout dirs ---')
        if os.path.exists(ds.tmpdata):
            shutil.rmtree(ds.tmpdata)
        os.mkdir(ds.tmpdata)
        print('--- Renaming files to Trappist format and copying to tmpdata ---')
        trap_reduction.pythrename(raw_dir, ds.tmpdata)
        if os.path.exists(ds.tmpout):
            shutil.rmtree(ds.tmpout)
        os.mkdir(ds.tmpout)
        print('--- Image reduction ---')
        trap_reduction.clreduce(ds.iraf)
        dark_warning = trap_reduction.check_darks(ds.iraf, ds.tmpout)
        if dark_warning == True:
            input("press enter")
        print('--- generating ephemeris ---')
        if ephemeris == None:
            ephemeris = get_ephem.ephemeris()
        ephemeris.retrieve_param_from_fits(ds.tmpout)
        ephemeris.query_input(unique_target=True)
        ephemeris.generate_ephem_files(ds.tmpout)
        
        print('--- generating calib file ---')
        ZP_warning = trap_reduction.generate_ZP(ds.calib, ephemeris, fitstable)
        if ZP_warning == True:
            input("press enter")
        pixsize = trap_reduction.set_pixsize_in_clafrhocalcext(fitstable)
        print('--- launching afrhocalcext ---')
        trap_reduction.clafrhocalcext(ds.iraf, pixsize[1], str(0), str(0), str(0), str(0)) #launch a first reduction of all the files by default
        trap_plot.plot_centering_profile(ds.tmpout)
        while True:
            centerlist = pd.read_csv(os.path.join(ds.tmpout, 'centerlist'),header=None, sep=' ',usecols=[0,2,3,5,10], names=['file', 'xcent', 'ycent', 'filt', 'ctnmethod'])
            print(centerlist)
            while True:
                solocomete = False
                inp = input('''check individual images for centering.
                            - Relaunch afrhocalcext for all images (r)
                            - Relaunch afrhocalcext for an individual file (IMINDEX XCENTER YCENTER BOXSIZE)
                            - Bypass (b)
                            :''')
                if inp == 'r' or inp == 'R' or inp == 'b' or inp == 'B':
                    break
                elif len(inp.split(' ')) == 4: #check of good format for a one file reduction
                    solocomete = True
                    try:
                        FILE = centerlist.iloc[[inp.split(' ')[0]]].file.values[0]
                    except:
                        solocomete = False
                        print('IMINDEX wrong format')
                    try:
                        XCENTER = int(inp.split(' ')[1])
                    except:
                        solocomete = False
                        print('XCENTER wrong format')
                    try:
                        YCENTER = int(inp.split(' ')[2])
                    except:
                        solocomete = False
                        print('YCENTER wrong format')
                    try:
                        BOXSIZE = int(inp.split(' ')[3])
                    except:
                        solocomete = False
                        print('BOXSIZE wrong format')
                    if solocomete == True: # break if good format
                        break
                else:
                    print('wrong input')
            if solocomete == True:
                trap_reduction.clafrhocalcext(ds.iraf, pixsize[1], FILE, str(XCENTER), str(YCENTER), str(BOXSIZE))
                trap_plot.plot_centering_profile(ds.tmpout, solocomet=True)
                continue
            if inp == 'r' or inp == 'R':
                print('relaunching afrhocalcext')
                trap_reduction.clafrhocalcext(ds.iraf, pixsize[1], str(0), str(0), str(0), str(0))
                trap_plot.plot_centering_profile(ds.tmpout)
                continue
            elif inp == 'b' or inp == 'B':
                print('continuing script')
                break
        trap_reduction.clean_afrhotot(ds.tmpout)
            
        print('--- initiating haserinput ---')
        while True:
            haserinput_warning = trap_reduction.generate_haserinput(ds.tmpout)
            if haserinput_warning == True:
                while True:
                    inp = input('regenerating haserinput (r) or bypass (b)? [r/b]')
                    if inp == 'r' or inp == 'R' or inp == 'b' or inp == 'B':
                        break
                if inp == 'r' or inp == 'R':
                    print('relaunching hasercalctest')
                    continue
                elif inp == 'b' or inp == 'B':
                    print('continuing script')
                    break
            else:
                break
            
        print('--- launching hasercalctest ---')
        while True:
            trap_reduction.clhasercalctest(ds.iraf, 'yes')
            while True:
                inp = input('relaunch hasercalctext (r) or bypass (b)? [r/b]')
                if inp == 'r' or inp == 'R' or inp == 'b' or inp == 'B':
                    break
            if inp == 'r' or inp == 'R':
                print('relaunching hasercalctest')
                continue
            elif inp == 'b' or inp == 'B':
                print('continuing script')
                break
            
            
        # save the files in the reduced folder
        print('--- initiating reduced directory structure ---')
        reduced_dir = os.path.join(ds.reduced, path.split('/')[-1] + pixsize[0])
        images_dir = os.path.join(reduced_dir, "images")
        profiles_dir = os.path.join(reduced_dir, "profiles")
        afrho_dir = os.path.join(reduced_dir, "afrho")
        haser_dir = os.path.join(reduced_dir, "haser")
        garbage_dir = os.path.join(reduced_dir, "probably_garbage")
        centering_dir = os.path.join(reduced_dir, "centering")
        
        if not os.path.exists(reduced_dir):
            os.makedirs(reduced_dir)
            print("created", reduced_dir)
        if not os.path.exists(profiles_dir):
            os.makedirs(profiles_dir)
            print("created", profiles_dir)
        if not os.path.exists(images_dir):
            os.makedirs(images_dir)
            print("created", images_dir)
        if not os.path.exists(afrho_dir):
            os.makedirs(afrho_dir)
            print("created", afrho_dir)
        if not os.path.exists(haser_dir):
            os.makedirs(haser_dir)
            print("created", haser_dir)
        if not os.path.exists(garbage_dir):
            os.makedirs(garbage_dir)
            print("created", garbage_dir)
        if not os.path.exists(centering_dir):
            os.makedirs(centering_dir)
            print("created", centering_dir)
        
        print('--- copying files in reduced directory ---')
        for path2, subdirs2, files2 in os.walk(ds.tmpout):
            for file in files2:
                if file.endswith('.fits') and 'tmp' not in file and 'coms' not in file:
                    shutil.copy(os.path.join(path2, file), os.path.join(images_dir, file))
                    print('copied', file, "in reduced dir")
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
                    shutil.copy(os.path.join(path2, file), os.path.join(centering_dir, file))
                    print('copied', file, "in reduced dir")
                # elif 'tmp' not in file:
                else:
                    shutil.copy(os.path.join(path2, file), os.path.join(garbage_dir, file))
                    print('copied', file, "in reduced dir")
                    
print('Executed in ', datetime.now() - dt)

# =============================================================================
# trap_plot.plot_radprof(input_dir=ds.reduced, saveplot=os.path.join(ds.reduced,'plots'))
# 
# print("launch haser for TS and TN!!")
# tel = 'TS'
# 
# # if os.path.exists(ds.hasercalc):
# #     shutil.rmtree(ds.hasercalc)
# # os.mkdir(ds.hasercalc)
# # trap_reduction.generate_haserinput_from_reduced((rad_dir = ds.reduced, haser_dir=ds.hasercalc, copy_rad = True, tel = tel)
# 
# trap_reduction.clhasercalctest(ds.iraf, 'yes')
# print("check if there is an error, if needed comment haserimput and relaunch")
# appendres = input('append results?(Y,y)')
# if appendres == 'Y' or appendres =="y":
#     timenow = datetime.now()
#     with open(os.path.join(ds.hasercalc, 'outputhasertestall-BC'),'r') as firstfile, open(os.path.join(ds.reduced, 'outputhasertestall-BC'),'a') as secondfile:
#         secondfile.write("########## haser with " + tel + ' #### ' + str(timenow) + '\n')
#         for line in firstfile:   
#             secondfile.write(line[:-1] + ' ' + tel + '\n')
#     with open(os.path.join(ds.hasercalc, 'Haserimput.testall-BC'),'r') as firstfile, open(os.path.join(ds.reduced, 'Haserimput.history'),'a') as secondfile:
#         secondfile.write("########## haser with " + tel + ' #### ' + str(timenow) + '\n')
#         for line in firstfile:   
#             secondfile.write(line)
# =============================================================================
            

