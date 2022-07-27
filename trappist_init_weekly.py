#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  8 15:33:43 2022

@author: Mathieu Vander Donckt
mathieu.vanderdonckt@uliege.be
"""

import query_NAS
import datetime
import pandas as pd
from trapconfig import param
import os
import trap_reduction
import shutil
import get_ephem
import trap_plot
import phase_angle

########################
# INPUT PARAMETERS
startdate = '2020-05-25' #the night starting
enddate = '2022-07-15' #starting that night not included
obs = 'TN'
comets = ['0067P'] # list of comets to take into account. set empty to take all 
skip = True # skip without asking raw data directory donwload if data already in raw_data.
# skip reduction if there is already a set of reduced data
# If set to False, will ask what to do in both cases

########################

#started mid-March? TS
# 2022 04 28
startdate = datetime.datetime.strptime(startdate + 'T12:00:00', '%Y-%m-%dT%H:%M:%S') #the night starting
enddate = datetime.datetime.strptime(enddate + 'T12:00:00', '%Y-%m-%dT%H:%M:%S') #starting that night not included

dt = datetime.datetime.now()

conda = True if param['conda'] == 'True' else False #wether to use 'source activate to launch cl or not'

def import_perihelions(file_path):
    file = pd.read_csv(file_path)
    perihelions = {}
    peri_serie = file['comet']
    for index, row in file.iterrows():
        perihelions[row['comet']] = pd.Timestamp(year=row['year'], month=row['month'], day=row['day'])
    return perihelions, peri_serie

if obs == 'TS':
    NASfitstable = query_NAS.loadcsvtable(param['TS_qNAS'])
elif obs == 'TN':
    NASfitstable = query_NAS.loadcsvtable(param['TN_qNAS'])

objects_names = query_NAS.check_objects_names(startdate,enddate,NASfitstable)
while True:
    perihelions, peri_serie = import_perihelions(param['perihelion'])
    # inobj = pd.merge(objects_names, peri_serie, how='inner', left_on='object', right_on='comet').drop(['comet'],axis=1)
    inperi_obj = objects_names.loc[objects_names['object'].isin(peri_serie)]
    outperi_obj = objects_names.loc[~objects_names['object'].isin(peri_serie)]
    if (len(comets) == 0):
        selected_obj = inperi_obj
        unselected_obj = []
    else:
        selected_obj = inperi_obj.loc[inperi_obj['object'].isin(comets)]
        unselected_obj = inperi_obj.loc[~inperi_obj['object'].isin(comets)]
        notinperi_selected = outperi_obj.loc[outperi_obj['object'].isin(comets)]
    
    inlist = selected_obj['object'].to_list()
    
    print('\nUnselected objects outside the perihelion file:')
    print('-------------------------------------------------')
    print(outperi_obj)
    print('\nUnselected objects in the perihelion file:')
    print('--------------------------------------------')
    print(unselected_obj)
    print('\nSelected objects:')
    print('-------------------')
    print(selected_obj)
    if len(notinperi_selected) > 0:
        print('\nTo add to the perihelion file (unselected):')
        print('-------------------')
        print(notinperi_selected)
    
    inp = input('\nbreak and download (d) or recheck (r)?')
    if inp == 'd' or inp =="D":
        break

# download the raw files
for comet in inlist:
    print('Downloading ' + comet)
    output_path = os.path.join(param['raw'], comet, obs)
    query_NAS.get_files(comet, NASfitstable, output_path, dayinterval=7, dateinterval=(startdate, enddate), skip_existing=skip)
        
# makes list of folders to reduce
list_to_reduce = []
for comet in inlist:
    output_path = os.path.join(param['raw'], comet, obs)
    for path, subdirs, files in os.walk(output_path):
        try:
            subdir = path.split('/')[-1]
            subdir_date = pd.Timestamp(year=int(subdir[:4]), month=int(subdir[4:6]), day=int(subdir[6:8]), hour=23, minute=59)
            if (subdir_date >= startdate) and (subdir_date < enddate):
                list_to_reduce.append(path)
        except:
            continue
list_to_reduce.sort()

# starting the reduction
for path in list_to_reduce:
    raw_dir = path
    comet = path.split('/')[-3]
    night = raw_dir.split('/')[-1]
    reduced_dir = os.path.join(param['reduced'], comet, night + obs)
    print("working with", raw_dir)
    print("reduced dir set as", reduced_dir)
    
    #check if data already reduced
    while True:
        if os.path.exists(reduced_dir) and skip == False:
            inp = input('reduced data detected in ' + reduced_dir + '. Delete old and reduce (d) or skip night (s)?')
            if inp == 'd' or inp =="D":
                # shutil.rmtree(reduced_dir)
                reduce_flag = True
                break
            elif inp == 's' or inp =="S":
                print('skipping dataset')
                reduce_flag = False
                break
            else:
                continue
        elif os.path.exists(reduced_dir) and skip == True:
            print('reduced data detected in ' + reduced_dir + '. skip = True, skipping dataset')
            reduce_flag = False
            break
        else:
            reduce_flag = True
            break
    
    if reduce_flag == True:
        while True: 
            print('--- Renaming .fts to .fits ---')
            trap_reduction.renameftsfits(raw_dir)
            print('--- checking for calib files ---')            
            print('\n----------------------------------------------------------------')
            fitstable = trap_reduction.get_fitstable(raw_dir)
            no_file = False
            if len(fitstable) == 0:
                inp = input('no fits file detected. press (s) for skip this night or enter to continue [s/enter]')
                if inp == 's' or inp == 'S':
                    no_file = True
                    break
                else:
                    no_file = False
            check_calib_warning, lighttable = trap_reduction.check_calib(fitstable)
            if check_calib_warning == True:
                inp = input('''Some calibration files are missing!
                            - Press enter to reload table
                            - Query for calibration files older than a week (give line index)
                            - bypass (b)
                            :''')
                            
                if inp == 'b' or inp == 'B':
                    break
                else:
                    try:
                        ind = int(inp)
                    except:
                        continue
                    if ind >= 0 and ind <= (len(lighttable) -1):
                        imline = lighttable.iloc[ind]
                        if imline['nb_bias'] == 0:
                            query_NAS.lookforcalib(NASfitstable, 'bias', raw_dir[:-8], night)
                        if imline['nb_flat'] == 0:
                            query_NAS.lookforcalib(NASfitstable, 'flat', raw_dir[:-8], night, filt=imline['filt'])
                        if imline['nb_dark'] == 0:
                            query_NAS.lookforcalib(NASfitstable, 'dark', raw_dir[:-8], night, exptime=int(imline['exptime']))
                        # query_NAS.lookforcalib(NASfitstable, imtype, output_fold, night, obj='', exptime=15, filt='R', dayinterval=0, copy=True)
                    else:
                        print('wrong index')
            else:
                break
        if no_file == True:
            continue
        # continue
        print('--- Initating tmpdata and tmpout dirs ---')
        if os.path.exists(param['tmpdata']):
            shutil.rmtree(param['tmpdata'])
        os.mkdir(param['tmpdata'])
        print('--- Renaming files to Trappist format and copying to tmpdata ---')
        trap_reduction.pythrename(raw_dir, param['tmpdata'])
        if os.path.exists(param['tmpout']):
            shutil.rmtree(param['tmpout'])
        os.mkdir(param['tmpout'])
        print('--- Image reduction ---')
        trap_reduction.clreduce(param['iraf'], conda=conda)
        dark_warning = trap_reduction.check_darks(param['iraf'], param['tmpout'])
        if dark_warning == True:
            input("press enter")
        print('--- generating ephemeris ---')
        ephemeris = get_ephem.ephemeris()
        ephemeris.retrieve_param_from_fits(param['tmpout'])
        ephemeris.query_input(unique_target=False, target=comet, convert_MPC_Horizon=True)
        ephemeris.generate_ephem_files(param['tmpout'])
        
        print('--- generating calib file ---')
        ZP_warning = trap_reduction.generate_ZP(param['calib'], ephemeris, fitstable, output_dir=param['tmpout'])
        if ZP_warning == True:
            input("press enter")
        pixsize = trap_reduction.set_pixsize_in_clafrhocalcext(fitstable)
        print('--- launching afrhocalcext ---')
        trap_reduction.clafrhocalcext(param['iraf'], pixsize[1], str(0), str(0), str(0), str(0), conda=conda) #launch a first reduction of all the files by default
        trap_plot.plot_centering_profile(param['tmpout'], comet_name=comet)
        while True:
            centerlist = pd.read_csv(os.path.join(param['tmpout'], 'centerlist'),header=None, sep=' ',usecols=[0,2,3,5,10], names=['file', 'xcent', 'ycent', 'filt', 'ctnmethod'])
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
                        XCENTER = float(inp.split(' ')[1])
                    except:
                        solocomete = False
                        print('XCENTER wrong format')
                    try:
                        YCENTER = float(inp.split(' ')[2])
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
                    if solocomete == True: # break if good format
                        break
                else:
                    print('wrong input')
            if solocomete == True:
                trap_reduction.clafrhocalcext(param['iraf'], pixsize[1], FILE, str(XCENTER), str(YCENTER), str(BOXSIZE), conda=conda)
                trap_plot.plot_centering_profile(param['tmpout'], solocomet=True, comet_name=comet)
                continue
            if inp == 'r' or inp == 'R':
                print('relaunching afrhocalcext')
                trap_reduction.clafrhocalcext(param['iraf'], pixsize[1], str(0), str(0), str(0), str(0), conda=conda)
                trap_plot.plot_centering_profile(param['tmpout'], comet_name=comet)
                continue
            elif inp == 'b' or inp == 'B':
                print('continuing script')
                break
        trap_reduction.clean_afrhotot(param['tmpout'])
        phase_angle.generate_palist_tmpout(comet, obs, param['tmpout'])
        print(len(fitstable.loc[fitstable['filt'].isin(['OH','CN','NH','C3','C2']) & fitstable['type'].isin(['LIGHT', 'Light Frame'])]))
        if len(fitstable.loc[fitstable['filt'].isin(['OH','CN','NH','C3','C2']) & fitstable['type'].isin(['LIGHT', 'Light Frame'])]) > 0:            
       	    print('--- initiating haserinput ---')
       	    while True:
                haserinput_warning = trap_reduction.generate_haserinput(param['tmpout'])
                if haserinput_warning == True:
               	    print('Problem generatin haserinput')
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
                trap_reduction.clhasercalctest(param['iraf'], 'yes', conda=conda)
                trap_plot.plot_haserprofile(param['tmpout'],comet_name=comet)
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
        else:
            print('No NB filters found, skipping Haser')
        # save the files in the reduced folder
        print('--- initiating reduced directory structure ---')
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
        for path2, subdirs2, files2 in os.walk(param['tmpout']):
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
                elif 'palist' in file and 'tmp' not in file:
                    shutil.copy(os.path.join(path2, file), os.path.join(centering_dir, file))
                    print('copied', file, "in reduced dir")
                # elif 'tmp' not in file:
                else:
                    shutil.copy(os.path.join(path2, file), os.path.join(garbage_dir, file))
                    print('copied', file, "in reduced dir")
                
print('Executed in ', datetime.datetime.now() - dt)
    
    
