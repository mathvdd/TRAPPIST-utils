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
import directory_structure
import os
import trap_reduction
import shutil
import get_ephem
import trap_plot

#started 2022-03-21 TS
# 2022 04 28
startdate = datetime.datetime.strptime('2022-01-05' + 'T12:00:00', '%Y-%m-%dT%H:%M:%S') #the night starting
enddate = datetime.datetime.strptime('2022-05-15' + 'T12:00:00', '%Y-%m-%dT%H:%M:%S') #starting that night not included
obs = 'TN'
comets = ['0019P', 'CK19L030'] # list of comets to take into account. set empty to take all 

########################

dt = datetime.datetime.now()
ds = directory_structure.directory_structure()
ephemeris = None

def import_perihelions(file_path):
    file = pd.read_csv(file_path)
    perihelions = {}
    for index, row in file.iterrows():
        perihelions[row['comet']] = pd.Timestamp(year=row['year'], month=row['month'], day=row['day'])
    return perihelions

if obs == 'TS':
    NASfitstable = query_NAS.loadcsvtable("/home/Mathieu/Documents/TRAPPIST/raw_data/TS_query.txt")
elif obs == 'TN':
    NASfitstable = query_NAS.loadcsvtable("/home/Mathieu/Documents/TRAPPIST/raw_data/TN_query.txt")

objects_names = query_NAS.check_objects_names(startdate,enddate,NASfitstable)
while True:
    perihelions = import_perihelions('/home/Mathieu/Documents/TRAPPIST/reduced_data/perihelions')
    inlist = []
    outlist = []
    for obj in objects_names:
        if obj in perihelions:
            if (len(comets) == 0):
                inlist.append(obj)
            elif obj in comets:
                inlist.append(obj)
            else:
                print(obj + ' in perihelion but not in the comets list')
                outlist.append(obj) 
        else:
            outlist.append(obj)
    # print('objects observed in the timeframe: ')
    # print(objects_names)
    print('already in the perihelions file:')
    print(inlist)
    print('out of the perihelions file:')
    print(outlist)
    inp = input('break and download (d) or recheck (r)?')
    if inp == 'd' or inp =="D":
        break

# download the raw files
for comet in inlist:
    print('Downloading ' + comet)
    output_path = os.path.join(ds.raw, comet, obs)
    query_NAS.get_files(comet, NASfitstable, output_path, dayinterval=7, dateinterval=(startdate, enddate))

        
# makes list of folders to reduce
list_to_reduce = []
for comet in inlist:
    output_path = os.path.join(ds.raw, comet, obs)
    for path, subdirs, files in os.walk(output_path):
        try:
            subdir = path.split('/')[-1]
            subdir_date = pd.Timestamp(year=int(subdir[:4]), month=int(subdir[4:6]), day=int(subdir[6:8]), hour=23, minute=59)
            if (subdir_date >= startdate) and (subdir_date < enddate):
                list_to_reduce.append(path)
        except:
            continue

# starting the reduction
for path in list_to_reduce:
    raw_dir = path
    comet = path.split('/')[-3]
    reduced_dir = os.path.join(ds.reduced, comet, raw_dir.split('/')[-1] + obs)
    print("working with", raw_dir)
    
    #check if data already reduced
    while True:
        if os.path.exists(reduced_dir):
            inp = input('reduced data detected in ' + reduced_dir + '. Delete old and reduce (d) or skip night (s)?')
            if inp == 'd' or inp =="D":
                shutil.rmtree(reduced_dir)
                reduce_flag = True
                break
            elif inp == 's' or inp =="S":
                print('skipping dataset')
                reduce_flag = False
                break
            else:
                continue
        else:
            reduce_flag = True
            break
    
    if reduce_flag == True:
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
        # changing the comet name to make it a good input for Horizons
        if (len(comet) == 8) and (comet[0:2] + comet[4] + comet[6] == 'CK00'):
            comet_name = '20' + comet[2:4] + ' ' + comet[5] + comet[7]
        else:
            comet_name = comet
        ephemeris = get_ephem.ephemeris()
        ephemeris.retrieve_param_from_fits(ds.tmpout)
        ephemeris.query_input(unique_target=False, target=comet_name)
        ephemeris.generate_ephem_files(ds.tmpout)
        
        print('--- generating calib file ---')
        ZP_warning = trap_reduction.generate_ZP(ds.calib, ephemeris, fitstable)
        if ZP_warning == True:
            input("press enter")
        pixsize = trap_reduction.set_pixsize_in_clafrhocalcext(fitstable)
        print('--- launching afrhocalcext ---')
        trap_reduction.clafrhocalcext(ds.iraf, pixsize[1], str(0), str(0), str(0), str(0)) #launch a first reduction of all the files by default
        while True:
            trap_plot.plot_centering_profile(ds.tmpout)
            centerlist = pd.read_csv(os.path.join(ds.tmpout, 'centerlist'),header=None, sep=' ',usecols=[0,2,3,5,10], names=['file', 'xcent', 'ycent', 'filt', 'ctnmethod'])
            print(centerlist)
            while True:
                print(comet)
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
                continue
            if inp == 'r' or inp == 'R':
                print('relaunching afrhocalcext')
                trap_reduction.clafrhocalcext(ds.iraf, pixsize[1], str(0), str(0), str(0), str(0))
                continue
            elif inp == 'b' or inp == 'B':
                print('continuing script')
                break
        trap_reduction.clean_afrhotot(ds.tmpout)
            
        print('--- initiating haserinput ---')
        while True:
            haserinput_warning, narrowcontlist = trap_reduction.check_haser_continuum(ds.tmpout)
            if haserinput_warning == True:
                print("no or more than one BC filter")
                print("narrowband continuum filters:")
                print(str(narrowcontlist))
                # with open(os.path.join(ds.tmpout, 'nohaser'),'w') as nohaserfile:
                #     nohaserfile.write("no BC filter\n")
                #     nohaserfile.write("narrowband continuum filters:\n")
                #     nohaserfile.write(str(narrowcontlist))
                # while True:
                #     skip_haser = input('Skip hasercalcext? [y/n]')
                #     if skip_haser in ['y','Y']:
                #         print('skipping hasercalctest')
                #         break
            else:
                trap_reduction.generate_haserinput(ds.tmpout)
                print('--- launching hasercalctest ---')
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
                
print('Executed in ', datetime.datetime.now() - dt)
    
    
