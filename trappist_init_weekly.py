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
startdate = '2010-12-07' #the night starting
enddate = '2024-12-08' #starting that night not included
obs = 'TS'
comets = ['0067P']
# comets = ['0081P','CK19L030','CK17K020','CK22P010','CK21Y010','CK20K010',
#           '0096P','CK22U020','CK20V020','CK22E030','CK19U050'] #'0364P' list of comets to take into account. set empty to take all 
skip = True # skip without asking raw data directory donwload if data already in raw_data.
# skip reduction if there is already a set of reduced data
# If set to False, will ask what to do in both cases
# Qfitlim = (4.0, 4.5) # limit in log10 km for the range over which Q is fitted
Qfitlim = (3.6, 4.1)
only_BVRI = False
no_CLEAR_Z=True
kitty = False
# fc = {'OH':15, #2022E3
#       'NH':21,
#       'CN':27,
#       'C3':231,
#       'C2':193,
#       'CO+':56,
#       'H2O':129}
fc = {'OH':19,
      'NH':24,
      'CN':30,
      'C3':248,
      'C2':170,
      'CO+':56,
      'H2O':129}
########################

dt = datetime.datetime.now()

startdate = datetime.datetime.strptime(startdate + 'T12:00:00', '%Y-%m-%dT%H:%M:%S') #the night starting
enddate = datetime.datetime.strptime(enddate + 'T12:00:00', '%Y-%m-%dT%H:%M:%S') #starting that night not included

if only_BVRI:
    if only_BVRI:
        print('only_BVRI is on')

conda = True if param['conda'] == 'True' else False #wether to use 'source activate to launch cl or not'

#def import_perihelions(file_path):
#    file = pd.read_csv(file_path)
#    perihelions = {}
#    peri_serie = file['comet']
#    for index, row in file.iterrows():
#        perihelions[row['comet']] = pd.Timestamp(year=row['year'], month=row['month'], day=row['day'])
#    return perihelions, peri_serie

if obs == 'TS':
    NASfitstable = query_NAS.loadcsvtable(param['TS_qNAS'])
    NASfitstable = NASfitstable.loc[NASfitstable['readmode'] == '1MHz 1CH']
elif obs == 'TN':
    NASfitstable = query_NAS.loadcsvtable(param['TN_qNAS'])

objects_names = query_NAS.check_objects_names(startdate,enddate,NASfitstable,only_BVRI = only_BVRI)
while True:
    peri_serie = trap_reduction.import_perihelion(param['perihelion'], update=False)['id']
    # inobj = pd.merge(objects_names, peri_serie, how='inner', left_on='object', right_on='comet').drop(['comet'],axis=1)
    inperi_obj = objects_names.loc[objects_names['object'].isin(peri_serie)]
    outperi_obj = objects_names.loc[~objects_names['object'].isin(peri_serie)]
    if (len(comets) == 0):
        selected_obj = inperi_obj
        unselected_obj = []
        notinperi_selected = []
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
    query_NAS.get_files(comet, NASfitstable, output_path, dayinterval=7, dateinterval=(startdate, enddate),
                        skip_existing=skip, only_BVRI = only_BVRI)
# input('download finished')      
# makes list of folders to reduce
list_to_reduce = []
# list_to_reduce = ['/home/Mathieu/Documents/TRAPPIST/raw_data/0073P/TN/20221221',
#                   '/home/Mathieu/Documents/TRAPPIST/raw_data/0081P/TN/20221221',
#                   '/home/Mathieu/Documents/TRAPPIST/raw_data/CK22A020/TN/20221221']
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
            fitstable = trap_reduction.get_fitstable(raw_dir)
            no_file = False
            if len(fitstable) == 0:
                inp = input('no fits file detected. press (s) for skip this night or enter to continue [s/enter]')
                if inp == 's' or inp == 'S':
                    no_file = True
                    break
                else:
                    no_file = False
            if len(fitstable.loc[(fitstable['type'].isin(['LIGHT', 'Light Frame'])) & ~(fitstable['filt'] == 'Clear')]) == 0:
                inp = input('Only "Clear" filter detected. press (s) for skip this night or enter to continue [s/enter]')
                if inp == 's' or inp == 'S':
                    no_file = True
                    break
                else:
                    no_file = False
            if only_BVRI == True:
                check_calib_warning, lighttable = trap_reduction.check_calib(fitstable,filt_list=['B','V','R','I','Rc','Ic'])
            else:
                check_calib_warning, lighttable = trap_reduction.check_calib(fitstable)
            if check_calib_warning == True:
                inp = input("\nSome calibration files are missing!\
                \n   - Press enter to reload table\
                \n   - Query for closest BC continuum file (BC)\
                \n   - Query for missing 15s darks for flats (df), 10s darks for old TS data (df10) (df5,df30)\
                \n   - Query for calibration files older than a week (give line index)\
                \n   - bypass (b)\
                \n   :")

                if inp == 'b' or inp == 'B':
                    break
                elif inp == 'bc' or inp == 'BC':
                    query_NAS.lookforcalib(NASfitstable, 'light', raw_dir[:-8], night,obj=comet, filt='BC')
                elif inp.startswith('df') or inp.startswith('df'):
                    if len(inp) == 2:
                        query_NAS.lookforcalib(NASfitstable, 'dark', raw_dir[:-8], night, exptime=15)
                    else:
                        print(int(inp[2:]))
                        query_NAS.lookforcalib(NASfitstable, 'dark', raw_dir[:-8], night, exptime=int(inp[2:]))
                #elif inp == 'df' or inp == 'DF':
            #        query_NAS.lookforcalib(NASfitstable, 'dark', raw_dir[:-8], night, exptime=15)
        #        elif inp == 'df10' or inp == 'DF10':
    #                query_NAS.lookforcalib(NASfitstable, 'dark', raw_dir[:-8], night, exptime=10)
#                elif inp == 'df5' or inp == 'DF5':
#                    query_NAS.lookforcalib(NASfitstable, 'dark', raw_dir[:-8], night, exptime=5)
#                elif inp == 'df30' or inp == 'DF30':
#                    query_NAS.lookforcalib(NASfitstable, 'dark', raw_dir[:-8], night, exptime=30)
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
        trap_reduction.pythrename(raw_dir, param['tmpdata'],only_BVRI=only_BVRI, no_CLEAR_Z=no_CLEAR_Z)
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
        # ZP_warning = trap_reduction.generate_ZP(param['calib'], ephemeris, fitstable, output_dir=param['tmpout'])
        filtlist = fitstable.loc[fitstable['type'].isin(['LIGHT', 'Light Frame']), 'filt'].drop_duplicates().values.tolist()
        if no_CLEAR_Z and 'Clear' in filtlist:
            filtlist.remove('Clear')
        if no_CLEAR_Z and 'z' in filtlist:
            filtlist.remove('z')
        ZP_warning = trap_reduction.generate_ZP_new_format(param['calib'], obs, night, filtlist, output_dir=param['tmpout'])
        if ZP_warning == True:
            input("press enter")
        pixsize = trap_reduction.set_pixsize_in_clafrhocalcext(fitstable)
        print('--- launching afrhocalcext ---')
        trap_reduction.clafrhocalcext(param['iraf'], pixsize[1], str(0), str(0), str(0), str(0), conda=conda) #launch a first reduction of all the files by default
        trap_plot.plot_centering_profile(param['tmpout'], comet_name=comet,kitty=kitty)
        trap_reduction.generate_center_comment(param['tmpout'])
        if kitty == True:
            try:
                os.system(f'for f in {param["tmpout"]}/*_centering.png ; do kitty +kitten icat "$f" ; done')
            except:
                print('kitty command failed')
        def add_comment(row, comment):
            row_comment = comment.loc[comment['file'] == row['file'], 'comment']
            return row_comment.to_string(index = False)
        
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
                trap_reduction.clafrhocalcext(param['iraf'], pixsize[1], FILE, str(XCENTER), str(YCENTER), str(BOXSIZE), conda=conda)
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
        # print(len(fitstable.loc[fitstable['filt'].isin(['OH','CN','NH','C3','C2','CO+','H2O']) & fitstable['type'].isin(['LIGHT', 'Light Frame'])]))
        if len(fitstable.loc[fitstable['filt'].isin(['CO+','H2O']) & fitstable['type'].isin(['LIGHT', 'Light Frame'])]) > 0:
            print(fitstable)
            input('CO+ or H2O filter detected')
                             
        if len(fitstable.loc[fitstable['filt'].isin(['OH','CN','NH','C3','C2','CO+','H2O']) & fitstable['type'].isin(['LIGHT', 'Light Frame'])]) > 0:            
       	    print('--- initiating haserinput ---')
            import gfactor
            gfactor.generate_tmptxt(ephemeris.v, ephemeris.rh)
            
       	    while True:
                haserinput_warning = trap_reduction.generate_haserinput(param['tmpout'], fc)
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
            trap_reduction.clhasercalctest(param['iraf'], arg='yes', Qproflow=Qfitlim[0], Qprofhigh=Qfitlim[1], conda=conda)
            trap_plot.plot_haserprofile(param['tmpout'],comet_name=comet,kitty=kitty)
            if kitty == True:
                try:
                    os.system(f'for f in {param["tmpout"]}/*_haserprofile.png ; do kitty +kitten icat "$f" ; done')
                except:
                    print('kitty command failed')
            # while True:
            #     trap_reduction.clhasercalctest(param['iraf'], arg='yes', Qproflow=Qfitlim[0], Qprofhigh=Qfitlim[1], conda=conda)
            #     trap_plot.plot_haserprofile(param['tmpout'],comet_name=comet)
            #     while True:
            #         inp = input('relaunch hasercalctext (r) or bypass (b)? [r/b]')
            #         if inp == 'r' or inp == 'R' or inp == 'b' or inp == 'B':
            #             break
            #     if inp == 'r' or inp == 'R':
            #         print('relaunching hasercalctest')
            #         continue
            #     elif inp == 'b' or inp == 'B':
            #         break
        else:
            print('No NB filters found, skipping Haser')
        phase_angle.generate_palist_tmpout(comet, obs, param['tmpout'])
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
    
    
