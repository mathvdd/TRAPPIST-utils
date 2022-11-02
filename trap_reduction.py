#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 25 15:51:19 2021

@author: Mathieu Vander Donckt
mathieu.vanderdonckt@uliege.be
"""

import os
import csv
import pandas as pd
from astropy.io import fits
import shutil
import datetime
from trapconfig import param
# Add a rename fts to fits for TS images
#Add check if all flats are there

fc = {'OH':5,
      'NH':20,
      'CN':25,
      'C3':190,
      'C2':170,
      'CO+':56,
      'H2O':129}

#see calib08140914.dat
ZP = {'OH': [3090, 10.56e-9,   1.791,  1.698e-2,  0.98,   1,  1.60],
      'NH': [3362,  8.42e-9,   1.188,  1.907e-2,  0.99,   2,  0.65],
      'UC': [3448,  7.802e-9,  1.101,  1.,        1.,     3,  0.59],
      'CN': [3870,  8.6e-9,    1.031,  1.812e-2,  0.99,   4,  0.36],
      'C3': [4062,  8.16e-9,   0.497,  3.352e-3,  0.19,   5,  0.29],
      'BC': [4450,  6.21e-9,   0.0,    1.,        1.,     7,  0.25],
      'C2': [5141,  3.887e-9, -0.423,  5.433e-3,  0.66,   8,  0.15],
      'GC': [5260,  3.616e-9, -0.507,  1.,        1.,     9,  0.14],
      'RC': [7128,  1.316e-9, -1.276,  1.,        1.,    11,  0.05],
      'B': [4440,  6.4e-9,    0.0,    1.,        1,     12,  0.25],
      'V': [5483,  3.67e-9,  -0.649,  1.,        1,     13,  0.14],
      'R': [6855,  1.92e-9,  -1.019,  1.,        1,     14,  0.098],
      'I': [8637,  9.39e-10, -1.375,  1.,        1,     15,  0.043],
      'CO+': [4266, 7.323e-9, 0.338, 1.549e-2, 0.99, 6, 0.25],
      'H2O' : [7020, 1.38e-9, -1.249, 5.424e-3, 1., 10, 0.07]
      }

filt_list = ['OH','CN','C2','C3','NH','CO+','H2O','UC','BC','RC','GC','R','I', 'B', 'V'] #to be coherent with the subsets file for progtrap2.cl

def renameftsfits(raw_path):
    """
    Rename 'fts' into 'fits'

    Parameters:
        raw_path (str): path to the folder containing the fits files (subfolders will be affected as well)
    """
    print("rename fts into fits")
    count = 0
    for path, subdirs, files in os.walk(raw_path):
        for file in files:
            if os.path.join(path, file).endswith(".fts"):
                os.rename(os.path.join(path, file), os.path.join(path, file)[:-3] + 'fits')
                count += 1
    print("renamed", count, "fts files")

def pythrename(raw_path, tmpdata_dir,only_BVRI=False):
    """
    Rename fits files to the trappist format and copy them to the tmpdata folder

    Parameters:
        raw_path (str): path to the raw file folder
        tmpdata_dir (str): path to the directory where the files are to be copied
    """
    count = 0
    for path, subdirs, files in os.walk(raw_path):
        for file in files:
            if os.path.join(path, file).endswith(".fits"):
                print(os.path.join(path, file))
                with fits.open(os.path.join(path, file)) as hdul:
                    obsdate = hdul[0].header['DATE-OBS']
                    if only_BVRI == True and hdul[0].header['IMAGETYP'] in ['LIGHT', 'Light Frame']:
                        obsfilt = hdul[0].header['FILTER']
                        if obsfilt in ['B','V','R','I','Rc','Ic']:
                            shutil.copy(os.path.join(path, file), os.path.join(tmpdata_dir, "TRAP." + obsdate[:19] + ".fits"))
                    else:
                        shutil.copy(os.path.join(path, file), os.path.join(tmpdata_dir, "TRAP." + obsdate[:19] + ".fits"))
                count += 1
            elif os.path.join(path, file).endswith(".fts"):
                print("WARNING: convert fts to fits")
    print("renamed and copied", count, "files in", tmpdata_dir)

def clrename(iraf_dir, raw_path, conda=True):
    print('rename and set the data in tmpdata')
    with open(os.path.join(iraf_dir, 'wrapper_rename.cl'), 'w') as f:
        f.write('\n'.join(['renamefits ' + raw_path + ' yes',
                      'logout']))
    if conda == True:
        os.system('\n'.join(['cd ' + iraf_dir,
                          'source activate iraf27', #conda activate not working
                          'cl < wrapper_rename.cl']))
    else:
        os.system('\n'.join(['cd ' + iraf_dir,
                          'cl < wrapper_rename.cl']))

def check_calib(fitstable, filt_list=filt_list):
    """
    Check whether all necessary calibration are available in the directory before recduction
    Return False if there is a warning flag

    Parameters:
        fitstable: table with fits file parameters (see get_fitstable())
        filt_list, optional: list of filters to take into consideration. Default is ['OH','CN','C2','C3','NH','UC','BC','RC','GC','R','I', 'B', 'V']
    """
    pd.set_option('display.max_rows', 500)
    pd.set_option('display.max_columns', 500)
    pd.set_option('display.width', 1000)
    # print("filter list set to", filt_list)

    # light_filters = fitstable.loc[fitstable['type'].isin(['LIGHT', 'Light Frame']), 'filter'].drop_duplicates().values.tolist()
    lighttable = fitstable.loc[(fitstable['type'].isin(['LIGHT', 'Light Frame'])) & (fitstable['filt'].isin(filt_list)), ['file', 'filt', 'exptime']]
    lighttable['nb_flat'] = None
    lighttable['nb_dark'] = None
    lighttable['nb_bias'] = (fitstable['type'].isin(['BIAS', 'Bias Frame'])).sum()

    for index, row in lighttable.iterrows() :
        lighttable.at[index, 'nb_flat'] = ((fitstable['type'].isin(['FLAT', 'Flat Frame'])) & (fitstable.filt == row['filt'])).sum()
        lighttable.at[index, 'nb_dark'] = ((fitstable['type'].isin(['DARK', 'Dark Frame'])) & (fitstable.exptime == row['exptime'])).sum()
    print('\n-------------------------------------------------')
    print('CHECK CALIB')
    print('-------------------------------------------------\n')
    print('NOTE: If no right exposure time dark is found, a linear extrapolation from another master dark will be used \
          \nThe extrapolated darks still need to be in the data folder\nCurrent configuration:')
    print(pd.read_csv(os.path.join(param['calib'], 'dark_substitution'), names=['exptime','scaled from'], sep=' ').transpose().to_string(header=False))
    print('\n')

    print(lighttable)

    nb_flat_dark15 = ((fitstable['type'].isin(['DARK', 'Dark Frame'])) & (fitstable.exptime == 15)).sum()
    nb_flat_dark10 = ((fitstable['type'].isin(['DARK', 'Dark Frame'])) & (fitstable.exptime == 10)).sum()
    if (nb_flat_dark15 == 0 ) and (nb_flat_dark10 > 0 ):
        print(f"\nNumber of dark frames for flat correction (exptime = 10s): {nb_flat_dark10}\n")
    else:
        print(f"\nNumber of dark frames for flat correction (exptime = 15s): {nb_flat_dark15}\n")

    warning_flag = False
    for index, row in lighttable.iterrows() :
        if row['nb_flat'] == 0:
            print("WARNING: no flat for", row['file'])
            warning_flag = True
        elif row['nb_flat'] < 5:
            print("WARNING: less than 5 flats for", row['file'])
            warning_flag = True
        if row['nb_dark'] == 0:
            print("WARNING: no dark (",row['exptime'],") for", row['file'])
            warning_flag = True
        elif row['nb_dark'] < 5:
            print("WARNING: less than 5 darks (",row['exptime'],") for", row['file'])
            warning_flag = True
        if row['nb_bias'] == 0:
            print("WARNING: no bias for", row['file'])
            warning_flag = True
        elif row['nb_bias'] < 5:
            print("WARNING: less than 5 bias for", row['file'])
            warning_flag = True
    if (nb_flat_dark15 == 0) and (nb_flat_dark10 == 0):
        print("WARNING: no darks (15s or 10s) for flats")
        warning_flag = True
    elif (nb_flat_dark15 != 0) and (nb_flat_dark15 < 5):
        print("WARNING: less than 5 darks (15s) for flats")
        warning_flag = True
    elif (nb_flat_dark10 != 0) and (nb_flat_dark10 < 5):
        print("WARNING: less than 5 darks (10s) for flats")
        warning_flag = True


    # check if there is a BC filter for narrow bands
    if lighttable['filt'].isin(['OH', 'C2', 'C3', 'CN', 'NH','CO+','H2O']).any() and (lighttable['filt'].isin(['BC']).any() == False):
        print("WARNING: Narrow band images while no continuum BC image")
        warning_flag = True


    return warning_flag, lighttable #can be used to stop the script only if there is a warning

def get_fitstable(raw_dir):
    """
    Scan a directory to get a table with relevant infos from the fits headers
    Return the table

    Parameters:
        raw_dir (str): path to the directory to scan
    """
    fitstable = pd.DataFrame(columns=('file','type', 'filt', 'exptime', 'observatory'))
    dir_list = [raw_dir,
                os.path.join(raw_dir, "Calibration"),
                os.path.join(raw_dir, "AutoFlat"),
                os.path.join(raw_dir, "Flat"),
                os.path.join(raw_dir, "Autoflat")]
    for directory in dir_list:
        if os.path.exists(directory):
            for file in os.listdir(directory):
                if 'tmp' not in file and file.endswith(".fits") or file.endswith(".fts"):
                    with fits.open(os.path.join(directory, file)) as hdul:
                        imfile = os.path.join(directory, file)[len(raw_dir)+1:]
                        imtype = hdul[0].header['IMAGETYP']
                        try:
                            imfilter = hdul[0].header['FILTER']
                        except:
                            imfilter = None
                        imexptime = hdul[0].header['EXPTIME']
                        try:
                            imobservatory = hdul[0].header['OBSERVAT']
                        except:
                            imobservatory = None
                        fitstable.loc[fitstable.shape[0]] = [imfile, imtype, imfilter, imexptime, imobservatory]
                        
    return fitstable

def generate_ZP(calib_dir, ephemeris, fitstable, ZPparams=ZP, output_dir=None):
    """
    Generate calib.dat file with closest zero points

    Parameters:
        calib_dir (str): path to the directory containing zero point files. calib.dat will be copied there
        ephemeris: ephemeris object (see get_ephem.py)
        fitstable: table containing the fits files info (see get_fitstable())
        ZPparams, optional: table containing constant values in the calib.dat file. default should be fine at all time
        output_dir (str), optional, default=None: if not given, output in calib dir
    """
    import datetime
    import jdcal
    warning_flag = False
    if ephemeris.observatory == 'Oukaimeden':
        ZPtable = pd.read_csv(os.path.join(calib_dir, "median_ZPC_North.log"), header=None, sep="\s+")
        observ = 'TN'
    elif ephemeris.observatory == 'TRAPPIST':
        ZPtable = pd.read_csv(os.path.join(calib_dir, "median_ZPC_South.log"), header=None, delim_whitespace=True)
        observ = 'TS'
    else:
        print('WARNING: observatory not well defined')
        warning_flag = True
    ZPtable = ZPtable.drop([1,2,4,5,8,9,10,12,14,16], axis=1)
    ZPtable.columns = ['filt','lambda','JDstart','JDend','ZPmed','ZPav','ZPsig','n']

    # get the closest date
    startdate = datetime.datetime.strptime(ephemeris.parameters['START_TIME'], '%Y-%m-%d %H:%M')
    # enddate = datetime.datetime.strptime(eph.parameters['STOP_TIME'], '%Y-%m-%d %H:%M')
    startdateJD = sum(jdcal.gcal2jd(startdate.year, startdate.month, startdate.day))
    if startdateJD > 2459950:
        print("WARNING: need tu update the code on February 24, 2023")
        warning_flag = True
    datediff = 10000000000
    for index, row in ZPtable.iterrows():
        midpointJD = 2450000 + (row['JDstart'] + (row['JDend'] - row['JDstart'])/2)
        ddiff = startdateJD - midpointJD
        if abs(datediff) > abs(ddiff):
            datediff = ddiff
            periodJD = (row['JDstart'], row['JDend'])
    if datediff > 50:
        print("WARNING: update the ZP file, closest calibration midpoint at", datediff, "days")
        warning_flag = True
    ZPtable_closest = ZPtable.loc[(ZPtable['JDstart'] == periodJD[0]) & (ZPtable['JDend'] == periodJD[1])]

    # get the list of matching filters and print the calib file
    filtlist = fitstable.loc[fitstable['type'].isin(['LIGHT', 'Light Frame']), 'filt'].drop_duplicates().values.tolist()
    # print('filtlist', fitstable)
    if output_dir == None:
        output_path = os.path.join(calib_dir, "calib.dat")
    else:
        output_path = os.path.join(output_dir, "calib.dat")
    with open(output_path, 'w') as f:
        # the Rc filter is refered as R in the fits headers and the calib.dat file # ! only for TS data after 2012-10-17
        ZPtable_closest['filt'] = ZPtable_closest['filt'].replace(['Rc'],'R')
        ZPtable_closest['filt'] = ZPtable_closest['filt'].replace(['Ic'],'I')
        ZPtable_closest['filt'] = ZPtable_closest['filt'].replace(['CO'],'CO+')
        # ZPtable_closest.loc[ZPtable_closest['filt'] == "Rc", 'filt']  = 'R'
        # ZPtable_closest['filt'] = ZPtable_closest['filt'].replace({'Rc':'R'})
        for filt in filtlist:
            if filt in ZPtable_closest['filt'].values.tolist():
                ZPval = ZPtable_closest.loc[ZPtable_closest['filt'] == filt, 'ZPmed']
                ZPsig = ZPtable_closest.loc[ZPtable_closest['filt'] == filt, 'ZPsig']
                line = filt + ' ' + ' '.join([str(int) for int in ZPparams.get(filt)]) + ' ' + str(ZPval.values)[1:-1] + ' ' + str(ZPsig.values)[1:-1]
                f.write(line + '\n')
            else:
                print("WARNING:", filt, "filter not found in the ZP table for the closest date", periodJD)
                warning_flag = True
        f.write('\n# Generated by trap_reduction.generate_ZP() for ' + observ + ' for the ' + startdate.strftime('%m/%d/%Y') + " JD " + str(startdateJD) + '.\n# see calib08140914.dat\n')

    return warning_flag

# import directory_structure
# import get_ephem
# ds = directory_structure.directory_structure()
# fitslist = get_fitstable(os.path.join(ds.raw, "20210403"))
# eph = get_ephem.ephemeris()
# eph.retrieve_param_from_fits(ds.tmpout)
# generate_ZP(ds.calib, eph, fitslist)


def clreduce(iraf_dir, conda=True):
    """
    Wrapper to launch iraf progtrap3 in python

    Parameters:
        iraf_dir (str): path to the home iraf directory
        conda (boolean, default=True): if False, will not use the source activate line when launching cl
    """
    print('reduceand calibrate images from the tmpdata into the tmpout folder')
    with open(os.path.join(iraf_dir, 'wrapper_reduce.cl'), 'w') as f:
        f.write('\n'.join(['progtrap3',
                      'logout']))
    if conda == True:
        os.system('\n'.join(['cd ' + iraf_dir,
                          'source activate iraf27', #conda activate not working
                          'cl < wrapper_reduce.cl']))
    else:
        os.system('\n'.join(['cd ' + iraf_dir,
                          'cl < wrapper_reduce.cl']))

def set_pixsize_in_clafrhocalcext(fitstable):
    obs = fitstable.loc[fitstable['type'].isin(['LIGHT','Light Frame']), 'observatory'][0]
    if obs == 'Oukaimeden':
        pixsize = ('TN', '0.60')
        print('Observatory set to TRAPPIST-North, Oukaimeden (observatory). pixsize = 0.60')
    elif obs == 'TRAPPIST':
        pixsize = ('TS', '0.65')
        print('Observatory set to La Silla--TRAPPIST (observatory). pixsize = 0.65')
    else:
        print('problem defining observatory and pixsize')
    return pixsize

def clafrhocalcext(iraf_dir, pixsize, solocomete, soloinitx, soloinity, soloinitcboxsize, conda=True):
    """
    Wrapper to launch iraf afrhocalcext in python

    Parameters:
        iraf_dir (str): path to the home iraf directory
        pixsize (float): size of the image pixel
        solocomete ()
        soloinitx
        soloinity
        soloinitcboxsize
        conda (boolean, default=True): if False, will not use the source activate line when launching cl
    """
    print('launch afrhocalext.cl')
    with open(os.path.join(iraf_dir, 'wrapper_afrhocalcext.cl'), 'w') as f:
        f.write('\n'.join(['afrhocalcext ' + pixsize + ' ' + solocomete + ' ' + soloinitx + ' ' + soloinity + ' ' + soloinitcboxsize,
                      'logout']))
    if conda == True:
        os.system('\n'.join(['cd ' + iraf_dir,
                          'source activate iraf27', #conda activate not working
                          'cl < wrapper_afrhocalcext.cl']))
    else:
        os.system('\n'.join(['cd ' + iraf_dir,
                          'cl < wrapper_afrhocalcext.cl']))

def check_darks(iraf_dir,tmpout_dir):
    """
    Checks if all darks are made for every line of list_D_exptime
    return False if not

    Parameters:
        iraf_dir (str): path to the home iraf directory
        tmpout_dir (str): path to the directory to check
    """
    warning_flag = False
    with open(os.path.join(iraf_dir, 'list_D_exptime'), newline='') as f:
        reader = csv.reader(f)
        exptimes = []
        for row in reader:
          exptimes.append(row[0])
        print('checking if darks exists for exptimes', exptimes)
        for exptime in exptimes:
            dark_file = os.path.join(tmpout_dir, 'Dark_' + exptime + '.fits')
            if os.path.exists(dark_file) == False:
                print("WARNING: " + dark_file + " not found")
                warning_flag = True
    return warning_flag



# Old version.
# Will generate them all in one file instead
# def generate_haserinput(folder_path, fc=fc, fz=0):
#     #add an option if there is no BC filter
#     print('generate haser input. Check fc and fz')
#     filelist = os.listdir(folder_path)
#     filelist_filter = []
#     for file in filelist:
#             if file[0:4] == 'rad_':
#                 with open(os.path.join(folder_path, file), newline='') as f:
#                     reader = csv.reader(f)
#                     filter_type = next(reader)[0].split()[14]
#                     filelist_filter.append((file, filter_type))
#     #check if it will run with the current folder configuration. May need to change this function later for bathching
#     i = 0
#     for item in filelist_filter:
#         if item[1] == 'BC':
#             BC_file = item
#             i +=1
#     if i != 1:
#         input(i, "BC file(s) found. Check BC files date and rerun/modify the script")
#     else:
#         print("BC file found")
#     for filt in fc:
#         output_list = []
#         for item in filelist_filter:
#             if filt == item[1]:
#                 output_list.append(item[0][4:] + ' ' + BC_file[0] + ' ' + str(fc.get(filt)) + ' ' + str(fz))
#         if len(output_list) > 0:
#             file_path = os.path.join(folder_path, 'Haserimput.test' + filt + '-BC')
#             print("creating " + file_path)
#             with open(file_path, 'w') as f:
#                 f.write('\n'.join(output_list))
#                 f.write('\n') #need a blank end line for hasercalcext to work

def generate_haserinput_from_reduced(rad_dir, haser_dir, copy_rad = False, fc=fc, fz=0, tel =''):
    '''obselete. was made to be run after files from other recipes were ordered in the reduced file
    see generate_haserinput for calculations in the tmpout directory'''
    print('generate haser input. Check fc and fz. will just do the filters documented in fc')
    #creates a table with all the profiles in the rad_dir
    radtable = pd.DataFrame(columns=('path','file', 'filt'))
    for path, subdirs, files in os.walk(rad_dir):
        if tel == '':
            input('need to specify telescope')
        elif tel in path: #need to distinguish between TN and TS for ephemerids
            for file in files:
                if 'rad_' in file and ".png" not in file:
                    with open(os.path.join(path, file), newline='') as f:
                        reader = csv.reader(f)
                        filter_type = next(reader)[0].split()[14]
                        radtable.loc[radtable.shape[0]] = [path, file, filter_type]
                #creates ephem.brol for the haser recipe
                if file == "ephem.brol":
                    with open(os.path.join(path, file),'r') as firstfile, open(os.path.join(haser_dir, 'ephem.brol'),'a') as secondfile:
                        for line in firstfile:
                            secondfile.write(line)
    #looks if the filer is in the filter list, search for BC filter and write the line in the haserinput list
    output_list = []
    for index, row in radtable.iterrows():
        if row['filt'] in fc:
            BCtable = radtable.loc[(radtable['path'] == row['path']) & (radtable['filt'] == 'BC')]
            if len(BCtable.index) != 1:
                print(BCtable)
                print(path)
                print("table length:", len(BCtable.index))
                input('less or more than one BC radial profile in the folder')
            output_list.append(row['file'][4:]
                               + ' '
                               + BCtable.iloc[0]['file']
                               + ' ' +
                               str(fc.get(row['filt']))
                               + ' ' +
                               str(fz))
            if copy_rad == True:
                shutil.copy(os.path.join(row['path'], row['file']), os.path.join(haser_dir, row['file']))
                shutil.copy(os.path.join(row['path'], 'radplus_' + row['file'][4:]), os.path.join(haser_dir, 'radplus_' + row['file'][4:]))
                shutil.copy(os.path.join(row['path'], 'radeplus_' + row['file'][4:]), os.path.join(haser_dir, 'radeplus_' + row['file'][4:]))
                shutil.copy(os.path.join(row['path'], 'radmoins_' + row['file'][4:]), os.path.join(haser_dir, 'radmoins_' + row['file'][4:]))
                shutil.copy(os.path.join(row['path'], 'rademoins_' + row['file'][4:]), os.path.join(haser_dir, 'rademoins_' + row['file'][4:]))
                shutil.copy(os.path.join(BCtable.iloc[0]['path'], BCtable.iloc[0]['file']), os.path.join(haser_dir, BCtable.iloc[0]['file']))
    # print hasercalcinput
    if len(output_list) > 0:
        file_path = os.path.join(haser_dir, 'Haserimput.testall-BC')
        print("creating " + file_path)
        with open(file_path, 'w') as f:
            f.write('\n'.join(output_list))
            f.write('\n') #need a blank end line for hasercalcext to work

def check_haser_continuum(tmpout):
    """
    Obscelete
    Check if there is a BC image for duct continuum correction

    Parameters:
        tmpout (str): directory to check
    """
    # check if there is a continuum image to correct for the dust continuum with the haser script
    # for now only BC is considered but also look for the others
    warning_flag = False
    fitslist = get_fitstable(tmpout)
    # narrowcontlist = fitslist.loc[(fitslist.type == 'Light Frame') & fitslist.filt.isin(['BC', 'RC', 'GC', 'UC'])]
    narrowcontlist = fitslist.loc[fitslist.type.isin(['Light Frame', 'LIGHT'])
                                  & fitslist.filt.isin(['BC', 'RC', 'GC', 'UC'])
                                  & fitslist.file.str.contains('TRAP')]
    nb_BC = len(narrowcontlist.loc[fitslist.filt == 'BC'])
    if nb_BC == 1:
        print('BC image found')
    elif nb_BC > 1:
        print("WARNING: more than one BC image found (" + str(nb_BC) + ')')
        warning_flag = True
    elif nb_BC == 0:
        print("WARNING: no BC images found for dust continuum correction")
        print('Narrow continuum images found:')
        print(narrowcontlist)
        warning_flag = True
    return warning_flag, narrowcontlist

def generate_haserinput(tmpout, fc=fc, fz=0):
    """
    generate the haserinput file for the iraf script. Select the BC image. Gives a warning on error

    Parameters:
        tmpout (str): working directory
        fc (dic, optional): dictionnaru containing the fc for each filter. See file for default
        fz (int, optional): value of fz. default =0
    """
    warning_flag = False

    output_list = []
    fitstable = get_fitstable(tmpout)
    filelist = fitstable.loc[fitstable['type'].isin(['LIGHT', 'Light Frame']) &
                             fitstable.file.str.contains('TRAP')]
    BCtable = filelist.loc[(filelist['filt'] == 'BC')]
    BCtable.reset_index(drop=True, inplace=True)

    #behaviour depending on the number of BC images in the folder
    if len(BCtable) == 1:
        print(BCtable.iloc[0])
        print('selecting the only BC file: ', BCtable.iloc[[0]]['file'])
        cont_file = BCtable.iloc[0]['file']
    elif len(BCtable) == 0:
        print("No BC images found.")
        NBcont = filelist.loc[filelist.filt.isin(['BC', 'RC', 'GC', 'UC'])]
        NBcont.reset_index(drop=True, inplace=True)
        if len(NBcont) > 0:
            print("other NB continuum images available:")
            print(NBcont)
            while True:
                inp = input('enter the index of the NB continuum image to use for dust subtraction, \
                            or b to bypass: ')
                try:
                    inp = int(inp)
                except:
                    inp = inp
                if (inp == 'b') or (inp == 'B'):
                    warning_flag = True
                    return warning_flag
                elif isinstance(inp, int) and (inp <= len(NBcont)-1) and (inp >=0):
                    cont_file = NBcont.iloc[inp]['file']
                    break
                else:
                    print('wrong input')
                    continue
        else:
            print('No other NB images found')
            warning_flag = True
            return warning_flag
    elif len(BCtable) > 1:
        print('More than one BC image found')
        print(BCtable)
        while True:
            inp = input("Input 'c' to use the closest BC image in time for each image.\
                        \nAlternatively, select the index of the BC image to use: ")
            try:
                inp = int(inp)
            except:
                inp = inp
            if (inp == 'c') or (inp == 'C'):
                cont_file = "BC_CLOSEST_DATE"
                break
            elif isinstance(inp, int) and (inp <= len(BCtable)-1) and (inp >=0):
                cont_file = BCtable.iloc[inp]['file']
                break
            else:
                print('wrong input')
                continue

    for index, row in filelist.iterrows():
        if row['filt'] in fc:
            if cont_file != "BC_CLOSEST_DATE":
                output_list.append(row['file'] + '.txt'
                    + ' '
                    + 'rad_' + cont_file + '.txt'
                    + ' ' +
                    str(fc.get(row['filt']))
                    + ' ' +
                    str(fz))
            else:
                delta = datetime.timedelta(days=9999)
                date_NB = datetime.datetime.strptime(row['file'], 'TRAP.%Y-%m-%dT%H:%M:%S.fits')
                for indexb, rowb in BCtable.iterrows():
                    date_BC = datetime.datetime.strptime(rowb['file'], 'TRAP.%Y-%m-%dT%H:%M:%S.fits')
                    if abs(date_NB - date_BC) < delta:
                        delta = abs(date_NB - date_BC)
                        BC_file = rowb['file']
                output_list.append(row['file'] + '.txt'
                    + ' '
                    + 'rad_' + BC_file + '.txt'
                    + ' ' +
                    str(fc.get(row['filt']))
                    + ' ' +
                    str(fz))

    if len(output_list) > 0:
        file_path = os.path.join(tmpout, 'inputhaser-BC')
        print("creating " + file_path)
        with open(file_path, 'w') as f:
            f.write('\n'.join(output_list))
            f.write('\n') #need a blank end line for hasercalcext to work
    else:
        print('No content to put in inputhaser file')
        warning_flag = True
        return warning_flag



# import directory_structure
# ds = directory_structure.directory_structure()
# check_haser_continuum(ds.tmpout)
# generate_haserinput(ds.tmpout)


    #             filelist_filter = []
    #             with open(os.path.join(path, file), newline='') as f:
    #                 reader = csv.reader(f)
    #                 filter_type = next(reader)[0].split()[14]
    #                 filelist_filter.append((file, filter_type))
    # #check if it will run with the current folder configuration. May need to change this function later for bathching
    # i = 0
    # for item in filelist_filter:
    #     if item[1] == 'BC':
    #         BC_file = item
    #         i +=1
    # if i != 1:
    #     input(i, "BC file(s) found. Check BC files date and rerun/modify the script")
    # else:
    #     print("BC file found")
    # output_list = []
    # inc_filt = []
    # for item in filelist_filter:
    #     if item[1] in fc:
    #         output_list.append(item[0][4:] + ' ' + BC_file[0] + ' ' + str(fc.get(item[1])) + ' ' + str(fz))
    #         if item[1] not in inc_filt:
    #             inc_filt.append(item[1])
    # if len(output_list) > 0:
    #     file_path = os.path.join(folder_path, 'Haserimput.testall-BC')
    #     print("creating " + file_path)
    #     print("include filters:", inc_filt)
    #     with open(file_path, 'w') as f:
    #         f.write('\n'.join(output_list))
    #         f.write('\n') #need a blank end line for hasercalcext to work
    # else:
    #     print("no radial profile with associated filter found. Please check fc or files. No haserinput written")


# =============================================================================
# def generate_haserinput(folder_path, fc=fc, fz=0):
#     #add an option if there is no BC filter
#     print('generate haser input. Check fc and fz. will just do the filters documented in fc')
#     filelist = os.listdir(folder_path)
#     filelist_filter = []
#     for file in filelist:
#             if file[0:4] == 'rad_':
#                 with open(os.path.join(folder_path, file), newline='') as f:
#                     reader = csv.reader(f)
#                     filter_type = next(reader)[0].split()[14]
#                     filelist_filter.append((file, filter_type))
#     #check if it will run with the current folder configuration. May need to change this function later for bathching
#     i = 0
#     for item in filelist_filter:
#         if item[1] == 'BC':
#             BC_file = item
#             i +=1
#     if i != 1:
#         input(i, "BC file(s) found. Check BC files date and rerun/modify the script")
#     else:
#         print("BC file found")
#     output_list = []
#     inc_filt = []
#     for item in filelist_filter:
#         if item[1] in fc:
#             output_list.append(item[0][4:] + ' ' + BC_file[0] + ' ' + str(fc.get(item[1])) + ' ' + str(fz))
#             if item[1] not in inc_filt:
#                 inc_filt.append(item[1])
#     if len(output_list) > 0:
#         file_path = os.path.join(folder_path, 'Haserimput.testall-BC')
#         print("creating " + file_path)
#         print("include filters:", inc_filt)
#         with open(file_path, 'w') as f:
#             f.write('\n'.join(output_list))
#             f.write('\n') #need a blank end line for hasercalcext to work
#     else:
#         print("no radial profile with associated filter found. Please check fc or files. No haserinput written")
#
# =============================================================================
def clhasercalctest(iraf_dir, arg='no', Qproflow=3.6, Qprofhigh=4.1, conda=True):
    """
    Wrapper to launch iraf hasercalctest in python

    Parameters:
        iraf_dir (str): path to the home iraf directory
        conda (boolean, default=True): if False, will not use the source activate line when launching cl
    """
    print('launch hasercalctest.cl')
    with open(os.path.join(iraf_dir, 'wrapper_hasercalctest.cl'), 'w') as f:
        f.write('\n'.join([f'hasercalctest {arg} {Qproflow} {Qprofhigh}',
                      'logout']))
    if conda == True:
        os.system('\n'.join(['cd ' + iraf_dir,
                          'source activate iraf27', #conda activate not working
                          'cl < wrapper_hasercalctest.cl']))
    else:
        os.system('\n'.join(['cd ' + iraf_dir,
                          'cl < wrapper_hasercalctest.cl']))

def clean_afrhotot(direc):
    """
    Clean the afrhotot files to keep only the last computed value for each image

    Parameters:
        direc (str): path to the working directory
    """
    for path, subdirs, files in os.walk(direc):
        for file in files:
            if 'afrho' in file and 'tot' in file and ".png" not in file:
                afrhofile = os.path.join(path, file)
                print(afrhofile)
                df = pd.read_csv(afrhofile, header=None, sep="\s+")
                df.drop_duplicates(subset =0, keep = 'last', inplace = True)
                df.to_csv(afrhofile, index=False, sep=" ", header=False)

def generate_center_comment(fold):
    centerlist = pd.read_csv(os.path.join(fold, 'centerlist'),header=None, sep=' '
                             ,usecols=[0,2,3,5,10],
                             names=['file', 'xcent', 'ycent', 'filt', 'ctnmethod'])
    
    center_comment = centerlist[['file', 'filt']]
    center_comment.drop_duplicates(subset ='file', keep = 'last', inplace = True)
    center_comment['comment'] = ''
    center_comment.to_csv(os.path.join(fold, 'center_comment'), index=False, sep=",", header=False)
    
def import_perihelion(file_path, update=False):
    #columns designation at https://www.minorplanetcenter.net/iau/info/CometOrbitFormat.html
    columns = ['id','pyear','pmonth','pday','pdist','e','param1','param2','param3',
                'param4','absmag','slope','name']
    if update == True:
        #url='132456789'
        url="https://www.minorplanetcenter.net/iau/MPCORB/CometEls.txt"
        tab=pd.read_fwf(url, header=None, widths=[13,18-13,22-18,30-22,40-30,50-40,60-50,70-60,
                        80-70,90-80,96-90,101-96,158-101])
        # except:
        #     print('could not download comet MPC database')
        tab.columns=columns
        tab.drop(['param1','param2','param3','param4'], inplace=True, axis=1)
        tab.to_csv(file_path, index=False)
    else:
        tab=pd.read_csv(file_path)
    return tab

#import_perihelion('/home/Mathieu/Documents/TRAPPIST/scripts/utils/test', update=False)
# generate_center_comment('/home/Mathieu/Documents/TRAPPIST/reduced_data/CK20R070/20220906TS/centering')
# import directory_structure
# ds = directory_structure.directory_structure()
# clean_afrhotot("/home/Mathieu/Documents/TRAPPIST/reduced_data/2020T2/20211004TS")
