#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  6 13:45:12 2021

@author: Mathieu Vander Donckt
mathieu.vanderdonckt@uliege.be

"""

import os
import pandas as pd
from astropy.io import fits
import datetime
import shutil

pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)

# def NAS_build(NAS_path, export_path, keyword=''):
#     """
#     Generate or update an indexed database of .fits and .fts file in the NAS of binning 2, containing keywords from the fits headers
#     An indexed database can be queryed much faster than going through all the files each time
    
#     Parameters
#     ----------
#     NAS_path : str
#         path the NAS to query
#     export_path : str
#         path for export of the database (with the filename)
#     keyword : str, optional
#         if different than '', needs to be present in the folder path for it to be considered. The default is ''.

#     Returns
#     -------
#     None.

#     """
    
#     dt = datetime.datetime.now()
#     fitstable = pd.DataFrame(columns=('file','object', 'type', 'filter', 'date', 'exptime', 'binning'))
    
#     # count = 0
#     for path, subdirs, files in sorted(os.walk(NAS_path)):
#         # count += 1
#         if keyword in path:
#             print(path)
#             for name in files:
#                 if name.endswith(".fits") or name.endswith(".fts"):
#                     try:
#                         with fits.open(os.path.join(path, name)) as hdul:
#                             try:
#                                 imobject = hdul[0].header['OBJECT']
#                             except:
#                                 imobject = None
#                             try:
#                                 imtype = hdul[0].header['IMAGETYP']
#                             except:
#                                 imtype = None
#                             try:
#                                 imfilter = hdul[0].header['FILTER']
#                             except:
#                                 imfilter = None
#                             try:
#                                 imdate = hdul[0].header['DATE-OBS']
#                             except:
#                                 imdate = None
#                             try:
#                                 imexptime = hdul[0].header['EXPTIME']
#                             except:
#                                 imexptime = None
#                             try:
#                                 imbinning = hdul[0].header['XBINNING']
#                             except:
#                                 imbinning = None
#                     except:
#                         print("error with", os.path.join(path, name))
#                         continue
#                 # print(os.path.join(path, name))
#                     if imbinning == 2:
#                         fitstable.loc[fitstable.shape[0]] = [os.path.join(path, name), imobject, imtype, imfilter, imdate, imexptime, imbinning]
#             # print(path)L
#         # if count % 10 == 0:
#             # print(path)
#             # print(len(fitstable.index))
    
#     fitstable.sort_values(by=['date']).to_csv(export_path, index=False)
    
#     print('Executed in ', datetime.datetime.now() - dt)
    

# NAS_build("/NASTN/Data_TrappistNord/ACP Astronomy/Images", "/home/Mathieu/Documents/TRAPPIST/raw_data/TN_query_update.txt", "202202")
# if __name__ == "__main__":
#     NAS_build("/NASTS2/Data_Trappist/Data_Trappist/ACP Astronomy/Images/2022", "/home/Mathieu/Documents/TRAPPIST/raw_data/TS2_query_update.txt", '202210')
#     NAS_build("/NASTN/Data_TrappistNord/ACP Astronomy/Images", "/home/Mathieu/Documents/TRAPPIST/raw_data/TN_query_update.txt", '202210')

def NAS_update(NAS_path, export_path, keyword='', cometlist = 'current'):
    
    dt = datetime.datetime.now()
    discarded_obj = []

    if os.path.isfile(export_path):
        print("importing existing table")
        fitstable = pd.read_csv(export_path)
        print("done importing")
    else:
        fitstable = pd.DataFrame(columns=('file','object', 'type', 'filter', 'date', 'exptime', 'binning', 'readmode'))
    
    blacklist_path = os.path.join('/'.join(export_path.split('/')[:-1]), 'blacklist_' +export_path.split('/')[-1] )
    if os.path.isfile(blacklist_path):
        blacklisttable = pd.read_csv(blacklist_path, header=None,names=['file', 'object'], sep=',')
        blacklist = pd.concat([fitstable['file'], blacklisttable['file']], ignore_index=True)
    else:
        blacklisttable = pd.DataFrame(columns=['file', 'object'])
        blacklist = fitstable['file']
    # print(blacklist.str.contains('/NASTN/Data_TrappistNord/ACP Astronomy/Images/20221024/Calibration/Bias-S001-R001-C008-B1.fts').any())
    # input()

    add_to_list = []
    # ['CK15E61R', 'CK15V020','CK17U010','CK17E030','CK14B010','CK14R030', '19K7',
    #                'CK19Q040','AK19U060','K14UR1N','CK09P010','CK10X010','CK10B010', 'CK10F87B',
    #                'CK11C010','CK09F040','CK09K030','CK10M010','CK06W030','CK11C030','PK11A020',
    #                '0073P      c',]
    if cometlist == 'current':
        from trapconfig import param
        from trap_reduction import import_perihelion
        objlist = import_perihelion(param['perihelion'], update=False)['id']
    elif cometlist == 'update':
        from trapconfig import param
        from trap_reduction import import_perihelion
        objlist = import_perihelion(param['perihelion'], update=True)['id']
    else:
        objlist = pd.read_csv(cometlist)
    objlist = objlist.append(pd.Series(add_to_list), ignore_index=True)
    
    for path, subdirs, files in sorted(os.walk(NAS_path)):
        # count += 1
        if keyword in path:
            print(datetime.datetime.now())
            print(path)
            for name in files:
                if name.endswith(".fits") or name.endswith(".fts"):                    # if blacklist.str.contains(os.path.join(path,name)).any() == False:
                    if (os.path.join(path,name) in blacklist.values) == False:
                        try:
                            with fits.open(os.path.join(path, name)) as hdul:
                                try:
                                    imobject = hdul[0].header['OBJECT']
                                except:
                                    imobject = None
                                try:
                                    imtype = hdul[0].header['IMAGETYP']
                                except:
                                    imtype = None
                                try:
                                    imfilter = hdul[0].header['FILTER']
                                except:
                                    imfilter = None
                                try:
                                    imdate = hdul[0].header['DATE-OBS']
                                except:
                                    imdate = None
                                try:
                                    imexptime = hdul[0].header['EXPTIME']
                                except:
                                    imexptime = None
                                try:
                                    imbinning = hdul[0].header['XBINNING']
                                except:
                                    imbinning = None
                                try:
                                    imreadmode = hdul[0].header['READOUTM']
                                except:
                                    imreadmode = None
                        except:
                            print("error with", os.path.join(path, name))
                            continue
                        # print(os.path.join(path, name), imobject, imtype, imfilter, imdate, imexptime, imbinning) 

                        fitstable.loc[fitstable.shape[0]] = [os.path.join(path, name), imobject, imtype, imfilter, imdate, imexptime, imbinning,imreadmode]
# =============================================================================
#                         if objlist.isin([imobject]).any():
#                             fitstable.loc[fitstable.shape[0]] = [os.path.join(path, name), imobject, imtype, imfilter, imdate, imexptime, imbinning,imreadmode]
#                             if imbinning != 2:
#                                 print(f"WARNING: found {name} ({imobject}) to be binning {imbinning}")
#                         elif imtype in (['DARK', 'Dark Frame','FLAT', 'Flat Frame','BIAS', 'Bias Frame']) and (imbinning == 2):
#                             fitstable.loc[fitstable.shape[0]] = [os.path.join(path, name), imobject, imtype, imfilter, imdate, imexptime, imbinning,imreadmode]
#                         else:
#                             blacklisttable.loc[blacklisttable.shape[0]] = [os.path.join(path, name), imobject]
#                             if (len(imobject) != 0) and (imobject in discarded_obj) == False:
#                                 discarded_obj.append(imobject)
# =============================================================================
                    # else:
                    #     print('rejected', os.path.join(path, name))
    fitstable.sort_values(by=['date']).to_csv(export_path, index=False)
    # blacklisttable.to_csv(blacklist_path, index=False, sep=",", header=False)
    # print('Discarded objects:', discarded_obj)
    print('Executed in ', datetime.datetime.now() - dt)
    
# if __name__ == "__main__":
#     NAS_update("/NASTN/Data_TrappistNord/ACP Astronomy/Images", "/home/Mathieu/Documents/TRAPPIST/raw_data/TN_query_all.db", '202311')
#     NAS_update("/NASTS2/Data_Trappist/Data_Trappist/ACP Astronomy/Images", "/home/Mathieu/Documents/TRAPPIST/raw_data/TS_query_all.db", '202311')
#     NAS_update("/NASTN/Data_TrappistNord/ACP Astronomy/Images", "/home/Mathieu/Documents/TRAPPIST/raw_data/TN_query_all.db", '202312')
#     NAS_update("/NASTS2/Data_Trappist/Data_Trappist/ACP Astronomy/Images", "/home/Mathieu/Documents/TRAPPIST/raw_data/TS_query_all.db", '202312')
    # NAS_update("/home/Mathieu/Documents/TRAPPIST/raw_data_TRAPPIST/CK17K020", "/home/Mathieu/Documents/TRAPPIST/prosetest_22E3_20230216/test.db", '')
    
    # NAS_update("/NASTS/Data_Trappist/ACP Astronomy/Images", "/home/Mathieu/Documents/TRAPPIST/raw_data/TS_query_all.db", '')
  
  
def queryZ(NAS_path):

    year = ("2017","2018","2019")
    
    count = 0
    for path, subdirs, files in os.walk(NAS_path):
        count += 1
        if year[0] in path or year[1] in path or year [2] in path:
            for name in files:
                if name.endswith(".fits.Z") or name.endswith(".fts.Z"):
                    count += 1
                    print(count, os.path.join(path, name))

def convertdate(date):
    """
    Convert date from the format in the fits into a datetime object.
    Used by the other scripts
    
    Parameters:
        date (str): date in the fits header format
    """
    if len(date) == 23:
        converted_date  = datetime.datetime.strptime(date, '%Y-%m-%dT%H:%M:%S.%f')
    elif len(date) == 19:
        converted_date = datetime.datetime.strptime(date, '%Y-%m-%dT%H:%M:%S')
    return converted_date

def date_to_startnight(date):
    """
    Take a date in the datetime format and return the night of the observation at 23:59 (e.g. the previous night if the observation is taken after midnight)
    Used by the other functions
    Used by other scripts
    
    Parameters:
        date (datetime): date in datetime format
    """
    if date.hour < 12:
        return pd.Timestamp(date.year, date.month, date.day, 23, 59) - datetime.timedelta(days=1)
    else:
        return pd.Timestamp(date.year, date.month, date.day, 23, 59)

def loadcsvtable(path):
    """
    Load an already existing indexed database and return it in the form of a pandas table
    Used by the other functions
    
    Parameters
    ----------
    path : str
        path to the database
    
    Returns
    -------
    None.

    """
      
    NASfitstable = pd.read_csv(path)
    # fitstable['date'] =  pd.to_datetime(fitstable['date'], format='%%Y-%m-%dT%H:%M:%S.%f')
    
    NASfitstable['date'] =  pd.to_datetime(NASfitstable['date'], format='%Y-%m-%dT%H:%M:%S.%f')
    NASfitstable['start_night'] = NASfitstable.apply(lambda row: date_to_startnight(row['date']), axis=1)
    return NASfitstable

# =============================================================================
# def datelist(fitstable):
#     unformated_nights = fitstable['start_date']
#     nights = []
#     for index, date in unformated_nights.items():
#         night = date_to_startnight(date)
#         if night not in nights:
#             nights.append(night)
#     #         print(night)
#     # print(len(nights))
#     return nights            
# =============================================================================
  
# =============================================================================
# def datelist(fitstable):
#     nightslist = []
#     timelist = fitstable.loc[fitstable['object'].isin(['CK20T020']), 'date'].drop_duplicates().values.tolist()
#     for date in timelist:
#         converted_date = convertdate(date)
#         # check if before or after midnight and set the night to the day the observation start
#         if converted_date.hour < 12:
#             night = datetime.datetime(converted_date.year, converted_date.month, converted_date.day -1)
#         else:
#             night = datetime.datetime(converted_date.year, converted_date.month, converted_date.day)
#         print(night)
#         if night not in nightslist:
#             nightslist.append(night)
#     print(nightslist)
#     return nightslist
# =============================================================================

def check_objects_names(startdate, enddate, NASfitstable, only_BVRI = False):
    """
    Makes a list of all the objects in a table between given dates
    if the dates are at 12:00, will include the startdate as starting night but not the enddate as starting night
    
    Parameters:
        startdate (datetime): starting time
        enddate (datetime): ending time
        NASfitstable (pandas dataframe): indexed database-style table
    """
    
       
    objects_table = NASfitstable.loc[(NASfitstable['start_night'] >= startdate)
                                & (NASfitstable['start_night'] <= enddate)
                                & (NASfitstable['type'].isin(['LIGHT', 'Light Frame']))]
    
    if only_BVRI == True:
        objects_table = objects_table.loc[objects_table['filter'].isin(['B','V','R','I','Rc','Ic'])]
    # objects_list = objects_table['object'].drop_duplicates().values.tolist() #old way to do it
    objects_names = objects_table['object'].drop_duplicates().reset_index(drop=True)
    objects = pd.DataFrame({'object':objects_names})
    for index, row in objects.iterrows():
        objects.loc[objects['object'] == row['object'], 'nb_nights'] = len(objects_table.loc[objects_table['object'] == row['object']]['start_night'].drop_duplicates())
        objects.loc[objects['object'] == row['object'], 'nb_images'] = len(objects_table.loc[objects_table['object'] == row['object']])
        objects.loc[objects['object'] == row['object'], 'nb_NB'] = len(objects_table.loc[(objects_table['object'] == row['object'])
                                                                                         & objects_table['filter'].isin(['OH','CN','C2','C3','NH','RC','BC','GC','UC'])])
    # objects.astype({'nb_nights': 'int32','nb_images': 'int32'}).dtypes
    return objects
    
def get_files(obj_name, NASfitstable, output_path, dayinterval, dateinterval=['',''], skip_existing=False,only_BVRI=False):
    '''Look for and prompt for downloading files in the NAS
    It is advised to not mix TN and TS queries
    look on the TN or TS NAS for all the nights containing the object (obj_name).
    Night by night, look for the calibration files for a given dayinterval.
    Download the lights and found calibs and put them in a date folder.
    Date interval can be given to select data betwenn two starting nights intervals
    
    Parameters:
        obj_name (str): object name in the NASfitstable
        NASfitstable (pandas dataframe): indexed table-style table
        output_path (str): path to a folder to save the files
        dayinterval (int): day interval to look for calib files before and after the observation night
        date_interval (list [datetime, datetime], optional, default=['','']): of 2 terms: the starting date for the search and second the ending date, in datetime format
            if '' given for one of the argument, consider no date limit on that end
        skip_existing(boolean, default = False): if set to True, will skip download if output_path already exists
            if False will ask if delete old dir and download or skip.
    '''
    
    # starting and ending date for the query
    if dateinterval[0] == '':
        startdate = pd.to_datetime('2000-01-01', format='%Y-%m-%d')
    else:
        startdate = dateinterval[0]
    if dateinterval[1] == '':
        enddate = pd.to_datetime('2200-01-01', format='%Y-%m-%d')
    else:
        enddate = dateinterval[1]

    LIGHTS = NASfitstable.loc[NASfitstable['object'].isin([obj_name])
                              & (NASfitstable['start_night'] >= startdate)
                              & (NASfitstable['start_night'] <= enddate)]
    if only_BVRI == True:
        LIGHTS = LIGHTS.loc[LIGHTS['filter'].isin(['B','V','R','I','Rc','Ic'])]
        
    nightslist = LIGHTS['start_night'].drop_duplicates().tolist()

    for night in nightslist:
        lower_interval = night - datetime.timedelta(days = dayinterval)
        upper_interval = night + datetime.timedelta(days = dayinterval)
        
        LIGHTtable = LIGHTS.loc[(LIGHTS['start_night'] == night)]
        
        filtlist = LIGHTtable['filter'].drop_duplicates().values.tolist()
        binnings = LIGHTtable['binning'].drop_duplicates().values.tolist()
        if len(binnings) != 1:
            input(f'WARNING found more than one binning ({binnings}) in light frames. ')
        binning = binnings[0]
        exptimelist = LIGHTtable['exptime'].drop_duplicates().values.tolist()
        exptimelist.append(10)
        exptimelist.append(15)
        
        FLATtable = NASfitstable.loc[NASfitstable['type'].isin(['FLAT', 'Flat Frame'])
                                     & (NASfitstable['start_night'] > lower_interval)
                                     & (NASfitstable['start_night'] <= upper_interval)
                                     & NASfitstable['filter'].isin(filtlist)
                                     & (NASfitstable['binning'] == binning)]
        BIAStable = NASfitstable.loc[NASfitstable['type'].isin(['BIAS', 'Bias Frame'])
                                      & (NASfitstable['start_night'] > lower_interval)
                                      & (NASfitstable['start_night'] <= upper_interval)
                                      & (NASfitstable['binning'] == binning)]
        DARKtable = NASfitstable.loc[NASfitstable['type'].isin(['DARK', 'Dark Frame'])
                                      & (NASfitstable['start_night'] > lower_interval)
                                      & (NASfitstable['start_night'] <= upper_interval)
                                      & NASfitstable['exptime'].isin(exptimelist)
                                      & (NASfitstable['binning'] == binning)]

        # continue to search for flats if less than 5
        for filt in filtlist:
            flat_li = lower_interval
            flat_ui = upper_interval
            while len(FLATtable.loc[FLATtable['filter'] == filt]) < 5:
                flat_li -= pd.Timedelta(days=7)
                flat_ui += pd.Timedelta(days=7)
                print(f'Looking for flats {filt} with bigger interval')
                new_FLAT = NASfitstable.loc[NASfitstable['type'].isin(['FLAT', 'Flat Frame'])
                                             & (NASfitstable['start_night'] > flat_li)
                                             & (NASfitstable['start_night'] <= flat_ui)
                                             & (NASfitstable['filter'] == filt)
                                             & (NASfitstable['binning'] == binning)]

                FLATtable = FLATtable.loc[FLATtable['filter'] != filt]
                FLATtable = pd.concat([FLATtable, new_FLAT])

        all_calibtable = pd.concat([FLATtable,BIAStable,DARKtable])
        print(night, "lights",len(LIGHTtable), "flats",len(FLATtable),"bias", len(BIAStable),"darks", len(DARKtable))

        # check if the directory already exists and ask if erase if already exists
        output_dir = os.path.join(output_path, str(night.strftime('%Y%m%d')))
        while True:
            if os.path.exists(output_dir) and skip_existing == False:
                inp = input(output_dir + ' already existing. Download data anyway and delete old dir (d) or skip dir (s)?')
                if inp == 'd' or inp =="D":
                    shutil.rmtree(output_dir)
                    os.makedirs(output_dir)
                    print("created", output_dir)
                    for index, row in LIGHTtable.iterrows():
                        shutil.copyfile(row['file'], os.path.join(output_dir, row['file'].split('/')[-1]))
                    # same for calib files
                    os.makedirs(os.path.join(output_dir, "Calibration"))
                    print("created", os.path.join(output_dir, "Calibration"))
                    count = 0
                    for index, row in all_calibtable.iterrows():
                        count +=1
                        shutil.copyfile(row['file'], os.path.join(output_dir, "Calibration", str(count) + row['file'].split('/')[-1]))
                    break
                elif inp == 's' or inp =="S":
                    print('skipping ' + output_dir)
                    break
                else:
                    continue
            elif os.path.exists(output_dir) and skip_existing == True:
                print(output_dir + ' already existing. skip_existing = True, skipping download')
                break
            else:
                os.makedirs(output_dir)
                print("created", output_dir)
                for index, row in LIGHTtable.iterrows():
                    shutil.copyfile(row['file'], os.path.join(output_dir, row['file'].split('/')[-1]))
                # same for calib files
                os.makedirs(os.path.join(output_dir, "Calibration"))
                print("created", os.path.join(output_dir, "Calibration"))
                count = 0
                for index, row in all_calibtable.iterrows():
                    count +=1
                    shutil.copyfile(row['file'], os.path.join(output_dir, "Calibration", str(count) + row['file'].split('/')[-1]))
                break
                break

        # if not os.path.exists(output_dir):
        #     os.makedirs(output_dir)
        #     print("created", output_dir)
        # for index, row in LIGHTtable.iterrows():
        #     shutil.copyfile(row['file'], os.path.join(output_dir, row['file'].split('/')[-1]))
            
        # count = 0
        # if not os.path.exists(os.path.join(output_dir, "Calibration")):
        #     os.makedirs(os.path.join(output_dir, "Calibration"))
        #     print("created", os.path.join(output_dir, "Calibration"))
        # for index, row in all_calibtable.iterrows():
        #     count +=1
        #     shutil.copyfile(row['file'], os.path.join(output_dir, "Calibration", str(count) + row['file'].split('/')[-1]))
            

# from trap_reduction import *

# NASfitstable = loadcsvtable("/home/Mathieu/Documents/TRAPPIST/raw_data/TS_query.txt")
# output_path = "/home/Mathieu/Documents/TRAPPIST/raw_data/CK17K020/TS/"
# obj_name = 'CK17K020'
# get_files(obj_name, NASfitstable, output_path, dayinterval=7, filt_list=[], dateinterval=('2022-04-01',''))


# for item in NASfitstable.loc[NASfitstable['object'].isin([obj_name]) & NASfitstable['file'].str.contains('20210609'), 'file']:
#     print(item)

###look for specific calib files around a given date
def lookforcalib(NASfitstable, imtype, output_fold, night, obj='', exptime=15, filt='R', dayinterval=0, binning = 2):
    """
    Look for a specific set of file the closest to the observation file in the database and download it. By default look for the files closest to the observation night
    
    Parameters:
        NASfitstable (pandas dataframe): indexed-type database, can be imported with loadcsvtable()
        imtype (str): image type. can take values of 'light', 'dark', 'flat' or 'bias'
        output_fold (str): path to the output directory.
            Inside this directory needs a subfolder with the observation night in YYYYMMDD format and a subsubfolder Calibration
        night (str YYYYMMDD): night of observation in an int (YYYYMMDD)
        obj (str, semi-optional, defaut=''): name of the object if imtype='light'
        exptime (int, semi-optional, defaut=15): exposure time if imtype='dark'
        filt (str, semi-optional, default='R'): name of the filter if imtype='flat'
        dayinterval (int, optional, default=0): the starting point in time to look for image from the observation night
    """
    ####################################
    ### PARAMETERS
    
    if imtype == 'light':
        imtype = ['LIGHT', 'Light Frame']
    elif imtype == 'dark':
        imtype = ['DARK', 'Dark Frame']
    elif imtype == 'bias':
        imtype = ['BIAS', 'Bias Frame']
    elif imtype == 'flat':
        imtype = ['FLAT', 'Flat Frame']
    # elif imtype == 'BC':
    #     print("Let's roll for the special BC round")
    else:
        print('wrong imtype')
        return
    
    ####################################
    
    night = (int(night[0:4]), int(night[4:6]), int(night[6:8]))
    
    obsnight = pd.Timestamp(year=night[0], month=night[1], day=night[2], hour=23, minute=59)

    while True:
        lower_interval = obsnight - datetime.timedelta(days = dayinterval, hours = 12)
        upper_interval = obsnight + datetime.timedelta(days = dayinterval, hours = 12)
        if 'DARK' in imtype:
            calibtable = NASfitstable.loc[(NASfitstable['date'] > lower_interval)
                                          & (NASfitstable['date'] <= upper_interval)
                                          # & (NASfitstable['binning'] == 2)
                                          & NASfitstable['type'].isin(imtype)
                                          & (NASfitstable['exptime'] == exptime)
                                          & (NASfitstable['binning'] == binning)]
        elif 'FLAT' in imtype:
            calibtable = NASfitstable.loc[(NASfitstable['date'] > lower_interval)
                                          & (NASfitstable['date'] <= upper_interval)
                                          # & (NASfitstable['binning'] == 2)
                                          & NASfitstable['type'].isin(imtype)
                                          & (NASfitstable['filter'] == filt)
                                          & (NASfitstable['binning'] == binning)]
        elif 'LIGHT' in imtype:
            calibtable = NASfitstable.loc[(NASfitstable['date'] > lower_interval)
                                          & (NASfitstable['date'] <= upper_interval)
                                          # & (NASfitstable['binning'] == 2)
                                          & NASfitstable['type'].isin(imtype)
                                          & (NASfitstable['object'] == obj)
                                          & (NASfitstable['filter'] == filt)
                                          ]
        elif 'BIAS' in imtype:
            calibtable = NASfitstable.loc[(NASfitstable['date'] > lower_interval)
                                          & (NASfitstable['date'] <= upper_interval)
                                          # & (NASfitstable['binning'] == 2)
                                          & NASfitstable['type'].isin(imtype)
                                          & (NASfitstable['binning'] == binning)]
        # elif imtype == 'BC':
        #     calibtable = NASfitstable.loc[(NASfitstable['date'] > lower_interval)
        #                                   & (NASfitstable['date'] <= upper_interval)
        #                                   & (NASfitstable['binning'] == 2)
        #                                   & NASfitstable['type'].isin(imtype)
        #                                   & (NASfitstable['object'] == obj)
        #                                   & (NASfitstable['filter'] == filt)
        #                                   ]
            
            # print(NASfitstable['object'][0])
        if len(calibtable.index) == 0:
            dayinterval += 1
        else:
            caliblist = []
            print('\nFound:')
            for item in sorted(calibtable['file']):
                print(item)
                caliblist.append(item)
            # print(calibtable)
            
                
            ### copy the files from the NAS into the output folder 
            # date for the subfolder
            subf_year = str(night[0])
            if night[1] > 9:
                subf_month = str(night[1])
            else:
                subf_month = '0' + str(night[1])
            if night[2] > 9:
                subf_day = str(night[2])
            else:
                subf_day = '0' + str(night[2])
            subfold_name = subf_year + subf_month + subf_day
            if 'LIGHT' in imtype:
                output_fold = os.path.join(output_fold, subfold_name)
            else:
                output_fold = os.path.join(output_fold, subfold_name, "Calibration")
            print("\noutput_fold set to", output_fold)
            
            print("\nDay interval of ", dayinterval)
            trans = input("Transfer files?: (y/n)")
            count = 0
            if trans == 'y' or trans == 'Y':
                for item in caliblist:
                    # print(os.path.join(output_fold, item.split('/')[-1]))
                    if 'DARK' in imtype:
                        shutil.copy(item, os.path.join(output_fold, 'extra' + str(count) + '_' + str(exptime) + item.split('/')[-1]))
                    elif 'FLAT' in imtype:
                        shutil.copy(item, os.path.join(output_fold, 'extra' + str(count) + '_' + filt + item.split('/')[-1]))
                    elif 'BIAS' in imtype:
                        shutil.copy(item, os.path.join(output_fold, 'extra' + str(count) + '_' + '0' + item.split('/')[-1]))
                    elif 'LIGHT' in imtype:
                        shutil.copy(item, os.path.join(output_fold, item.split('/')[-1]))
                    print("copied", item.split('/')[-1])
                    count +=1
            # if imtype == 'BC':
            #     BCexptime = caliblist.iloc[0, 'exptime']
            #     BCdate = caliblist.iloc[0, 'date']
                
            break

def lookforcalib_old():
    """
    old version
    """
    ####################################
    ### PARAMETERS
    
    ### uncomment line bellow if querying for a light image
    # imtype = ['LIGHT', 'Light Frame']
    obj = "0144P" #target name in the fits header. only for lights and for the output path
    
    ### uncomment line bellow if querying for dark frames
    imtype = ['DARK', 'Dark Frame']
    exptime = 900 #exposure time. only for darks
        
    ### uncomment line bellow if querying for flat frames
    # imtype = ['FLAT', 'Flat Frame']
    filt = 'I' #filter. only for flats
    
    ### uncomment line bellow if querying for bias frames
    # imtype = ['BIAS', 'Bias Frame']
    
    telescope = 'TN'
    night = (2024,2,5) ### set the observation night
    
    dayinterval = 10
    # starting point for the search
    
    binning = 2
    
    ####################################
    
    NASfitstable = loadcsvtable("/home/Mathieu/Documents/TRAPPIST/raw_data/" + telescope + "_query_all.db") ### path to the indexed database
    # print(NASfitstable.loc[ (NASfitstable['binning'] == 2)])
    # output_fold = "/home/Mathieu/Documents/TRAPPIST/raw_data/2020T2/TS/20210703/Calibration" ### path to the output folder

    subf_year = str(night[0])
    if night[1] > 9:
        subf_month = str(night[1])
    else:
        subf_month = '0' + str(night[1])
    if night[2] > 9:
        subf_day = str(night[2])
    else:
        subf_day = '0' + str(night[2])
    subfold_name = subf_year + subf_month + subf_day
    output_fold = "/home/Mathieu/Documents/TRAPPIST/raw_data_TRAPPIST/" + obj + "/" + telescope + "/" +  subfold_name + "/Calibration"
    
    obsnight = pd.Timestamp(year=night[0], month=night[1], day=night[2], hour=23, minute=59)
    
    while True:
        lower_interval = obsnight - datetime.timedelta(days = dayinterval, hours = 12)
        upper_interval = obsnight + datetime.timedelta(days = dayinterval, hours = 12)
        if 'DARK' in imtype:
            calibtable = NASfitstable.loc[(NASfitstable['date'] > lower_interval)
                                          & (NASfitstable['date'] <= upper_interval)
                                          # & (NASfitstable['binning'] == 2)
                                          & NASfitstable['type'].isin(imtype)
                                          & (NASfitstable['exptime'] == exptime)
                                          & (NASfitstable['binning'] == binning)]
        elif 'FLAT' in imtype:
            calibtable = NASfitstable.loc[(NASfitstable['date'] > lower_interval)
                                          & (NASfitstable['date'] <= upper_interval)
                                          # & (NASfitstable['binning'] == 2)
                                          & NASfitstable['type'].isin(imtype)
                                          & (NASfitstable['filter'] == filt)
                                          & (NASfitstable['binning'] == binning)]
        elif 'LIGHT' in imtype:
            calibtable = NASfitstable.loc[(NASfitstable['date'] > lower_interval)
                                          & (NASfitstable['date'] <= upper_interval)
                                          # & (NASfitstable['binning'] == 2)
                                          & NASfitstable['type'].isin(imtype)
                                          & (NASfitstable['object'] == obj)
                                            & (NASfitstable['filter'] == filt)
                                          ]
        elif 'BIAS' in imtype:
            calibtable = NASfitstable.loc[(NASfitstable['date'] > lower_interval)
                                          & (NASfitstable['date'] <= upper_interval)
                                          # & (NASfitstable['binning'] == 2)
                                          & NASfitstable['type'].isin(imtype)
                                          & (NASfitstable['binning'] == binning)]
            # print(NASfitstable['object'][0])
        if len(calibtable.index) == 0:
            dayinterval += 1
        else:
            # print(calibtable)
            caliblist = []
            print('\n')
            for item in sorted(calibtable['file']):
                print(item)
                caliblist.append(item)
                
            ### copy the files from the NAS into the output folder 
            print("\noutput_fold set to", output_fold)
            
            print("\nday interval of", dayinterval)
            
            trans = input("\ntransfer files?: (y/n)")
            count = 0
            if trans == 'y' or trans == 'Y':
                for item in caliblist:
                    # print(os.path.join(output_fold, item.split('/')[-1]))
                    if 'DARK' in imtype:
                        shutil.copy(item, os.path.join(output_fold, 'extra' + str(count) + '_' + str(exptime) + item.split('/')[-1]))
                    elif 'FLAT' in imtype:
                        shutil.copy(item, os.path.join(output_fold, 'extra'  + str(count) + '_'+ filt + item.split('/')[-1]))
                    elif 'BIAS' in imtype:
                        shutil.copy(item, os.path.join(output_fold, 'extra' + str(count) + '_' + '0' + item.split('/')[-1]))
                    elif 'LIGHT' in imtype:
                        shutil.copy(item, os.path.join(output_fold[:-12], item.split('/')[-1]))
                    print("copied", item.split('/')[-1])
                    count +=1
            break

#if __name__ == "__main__":
#    lookforcalib_old()

if __name__ == "__main__":
    # NAS_update("/NASTN/Data_TrappistNord/ACP Astronomy/Images", "/home/Mathieu/Documents/TRAPPIST/raw_data/TN_query_all.db", '202402')
    NAS_update("/NASTN/Data_TrappistNord/ACP Astronomy/Images", "/home/Mathieu/Documents/TRAPPIST/raw_data/TN_query_all.db", '202403')
    # NAS_update("/NASTS2/Data_Trappist/Data_Trappist/ACP Astronomy/Images", "/home/Mathieu/Documents/TRAPPIST/raw_data/TS_query_all.db", '202402')
    NAS_update("/NASTS2/Data_Trappist/Data_Trappist/ACP Astronomy/Images", "/home/Mathieu/Documents/TRAPPIST/raw_data/TS_query_all.db", '202403')
 
