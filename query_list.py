#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  9 12:34:33 2022

@author: Mathieu Vander Donckt
"""

import query_NAS
from trapconfig import param
import pandas as pd
import datetime
import os
import shutil

copyfile = True #for test purpose
make_mapping = True

filt_list = ['B','V','R','I']
comet_list = ['0009P',
'0046P',
'0088P',
'0103P',
'0246P',
'CK09F040',
'CK09P010',
'CK11L040',
'CK12F060',
'CK12K010',
'CK13A010',
'CK13R010',
'CK13U10S',
'CK15E61R',
'CK19Q040']

TStab = (query_NAS.loadcsvtable(param['TS_qNAS']), 'TS')
TNtab = (query_NAS.loadcsvtable(param['TN_qNAS']), 'TN')
dt = pd.DataFrame()
dayinterval = 15
output_dir = '/NASTN/BVRI_comets_backup'
if os.path.exists(output_dir):
    shutil.rmtree(output_dir)
os.mkdir(output_dir)
TS_calib_dir = os.path.join(output_dir, 'TS_calib')
TN_calib_dir = os.path.join(output_dir, 'TN_calib')
os.mkdir(TS_calib_dir)
os.mkdir(TN_calib_dir)

def file_path_renaming(row):
    # print(row.str.split('/')['file'][-1])
    # print(type(row.str.split('/')))
    print(row)
    return '1'
    # return os.path.join(f'{obs}_calib',row.str.split('/')[-1])

for comet in comet_list:
    comet_dir = os.path.join(output_dir, comet)
    os.mkdir(comet_dir)
    for tabandobs in [TStab, TNtab]:
        tab = tabandobs[0]
        obs = tabandobs[1]
        
        lighttab = tab.loc[(tab['object'] == comet)
                             & (tab['filter'].isin(filt_list))]
        nightslist = lighttab['start_night'].drop_duplicates().tolist()
        print(comet, obs, len(nightslist))
        if len(nightslist) > 0:
            for night in nightslist:
                night_dir = os.path.join(comet_dir, str(night)[0:4] + str(night)[5:7] + str(night)[8:10] + obs)
                os.mkdir(night_dir)
                #takes the lights
                lights = lighttab.loc[lighttab['start_night'] == night]
                
                #getting the calib
                lower_interval = night - datetime.timedelta(days = dayinterval)
                upper_interval = night + datetime.timedelta(days = dayinterval)
                filtlist = lights['filter'].drop_duplicates().values.tolist()
                exptimelist = lights['exptime'].drop_duplicates().values.tolist()
                # exptimelist.append(10)
                exptimelist.append(15)
                
                flats = tab.loc[tab['type'].isin(['FLAT', 'Flat Frame'])
                                             & (tab['start_night'] > lower_interval)
                                             & (tab['start_night'] <= upper_interval)
                                             & tab['filter'].isin(filtlist)]
                bias = tab.loc[tab['type'].isin(['BIAS', 'Bias Frame'])
                                              & (tab['start_night'] > lower_interval)
                                              & (tab['start_night'] <= upper_interval)]
                darks = tab.loc[tab['type'].isin(['DARK', 'Dark Frame'])
                                              & (tab['start_night'] > lower_interval)
                                              & (tab['start_night'] <= upper_interval)
                                              & tab['exptime'].isin(exptimelist)]
                nighttable = pd.concat([lights, flats, bias, darks])
                #check the calib. More complicated than should be but took the code from trap_reduction
                checktable = lights.reindex(columns = ['file', 'filter', 'exptime'])
                checktable['nb_flat'] = None
                checktable['nb_dark'] = None
                checktable['nb_bias'] = len(bias)
                for index, row in checktable.iterrows() :
                    warning_flag = False
                    map_flats = flats.loc[flats['filter'] == row['filter']]
                    map_darks = darks.loc[darks['exptime'] == row['exptime']]
                    checktable.at[index, 'nb_flat'] = len(map_flats)
                    checktable.at[index, 'nb_dark'] = len(map_darks)              
                    if checktable.at[index, 'nb_flat'] == 0:
                        print(night)
                        # print("WARNING: no flat for", row['file'])
                        warning_flag = True
                        allflats = tab.loc[tab['type'].isin(['FLAT', 'Flat Frame']) & tab['filter'].isin([row['filter']])]
                        nextdate = allflats.iloc[(allflats['start_night'] -night).abs().argsort(),:].iloc[0]['start_night']
                        sup_flat = allflats.loc[tab['start_night'] == nextdate]
                        # print(row['filter'], nextdate, sup_flat)
                        map_flats = sup_flat
                        nighttable = pd.concat([nighttable, sup_flat])
                    elif checktable.at[index, 'nb_flat'] < 4:
                        # print(night)
                        print(row['nb_flat'])
                        print("WARNING: less than 5 flats for", row['file'])
                        warning_flag = True
                    if checktable.at[index, 'nb_dark'] == 0:
                        # print("WARNING: no dark (",row['exptime'],") for", row['file'])
                        warning_flag = True
                        alldarks = tab.loc[tab['type'].isin(['DARK', 'Dark Frame']) & (tab['start_night'] > lower_interval) & (tab['start_night'] <= upper_interval)]
                        closest_exptime = alldarks.iloc[(alldarks['exptime'] -row['exptime']).abs().argsort(),:].iloc[0]['exptime']
                        sup_dark = alldarks.loc[alldarks['exptime'] == closest_exptime]
                        map_darks = sup_dark
                        # print(row['exptime'], closest_exptime, sup_dark)
                        # input()
                        nighttable = pd.concat([nighttable, sup_dark])
                    elif checktable.at[index, 'nb_dark'] < 5:
                        print("WARNING: less than 5 darks (",row['exptime'],") for", row['file'])
                        warning_flag = True
                    if row['nb_bias'] == 0:
                        print("WARNING: no bias for", row['file'])
                        warning_flag = True
                    elif row['nb_bias'] < 5:
                        print("WARNING: less than 5 bias for", row['file'])
                        warning_flag = True
                    
                    if copyfile == True:
                        shutil.copy(row['file'], os.path.join(night_dir, row['file'].split('/')[-1]))

                    map_content = pd.concat([map_flats, map_darks, bias])
                    # map_content.loc[map_content['file'],'file'] = map_content['file'].split('/')
                    if make_mapping == True:
                        for index2, row2 in map_content.iterrows():
                            date = str(row2['start_night'])[0:4] + str(row2['start_night'])[5:7] + str(row2['start_night'])[8:10]
                            map_content.loc[index2, 'file'] = obs + '_calib/' + date + row2['file'].split('/')[-1]
                        # map_content['file'] = [item.split('/')[-1] for item in map_content['file']]
                        map_content.to_csv(os.path.join(night_dir, row['file'].split('/')[-1] + '_calibmapping.txt'), index=False)
                    
                map_dark15 = darks.loc[darks['exptime'] == 15]
                nb_flat_dark15 = len(map_dark15)
                if (nb_flat_dark15 == 0 ):
                    # print("WARNING: less than 5 darks (15s) for flats")
                    alldarks = tab.loc[tab['type'].isin(['DARK', 'Dark Frame']) & (tab['start_night'] > lower_interval) & (tab['start_night'] <= upper_interval)]
                    closest_exptime = alldarks.iloc[(alldarks['exptime'] -15).abs().argsort(),:].iloc[0]['exptime']
                    sup_dark15 = alldarks.loc[alldarks['exptime'] == closest_exptime]
                    map_dark15 = sup_dark15
                    # print('15', closest_exptime, sup_dark15)
                    # input()
                    nighttable = pd.concat([nighttable, sup_dark15])
                elif (nb_flat_dark15 < 5):
                    print("WARNING: less than 5 darks (15s) for flats") 
                # input()
                
                if make_mapping == True:
                    for index2, row2 in map_dark15.iterrows():
                        date = str(row2['start_night'])[0:4] + str(row2['start_night'])[5:7] + str(row2['start_night'])[8:10]
                        map_dark15.loc[index2, 'file'] = obs + '_calib/' + date + row2['file'].split('/')[-1]
                    map_dark15.to_csv(os.path.join(night_dir, 'flats_calibmapping.txt'), index=False)


                dt = pd.concat([dt, nighttable])
                dt_clean = dt.drop_duplicates(subset='file', keep='first', inplace=False)
                
obs = None
print(len(dt),len(dt_clean))
calib_table = dt_clean.loc[~dt_clean['type'].isin(['LIGHT', 'Light Frame'])]
print(len(calib_table))

if copyfile == True:
    for index, row in calib_table.iterrows():
        date = str(row['start_night'])[0:4] + str(row['start_night'])[5:7] + str(row['start_night'])[8:10]
        filename = date + row['file'].split('/')[-1]
        if 'TS' in row['file'].split('/')[1]:
            calib_dir = TS_calib_dir
        if 'TN' in row['file'].split('/')[1]:
            calib_dir = TN_calib_dir
            
        shutil.copy(row['file'], os.path.join(calib_dir, date + row['file'].split('/')[-1]))
                    