#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 15 13:14:23 2022

@author: Mathieu

Phase angle correction
"""

from trapconfig import param
import os
import pandas as pd
import get_ephem
import datetime

# available on https://asteroid.lowell.edu/comet/dustphaseHM_table.txt
pf_schleicher = pd.read_csv(os.path.join(param['calib'], 'pf_schleicher.dat'),
                            skiprows=7, names=['pa','0deg','90deg'], sep='\s+')

def import_pa_from_eph(imname, target, observatory): #function to import and format the ephemeris
    eph = get_ephem.ephemeris()
    if (len(target) == 8) and (target[0:2] + target[5] + target[7] == 'CK00'):
        eph.parameters['COMMAND'] = '20' + target[2:4] + ' ' + target[4] + target[6]
    elif (len(target) == 5) and (target[0:3] + target[-1] == '000P'):
        eph.parameters['COMMAND'] = target[3:]
    elif (len(target) == 5) and (target[0:2] + target[-1] == '00P'):
        eph.parameters['COMMAND'] = target[2:]
    elif (len(target) == 5) and (target[0:1] + target[-1] == '0P'):
        eph.parameters['COMMAND'] = target[1:]
    else:
        eph.parameters['COMMAND'] = target
    eph.parameters['QUANTITIES'] = '24'
    if observatory == 'TN':
        eph.parameters['CENTER'] = 'Z53@399'
    elif observatory == 'TS':
        eph.parameters['CENTER'] = 'I40@399'
    eph.parameters['START_TIME'] = datetime.datetime.strptime(imname, 'TRAP.%Y-%m-%dT%H:%M:%S.fits')
    eph.parameters['STOP_TIME'] = eph.parameters['START_TIME'] + datetime.timedelta(minutes=1)
    eph.parameters['STEP_SIZE'] = '1 m'
    eph.query_horizons()
    if eph.query_result[-4][:10] == '    Author':
        pass
    elif "To SELECT, enter record # (integer), followed by semi-colon.)" in eph.query_result[-3]:
        last_line = eph.query_result[-5]
        eph.record = [ elem for elem in last_line.split(" ") if elem != ''][0]
        #check if the record is the last epoch
        year_record = [ elem for elem in last_line.split(" ") if elem != ''][1]
        previous_year = [ elem for elem in eph.query_result[-6].split(" ") if elem != ''][1]
        if int(previous_year) > int(year_record):
            print('selected record: ', eph.record)
            input('WARNING: check the epoch selected is the las one')
        eph.parameters['COMMAND'] = eph.record
        eph.query_horizons()
    else:
        print(eph.query_result)
        input('debug this')
    data = []
    for line in eph.query_result[eph.query_result.index('$$SOE')+1:eph.query_result.index('$$EOE')]:
        data.append((imname +',' + line[:-1]).split(','))
        break #only need the first line here
    df = pd.DataFrame(data, columns=('imname', 'date','dateJD', 'sol_marq', 'lun_marq', 'pa')).drop(['date','dateJD', 'sol_marq', 'lun_marq'],axis=1)
    df = df.astype({'pa':'float'})
    df.reset_index(drop=True, inplace=True)
    return df

def generate_palist_reddir(dir_to_scan):
    # to add the files in already existing reduced dir
    for path, subdirs, files in os.walk(dir_to_scan):
        if path.endswith('TS') or path.endswith('TN'):
            print(path)
            comet = path.split('/')[-2]
            observatory = path.split('/')[-1][-2:]
            centerfile = os.path.join(path, 'centering', 'centerlist')
            tab = pd.read_csv(centerfile, header=None, sep="\s+")
            tab.drop_duplicates(subset =0, keep = 'last', inplace = True)
            pa_table = pd.DataFrame(columns=('imname','pa'))
            for index, row in tab.iterrows():
                pa = import_pa_from_eph(row[0], comet, observatory)
                pa_table = pd.concat([pa_table, pa])
            export_path = os.path.join(path, 'centering', 'palist')
            pa_table.sort_values(by=['imname']).to_csv(export_path, index=False)
            
# generate_palist_reddir(param['reduced'])

def generate_palist_tmpout(comet, observatory, working_folder=param['tmpout']):
    print('generate phase angle file')
    pa_table = pd.DataFrame(columns=('imname','pa'))
    for path, subdirs, files in os.walk(working_folder):
        for file in files:
            if file.startswith('TRAP.') and file.endswith(".fits"):
                pa = import_pa_from_eph(file, comet, observatory)
                pa_table = pd.concat([pa_table, pa])
        path = os.path.join(working_folder, 'palist')
        pa_table.sort_values(by=['imname']).to_csv(path, index=False)

def schleicher_0deg(afrho, pa):
    u_pa = pf_schleicher.loc[pf_schleicher['pa'] > pa][:1]['pa'].values[0]
    l_pa = pf_schleicher.loc[pf_schleicher['pa'] <= pa][-1:]['pa'].values[0]
    u_pf = pf_schleicher.loc[pf_schleicher['pa'] > pa][:1]['0deg'].values[0]
    l_pf = pf_schleicher.loc[pf_schleicher['pa'] <= pa][-1:]['0deg'].values[0]
    pf =  l_pf + (u_pf-l_pf)*((pa-l_pa)/(u_pa-l_pa))
    
    afrho_corr = afrho / pf
    return afrho_corr
    
schleicher_0deg(10, 5.9)