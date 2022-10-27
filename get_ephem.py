#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 10 16:29:43 2021

@author: Mathieu Vander Donckt
mathieu.vanderdonckt@uliege.be
"""

import requests
from astropy.io import fits
import os
import datetime
import pandas as pd

# home_dir = os.path.expanduser('~')
# fits_dir = os.path.join(home_dir, "Documents/TRAPPIST/tmpout/")
# output_dir = os.path.join(home_dir, "Documents/TRAPPIST/")

# print('Fits input directory set to: ' + fits_dir)
# print('Output directory set to: ' + output_dir)


class ephemeris:
    """
    Initiate an ephemeris class with defaults parameters to query NASA Horizons
    """
    
    def __init__(self):
        self.parameters = {
            'COMMAND' : None,
            'OBJ_DATA' : 'YES',
            'MAKE_EPHEM' : 'YES',
            'TABLE_TYPE' : 'OBSERVER',
            'CENTER' : None,
            'START_TIME' : None,
            'STOP_TIME' : None,
            'STEP_SIZE' : '10 m',
            'QUANTITIES' : '1,9,19,20',
            'CAL_FORMAT' : 'BOTH',
            'CSV_FORMAT' : 'YES',
            'REF_SYSTEM' : 'ICRF', #here and bellow are defaults values, not mandatory
            'TIME_DIGITS' : 'MINUTES',
            'ANG_FORMAT' : 'DEG',
            'APPARENT' : 'AIRLESS',
            'RANGE_UNITS' : 'AU',
            'SUPPRESS_RANGE_RATE' : 'NO',
            'SKIP_DAYLT' : 'NO',
            'SOLAR_ELONG' : '0,180',
            'EXTRA_PREC' : 'NO',
            'R_T_S_ONLY' : 'NO'
            }
        self.already_quered = False #Control not to ask for target name every time if query for multiple obs
        self.good_query = False
        self.record = None
        
    def retrieve_param_from_fits(self, fits_dir):
        """
        Open the first fits files in 'fits_dir' starting with TRAP to set the observatory in the ephemeris class
        Open all the fits files in 'fits_dir' to get the dates and set the ephemeris range in the ephemeris class
        
        Parameters:
            fits_dir (str): path to the directory
        """
        
        # should look inside the fits header for image type rather than looking at the file name
        # would allow to be used indepedently of the rename script
        
        # gets useful content from the fits headers
        
        self.fits_dir = fits_dir
        self.fits_filelist = os.listdir(self.fits_dir)
        
        # define observatory from the fits headers
        
        for item in self.fits_filelist:
            if item[0:4] == 'TRAP' and item.endswith('.fits'):
                print(item)
                with fits.open(os.path.join(self.fits_dir, item)) as hdul:
                    print(os.path.join(self.fits_dir, item))
                    self.observatory = hdul[0].header['OBSERVAT']
                break
        if self.observatory == 'Oukaimeden':
            self.parameters['CENTER'] = 'Z53@399'
            obs_print = 'Observatory set to TRAPPIST-North, Oukaimeden (observatory)'
        elif self.observatory == 'TRAPPIST':
            self.parameters['CENTER'] = 'I40@399'
            obs_print = 'Observatory set to La Silla--TRAPPIST (observatory)'
        else:
            obs_print = 'Observatory could not be set (fits header returns' + self.observatory + ')'
        print(obs_print)
        
        # defines the date range from the fits headers
        
        self.fits_datelist=[]
        for item in self.fits_filelist:
            if (item[0:4] == 'TRAP') and (item.endswith('.fits')):
                
                print(os.path.join(self.fits_dir, item))
                with fits.open(os.path.join(self.fits_dir, item)) as hdul:
                    date = hdul[0].header['DATE-OBS']
                    if len(date) == 23:
                        self.fits_datelist.append(datetime.datetime.strptime(date, '%Y-%m-%dT%H:%M:%S.%f'))
                    elif len(date) == 19:
                        self.fits_datelist.append(datetime.datetime.strptime(date, '%Y-%m-%dT%H:%M:%S'))
                    else:
                        print("problem with date", date)
        self.parameters['START_TIME'] = (min(self.fits_datelist) - datetime.timedelta(hours=1)).strftime('%Y-%m-%d %H:%M')
        self.parameters['STOP_TIME'] = (max(self.fits_datelist) + datetime.timedelta(hours=1)).strftime('%Y-%m-%d %H:%M')
    
        print('Ephemeris range set from ' + self.parameters['START_TIME'] + ' to ' + self.parameters['STOP_TIME'])
    
    # sets the object and query
    
    def query_horizons(self):
        """
        Query NASA Horizons and format the result
        """
        # makes the url
        self.query_url = "https://ssd.jpl.nasa.gov/api/horizons.api?format=text"
        for param in self.parameters:
            self.query_url += '&' + param + "='" + str(self.parameters[param]) + "'"
            
        # fetch the data
        self.query_result = requests.get(self.query_url).text.split('\n')
        
    def query_input(self, unique_target=False, target=None, convert_MPC_Horizon = False, last_epoch = True):
        """
        Query for the object name and launch query_horizons()
        Loop until a query is successful, otherwise ask for a new input
        
        Parameters:
            unique_target (str, optional, default=False): set to True to remember the name of the target of the first successful query.
                Useful for reducing different night with the same object.
                Set to False if using nights from different objects.
            target (str, optional, default=None): if different from None, use as initial object name input for query_horizons()
            convert_MPC_Horizon (boolean, optional, default=False): if is True and target different than None, covert target_name from MPC to NASA Horizon format
            last_epoch (boolean, optional, default=True): if True relauch a query with the last epoch if there is multiple epochs returned after a successful query'
        """
        while True:
            if (self.record != None) and (last_epoch == True):
                print('Making a second query with the record #', self.record)
                self.parameters['COMMAND'] = self.record
                self.record = None
            elif (self.good_query == True) and (unique_target == True):
                print('Target already defined as ' + self.obj_fullname)
            elif (self.already_quered == False) and (target != None):
                if convert_MPC_Horizon == True:
                    if (len(target) == 8) and (target[0:2] + target[5] + target[7] == 'CK00'):
                        self.parameters['COMMAND'] = '20' + target[2:4] + ' ' + target[4] + target[6]
                    elif (len(target) == 5) and (target[0:3] + target[-1] == '000P'):
                        self.parameters['COMMAND'] = target[3:]
                    elif (len(target) == 5) and (target[0:2] + target[-1] == '00P'):
                        self.parameters['COMMAND'] = target[2:]
                    elif (len(target) == 5) and (target[0:1] + target[-1] == '0P'):
                        self.parameters['COMMAND'] = target[1:]
                    else:
                        self.parameters['COMMAND'] = target
                else:
                    self.parameters['COMMAND'] = target
            else:
                self.parameters['COMMAND'] = input('Object name: ')
            
            self.query_horizons()
            
            for line in self.query_result:
                print(line)
            
            if self.query_result[-4][:10] == '    Author':
                self.already_quered = True
                self.good_query = True
                break
            elif (last_epoch == True) and ("To SELECT, enter record # (integer), followed by semi-colon.)" in self.query_result[-3]):
                # added an exception for 73P because fragments
                if self.parameters['COMMAND'] == '73P':
                    i = -5
                    while True:
                        if [ elem for elem in self.query_result[i].split(" ") if elem != ''][2] == '73P':
                            self.record = [ elem for elem in self.query_result[i].split(" ") if elem != ''][0]
                            break
                        else:
                            i -=1
                else:
                    last_line = self.query_result[-5]
                    # record = last_line.split(" ")
                    self.record = [ elem for elem in last_line.split(" ") if elem != ''][0]
                #check if the record is the last epoch
                    year_record = [ elem for elem in last_line.split(" ") if elem != ''][1]
                    previous_year = [ elem for elem in self.query_result[-6].split(" ") if elem != ''][1]
                    if [ elem for elem in self.query_result[-5].split(" ") if elem != ''][2] != self.parameters['COMMAND']:
                        print(self.parameters['COMMAND'])
                        print(self.query_result[-5])
                        input("WARNING: last line does not correspond to object?")
                    if int(previous_year) > int(year_record):
                        print('selected record: ', self.record)
                        input('WARNING: check the epoch selected is the las one')
            else:
                self.already_quered = True
                print('previous attempt: ', self.parameters['COMMAND'])
    
    def generate_ephem_files(self, output_dir):
        """
        Generates ephem.brol and eph.dat file in 'output_dir'.
        Dates are expressed in MJD at the moment but can be changed in the script (then also need to be changed in afrhocalcext)
        
        Parameters:
            output_dir (str): path to the output directory
        """
        
        self.output_dir = output_dir
        
        # table the data from the query
        
        dataline = False
        SOE = False # used for not printing the SOE line
        data = []
        
        for line in self.query_result:
            if line == '$$SOE':
                SOE = True
                dataline = True
            elif line == '$$EOE':
                dataline = False
            elif line[:18] == " Date__(UT)__HR:MN":
                column_head = line
            elif line[:17] == "Target body name:":
                self.obj_fullname = line[18:].split('  ')[0]
                
            if dataline == True and SOE == True:
                SOE = False
            elif dataline == True and SOE == False:
                data.append(line.split(','))
        
        
        data_table = pd.DataFrame()
        data_table = data_table.append(data)
        data_table.columns=column_head.split(',')
        
        #formating for the trappist pipeline
        # data_table[' Date_________JDUT'] = data_table[' Date_________JDUT'].astype(float)
        self.ephem_brol = data_table[[' Date_________JDUT', ' R.A._(ICRF)', ' DEC_(ICRF)','                r',
                                  '       rdot', '             delta', '     deldot']]
        self.eph_dat = data_table[[' Date_________JDUT', ' R.A._(ICRF)', ' DEC_(ICRF)','                r', '             delta']]
        
        self.ephem_brol = self.ephem_brol.astype(float)
        self.eph_dat = self.eph_dat.astype(float)
        self.ephem_brol[' Date_________JDUT'] = self.ephem_brol[' Date_________JDUT'] - 2400000.5
        self.eph_dat[' Date_________JDUT'] = self.eph_dat[' Date_________JDUT'] - 2400000.5
        
        # saving the files
        self.ephem_brol.to_csv(os.path.join(self.output_dir, 'ephem.brol'), header=None, index=None, sep=' ', mode='w')
        self.eph_dat.to_csv(os.path.join(self.output_dir, 'eph.dat'), header=None, index=None, sep=' ', mode='w')
        
        print('Image directory: ' + self.fits_dir)
        print('Observatory: ' + self.observatory)
        print('Target: ' + self.obj_fullname)
        # print('Ephemeris range set from ' + parameters['START_TIME'] + ' to ' + parameters['STOP_TIME'])
        print('\n')
        
        print('ephem.brol printed in ' + self.output_dir)
        print('eph.dat printed in ' + self.output_dir)
        
        
# ephemeris = ephemeris()
# ephemeris.retrieve_param_from_fits(fits_dir)
# ephemeris.query_input()
# ephemeris.generate_ephem_files(output_dir)

if __name__ == "__main__":
    import matplotlib.pyplot as plt
    import matplotlib.dates as mdates

    ### PARAMETERS ###
    comets = ['2017 K2', '2022 L1'
                , '2022 E3', '2019 E3','2021 E3'
                # ,'2020 R7'
                 # ,'2021 E3','2020 R7', '2021 G2', '2020 K1', '2021 C5', '2022 A2'
                 # ,'2019 E3', '2020 V2', '2021 P4', '2021 T2'
                    ,'29P','22P','2022 P1', '29P', '73P'
                  #,'408P', '117P', '61P', '81P', '118P', '116P', '100P', '327P', '71P', '73P'
                    # ,'107P', '169P'
              ]
    observatory = 'TN'
    night = '2022-10-26'
    save_path = f'/home/Mathieu/visibility_plot{observatory}.png'
    ##################
    
    dt = datetime.datetime.now()
    def import_eph(target): #function to import and format the ephemeris
        print(target)
        eph.parameters['STEP_SIZE'] = '1 m'
        eph.parameters['COMMAND'] = target
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
            input('debug this')
        data = []
        for line in eph.query_result[eph.query_result.index('$$SOE')+1:eph.query_result.index('$$EOE')]:
            data.append(line[:-1].split(','))
        df = pd.DataFrame(data, columns=('date','dateJD', 'sol_marq', 'lun_marq', 'azi', 'elev', 'airmass', 'magex', 'lun_ang', 'lun_ill'))
        df['date'] =  pd.to_datetime(df['date'], format=' %Y-%b-%d %H:%M')
        df = df.astype({'elev':'float'})
        df.reset_index(drop=True, inplace=True)
        return df
    
    #initiate ephemeris
    eph = ephemeris()
    eph.parameters['QUANTITIES'] = '4,8,25'
    if observatory == 'TN':
        eph.parameters['CENTER'] = 'Z53@399'
    elif observatory == 'TS':
        eph.parameters['CENTER'] = 'I40@399'
    # datetime.datetime.strptime(date, '%Y-%m-%dT%H:%M:%S.%f')
    # - datetime.timedelta(hours=1)).strftime('%Y-%m-%d %H:%M')
    eph.parameters['START_TIME'] = datetime.datetime.strptime(night + 'T12:00', '%Y-%m-%dT%H:%M')
    eph.parameters['STOP_TIME'] = datetime.datetime.strptime(night + 'T12:00', '%Y-%m-%dT%H:%M') + datetime.timedelta(days=1)
    
    #initiate the plot
    fig = plt.figure(figsize=(12,8))
    ax = fig.gca()
    
    # for the moon
     
    df_moon = import_eph('301')
    
    naut_start = df_moon.loc[df_moon['sol_marq'] == "A"][:1].index[0]
    naut_end = df_moon.loc[df_moon['sol_marq'] == "A"][-1:].index[0]
    df_moon = df_moon.loc[naut_start:naut_end]
    
    if len(df_moon.loc[df_moon['elev']>0]) > 0:
        # print(df_moon.loc[df_moon['elev']>0].to_string())
        plt.plot(df_moon['date'], df_moon['elev'], ls='dashed', color='black')
        ax.text(df_moon.loc[df_moon['elev'].idxmax(), 'date'],
                df_moon.loc[df_moon['elev'].idxmax(), 'elev'] + 1,
                'Moon', ha='center')
    
    #get and plot the comets
    for comet in comets:
        
        df = import_eph(comet)
        df = df.loc[naut_start:naut_end]
        
        if len(df.loc[df['elev']>0]) > 0:
            plt.plot(df['date'], df['elev'])
            lunang_start = int(float(df.loc[naut_start,'lun_ang']))
            lunang_end = int(float(df.loc[naut_end,'lun_ang']))
            ax.text(df.loc[df['elev'].idxmax(), 'date'],
                    df.loc[df['elev'].idxmax(), 'elev'] + 1,
                    f"{comet}\n({str(lunang_start)}-{str(lunang_end)}°)", ha='center')
    
    
    
    # #pretty plot
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
    plt.xlim(df_moon.loc[df_moon['sol_marq'] == "A"][:1]['date'],
              df_moon.loc[df_moon['sol_marq'] == "A"][-1:]['date'])
    plt.ylim(0,90)
    plt.axvline(x=df_moon.loc[df_moon['sol_marq'] == " "][:1]['date'],
                color='black', ls='dotted', alpha=0.5)
    plt.axvline(x=df_moon.loc[df_moon['sol_marq'] == " "][-1:]['date'],
                color='black', ls='dotted', alpha=0.5)
    plt.axhline(y=20, color='black', ls='dotted')
    moon_illu_min = int(float(df_moon.loc[naut_start, 'lun_ill']))
    moon_illu_max = int(float(df_moon.loc[naut_end, 'lun_ill']))
    title = f'{night}\nMoon illumination {str(moon_illu_min)}-{str(moon_illu_max)}%'
    xmiddle = df_moon.loc[df_moon['sol_marq'] == "A"][:1]['date']+ (df_moon.loc[df_moon['sol_marq'] == "A"][-1:]['date'].values[0] - df_moon.loc[df_moon['sol_marq'] == "A"][:1]['date'].values[0])/2
    ax.text(xmiddle, 80, title, ha='center', backgroundcolor='white', fontsize='large')
    ax.set_xlabel('Time UT')
    ax.set_ylabel('Elevation (°)')
    # airmass from “Revised Optical Air Mass Tables and Approximation Formula”, Kasten F., Young A., Applied Optics, vol 28, no. 22, p. 4735-4738, Nov. 15, 1989.
    ax2 = ax.twinx()
    ax2.set_ylim(0,90)
    plt.yticks(ticks=[10,20,30,40,50,60], labels=[5.6, 2.9, 2.0, 1.6, 1.6, 1.1])
    ax2.set_ylabel('Airmass')
    # ax2.set_xticks([10,20,30,40,50])
    # ax2.set_ylabels([38, 5.6, 2.9, 2, 1.6, 1.6, 1.1])
    print('Executed in', datetime.datetime.now() - dt)
    plt.tight_layout()
    # plt.show()
    plt.savefig(save_path, bbox_inches='tight')

    
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
    
        
