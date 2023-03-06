#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 22 11:41:28 2023

@author: Mathieu
"""

from prose import FITSImage, Sequence, blocks
from prose import FitsManager
from prose import Image
import pandas as pd
import os
import warnings
warnings.filterwarnings('ignore', category=UserWarning, append=True)
pd.set_option('display.max_columns', None)

targets = ['CK22E030','hkj']
startdate = "2023-02-11"
enddate = "2023-02-20"

fm = FitsManager('/home/Mathieu/Documents/TRAPPIST/prosetest_22E3_20230216/20230216',
                 depth = 10,
                 file = '/home/Mathieu/Documents/TRAPPIST/prosetest_22E3_20230216/fitsmanager.sqlite',
                 batch_size = False)

def import_fmcalib():
    fmcalib = FitsManager('/home/Mathieu/Documents/TRAPPIST/reduced_data/masters_calib',
                          depth = 10,
                          file = '/home/Mathieu/Documents/TRAPPIST/prosetest_22E3_20230216/fm_calib.sqlite',
                          batch_size = False)
    return fmcalib 

for target in targets:
    try:
        science = pd.concat([science, fm.files(type = 'light', target = target, path=True)], ignore_index=True)
    except:
        science = fm.files(type = 'light', target = target, path=True)
science.sort_values(by='date', axis=0, ascending=True, inplace=True)
# science.sort_values(by='target', axis=0, ascending=True, inplace=True)
science = science.loc[(science.date > startdate) & (science.date < enddate)]
# print(science)

def custom_sql(date, dayinterval, telescope, width):
    return f"AND telescope LIKE '{telescope}' \
    AND width LIKE {width} \
        AND date >= date('{date}', '-{dayinterval:.0f} days') \
            AND date <= date('{date}', '+{dayinterval:.0f} days')"

def fetch_dark_files(fm, fmcalib, date, exposure, telescope, width):
    tolerance = 0
    nodark = True
    while nodark:
        master_dark_path = fmcalib.files(type = 'dark', exposure = exposure,
                                         tolerance=tolerance, path=True,
                     date = date, width = width, telescope = telescope).path.tolist()
        
        dayinterval = 7
        if len(master_dark_path) > 0:
            dark_path = master_dark_path[0]
            darks = FITSImage(dark_path).data
            print('found master dark')
            meantime = exposure
            dark_name = dark_path.split('/')[-1]
            save = False
            nodark = False
        
        else:
            print('create dark list')
            dark_name = None
            # create the dark list
            darks = []
            # tolerance = 0
            while len(darks) < 5:
                darks_tb = fm.files(type = 'dark', exposure = exposure, tolerance=tolerance, path=True,
                             custom=custom_sql(date, dayinterval, telescope, width))
                meantime = darks_tb['exposure'].mean()
                darks = darks_tb.path.tolist()
                nodark = False
                if dayinterval < 21:
                    dayinterval += 7
                else:
                    tolerance += exposure*0.2
                    print(f"Could not find {exposure}s darks in {dayinterval}. \
                          Try with tolerance of {tolerance}")
                    nodark = True
                    break
            
            save = True
            
    if dark_name == None:
        dark_name = f'mdark_{date}_{str(int(meantime))}.fits'
    return darks,save,meantime,dark_name

def fetch_darksflats_files(fm, fmcalib, date, telescope, width):
    tolerance = 0
    exposure = 15
    
    master_dark_path = fmcalib.files(type = 'dark', exposure = exposure,
                                      tolerance=tolerance, path=True,
                  date = date, width = width, telescope = telescope).path.tolist()
    
    dayinterval = 5
    if len(master_dark_path) > 0:
        dark_path = master_dark_path[0]
        darks = FITSImage(dark_path).data
        print('found master flatsdark')
        dark_name = dark_path.split('/')[-1]
        save = False
    
    else:
        print('create flatsdarks list')
        # create the dark list
        darks = []
        # tolerance = 0
        while len(darks) < 5:
            darks = fm.files(type = 'dark', exposure = exposure, tolerance=0, path=True,
                          custom=custom_sql(date, dayinterval, telescope, width)).path.tolist()
            if dayinterval > 21:
                print(f"Looking for darksflats older than {dayinterval} days.")
            dayinterval += 5
            
        dark_name = f'mdark_{date}_{str(int(exposure))}.fits'
        save = True
        
        
    return darks,save,dark_name

def fetch_flat_files(fm, fmcalib, date, filt, telescope, width):

    dayinterval = 7
        
    master_flat_path = fmcalib.files(type = 'flat', filter=filt, path=True,
                 date = date, width = width, telescope = telescope).path.tolist()
    
    if len(master_flat_path) > 0:
        flat_path = master_flat_path[0]
        flats = FITSImage(flat_path).data
        print('found master flat')
        flat_name = flat_path.split('/')[-1]
        save = False
    
    else:
        print('create flat list')
        flats = []
        while len(flats) < 5:
            flats = fm.files(type = 'flat', filter=filt, path=True,
                         custom=custom_sql(date, dayinterval, telescope, width)).path.tolist()
            if dayinterval > 21:
                print(f"Looking for {filt} flats older than {dayinterval} days.")
            dayinterval += 7

        flat_name = f'mflat_{date}_{filt}.fits'
        save = True            
            
    return flats,save,flat_name

def fetch_bias_files(fm, fmcalib, date, telescope, width):
    
    master_bias_path = fmcalib.files(type = 'bias', path=True,
                  date = date, width = width, telescope = telescope).path.tolist()
    
    dayinterval = 5
    if len(master_bias_path) > 0:
        bias_path = master_bias_path[0]
        bias = FITSImage(bias_path).data
        print('found master bias')
        bias_name = bias_path.split('/')[-1]
        save = False
    
    else:
        print('create bias list')
        # create the dark list
        bias = []
        # tolerance = 0
        while len(bias) < 5:
            bias = fm.files(type = 'bias', path=True,
                          custom=custom_sql(date, dayinterval, telescope, width)).path.tolist()
            if dayinterval > 21:
                print(f"Looking for darksflats older than {dayinterval} days.")
            dayinterval += 5
            
        bias_name = f'mbias_{date}.fits'
        save = True
    
    return bias,save,bias_name

def create_calib_list(fm, date, filt, exposure, telescope, width):
    
    
  
    # create the dark list
    darks = []
    tolerance = 0
    dayinterval = 5
    while len(darks) < 5:
        darks = fm.files(type = 'dark', exposure = exposure, tolerance=tolerance, path=True,
                     custom=custom_sql(date, dayinterval, telescope, width)).path.tolist()
        if dayinterval < 21:
            dayinterval += 5
        else:
            dayinterval = 5
            tolerance += exposure*0.2
            print(f"Could not find {exposure}s darks in {dayinterval}. \
                  Try with tolerance of {tolerance}")
        
    darksflats = []
    dayinterval = 5
    while len(darksflats) < 5:
        darksflats = fm.files(type = 'dark', exposure = 15, tolerance=0, path=True,
                     custom=custom_sql(date, dayinterval, telescope, width)).path.tolist()
        if dayinterval < 21:
            dayinterval += 5
        else:
            print(f"Looking for 15 darks flat older than {dayinterval} days.")
            dayinterval += 5   
    
    flats = []
    dayinterval = 5
    while len(flats) < 5:
        flats = fm.files(type = 'flat', filter=filt, path=True,
                     custom=custom_sql(date, dayinterval, telescope, width)).path.tolist()
        if dayinterval < 21:
            dayinterval += 5
        else:
            print(filt)
            print(f"Looking for flats older than {dayinterval} days.")
            dayinterval += 5
    
    bias = []
    dayinterval = 5
    while len(bias) < 5:
        bias = fm.files(type = 'bias', path=True,
                     custom=custom_sql(date, dayinterval, telescope, width)).path.tolist()
        if dayinterval < 20:
            dayinterval += 5
        else:
            print(filt)
            input(f"Looking for flats older than {dayinterval} days.")
            dayinterval += 5
    
    return darks,darksflats,flats,bias

for index, r in science.iterrows():
    im_path = r['path']
    image = FITSImage(im_path)
    print(f'Calibrating {image.label}')
    # darks,darksflats,flats,bias= create_calib_list(fm,
    #                                             r['date'], r['filter'], r['exposure'],
    #                                                 r['telescope'], r['width'])
    fmcalib = import_fmcalib()
    darks, save_dark, dark_meantime, mdark_name = fetch_dark_files(fm, fmcalib, r['date'],
                                                                   r['exposure'],
                                 r['telescope'], r['width'])
    flats,save_flat,mflat_name = fetch_flat_files(fm, fmcalib, r['date'], r['filter'],
                                 r['telescope'], r['width'])
    
    if isinstance(flats, list) and dark_meantime != 15:
        darksflats, save_darksflats, mdarksflats_name = fetch_darksflats_files(fm, fmcalib,
                                    r['date'], r['telescope'], r['width'])
    else:
        darksflats, save_darksflats, mdarksflats_name = None,False,None
    
    bias, save_bias, mbias_name = fetch_bias_files(fm, fmcalib,
                                r['date'], r['telescope'], r['width'])
    
    # image.show()
    calibration = Sequence([
        blocks.Calibration(darks=darks, darksflats=darksflats, bias=bias),
        # blocks.Trim(),
        # blocks.PointSourceDetection(n=20), # stars detection
        # blocks.Cutouts(21),                # stars cutouts
        # blocks.MedianEPSF(),               # building EPSF
        # blocks.psf.Moffat2D(),             # modeling EPSF
    ])

    calibration.run(image)
    
    image.fits_header['MBIAS'] = mbias_name
    image.fits_header['MFLAT'] = mflat_name
    image.fits_header['MDARK'] = mdark_name
    image.writeto(os.path.join('/home/Mathieu/Documents/TRAPPIST/prosetest_22E3_20230216',
                  image.fits_header['ORIGFILE']))
    
    if save_dark:
        master_data = calibration[type=='Calibration'].master_dark
        master_ref = FITSImage(darks[0])
        master_ref.data = master_data
        master_ref.fits_header['DATE-OBS'] = image.fits_header['DATE-OBS']
        master_ref.fits_header['IMAGETYP'] = 'MASTER_DARK'
        master_ref.fits_header['EXPTIME'] = dark_meantime
        master_ref.fits_header['EXPOSURE'] = dark_meantime
        master_ref.fits_header['MBIAS'] = mbias_name
        master_ref.writeto(os.path.join('/home/Mathieu/Documents/TRAPPIST/reduced_data/masters_calib/',mdark_name))
    
    if save_darksflats:
        master_data = calibration[type=='Calibration'].master_darksflats
        master_ref = FITSImage(darksflats[0])
        master_ref.data = master_data
        master_ref.fits_header['DATE-OBS'] = image.fits_header['DATE-OBS']
        master_ref.fits_header['IMAGETYP'] = 'MASTER_DARK'
        master_ref.fits_header['MBIAS'] = mbias_name
        master_ref.writeto(os.path.join('/home/Mathieu/Documents/TRAPPIST/reduced_data/masters_calib/',mdarksflats_name))
        
    if save_flat:
        master_data = calibration[type=='Calibration'].master_flat
        master_ref = FITSImage(flats[0])
        master_ref.data = master_data
        master_ref.fits_header['DATE-OBS'] = image.fits_header['DATE-OBS']
        master_ref.fits_header['IMAGETYP'] = 'MASTER_FLAT'
        master_ref.fits_header['MDARK'] = mdarksflats_name
        master_ref.fits_header['MBIAS'] = mbias_name
        master_ref.writeto(os.path.join('/home/Mathieu/Documents/TRAPPIST/reduced_data/masters_calib/',mflat_name))
        
    if save_bias:
        master_data = calibration[type=='Calibration'].master_bias
        master_ref = FITSImage(bias[0])
        master_ref.data = master_data
        master_ref.fits_header['DATE-OBS'] = image.fits_header['DATE-OBS']
        master_ref.fits_header['IMAGETYP'] = 'MASTER_BIAS'
        master_ref.writeto(os.path.join('/home/Mathieu/Documents/TRAPPIST/reduced_data/masters_calib/',mbias_name))
        
        # master_dark = Image(calibration[0].master_dark)
        # master_dark.metadata = FITSImage(darks[0]).metadata
        # master_dark.metadata['date'] = image.metadata['date']
        # master_dark.metadata['jd'] = image.metadata['jd']
        # master_dark.metadata['path'] = None
        # master_dark.metadata['wcs'] = None
        # master_dark.files = darks
        # master_dark.save(filepath='/home/Mathieu/Documents/TRAPPIST/reduced_data/masters_calib/master_dark.improse',
                         # low_data=False, image_dtype='float64')

    
    # input()
    # image.show()

