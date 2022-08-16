#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  8 13:29:26 2021

@author: Mathieu Vander Donckt
mathieu.vanderdonckt@uliege.be
"""

import os
import pandas as pd
import matplotlib.pyplot as plt  
plt.rcParams.update({'font.size': 22})
# from mpl_toolkits.mplot3d import axes3d
import numpy as np
from astropy.visualization import astropy_mpl_style
from astropy.io import fits
from matplotlib import colors
from matplotlib.patches import Circle
import numpy as np
import math 
from matplotlib.ticker import MultipleLocator

def list_radprof(input_dir, output_dir):
    radtable = pd.DataFrame()
    for path, subdirs, files in os.walk(input_dir):
        for file in files:
            if 'rad_' in file and ".png" not in file:
                radfile = os.path.join(path, file)
                # print(file)
                df = pd.read_csv(radfile, header=None, sep="\s+")
                if 'TN' in radfile:
                    df['observatory'] = 'TN'
                elif 'TS' in radfile:
                    df['observatory'] = 'TS'
                else:
                    input("problem with setting observatory")
                df['obsname'] = radfile.split('/')[-3] 
                radtable = radtable.append(df)
    radtable.sort_values(by=['date']).to_csv(export_path, index=False)

# =============================================================================
# def compil_haser(reduced_path):
#     hasertab = pd.DataFrame()
#     for path, subdirs, files in os.walk(reduced_path):
#         for file in files:
#             if file == 'outputhaser-BC':
#                 hasfile = os.path.join(path, file)
#                 # print(file)
#                 df = pd.read_csv(hasfile, header=None, sep="\s+")
#                 if 'TN' in hasfile:
#                     df['observatory'] = 'TN'
#                 elif 'TS' in hasfile:
#                     df['observatory'] = 'TS'
#                 else:
#                     input("problem with setting observatory")
#                 df['obsname'] = hasfile.split('/')[-3] 
#                 hasertab = hasertab.append(df)
#                 
#     filtlist = hasertab[15-1].drop_duplicates().values.tolist()
#     
#     for filt in filtlist:
#         
#         filteredtable = hasertab.loc[hasertab[15-1] == filt]
#         
#         fig2D = plt.figure(figsize=(12,9))
#         ax2D = fig2D.add_subplot()
#         ax2D.scatter(filteredtable[5-1], filteredtable[12-1])
#         # ax2D.axvline(perihelion)
#         fig2D.suptitle("haser_" + filt)
#         ax2D.set_ylabel('Q (s-1)')
#         ax2D.set_xlabel('MJD to perihelion')
#                 
#     # plt.plot(hasertab.loc[hasertab[14] == 'C2'].sort_values(by=4)[4], hasertab.loc[hasertab[14] == 'C2'].sort_values(by=4)[11])
#     # print(hasertab)
#                 
# 
# # compil_haser('/home/Mathieu/Documents/TRAPPIST/reduced_data/2020T2')
# =============================================================================
            
    
    
# =============================================================================
#     db_path = os.path.join(reduced_path, 'results.db')
#     if os.path.isfile(db_path) == False:
#         print('generating results.db')
#         results = pd.DataFrame(columns=('path','image','filter', 'magnitude', ))
#         results.to_csv(db_path, index=False)
#         # with open(db_path, 'a')
# =============================================================================
# generating_plotfiles('/home/Mathieu/Documents/TRAPPIST/reduced_data/2020T2')

def plot_radprof(input_dir, saveplot=''):
    radtable = pd.DataFrame()
    for path, subdirs, files in os.walk(input_dir):
        for file in files:
            if 'rad_' in file and ".png" not in file:
                radfile = os.path.join(path, file)
                # print(file)
                df = pd.read_csv(radfile, header=None, sep="\s+")
                if 'TN' in radfile:
                    df['observatory'] = 'TN'
                elif 'TS' in radfile:
                    df['observatory'] = 'TS'
                else:
                    input("problem with setting observatory")
                df['obsname'] = radfile.split('/')[-3] 
                radtable = radtable.append(df)
    # rad_columns = [
    #             "1. Name of the image",
    #             "2. Distance rx from the comet's nucleus in pixels",
    #             "3. Number of pixels forming the circle of radius rx around the nucleus",
    #             "4. Median flux at a distance rx (in ADU s −1 )",
    #             "5. Distance r from the nucleus (in arcsec)",
    #             "6. Total number of pixels in a disk of radius rx",
    #             "7. Total flux in a disk of radius rx (in ADU s −1 )",
    #             "8. Median flux at a distance r from the comet center (in ADU s −1 arcsec −2 )",
    #             "9. Median magnitude per arcsec 2 at a distance r from the comet center",
    #             "10. Flux at a distance r from the nucleus per unit wavelenght (in erg cm −2 s −1 Å arcsec −2 )",
    #             "11. Flux at a distance r from the nucleus in the entire filter's band (in erg cm −2 s −1 arcsec −2 )",
    #             "12. Total magnitude of a disk of radius r",
    #             "13. Integrated flux in a disk of radius r (in erg cm −2 s −1 Å )",
    #             "14. Integrated flux in a disk of radius r in the filter's band (in erg cm −2 s −1 )",
    #             "15. Name of the filter",
    #             "16. Time in Julian Days",
    #             "17. Heliocentric distance (in AU)",
    #             "18. Geocentric distance (in AU)"
    #             ]
    #               + observatory
    
    filtlist = radtable[15-1].drop_duplicates().values.tolist()
    
    for filt in filtlist:        
        filteredtable = radtable.loc[radtable[15-1] == filt]
        # print(filteredtable)
        MJDlist = filteredtable[16-1].drop_duplicates().values.tolist()
        Xlist = []
        Ylist = []
        Zlist = []
        Zlabel = []
        colorlist = []
        obsname = []
        filenamelist = []
        for index, MJD in enumerate(sorted(MJDlist)):
            # print(index, JD)
            MJDfilttab = filteredtable.loc[filteredtable[16-1] == MJD]
            Xlist.append(MJDfilttab[2-1][:40])
            Ylist.append(MJDfilttab[4-1][:40])
            # Zlist.append(np.full(len(Xlist[0]), index))
            # Zlist.append(np.full(len(Xlist[0]), MJD-perihelion))
            Zlist.append(np.full(len(Xlist[0]), MJD))
            Zlabel.append(MJDfilttab['obsname'][0][:-2])
            obsname.append(MJDfilttab['obsname'][0])
            filenamelist.append(MJDfilttab[0][0])
            if MJDfilttab['obsname'][0][-2:] == 'TS':
                colorlist.append('r')
            elif MJDfilttab['obsname'][0][-2:] == 'TN':
                colorlist.append('b')
            else:
                colorlist.append('g')
        
        fig3D = plt.figure(figsize=(10,10))
        ax3D = fig3D.add_subplot(111, projection='3d')
        for index, line in enumerate(Xlist):
            ax3D.plot(Xlist[index], Zlist[index], Ylist[index], color=colorlist[index])
            ax3D.set_ylabel('MJD to perihelion')
            ax3D.set_zlabel('Median flux (ADU/s)')
            ax3D.set_xlabel('Distance from the nucleus (pixel)')
            # ax.set_yticklabels(Zlabel)
            # ax3D.set_title(filt)
            # ax.view_init(elev=40., azim=-70)
            fig2D = plt.figure(figsize=(12,9))
            ax2D = fig2D.add_subplot()
            ax2D.plot(Xlist[index], Ylist[index])
            # fig2D.tight_layout()
            fig2D.suptitle(radfile.split('/')[-4] + "_rad_" + filt + "_"+ obsname[index])
            ax2D.set_title(filenamelist[index] + ' (MJD ' + str(int(Zlist[index][0])) + ')')
            ax2D.set_ylabel('Median flux (ADU/s)')
            ax2D.set_xlabel('Distance from the nucleus (pixel)')
            
            # fig3D.tight_layout()
            fig3D.suptitle(radfile.split('/')[-4] + "_rad3D_" + filt)
            
            if saveplot != '':
                if not os.path.exists(saveplot):
                    os.makedirs(saveplot)
                fig2D.savefig(os.path.join(saveplot, radfile.split('/')[-4] + "_rad_" + filt + "_"+ obsname[index] + ".png"))
                fig3D.savefig(os.path.join(saveplot, radfile.split('/')[-4] + "_rad3D_" + filt + ".png"))
        
def plot_mag(input_dir, saveplot=''):
    colordict = {
        'R':'tab:red',
        'RC':'tab:orange',
        'V':'tab:green',
        'B':'tab:blue',
        'BC':'tab:purple',
        'I':'tab:olive'
        }
    radtable = pd.DataFrame()
    for path, subdirs, files in os.walk(input_dir):
        for file in files:
            if 'rad_' in file and ".png" not in file:
                radfile = os.path.join(path, file)
                # print(file)
                df = pd.read_csv(radfile, header=None, sep="\s+")
                if 'TN' in radfile:
                    df['observatory'] = 'TN'
                elif 'TS' in radfile:
                    df['observatory'] = 'TS'
                else:
                    input("problem with setting observatory")
                df['obsname'] = radfile.split('/')[-3] 
                df['path'] = path
                radtable = radtable.append(df)
    
    filtlist = radtable[15-1].drop_duplicates().values.tolist()
    
    figtot = plt.figure(figsize=(12,9))
    figtot.suptitle(path.split('/')[-3] + "_mag_tot")
    axtot = figtot.add_subplot()  
    
    filtlist = (filt for filt in filtlist if filt in ['B','V','R','I','BC','RC','GC','UC'])
    for filt in filtlist:        
        filteredtable = radtable.loc[radtable[15-1] == filt]
        # print(filteredtable)
        MJDlist = filteredtable[16-1].drop_duplicates().values.tolist()
        Xlist = []
        Ylist = []
        uperrorlist = []
        downerrorlist = []
        colorlist = []
        r_aperture = 5
        for index, MJD in enumerate(sorted(MJDlist)):
            MJDfilttab = filteredtable.loc[filteredtable[16-1] == MJD]
            # print(filt,index, MJD, MJDfilttab.loc[MJDfilttab[5-1] >=5][:1].values[0][0])
            up_row = MJDfilttab.loc[MJDfilttab[5-1] >=r_aperture][:1]
            up_radius = up_row.values[0][5-1] 
            up_mag = up_row.values[0][12-1]
            up_error_file = os.path.join(up_row.values[0][21-1], 'radplus_' + up_row.values[0][1-1] + '.txt')
            up_error_row = pd.read_csv(up_error_file, header=None, sep="\s+").iloc[up_row.index]
            up_error = up_error_row.values[0][12-1] - up_mag
            
            down_row = MJDfilttab.loc[(MJDfilttab[5-1] < r_aperture) & (MJDfilttab[2-1] < 20)][-1:]
            down_radius = down_row.values[0][5-1]
            down_mag = down_row.values[0][12-1]
            down_error_file = os.path.join(down_row.values[0][21-1], 'radmoins_' + down_row.values[0][1-1] + '.txt')
            down_error_row = pd.read_csv(down_error_file, header=None, sep="\s+").iloc[down_row.index]
            down_error = down_mag - down_error_row.values[0][12-1]
            # print(upr, upmag)
            # print(downr, downmag)
            # linear interpolation y = y1 + ((x – x1) / (x2 – x1)) * (y2 – y1)
            mag5arcsec = down_mag + ((r_aperture - down_radius) / (up_radius - down_radius)) * (up_mag - down_mag)
            # error only considering error on magnitude (not aperture, and not for the )
            mag5_error = ( (-(r_aperture-down_radius)/(up_radius-down_radius)*down_error)**2 + ((r_aperture-down_radius)/(up_radius-down_radius)*up_error)**2 )**(1/2)
            
            Xlist.append(MJD-perihelion)
            Ylist.append(mag5arcsec)
            uperrorlist.append(mag5_error)
            downerrorlist.append(mag5_error)
            # print(type(mag5_error))
            # print(len(Ylist), len(uperrorlist))
            # print(mag5arcsec, mag5_error)
            if MJDfilttab[:1].values[0][18] == 'TS':
                colorlist.append('r')
            elif MJDfilttab[:1].values[0][18] == 'TN':
                colorlist.append('b')
            else:
                colorlist.append('g')
        
        model1 = np.poly1d(np.polyfit(Xlist, Ylist, 2))
        polyline = np.linspace(min(Xlist), max(Xlist), 50)
        axtot.plot(polyline, model1(polyline), color=colordict.get(filt))
        
        fig2D = plt.figure(figsize=(12,9))
        ax2D = fig2D.add_subplot()
        ax2D.plot(polyline, model1(polyline), color=colordict.get(filt))
        ax2D.scatter(Xlist, Ylist,color=colorlist, zorder=2)
        ax2D.errorbar(Xlist, Ylist, yerr=[uperrorlist, downerrorlist], fmt='.', color='black', capsize=4, zorder=1)
        ax2D.axvline(perihelion-perihelion)
        fig2D.suptitle(path.split('/')[-3] + "_mag_" + filt)
        ax2D.set_ylabel('mag (aperture radius = ' + str(r_aperture) + "'')")
        ax2D.set_xlabel('MJD to perihelion')
        
        axtot.scatter(Xlist, Ylist, zorder=2, label=filt, c=colordict.get(filt))
        axtot.errorbar(Xlist, Ylist, yerr=[uperrorlist, downerrorlist],fmt='.', capsize=4, zorder=1, c=colordict.get(filt))
        axtot.axvline(perihelion-perihelion)
        axtot.legend()
        axtot.set_ylabel('mag (aperture radius = ' + str(r_aperture) + "'')")
        axtot.set_xlabel('MJD to perihelion')
        
        if saveplot != '':
            if not os.path.exists(saveplot):
                os.makedirs(saveplot)
            # fig2D.tight_layout()
            fig2D.savefig(os.path.join(saveplot, path.split('/')[-3] + "_mag_" + filt + ".png"))
    
    axtot.legend()
    if saveplot != '':
        if not os.path.exists(saveplot):
            os.makedirs(saveplot)
        # fig2D.tight_layout()
        figtot.savefig(os.path.join(saveplot, path.split('/')[-3] + "_mag_tot.png"))
        

def plot_centering(input_dir, output_dir=None):
    """
    Create a png for each image with the centering given by the pipeline at a radius of 5'' and 10 000 km.

    Parameters
    ----------
    input_dir : str
        path of the folder containing the image and centerlist (typicaly tmpout)
    output_dir : str, optional
        path for outputting the images. If None, uses input_dir. The default is None.

    Returns
    -------
    None.

    """
    
    plt.style.use(astropy_mpl_style)
    for path, subdirs, files in os.walk(input_dir):
        for file in files:
            if 'centerlist' in file:
                centerfile = os.path.join(path, file)
                tab = pd.read_csv(centerfile, header=None, sep="\s+")
                for index, row in tab.iterrows():
                    fitsname = row[0]
                    if os.path.isfile(os.path.join(path, fitsname)):
                        fitspath = os.path.join(path, fitsname)
                    elif os.path.isfile(os.path.join(path[:-9],'images', fitsname)):
                        print('found in the reduced  images folder?')
                        fitspath = os.path.join(path[:-9],'images', fitsname)
                    else:
                        print('error finding fits file')
                    xcent = row[2]
                    ycent = row[3]
                    filt = row[5]
                    delta = row[8] #distance to the earth
                    pixsize = row[9]
                    ctnmethod = row[10]
                    image_data = fits.getdata(fitspath, ext=0)
                    fig = plt.figure()
                    ax = fig.gca()
                    # print(fitspath)
                    pixelcropping = 350 #remove border
                    # plt.colorbar()
                    # plt.axvline(x=xcent-pixelcropping,color='red', alpha=0.2)
                    # plt.axhline(y=ycent-pixelcropping,color='red', alpha = 0.2)
                    plt.axis('off')
                    # plt.scatter(xcent,ycent,color='r')
                    plt.suptitle(filt + ' ' + str(xcent) + ' ' + str(ycent) + ' (' + ctnmethod + ')')
                    ax.set_title(fitsname)
                    circ5arcsec = Circle((xcent-pixelcropping,ycent-pixelcropping),radius = 5/pixsize, alpha=0.5, fill=False, color='red', label='5arcsec')
                    arcsec10k = 206265*10000/(delta*1.5*100000000)
                    circ10k = Circle((xcent-pixelcropping,ycent-pixelcropping),radius = arcsec10k/pixsize, alpha=0.5, fill=False, color='blue', label='10k km')
                    ax.add_patch(circ5arcsec)
                    ax.add_patch(circ10k)
                    output_dir = path if output_dir is None else output_dir
                    ax.legend(loc='upper right')
                    plt.tight_layout()
                    # print(image_data[int(xcent), int(ycent)])
                    # print(np.median(image_data))
                    
                    # if raw_dir != None:
                    #     raw_fitspath = os.path.join(raw_dir, fitsname)
                    #     raw_image_data = fits.getdata(raw_fitspath, ext=0)
                        
                    
                    try:
                        plt.imshow(image_data[pixelcropping:-pixelcropping,pixelcropping:-pixelcropping], cmap='pink', norm=colors.LogNorm(vmin=np.median(image_data), vmax=image_data[int(xcent), int(ycent)]*2))
                    except:
                        plt.close()
                        print('ERROR logscale ', fitspath)
                        fig = plt.figure()
                        ax = fig.gca()
                        circ5arcsec = Circle((xcent-pixelcropping,ycent-pixelcropping),radius = 5/pixsize, alpha=0.5, fill=False, color='red')
                        arcsec10k = 206265*10000/(delta*1.5*100000000)
                        circ10k = Circle((xcent-pixelcropping,ycent-pixelcropping),radius = arcsec10k/pixsize, alpha=0.5, fill=False, color='blue')
                        ax.add_patch(circ5arcsec)
                        ax.add_patch(circ10k)
                        plt.suptitle(filt + ' ' + str(xcent) + ' ' + str(ycent) + ' (' + ctnmethod + ') error logscale')
                        ax.set_title(fitsname)
                        plt.imshow(image_data[pixelcropping:-pixelcropping,pixelcropping:-pixelcropping], cmap='binary')
                        
                    plt.savefig(os.path.join(output_dir, fitsname[:-5] + '_centering.png'))
                    plt.show()
                    # plt.savefig(os.path.join(path, fitsname[:-5] + '.png'))
                    # plt.close()
                    
# plot_centering('/home/Mathieu/Documents/TRAPPIST/tmpout')

def plot_centering_profile(input_dir, output_dir=None, solocomet=False, comet_name=''):
    """
    Same as plot_centering + a plot of the radial profile
    Create a png for each image with the centering given by the pipeline at a radius of 5'' and 10 000 km.

    Parameters
    ----------
    input_dir : str
        path of the folder containing the image and centerlist (typicaly tmpout)
    output_dir : str, optional
        path for outputting the images. If None, uses input_dir. The default is None.
    solocomet: boolean, optional
        if solocomet = True, only the last reduced image (the last in the centerlist) will be replotted
    comet_name : str, optional
        name of the comet for the plot tile. default is ''  
        
    Returns
    -------
    None.

    """
    
    plt.style.use(astropy_mpl_style)
    for path, subdirs, files in os.walk(input_dir):
        for file in files:
            if 'centerlist' in file:
                centerfile = os.path.join(path, file)
                tab = pd.read_csv(centerfile, header=None, sep="\s+")
                if solocomet == True:
                    tab = tab.iloc[[-1]]
                for index, row in tab.iterrows():
                    fitsname = row[0]
                    if os.path.isfile(os.path.join(path, fitsname)):
                        fitspath = os.path.join(path, fitsname)
                        radpath = os.path.join(path, 'rad_' + fitsname + '.txt')
                    elif os.path.isfile(os.path.join(path[:-9],'images', fitsname)):
                        # print('found in the reduced images folder')
                        fitspath = os.path.join(path[:-9],'images', fitsname)
                        radpath = os.path.join(path[:-9], 'profiles', 'rad_' + fitsname + '.txt')
                    else:
                        print('error finding fits file')
                        print(fitsname)
                        input('Acknowledge press enter to continue')
                    
                    
                    save_dir = path if output_dir is None else output_dir
                    # print(path, save_dir)
                    
                    filt = row[5]
                    ctnmethod = row[10]
                    
                    # preparing the centering image
                    xcent = row[2]
                    ycent = row[3]
                    delta = row[8] #distance to the earth
                    pixsize = row[9]

                    centerlimit = 200
                    # pixelcropping = 350 #remove border of the image
                    
                    image_data = fits.getdata(fitspath, ext=0)
                    
                    ximcent = len(image_data[1])/2
                    yimcent = len(image_data)/2
                    # circ5arcsec = Circle((xcent-pixelcropping,ycent-pixelcropping),radius = 5/pixsize, alpha=0.5, fill=False, color='red', label='5arcsec')
                    # arcsec10k = 206265*10000/(delta*1.5*100000000)
                    # circ10k = Circle((xcent-pixelcropping,ycent-pixelcropping),radius = arcsec10k/pixsize, alpha=0.5, fill=False, color='blue', label='10k km')
                    
                    
                    # plotting radial profile
                    
                    df = pd.read_csv(radpath, header=None, sep="\s+")
                
                    # add this because savefig can bug after imshow with lognorm

                    def plot_this_thing():
                        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14,6))
                        circ5arcsec = Circle((xcent,ycent),radius = 5/pixsize, alpha=0.5, fill=False, color='red', label='5arcsec')
                        arcsec10k = 206265*10000/(delta*1.5*100000000)
                        circ10k = Circle((xcent,ycent),radius = arcsec10k/pixsize, alpha=0.5, fill=False, color='blue', label='10k km')
                        
                        ax1.add_patch(circ5arcsec)
                        ax1.add_patch(circ10k)
                        # ax1.axis('off')
                        # ax1.grid()
                        ax1.xaxis.set_major_locator(MultipleLocator(100))
                        ax1.yaxis.set_major_locator(MultipleLocator(100))
                        ax1.grid(which='major', alpha=0.8, color='#CCCCCC', linestyle=':')                 
                        ax1.grid(True)
                        ax1.tick_params(direction='out',length=5)
                        ax1.set_ylim(yimcent - centerlimit, yimcent+centerlimit)
                        ax1.set_xlim(ximcent - centerlimit, ximcent+centerlimit)
                        ax1.legend()
                        centpixel_value = np.max(image_data[int(ycent)-2:int(ycent)+2, int(xcent)-2:int(xcent)+2])
                        t = ax1.text(x=0.95,y=0.05,s='Center pixel value: ' + str(int(centpixel_value)).rjust(5), transform=ax1.transAxes,
                                 horizontalalignment='right')
                        t.set_bbox(dict(facecolor='white', alpha=0.8))
                        ax2.plot(df[2-1][:40], df[4-1][:40])
                        ax2.axvline(x=5/pixsize,color='red', alpha=0.5, linestyle='--')
                        ax2.axvline(x=arcsec10k/pixsize,color='blue', alpha=0.5, linestyle='--')
                        ax2.set_ylabel('Median flux (ADU/s)')
                        ax2.set_xlabel('Distance from the nucleus (pixel)')
                        # add a point in 0,0 to visualize wider profile for gaz filter
                        ax2.plot(0,0, alpha = 0)
                        return fig, ax1
                          
                    try:
                        fig,ax1=plot_this_thing()
                        
                        plt.suptitle(comet_name + ' ' + filt + ' ' + fitsname +  '\n' + str(xcent) + ' ' + str(ycent) + ' (' + ctnmethod + ')')
                        ax1.imshow(image_data, cmap='pink', norm=colors.LogNorm(vmin=np.median(image_data), vmax=image_data[int(xcent), int(ycent)]*2))
                        
                        fig.savefig(os.path.join(save_dir, fitsname[:-5] + '_centering.png'), bbox_inches='tight')
                    except:
                        print('ERROR LOGSCALE')
                        plt.close()
                        fig,ax1=plot_this_thing()
                       
                        plt.suptitle(comet_name + ' ' + filt + ' ' + fitsname + '\n' + str(xcent) + ' ' + str(ycent) + ' (' + ctnmethod + ') error logscale')
                        ax1.imshow(image_data, cmap='binary')
                        plt.savefig(os.path.join(save_dir, fitsname[:-5] + '_centering.png'), bbox_inches='tight')
                    
                    
                    plt.tight_layout()
                    plt.show()
                    plt.close()
                    
# plot_centering_profile('/home/Mathieu/Documents/TRAPPIST/reduced_data/0007P/20210511TN/centering')


def plot_afrho(input_dir, saveplot=''):
    colordict = {
        'R':'tab:red',
        'RC':'tab:orange',
        'V':'tab:green',
        'B':'tab:blue',
        'BC':'tab:purple',
        'I':'tab:olive'
        }
    afrhotable = pd.DataFrame()
    for path, subdirs, files in os.walk(input_dir):
        for file in files:
            if 'afrho' in file and 'tot' in file and ".png" not in file:
                afrhofile = os.path.join(path, file)
                # print(file)
                df = pd.read_csv(afrhofile, header=None, sep="\s+")
                if 'TN' in afrhofile:
                    df['observatory'] = 'TN'
                elif 'TS' in afrhofile:
                    df['observatory'] = 'TS'
                else:
                    input("problem with setting observatory")
                df['obsname'] = afrhofile.split('/')[-3] 
                df['filt'] = file[5:-7]
                afrhotable= afrhotable.append(df)
    
    filtlist = afrhotable['filt'].drop_duplicates().values.tolist()
    
    figtot = plt.figure(figsize=(12,9))
    figtot.suptitle(path.split('/')[-3] + "_afrho_tot")
    axtot = figtot.add_subplot()   
    # axtot.set_yscale("log")
    # axtot.grid(True, which="major", axis='y', ls="-", linewidth=2)
    # axtot.grid(True, which="minor", axis='y', ls="-", linewidth=1)
    afrhotable = afrhotable.sort_values(by=4-1)
    for filt in filtlist:
        filteredtable = afrhotable.loc[afrhotable['filt'] == filt]
        
        model1 = np.poly1d(np.polyfit(filteredtable[4-1]-perihelion, filteredtable[9-1], 2))
        polyline = np.linspace((filteredtable[4-1]-perihelion).min(), (filteredtable[4-1]-perihelion).max(), 50)
        axtot.plot(polyline, model1(polyline), color=colordict.get(filt))
        axtot.scatter(filteredtable[4-1]-perihelion, filteredtable[9-1], zorder=2, label=filt, c=colordict.get(filt))
        axtot.errorbar(filteredtable[4-1]-perihelion, filteredtable[9-1], yerr=[filteredtable[10-1], filteredtable[11-1]],fmt='.', capsize=4, zorder=1, c=colordict.get(filt))
        axtot.axvline(perihelion-perihelion)
        axtot.set_ylabel('Afrho (cm)')
        axtot.set_xlabel('MJD to perihelion')

        MJDlist = []
        afrholist = []
        colorlist = []
        uperrorlist = []
        doerrorlist = []
        for index, row in filteredtable.iterrows() :
            MJDlist.append(row[4-1]-perihelion)
            # JDlist.append(row[4-1])
            afrholist.append(row[9-1])
            uperrorlist.append(row[10-1])
            doerrorlist.append(row[11-1])
            if row['observatory'] == 'TS':
                colorlist.append('r')
            elif row['observatory'] == 'TN':
                colorlist.append('b')
            else:
                colorlist.append('g')
        fig2D = plt.figure(figsize=(12,9))
        ax2D = fig2D.add_subplot()
        ax2D.scatter(MJDlist, afrholist,color=colorlist, zorder=2)
        ax2D.errorbar(MJDlist, afrholist, yerr=[uperrorlist, doerrorlist], fmt='.', color='black', capsize=4, zorder=1)
        # ax2D.scatter(MJDlist, afrholist,color=colorlist)
        # ax2D.axvline(perihelion)
        ax2D.plot(polyline, model1(polyline), color=colordict.get(filt))
        fig2D.suptitle(afrhofile.split('/')[-4] + "_afrho_" + filt)
        ax2D.set_ylabel('Afrho (cm)')
        ax2D.set_xlabel('MJD to perihelion')
        
        if saveplot != '':
            if not os.path.exists(saveplot):
                os.makedirs(saveplot)
            # fig2D.tight_layout()
            fig2D.savefig(os.path.join(saveplot, afrhofile.split('/')[-4] + "_afrho_" + filt + ".png"))

    axtot.legend()
    if saveplot != '':
        if not os.path.exists(saveplot):
            os.makedirs(saveplot)
         # fig2D.tight_layout()
        figtot.savefig(os.path.join(saveplot, path.split('/')[-3] + "_afrho_tot.png"))

# plot_afrho('/home/Mathieu/Documents/TRAPPIST/reduced_data/2020T2')        

def plot_haser(input_dir, saveplot=''):

    # for path, subdirs, files in os.walk(input_dir):
    #     for file in files:
    #         if file == 'outputhasertestall-BC':
    #             haserfile = os.path.join(path, file)
        #             hasertable = pd.read_csv(haserfile, header=None, sep="\s+", comment='#')
    
    hasertable = pd.DataFrame()
    for path, subdirs, files in os.walk(input_dir):
        for file in files:
            if file == 'outputhaser-BC':
                hasfile = os.path.join(path, file)
                # print(file)
                df = pd.read_csv(hasfile, header=None, sep="\s+")
                if 'TN' in hasfile:
                    df['observatory'] = 'TN'
                elif 'TS' in hasfile:
                    df['observatory'] = 'TS'
                else:
                    input("problem with setting observatory")
                df['obsname'] = hasfile.split('/')[-3] 
                hasertable = hasertable.append(df)
                
    filtlist = hasertable[15-1].drop_duplicates().values.tolist()

    figtot = plt.figure(figsize=(12,9))
    axtot = figtot.add_subplot()  
    figtot.suptitle(path.split('/')[-3] + "_haser_tot")
    axtot.set_yscale("log")
    axtot.grid(True, which="major", axis='y', ls="-", linewidth=2)
    axtot.grid(True, which="minor", axis='y', ls="-", linewidth=1)
    for filt in filtlist:
        axtot.scatter(hasertable.loc[(hasertable[15-1] == filt)][5-1]-perihelion, hasertable.loc[(hasertable[15-1] == filt)][12-1], zorder=2, label=filt)
        axtot.errorbar(hasertable.loc[(hasertable[15-1] == filt)][5-1]-perihelion, hasertable.loc[(hasertable[15-1] == filt)][12-1], yerr=[hasertable.loc[(hasertable[15-1] == filt)][13-1], hasertable.loc[(hasertable[15-1] == filt)][14-1]],fmt='.', capsize=4, zorder=1)
        axtot.axvline(perihelion-perihelion)
        axtot.set_ylabel('Q (s-1)')
        axtot.set_xlabel('MJD to perihelion')
        
        filteredtable = hasertable.loc[hasertable[15-1] == filt]
        MJDlist = []
        Qlist = []
        Quppererror = []
        Qlowererror = []
        colorlist = []
        hdist = []
        for index, row in filteredtable.iterrows():
            MJDlist.append(row[5-1]-perihelion)
            hdist.append(row[6-1])
            # JDlist.append(row[5-1])
            Qlist.append(row[12-1])
            Quppererror.append(row[13-1])
            Qlowererror.append(row[14-1])
            if row.observatory == 'TS':
                colorlist.append('r')
            elif row.observatory == 'TN':
                colorlist.append('b')
            else:
                colorlist.append('g')
        
        
        fig2D = plt.figure(figsize=(12,9))
        ax2D = fig2D.add_subplot()
        ax2D.scatter(MJDlist, Qlist,color=colorlist, zorder=2)
        ax2D.errorbar(MJDlist, Qlist, yerr=[Quppererror, Qlowererror], fmt='.', color='black', capsize=4, zorder=1)
        # ax2D.axvline(perihelion)
        fig2D.suptitle(path.split('/')[-3] + "_haser_" + filt)
        ax2D.set_ylabel('Q (s-1)')
        ax2D.set_xlabel('MJD to perihelion')
        
        if saveplot != '':
            if not os.path.exists(saveplot):
                os.makedirs(saveplot)
            # fig2D.tight_layout()
            fig2D.savefig(os.path.join(saveplot, path.split('/')[-3] + "_haser_" + filt + ".png"))
    
    axtot.legend()
    if saveplot != '':
        if not os.path.exists(saveplot):
            os.makedirs(saveplot)
        # fig2D.tight_layout()
        figtot.savefig(os.path.join(saveplot, path.split('/')[-3] + "_haser_tot.png"))

# import directory_structure
# ds = directory_structure.directory_structure()
# perihelion = 59406
# input_path = os.path.join(ds.reduced, "CK21E030")
# savepath = os.path.join(ds.reduced, "CK21E030",'plots')

# plot_afrho(input_dir=input_path, saveplot=savepath)
# plot_haser(input_dir=input_path, saveplot=savepath)
# plot_radprof(input_dir=input_path, saveplot=savepath)
# plot_mag(input_dir=input_path, saveplot='')

def plot_haserprofile(input_dir, output_dir=None, comet_name=''):
    """
    Plot the flux profile for the NB images, the continuum, their difference and the Haser model.
    Create a png for each NB image.

    Parameters
    ----------
    input_dir : str
        path of the folder containing the image and outputhaser-BC (typicaly tmpout)
    output_dir : str, optional
        path for outputting the images. If None, uses input_dir. The default is None.
    comet_name : str, optional
        name of the comet for the plot tile. default is ''
    
    Returns
    -------
    None.

    """
    
    plt.style.use(astropy_mpl_style)
    for path, subdirs, files in os.walk(input_dir):
        for file in files:
            if 'outputhaser-BC' in file:
                centerfile = os.path.join(path, file)
                tab = pd.read_csv(centerfile, header=None, sep="\s+")
                for index, row in tab.iterrows():
                    imname = row[0]
                    Q = row[11]
                    fc = row[2]
                    if len(tab.columns) == 17: #old version
                        filt = row[14]
                        Qproflow = row[15]
                        Qprofhigh = row[16]
                    elif len(tab.columns) == 18: #new version avec error details and only uperror
                        filt = row[15]    
                        Qproflow = row[16]
                        Qprofhigh = row[17]
                    error = row[12]
                    if filt != 'CO+' and filt != 'H2O':
                        haserprofile_path = os.path.join(path, 'haserprofile_' + imname)
                        haserprofilecont_path = os.path.join(path, 'haserprofilecont_' + imname)
                        hasermodel_path = os.path.join(path, 'hasermodel_' + imname)
                        if os.path.isfile(haserprofile_path) and os.path.isfile(haserprofilecont_path) and os.path.isfile(hasermodel_path):
                            obs = pd.read_csv(haserprofile_path, header=None, sep="\s+")
                            cont = pd.read_csv(haserprofilecont_path, header=None, sep="\s+")
                            model = pd.read_csv(hasermodel_path, header=None, sep="\s+")
                            obscont = obs-cont
                        else:
                            print(path)
                            print(imname)
                            print('error finding haser profile paths')
                            input('Acknowledge press enter to quit')
                            return
                        
                        # replace negative values by np.nan for the log scale
                        obscont[obscont < 0] = np.nan
                        obs[obs < 0] = np.nan
                        cont[cont < 0] = np.nan
                        # print(cont.loc[cont[1]<0])
                        
                        save_dir = path if output_dir is None else output_dir
                        
                        # little warning if scale is not the same but should not be a problem
                        if (obs[0][0] != cont[0][0]) or (obs[0][2] != cont[0][2]):
                            print(imname)
                            print('WARNING: x scale may be different for NB and continuum spectrum')
                            print('difference ratio:', ((cont[0]- obs[0])/obs[0])[0])
                        
                        # define the plot
                        fig = plt.figure(figsize=(12,9))
                        ax = fig.gca()
                        
                        # added this so don't see the extreme right of the continuum which is very noisy and mess the scale
                        xlim1 = 0
                        xlim2 = len(cont)-100
                        ax.set_xlim([np.log10(cont[0][xlim1]), np.log10(cont[0][xlim2])])
                        
                        # need to redefine the limit otherwise the plot does not rescale du to high continuum values outside the plotting limits
                        plt.plot(np.log10(obs[0][:xlim2]),np.log10(obs[1][:xlim2]),color='tab:red',ls=':',label='Observed profile')
                        plt.plot(np.log10(cont[0][:xlim2]),np.log10(cont[1][:xlim2]),color='tab:blue',ls=':', label=f'Continuum profile (fc={fc})')
                        plt.plot(np.log10(cont[0][:xlim2]),np.log10(obscont[1][:xlim2]),color='tab:red', label='Profile after dust subtraction')
                        # plt.plot(np.log10(cont[0][:lim]),np.log10(model[1][:lim]),label='tmpmodel.dat')
                        plt.plot(np.log10(model[0]),np.log10(model[1]),color='tab:green',label='Haser model')
                        plt.axvline(x=Qproflow,color='black', alpha=0.5, linestyle='--', label=f'limits for the fit ({Qproflow},{Qprofhigh})')
                        plt.axvline(x=Qprofhigh,color='black', alpha=0.5, linestyle='--')
                        ax.set_ylabel('Column density (unit?)')
                        ax.set_xlabel('Log rho (km)')
                        
                        plt.legend(loc='lower left')
                        haser_magorder = int(math.modf(math.log10(Q))[1])
                        haser_str = "{:.2f}".format(Q/10**haser_magorder)
                        error_str = "{:.2f}".format(error/10**haser_magorder)
                        plt.suptitle(f"{comet_name} {filt} {imname}\nQ({filt}) = {haser_str} ± {error_str} E{haser_magorder} s-1")
                        
                        plt.tight_layout()
                        plt.savefig(os.path.join(save_dir, imname[:-9] + '_haserprofile.png'), bbox_inches='tight')
                        
                        plt.show()
                        plt.close()
    
# plot_haserprofile('/home/Mathieu/Documents/TRAPPIST/tmpout', output_dir=None)

