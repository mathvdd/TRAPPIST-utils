#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 20 15:22:08 2023

@author: Mathieu
"""

from astropy.io import fits
import matplotlib.pyplot as plt
import requests
import math

# impath = '/home/Mathieu/Documents/TRAPPIST/reduced_data/0012P/20231118TN/images/TRAP.2023-11-18T18:43:42.fits'
# comet = '12P'
# crop = 200
# upshift =100
# rightshift=100
# vmin = 30
# vmax = 50

# impath = '/home/Mathieu/Documents/TRAPPIST/reduced_data/0012P/20231118TN/images/TRAP.2023-11-18T18:16:39.fits'
# comet = '12P'
# crop = 400
# upshift =50
# rightshift=10
# vmin = 14000
# vmax = 18000

impath = '/home/Mathieu/Documents/TRAPPIST/reduced_data/0012P/20231118TN/images/TRAP.2023-11-18T19:57:49.fits'
comet = '12P'
crop = 400
upshift =50
rightshift=10
vmin = 1500
vmax = 15000

obs= impath.split('/')[-3][-2:]

with fits.open(impath) as hdul:
    header = hdul[0].header
    image = hdul[0].data

if (obs == 'TN'):
    obs_code = 'I40'
    if (header['PIERSIDE'] == 'EAST'):    
        N = 0
    elif (header['PIERSIDE'] == 'WEST'):
        N = 180
elif (obs == 'TS'):
    obs_code = 'Z53'
    if (header['PIERSIDE'] == 'EAST'):
        N = 180
    elif (header['PIERSIDE'] == 'WEST'):
        N = 0
        
jd = header['JD']
filt = header['FILTER']
date_obs = header['DATE-OBS']

parameters = {
            'COMMAND' : 'DES='+comet+ '%3BCAP%3BNOFRAG%3B',
            'OBJ_DATA' : 'NO',
            'MAKE_EPHEM' : 'YES',
            'TABLE_TYPE' : 'OBSERVER',
            'CENTER' : obs_code,
            'START_TIME' : f'JD{jd}',
            'STOP_TIME' : f'JD{jd+(1/24/60/60)}',
            'STEP_SIZE' : '1',
            'QUANTITIES' : '27', #RA and DEC, rh, delta, phase angle 
            'CAL_FORMAT' : 'BOTH',
            'CSV_FORMAT' : 'YES',
            'REF_SYSTEM' : 'ICRF'}

query_url = "https://ssd.jpl.nasa.gov/api/horizons.api?format=text"
for param in parameters:
    query_url += '&' + param + "='" + str(parameters[param]) + "'"
query_result = requests.get(query_url).text.split('\n')

try:
    colnames = [i.strip(' ') for i in query_result[query_result.index('$$SOE') - 2].split(',')] 
    values = [i.strip(' ') for i in query_result[query_result.index('$$SOE') + 1].split(',')]
    PsAng = N+float(values[colnames.index('PsAng')])
    PsAMV = N+float(values[colnames.index('PsAMV')])
except:
    print(query_result)
    input('ERROR looking for JPL ephem')
    
"""
The position angles of the extended Sun-to-target radius vector ("PsAng")
and the negative of the targets' heliocentric velocity vector ("PsAMV"),
as seen in the observers' plane-of-sky, measured CCW (east) from reference
frame north-pole. Computed for small-bodies only (and primarily intended
for ACTIVE COMETS), "PsAng" is an indicator of the comets' gas-tail
orientation in the sky (being in the anti-sunward direction) while "PsAMV"
is an indicator of dust-tail orientation.

  Units: DEGREES

  Labels: PsAng PsAMV
"""
    

plt.imshow(image, cmap='binary', vmin=vmin, vmax=vmax)
xleft = crop+rightshift
xright = image.shape[0]-crop+rightshift
ybot = crop+upshift
ytop = image.shape[1]-crop+upshift
plt.xlim(left=xleft,right=xright)
plt.ylim(bottom=ybot,top=ytop)

# plt.arrow(400, 400, 100*math.cos(PsAng), 100*math.sin(PsAng), color='r')
# plt.arrow(400, 400, 100*math.cos(PsAng), 100*math.sin(PsAMV), color='r')

plt.annotate("", xy=((xleft+(xright-xleft)/5)+((xright-xleft)/6*math.cos(PsAng)),
                     (ytop-(ytop-ybot)/5)+((ytop-ybot)/6*math.sin(PsAng))),
             xytext=((xleft+(xright-xleft)/5),
                     (ytop-(ytop-ybot)/5)),
             arrowprops=dict(arrowstyle="->",color='r',lw=3),
             size=15)
plt.annotate("", xy=((xleft+(xright-xleft)/5)+((xright-xleft)/6*math.cos(PsAMV)),
                     (ytop-(ytop-ybot)/5)+((ytop-ybot)/6*math.sin(PsAMV))),
             xytext=((xleft+(xright-xleft)/5),
                     (ytop-(ytop-ybot)/5)),
             arrowprops=dict(arrowstyle="->",color='b',lw=3),
             size=15)
t1=plt.text(x=(xleft+(xright-xleft)/5)+((xright-xleft)/6*math.cos(PsAng)),
         y=(ytop-(ytop-ybot)/5)+((ytop-ybot)/6*math.sin(PsAng)),
         s= 'Anti-sunward', color='r',size=15)
t1.set_bbox(dict(facecolor='w', alpha=0.7))
t2=plt.text(x=(xleft+(xright-xleft)/5)+((xright-xleft)/6*math.cos(PsAMV)),
         y=(ytop-(ytop-ybot)/5)+((ytop-ybot)/6*math.sin(PsAMV)),
         s= 'Anti-velocity', color='b',size=15)
t2.set_bbox(dict(facecolor='w', alpha=0.7))

plt.title(f"{comet} {filt} {date_obs}")
plt.colorbar()
plt.tight_layout()
plt.savefig(impath[:-4] + 'png')

