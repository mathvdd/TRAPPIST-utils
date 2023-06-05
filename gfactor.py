#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  2 09:08:50 2023

@author: math
"""

import pandas as pd
import numpy as np
from trapconfig import param
import os

calib_path = param['calib']

jean = pd.read_csv(os.path.join(calib_path, 'Fluo_eff_Schleicher.dat'), sep = '\s+')
CN_tab = pd.read_csv(os.path.join(calib_path,'gCN_Schleicher_2010.dat'), sep = '\s+')
C2C3_tab = pd.read_csv(os.path.join(calib_path,'gC2C3_AHearn_1995.dat'), sep = '\s+')
OH_tab = pd.read_csv(os.path.join(calib_path,'gOH_Schleicher_1988.dat'), sep = '\s+')
# make NHtab in same format as the CN tab
NHtab = pd.DataFrame()
for i, row in jean.iterrows():
    NHtab.loc[row['v'],row['r']] = row['NH']*row['r']**2
NHtab.reset_index(inplace=True)

def intrapolation(x, a1, a2, b1, b2):
    # intrapolation of b1 and b2 between a1 and a2
    px = b1 + (b2-b1) * ( (x-a1) / (a2-a1) )
    return px

def g_CN(v, rh):
    #get the index of the lines with the colsest values to v
    vsub = CN_tab['v'].sub(v)
    v_min_index = abs(vsub.loc[vsub<=0]).idxmin()
    v_max_index = abs(vsub.loc[vsub>=0]).idxmin()  
    #compute the intrapolation of v if necessary and convert to pd series and return a line
    if v_min_index != v_max_index:
        low_vrow = CN_tab.loc[[v_min_index]]
        up_vrow = CN_tab.loc[[v_max_index]]
        intrav = intrapolation(v, low_vrow['v'].values[0], up_vrow['v'].values[0], low_vrow.squeeze(), up_vrow.squeeze())
    else:
        intrav = CN_tab.loc[[v_min_index]].squeeze()
    #some necessary transformation to the table/series
    rhtab = pd.DataFrame(np.vstack([intrav.to_frame().transpose().columns, intrav.to_frame().transpose()]))
    rhtab.drop(0, inplace=True, axis=1)
    rhtab = rhtab.transpose().astype('float')
    #compute the closest index to rh
    rhsub = rhtab[0].sub(rh)
    rh_min_index = abs(rhsub.loc[rhsub<=0]).idxmin()
    rh_max_index = abs(rhsub.loc[rhsub>=0]).idxmin()
    #makes the intrapolation on rh if necessary
    if rh_min_index != rh_max_index:
        g_rh2 = intrapolation(rh, rhtab.loc[[rh_min_index]][0].values[0], rhtab.loc[[rh_max_index]][0].values[0],
                        rhtab.loc[[rh_min_index]][1].values[0], rhtab.loc[[rh_max_index]][1].values[0])
    else:
        g_rh2 = rhtab.loc[rh_min_index][1]
    #the value in the table is multiplied by rh^2 and to 10*-13
    g = g_rh2/rh**2*10**(-13)
    return g

# g_CN = g_CN(-2,4)
# print(g_CN)

def g_C2(rh):
    gC2_1AU = C2C3_tab['g_C2'].values[0]
    # print(gC2_1AU)
    gC2 = gC2_1AU/rh**2
    return gC2

def g_C3(rh):
    gC3_1AU = C2C3_tab['g_C3'].values[0]
    # print(gC3_1AU)
    gC3 = gC3_1AU/rh**2
    return gC3

def g_OH(v, rh):
    vsub = OH_tab['v'].sub(v)
    v_min_index = abs(vsub.loc[vsub<=0]).idxmin()
    v_max_index = abs(vsub.loc[vsub>=0]).idxmin()
    #makes the intrapolation on rh if necessary
    if v_min_index != v_max_index:
        g_rh2 = intrapolation(v, OH_tab.loc[[v_min_index]]['v'].values[0], OH_tab.loc[[v_max_index]]['v'].values[0],
                        OH_tab.loc[[v_min_index]]['gOH(0-0)'].values[0], OH_tab.loc[[v_max_index]]['gOH(0-0)'].values[0])
    else:
        g_rh2 = OH_tab.loc[v_min_index][1]
    # scale by rh^2 and to 10*-15 (table values)    
    g = g_rh2/rh**2*10**(-15)
    return g
    
def g_NH(v, rh):
    
    #get the index of the lines with the colsest values to v
    vsub = NHtab['index'].sub(v)
    v_min_index = abs(vsub.loc[vsub<=0]).idxmin()
    v_max_index = abs(vsub.loc[vsub>=0]).idxmin()  
    #compute the intrapolation of v if necessary and convert to pd series and return a line
    if v_min_index != v_max_index:
        low_vrow = NHtab.loc[[v_min_index]]
        up_vrow = NHtab.loc[[v_max_index]]
        intrav = intrapolation(v, low_vrow['index'].values[0], up_vrow['index'].values[0], low_vrow.squeeze(), up_vrow.squeeze())
    else:
        intrav = NHtab.loc[[v_min_index]].squeeze()
    #some necessary transformation to the table/series
    rhtab = pd.DataFrame(np.vstack([intrav.to_frame().transpose().columns, intrav.to_frame().transpose()]))
    rhtab.drop(0, inplace=True, axis=1)
    rhtab = rhtab.transpose().astype('float')
    #compute the closest index to rh
    rhsub = rhtab[0].sub(rh)
    rh_min_index = abs(rhsub.loc[rhsub<=0]).idxmin()
    rh_max_index = abs(rhsub.loc[rhsub>=0]).idxmin()
    #makes the intrapolation on rh if necessary
    if rh_min_index != rh_max_index:
        g_rh2 = intrapolation(rh, rhtab.loc[[rh_min_index]][0].values[0], rhtab.loc[[rh_max_index]][0].values[0],
                        rhtab.loc[[rh_min_index]][1].values[0], rhtab.loc[[rh_max_index]][1].values[0])
    else:
        g_rh2 = rhtab.loc[rh_min_index][1]
    #the value in the table is multiplied by rh^2
    g = g_rh2/rh**2
    return g

def generate_tmptxt(v, rh):
    try:
        print('gfac')
        # tmpbrol = pd.read_csv(os.path.join(param['tmpout'],'tmpbrol'), sep='\s+', names=['jd','RA','DEC', 'rh','rhdot','Delta','Deltadot','phase'])
        # rh = tmpbrol['rh'][0]
        # v = tmpbrol['rhdot'][0]
        line= f"{'%.3f'%rh} {'%.3f'%v} {'%.2E'% g_OH(v=v,rh=rh)} {'%.2E'% g_NH(v=v,rh=rh)} {'%.2E'% g_CN(v=v,rh=rh)} {'%.2E'% g_C3(rh=rh)} {'%.2E'% g_C2(rh=rh)} 1"
        print(line)
        with open(os.path.join(param['tmpout'],'tmp.txt'), 'w') as f:
            f.write(line)
    except:
        input('ERROR generating g-factor file')

# v,rh=7.54,0.63
# print(g_OH(v=v,rh=rh),g_NH(v=v,rh=rh),g_CN(v=v,rh=rh),g_C3(rh=rh),g_C2(rh=rh))
# if __name__ == "__main__":
#     generate_tmptxt()