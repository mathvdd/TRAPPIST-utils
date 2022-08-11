#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 11 18:57:27 2022

@author: Mathieu
"""

from trapconfig import param
import os

f2 = []
with open(os.path.join(param['tmpout'], 'outf'), 'rb') as f:
    for line in f:
        try:
            # print(line.split())
            for item in line.split():
                # print(item)
                float(item)
            f2.append(line)
        except:
            print('removed from outf:', line)
            pass
            # print(line)
            
with open(os.path.join(param['tmpout'], 'outf2'), 'w') as f:
    for line in f2:
        f.write(line.decode('UTF-8'))