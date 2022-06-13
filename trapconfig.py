#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 13 11:41:56 2022

@author: Mathieu Vander Donckt
mathieu.vanderdonckt@uliege.be
"""

import os

# look if traputils.conf is located in the home dir, otherwise take the one in the current dir
if os.path.isfile(os.path.join(os.path.expanduser('~'), 'traputils.conf')):
    path_conf = os.path.join(os.path.expanduser('~'), 'traputils.conf')
elif os.path.isfile(os.path.join(os.path.realpath(os.path.dirname(__file__)), 'traputils.conf')):
    path_conf = os.path.join(os.path.realpath(os.path.dirname(__file__)), 'traputils.conf')
else:
    print('ERROR: could not find traputils.conf, see trapconfig.py')
print('Config file loaded from', path_conf)

#import parameters from the config file  
param = {} # will old the info given in trapconf
with open(path_conf, 'r') as f: 
    contents = f.readlines()
    for line in contents:
        if line[0] == '#':
            continue
        if line.endswith('\n'):
            line = line[:-1]
        param[line.split(':')[0]] = line.split(':')[1]

# print(param)
