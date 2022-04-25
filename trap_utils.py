#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 25 11:26:03 2022

@author: Mathieu Vander Donckt
mathieu.vanderdonckt@uliege.be

See the README.md
"""

import get_ephem
import query_NAS
import trap_plot
import trap_reduction

def ephemeris():
    return get_ephem.ephemeris()

def plot_centering(input_dir, output_dir=None):
    return trap_plot.plot_centering(input_dir, output_dir)

def NAS_build(NAS_path, export_path, keyword=''):
    return query_NAS.NAS_build(NAS_path, export_path, keyword)

# def check_objects_names(startdate, enddate, NASfitstable):
#     return