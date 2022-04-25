#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 25 12:31:42 2022

@author: Mathieu Vander Donckt
mathieu.vanderdonckt@uliege.be

see https://github.com/MichaelKim0407/tutorial-pip-package
"""

from setuptools import setup

setup(
    name='TRAPPIST-utils',
    version=0,

    url='https://github.com/mathvdd/TRAPPIST-utils',
    author='Mathieu Vander Donckt',
    author_email='mathieu.vanderdonckt@uliege.be',

    packages=['trap_reduction'],
    package_dir={'trap_reduction':'./'}
)
