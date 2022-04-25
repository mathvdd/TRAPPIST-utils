#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 28 15:07:56 2021

@author: Mathieu Vander Donckt
mathieu.vanderdonckt@uliege.be
"""
import os

# raw_directory = "/home/Mathieu/Documents/TRAPPIST/raw_data/"
# reduced_directory = "/home/Mathieu/Documents/TRAPPIST/reduction/"
# raw_path = "/home/Mathieu/Documents/TRAPPIST/raw_data/2020T2/TN/20210714"
# reduced_path = "/home/Mathieu/Documents/TRAPPIST_reduction/20210714"

class directory_structure:
    """
    Class containing path to different useful directories
    Can also be used to prompt for new directories
    """
    
    def __init__(self):
        self.home = os.path.expanduser('~')
        self.TRAPDIR = os.path.join(self.home, 'Documents/TRAPPIST')
        self.raw = os.path.join(self.TRAPDIR, 'raw_data')
        self.reduced = os.path.join(self.TRAPDIR, 'reduced_data')
        self.tmpout = os.path.join(self.TRAPDIR, 'tmpout')
        self.tmpdata = os.path.join(self.TRAPDIR, 'tmpdata')
        self.iraf = os.path.join(self.home, 'iraf')
        self.calib = os.path.join(self.TRAPDIR, 'calib')
        self.hasercalc = os.path.join(self.TRAPDIR, 'hasercalc')
        
    def set_dir_input(self, dir_name='new_directory', path='', create_new=True, complete_path=False):
        """Sets a directory path with user input.
        
        Parameters:
            dir_name (str, optional): name of the directory that will be created.
                Used for prompting the user and eventually suggest a path.
            path (str, optional): suggesed path (or above directory if complete_path = True)
            create_new (boolean): defines if the input directory needs should exists, otherwise will be created
            complete_path (boolean): set the path as the suggested path concatenated with the suggested path
        """
        
        if path == "":
            path = os.path.join(self.home, dir_name)

        while True:
            if complete_path == True:
                print("\nEnter ls for a list of directories")
                input_path = input('Set a path for ' + dir_name + " inside directory " + path + ": ")
            else:
                input_path = input('Set a path for ' + dir_name + " directory (" + path + "): ")
            
            if input_path == '':
                input_path = path
            
            if complete_path == True:
                if input_path == 'ls':
                    os.system('\n'.join(['cd ' + path,
                     'ls']))
                    continue
                else:
                    input_path = os.path.join(path, input_path)
            
            if os.path.isdir(input_path) == False:
                if create_new == False:
                    print(input_path + ' does not exist. Enter an existing directory')
                elif create_new == True:
                    print("Directory does not exists. Building directories.")
                    os.makedirs(input_path)
                    break
            else:
                break
        print(dir_name + " directory set to " + input_path)
        return input_path