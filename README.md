# TRAPPIST-utils
Utilities for the TRAPPIST image reduction for comets. Complementary to the iraf scripts. Main functions are highlighted bellow.

[trap_reduction](#trap_reductionpy)

[get_ephem](#get_ephempy)

[NAS_query](#query_NASpy)

[trap_plot](#trap_plotpy)

## trap_reduction.py
**renameftsfits**(raw_path):

    Rename 'fts' into 'fits'

    Parameters:
        raw_path (str): path to the folder containing the fits files (subfolders will be affected as well)

**pythrename**(raw_path, tmpdata_dir):


    Rename fits files to the trappist format and copy them to the tmpdata folder

    Parameters:
        raw_path (str): path to the raw file folder
        tmpdata_dir (str): path to the directory where the files are to be copied



**check_calib**(fitstable, filt_list=filt_list):


    Check whether all necessary calibration are available in the directory before recduction
    Return False if there is a warning flag

    Parameters:
        fitstable: table with fits file parameters (see get_fitstable())
        filt_list, optional: list of filters to take into consideration. Default is ['OH','CN','C2','C3','NH','BC','RC','GC','R','I', 'B', 'V']



**get_fitstable**(raw_dir):


    Scan a directory to get a table with relevant infos from the fits headers
    Return the table

    Parameters:
        raw_dir (str): path to the directory to scan



**generate_ZP**(calib_dir, ephemeris, fitstable, ZPparams=ZP):


    Generate calib.dat file with closest zero points

    Parameters:
        calib_dir (str): path to the directory containing zero point files. calib.dat will be copied there
        ephemeris: ephemeris object (see get_ephem.py)
        fitstable: table containing the fits files info (see get_fitstable())
        ZPparams, optional: table containing constant values in the calib.dat file. default should be fine at all time



**clreduce**(iraf_dir):


    Wrapper to launch iraf progtrap3 in python

    Parameters:
        iraf_dir (str): path to the home iraf directory



**clafrhocalcext**(iraf_dir, pixsize, solocomete, soloinitx, soloinity, soloinitcboxsize):

    Wrapper to launch iraf afrhocalcext in python

    Parameters:
        iraf_dir (str): path to the home iraf directory
        pixsize (float): size of the image pixel
        solocomete ()
        soloinitx
        soloinity
        soloinitcboxsize


**check_darks**(iraf_dir,tmpout_dir):


    Checks if all darks are made for every line of list_D_exptime
    return False if not

    Parameters:
        iraf_dir (str): path to the home iraf directory
        tmpout_dir (str): path to the directory to check



**check_haser_continuum**(tmpout):


    Obscelete
    Check if there is a BC image for duct continuum correction

    Parameters:
        tmpout (str): directory to check



**generate_haserinput**(tmpout, fc=fc, fz=0):


    generate the haserinput file for the iraf script. Select the BC image. Gives a warning on error

    Parameters:
        tmpout (str): working directory
        fc (dic, optional): dictionnaru containing the fc for each filter. See file for default
        fz (int, optional): value of fz. default =0



**clhasercalctest**(iraf_dir, arg='no'):


    Wrapper to launch iraf hasercalctest in python

    Parameters:
        iraf_dir (str): path to the home iraf directory



**clean_afrhotot**(direc):


    Clean the afrhotot files to keep only the last computed value for each image

    Parameters:
        direc (str): path to the working directory

## get_ephem.py
**ephemeris()**:


    Initiate an ephemeris class with defaults parameters to query NASA Horizons



**ephemeris.retrieve_param_from_fits**(self, fits_dir):


        Open the first fits files in 'fits_dir' starting with TRAP to set the observatory in the ephemeris class
        Open all the fits files in 'fits_dir' to get the dates and set the ephemeris range in the ephemeris class

        Parameters:
            fits_dir (str): path to the directory



**ephemeris.query_horizons**(self):


        Query NASA Horizons and format the result



**ephemeris.query_input**(self, unique_target=False, target=None, convert_MPC_Horizon=False):


        Query for the object name and launch query_horizons()
        Loop until a query is successful, otherwise ask for a new input

        Parameters:
            unique_target (str, optional, default=False): set to True to remember the name of the target of the first successful query.
                Useful for reducing different night with the same object.
                Set to False if using nights from different objects.
            target (str, optional, default=None): if different from None, use as initial object name input for query_horizons()
            convert_MPC_Horizon (boolean, optional, default=False): if is True and target different than None, covert target_name from MPC to NASA Horizon format



**ephemeris.generate_ephem_files**(self, output_dir):


        Generates ephem.brol and eph.dat file in 'output_dir'.
        Dates are expressed in MJD at the moment but can be changed in the script (then also need to be changed in afrhocalcext)

        Parameters:
            output_dir (str): path to the output directory


## query_NAS.py

**NAS_build**(NAS_path, export_path, keyword=''):


    Generate an indexed database of .fits and .fts file in the NAS of binning 2, containing keywords from the fits headers
    An indexed database can be queryed much faster than going through all the files each time

    Parameters
    ----------
    NAS_path : str
        path the NAS to query
    export_path : str
        path for export of the database (with the filename)
    keyword : str, optional
        if different than '', needs to be present in the folder path for it to be considered. The default is ''.

    Returns
    -------
    None.



**check_objects_names**(startdate, enddate, NASfitstable):


    Makes a list of all the objects in a table between given dates
    if the dates are at 12:00, will include the startdate as starting night but not the enddate as starting night

    Parameters:
        startdate (datetime): starting time
        enddate (datetime): ending time
        NASfitstable (pandas dataframe): indexed database-style table



**get_files**(obj_name, NASfitstable, output_path, dayinterval, dateinterval=['','']):

    Look for and prompt for downloading files in the NAS
    It is advised to not mix TN and TS queries
    look on the TN or TS NAS for all the nights containing the object (obj_name).
    Night by night, look for the calibration files for a given dayinterval.
    Download the lights and found calibs and put them in a date folder.
    Date interval can be given to select data betwenn two starting nights intervals

    Parameters:
        obj_name (str): object name in the NASfitstable
        NASfitstable (pandas dataframe): indexed table-style table
        output_path (str): path to a folder to save the files
        dayinterval (int): day interval to look for calib files before and after the observation night
        date_interval (list [datetime, datetime], optional, default=['','']): of 2 terms: the starting date for the search and second the ending date, in datetime format
            if '' given for one of the argument, consider no date limit on that end


**lookforcalib**(NASfitstable, imtype, output_fold, night, obj='', exptime=0, filt='R', dayinterval=0, copy=True):


    Look for a specific set of file the closest to the observation file in the database and download it. By default look for the files closest to the observation night

    Parameters:
        NASfitstable (pandas dataframe): indexed-type database, can be imported with loadcsvtable()
        imtype (str): image type. can take values of 'light', 'dark', 'flat' or 'bias'
        output_fold (str): path to the output directory.
            Inside this directory needs a subfolder with the observation night in YYYYMMDD format and a subsubfolder Calibration
        night ((int,int,int)): night of observation in a tuple of int (YYYY,MM,DD), to be changed
        obj (str, semi-optional, defaut=''): name of the object if imtype='light'
        exptime (int, semi-optional, defaut=0): exposure time if imtype='dark'
        filt (str, semi-optional, default='R'): name of the filter if imtype='flat'
        dayinterval (int, optional, default=0): the starting point in time to look for image from the observation night
        copy (boolean, optional, default=True): wheter to prompt for copying the files or not



## trap_plot.py

**plot_centering**(input_dir, output_dir=None):


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

**plot_centering_profile**(input_dir, output_dir=None):

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

    Returns
    -------
    None.

**plot_haserprofile**(input_dir, output_dir=None):

    Plot the flux profile for the NB images, the continuum, their difference and the Haser model.
    Create a png for each NB image.

    Parameters
    ----------
    input_dir : str
        path of the folder containing the image and outputhaser-BC (typicaly tmpout)
    output_dir : str, optional
        path for outputting the images. If None, uses input_dir. The default is None.

    Returns
    -------
    None.
