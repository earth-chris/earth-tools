#!/usr/bin/python
##
## reads and applies standard processing parameters to planet rgb imagery
##
####################
##
import aei as aei

# define global parameters
class default_params:
    def __init__(self):
        
        # set base file parameters
        self.imageFile = None
        self.metadataFile = None
        self.geoinfoFile = None
        self.inputPath = None
        self.basename = None
        
        # set some metadata parameters
        self.snr = None
        self.exposureTime = None
        self.altitude = None
        self.offNadir = None
        self.solarZenith = None
        self.solarAzimuth = None
        self.localTime = None
        self.cloudCover = None

# parse command line arguments
class parse_args:
    def __init__(self, arglist, params):
        
        # we'll be parsing the command line arguments passed
        #  via arglist, updating the default parameters,
        #  then returning the params list 
        
        # exit if no arguments passed
        if len(arglist) == 1:
            usage(exit=True)
            
        # get glob to search for files
        import glob
    
        # read arguments from command line
        i = 1
        print(arglist)
        
        while i < len(arglist):
            arg = arglist[i]
        
            # check input flag    
            if arg.lower() == '-i':
                i += 1
                arg = arglist[i]
                
                # if input is a path, find the input files
                if aei.params.os.path.isdir(arg):
                    params.inputPath = arg
                    search = glob.glob(arg + aei.params.pathsep + "*")
                    
                    # reject if no files found
                    if len(search) == 0:
                        usage()
                        print("[ ERROR ]: no files found in input directory")
                        sys.exit(1)
                        
                    for file in search:
                        if ".tif" in file:
                            params.imageFile = file
                        if "_metadata.json" in file:
                            params.metadataFile = file
                        if "_geoinfo.json" in file:
                            params.geoinfoFile = file
                
                # if input is a .tif file            
                else:
                    dirname = aei.params.os.path.dirname(arg)
                    params.inputPath = dirname
                    params.imageFile = arg
                    
                    # find json files
                    loc = arg.find("_analytic.tif")
                    
                    if loc == -1:
                        usage()
                        print("[ ERROR ]: input tif file does not end with _analytic.tif")
                        print("[ ERROR ]: input tif file: %s" % arg)
                        sys.exit(1)
                        
                    else:
                        params.metadataFile = arg[:loc] + "_metadata.json"
                        params.geoinfoFile = arg[:loc] + "_geoinfo.json"

# read the metadata json file        
def readMetadata(params):
    import json
    
    # read the json metadata info into memory
    jdict = json.load(params.metadataFile)
    j = jdict[u'properties']
    
    # extract useful metadata info
    self.exposureTime = j[u'camera'][u'exposure_time']
    self.cloudCover = j[u'cloud_cover'][u'estimated']
    self.snr = j[u'snr']
    self.altitude = j[u'sat'][u'alt']
    self.offNadir = j[u'sat'][u'off_nadir']
    self.solarZenith = 90-j[u'sun'][u'altitude']
    self.solarAzimuth = j[u'sun'][u'azimuth']
    self.localTime = j[u'sun'][u'local_time_of_day']

# make sure files exist    
def checkInputs(params):
    if not aei.fn.checkFile(params.imageFile):
        usage()
        aei.fn.checkFile(params.imageFile, False)
        sys.exit(1)
        
    if not aei.fn.checkFile(params.metadataFile):
        usage()
        aei.fn.checkFile(params.metadataFile, False)
        sys.exit(1)
           
    if not aei.fn.checkFile(params.geoinfoFile):
        usage()
        aei.fn.checkFile(params.geoinfoFile, False)
        sys.exit(1)
    
def usage(exit=False):
    """
    describes the aeiplanet.py procedure in case of incorrect parameter calls
    
    syntax: usage(exit=False)
    """

    print(
        """
$ aeiplanet.py -i input_files
    
    output band list is:
    1. Red  2. Green  3. Blue  4. Alpha  5. 
        """)
        

def main ():
    """
    the main program for aeiplanet.py

    syntax: main()
    """
    import gdal
    import numpy as np
    import cv2

    # first, read the default parameters to get the processing object
    params = default_params()
    
    # then parse the arguments passed via command line
    args = aei.params.sys.argv
    parse_args(args, params)
    
    # check inputs
    checkInputs(params)
    
    # parse the json data
    readMetadata(params)
    
    # read the rgb image
    gdalRef = gdal.Open(params.imageFile)
    rgba = gdalRef.ReadAsArray()
    
    