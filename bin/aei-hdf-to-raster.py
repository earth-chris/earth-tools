#!/usr/bin/python
#####
# aei-raster-to-hdf.py
#  
#  converts a raster file to an hdf5 file with a specific format to store
#   geotiff information
#
# c. 2016 Christopher Anderson
#####

import os
import sys
import aei
import gdal as gdal
import numpy as np
import h5py as hdf

class parse_args:
    def __init__(self, arglist):
        
        # set up main variables and defaults to parse
        self.inFile = None
        self.outFile = None
        self.noData = None
        self.useNoData = False
        self.compress = True
        self.ot = 'Float32'
        self.tempFile = aei.params.scratchdir + aei.params.pathsep + 'temp_hdf2raster.tif'

        # exit if no arguments passed
        if len(arglist) == 1:
            usage(exit=True)
    
        # read arguments from command line
        i = 1
        while i < len(arglist):
            arg = arglist[i]
        
            # check the input file   
            if arg.lower() == '-i':
                i += 1
                arg = arglist[i]
                
                if not aei.fn.checkFile(arg, quiet = True):
                    usage()
                    aei.fn.checkFile(arg)
                    aei.params.sys.exit(1)
                    
                self.inFile = arg
                    
            # check output file
            elif arg.lower() == '-o':
                i += 1
                arg = arglist[i]
                
                self.outFile = arg
                
            # check no-data flag
            elif arg.lower() == '-nodata':
                i += 1
                arg = arglist[i]
                
                try:
                    self.noData = float(arg)
                    self.useNoData = True
                except:
                    usage()
                    print("[ ERROR ]: Unable to set nodata value: %s" % arg)
                    sys.exit(1)
                    
            # check output data type
            elif arg.lower() == '-ot':
                i += 1
                arg = arglist[i]
                
                try:
                    self.ot = aei.read.ot(arg)
                except:
                    usage()
                    print("[ ERROR ]: Unable to read -ot value: %s" % arg)
                    sys.exit(1)
                    
            # check if compression is not set
            elif arg.lower() == '-nocompress':
                self.compress = False
            
            # set up catch-all for incorrect parameter call
            else:
                usage()
                print("[ ERROR ]: Unrecognized argument: %s" % arg)
                print("%s" % i)
                aei.params.sys.exit(1)
            
            # increment the counter for next argument
            i += 1

# set up function to check that arguments have been properly specified
def check_args(args):
    
    # check input file
    if len(args.inFile) == 0:
        usage()
        print("[ ERROR ]: No input file specified.")
        sys.exit(1)
        
    # check output file
    if len(args.inFile) == 0:
        usage()
        print("[ ERROR ]: No output file specified.")
        sys.exit(1)
    
    # return fixed arguments
    return(args)
            
def usage(exit=False):
    """
    describes the aei-hdf-to-raster.py procedure in case of incorrect parameter calls
    
    syntax: usage()
    """
    print(
        """
$ aei-hdf-to-raster.py -i inputFile -o outputFile [-nodata noDataVal]
    [-ot outputDataType] [-nocompress]
        """
        )
    if exit:
        sys.exit(1)

def main():
    """
    the main program for aei-hdf-to-raster.py
    
    the order of operations for this procedure is to
    1) read an input raster file into memory
    2) find and mask the no-data values
    3) write the output, with x-y indices, to an hdf5 file
    
    syntax: main()
    """
    # parse the argument list
    args = parse_args(sys.argv)
    
    # check the argument list to ensure consistent arguments are set
    args = check_args(args)
    
    # if compression is set, we need to create a temp file to write to initially
    if args.compress:
        firstFile = args.tempFile
    else:
        firstFile = args.outFile
    
    # report starting
    print("[ STATUS ]: Starting aei-hdf-to-raster.py")
    print("[ STATUS ]: Input file  : %s" % args.inFile)
    print("[ STATUS ]: Output file : %s" % args.outFile)
    print("[ STATUS ]: ----------")
    print("[ STATUS ]: Reading input data")
    
    # open the hdf  reference file
    with hdf.File(args.inFile, 'r') as inf:
        
        # get the geo and raster data stored in the hdf file
        GeoData = inf.get(inf.keys()[0])
        RasterData = inf.get(inf.keys()[1])
        
        # set a key list to loop through and get data from
        rasterKeys = RasterData.keys()
        
        # remove x and y indices to not read them
        rasterKeys.remove('X-Indices')
        rasterKeys.remove('Y-Indices')
        
        # the number of RasterData keys should be the number out output bands
        onb = len(rasterKeys)
    
        # read the GeoTransform data, Projection, and file dimensions
        GeoTransform = np.array(GeoData.get('GeoTransform'))
        Projection = str(np.array(GeoData.get('Projection')))
        nx, ny, nb = np.array(GeoData.get('[nx, ny, nb]'))
        
        # create the output file and set the parameters
        print("[ STATUS ]: Creating output file")
        outRef = gdal.GetDriverByName("GTiff").Create(firstFile, nx, ny, onb, gdal.GDT_Float32)
        outRef.SetGeoTransform(GeoTransform)
        outRef.SetProjection(Projection)
        
        # read the x and y indices
        print("[ STATUS ]: Reading X-Y Indices")
        x = np.array(RasterData.get('X-Indices'))
        y = np.array(RasterData.get('Y-Indices'))
        
        # create an output array to store outputs
        outArr = np.zeros((ny, nx), dtype = np.float32)
        
        # if no-data included, set it here
        if args.useNoData:
            outArr += args.noData
        
        # loop through each band, read the band info from the hdf file, and write to raster
        for i in range(onb):
            
            # report
            print("[ STATUS ]: Reading data set: %s" % (rasterKeys[i]))
            
            # get the 1-d data
            bandData = np.array(RasterData.get(rasterKeys[i]), dtype = np.float32)
            
            # put it in the 2-d array
            outArr[y, x] = bandData
            
            # get the band reference
            bandRef = outRef.GetRasterBand(i + 1)
            
            # write the array
            bandRef.WriteArray(outArr)
            
            # set no-data if necessary
            if args.useNoData:
                bandRef.SetNoDataValue(args.noData)
                
            # destroy the band reference and free memory
            bandRef = None
            bandData = None
        
        # once finished looping through each band, kill the gdal reference
        outArr = None
        outRef = None
        
    # compress the output file
    if args.compress:
        print("[ STATUS ]: Compressing output file")
        aei.cmd.gdal_translate(args.tempFile, args.outFile, 
            etc=['-co','COMPRESS=LZW'], ot=args.ot)
        aei.params.os.remove(args.tempFile)
            
    # report finished
    print("[ STATUS ]: ----------")
    print("[ STATUS ]: Finished writing raster data!")
    print("[ STATUS ]: Please see output file : %s" % args.outFile)
    
# call the aain routine when run from command lne
if __name__ == "__main__":
    main()