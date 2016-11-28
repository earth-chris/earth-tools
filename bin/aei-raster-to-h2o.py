#!/usr/bin/python
#####
# aei-raster-to-hdf.py
#  
#  
#
# c. 2016 Christopher Anderson
#####

import os
import sys
import aei
import gdal as gdal
import numpy as np
import h2o

class parse_args:
    def __init__(self, arglist):
        
        # set up main variables and defaults to parse
        self.inFile = None
        self.outFile = None
        self.noData = None
        self.useNoData = False
        self.mem = str("%01dG" % (0.6 * aei.params.mem / (1024 * 1024 * 1024.)))
        self.cores = aei.params.cores - 1

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
                    print("[ ERROR ]: Unable to set nodata value: %s" % arg)
            
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
    describes the aei-raster-to-h2o.py procedure in case of incorrect parameter calls
    
    syntax: usage()
    """
    print(
        """
$ aei-raster-to-h2o.py -i inputFile -o outputFile -nodata noDataVal
        """
        )
    if exit:
        sys.exit(1)

def main():
    """
    the main program for aei-raster-to-h2o.py
    
    the order of operations for this procedure is to
    1) read an input raster file into memory
    2) find and mask the no-data values
    3) write the output, with x-y indices, to an h2o file
    
    syntax: main()
    """
    # parse the argument list
    args = parse_args(sys.argv)
    
    # check the argument list to ensure consistent arguments are set
    args = check_args(args)
    
    # report starting
    print("[ STATUS ]: Starting aei-raster-to-h2o.py")
    print("[ STATUS ]: Input file  : %s" % args.inFile)
    print("[ STATUS ]: Output file : %s" % args.outFile)
    print("[ STATUS ]: ----------")
    
    # initialize the h2o cluster
    print("[ STATUS ]: Initializing h2o")
    h2o.init(ip='localhost', nthreads = args.cores, max_mem_size = args.mem)
    
    # report starting
    print("[ STATUS ]: Reading input data")
    
    # open the reference file
    inRef = gdal.Open(args.inFile)
    
    # get metadata info
    proj = inRef.GetProjection()
    geot = inRef.GetGeoTransform()
    innx = inRef.RasterXSize
    inny = inRef.RasterYSize
    innb = inRef.RasterCount
    
    # handle band 1 first to get indices to store
    bandRef = inRef.GetRasterBand(1)
        
    # check if no-data value set by user, otherwise read from file
    if not args.useNoData:
        args.noData = bandRef.GetNoDataValue()
        
        # check that it exists
        if args.noData is None:
            print("[ ERROR ]: No no-data value set for input file")
            print("[ ERROR ]: Please set -nodata value, or set in input file")
            sys.exit(1)
    
    # read into an array
    bandArr = bandRef.ReadAsArray()
    
    # find good data values
    print("[ STATUS ]: Reading good-data indices")
    gd = np.where(bandArr != args.noData)
    
    # error check here
    if not gd[0].any():
        print("[ ERROR ]: No good-data values found.")
        sys.exit(1)
    
    # create a list to store output names
    names = []
    names.append('X-Indices')
    names.append('Y-Indices')
    
    # create h2o array and save the x and y indices
    dframe = h2o.H2OFrame(gd[1])
    dframe = dframe.concat(h2o.H2OFrame(gd[0]))
    
    # save the band 1 data
    print("[ STATUS ]: Reading band 001")
    dframe = dframe.concat(h2o.H2OFrame(bandArr[gd[0], gd[1]]))
    names.append('Band-001')
    
    # kill references
    bandArr = None
    bandRef = None
    
    # loop through each other band
    if innb > 1: 
        for i in range(2, innb+1):
            print("[ STATUS ]: Reading band %03d" % i)
            bandRef = inRef.GetRasterBand(i)
            bandArr = bandRef.ReadAsArray()
            dframe = dframe.concat(h2o.H2OFrame(bandArr[gd[0], gd[1]]))
            bandArr = None
            bandRef = None
            names.append('Band-%03d' % i)
            
    # that's it for setting up the data frame, so kill reference
    inRef = None
    
    # set the data frame names
    dframe.names = names
    
    # now write the output file
    h2o.export_file(dframe, args.outFile, force=True)
            
    # report finished
    print("[ STATUS ]: ----------")
    print("[ STATUS ]: Finished writing h2o data!")
    print("[ STATUS ]: Please see output file : %s" % args.outFile)
    
# call the aain routine when run from command lne
if __name__ == "__main__":
    main()