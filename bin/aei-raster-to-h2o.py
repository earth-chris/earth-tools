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
    
    # open the reference file
    inRef = gdal.Open(args.inFile)
    
    # get metadata info
    proj = inRef.GetProjection()
    geot = inRef.GetGeoTransform()
    innx = inRef.RasterXSize
    inny = inRef.RasterYSize
    innb = inRef.RasterCount
    
    # open the output file
    with hdf.File(args.outFile, 'w') as outf:
        
        # create a group of data with georeferencing info
        g1 = outf.create_group('GeoData')
        g1.create_dataset('Projection', data = proj)
        g1.create_dataset('GeoTransform', data = geot)
        g1.create_dataset('[nx, ny, nb]', data = [innx, inny, innb])
        
        # create a new group to store raster data
        g2 = outf.create_group('RasterData')
        
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
        
        # save the x and y indices
        g2.create_dataset("Y-Indices", data = gd[0])
        g2.create_dataset("X-Indices", data = gd[1])
        
        # save the band 1 data
        print("[ STATUS ]: Writing band 001")
        g2.create_dataset("Band-001", data = bandArr[gd[0], gd[1]])
        
        # kill references
        bandArr = None
        bandRef = None
        
    # then read the full data set into an array
    fileArr = inRef.ReadAsArray
    
    # subset the data
    fileArr = fileArr[:, gd[0], gd[1]]
    
    # create a header for the file
    header = []
    for i in range(1, innb+1):
        header.append("Band-%03d" % i)
    
    # then transpose it and write as a csv
    np.savetxt(args.outFile, fileArr.transpose(), delimiter = ',', 
        header = aei.fn.strJoin(header))
                
    # that's it for writing to the hdf file, so kill reference
    inRef = None
            
    # report finished
    print("[ STATUS ]: ----------")
    print("[ STATUS ]: Finished writing hdf5 data!")
    print("[ STATUS ]: Please see output file : %s" % args.outFile)
    
# call the main routine when run from command lne
if __name__ == "__main__":
    main()