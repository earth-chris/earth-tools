#!/usr/bin/python
#####
# aei-align.py
#  
#  aligns and optionally stacks rasters to the same extent and resolution
#  as a reference raster file. can set the resampling options to be unique for
#  each raster to be resampled.
#
# c. 2016 Christopher Anderson
#####

import os
import sys
import aei
import gdal as gdal
import osr as osr
import numpy as np
srs = osr.SpatialReference()

class parse_args:
    def __init__(self, arglist):
        
        # set up main variables and defaults to parse
        self.refFile = ''
        self.inputFiles = []
        self.nInputFiles = 0
        self.outputFiles = []
        self.nOutputFiles = 0
        self.overwrite = False
        self.stack = False
        self.tempResampleFiles = []
        self.tempResampleBase = 'temp_resample_file_'
        self.resampling = []
        self.nResamplingOptions = 0
        self.defaultResampling = 'nearest'
        resamplingOptions = ['nearest', 'bilinear', 'cubic', 'average', 'mode']

        # exit if no arguments passed
        if len(arglist) == 1:
            usage(exit=True)
    
        # read arguments from command line
        i = 1
        while i < len(arglist):
            arg = arglist[i]
        
            # check input data paths            
            if arg.lower() == '-i':
                i += 1
                arg = arglist[i]
                
                # loop through inputs until we find a new parameter
                newArg = False
                while not newArg:
                    if arg[0] == '-':
                        newArg = True
                        i -= 1
                        continue
                    if not aei.fn.checkFile(arg, quiet=True):
                        usage()
                        aei.fn.checkFile(arg)
                        aei.params.sys.exit()
                        
                    self.inputFiles.append(arg)
                    
                    # increment counter and move on
                    i += 1
                    if i >= len(arglist): 
                        newArg = True
                        i -= 1
                        continue
                    arg = arglist[i]
                    
            # parse the base reference file flag    
            elif arg.lower() == '-ref':
                i += 1
                arg = arglist[i]
                
                if not aei.fn.checkFile(arg, quiet = True):
                    usage()
                    aei.fn.checkFile(arg)
                    aei.params.sys.exit(1)
                    
                self.refFile = arg
                    
            # check output flag
            elif arg.lower() == '-o':
                i += 1
                arg = arglist[i]
                
                # loop through outputs until we find a new parameter
                newArg = False
                while not newArg:
                    if arg[0] == '-':
                        newArg = True
                        i -= 1
                        continue
                    
                    self.outputFiles.append(arg)
                    
                    # increment counter and move on
                    i += 1
                    if i >= len(arglist): 
                        newArg = True
                        i -= 1
                        continue
                    arg = arglist[i]
            
            # check resampling options
            elif arg.lower() == "-r":
                i += 1
                arg = arglist[i]
                
                # loop through resampling args until we find a new parameter
                newArg = False
                while not newArg:
                    if arg[0] == '-':
                        newArg = True
                        i -= 1
                        continue
                    if not arg in resamplingOptions:
                        print("[ ERROR ]: Unrecognized resampling option specified: %s" % arg)
                        print("[ ERROR ]: Using default: %s" % self.defaultResampling)
                        arg = self.defaultResampling
                        
                    self.resampling.append(arg)
                    
                    # increment counter and move on
                    i += 1
                    if i >= len(arglist): 
                        newArg = True
                        i -= 1
                        continue
                    arg = arglist[i]
                    
            elif arg.lower() == '-stack':
                self.stack = True
                
            elif arg.lower() == '-overwrite':
                self.overwrite = True
            
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
    
    # count # of input files
    args.nInputFiles = len(args.inputFiles)
    
    # check that its not 0
    if args.nInputFiles == 0:
        usage()
        print("[ ERROR ]: No input files specified.")
        sys.exit(1)
        
    # check there is a reference file
    if len(args.refFile) == 0:
        usage()
        print("[ ERROR ]: No base reference file specified.")
        sys.exit(1)
        
    # make sure enough resampling arguments are passed, otherwise default to nearest neighbor
    args.nResamplingOptions = len(args.resampling)
    
    if args.nResamplingOptions != args.nInputFiles:
        
        # if too few options set, fill in with default
        if args.nResamplingOptions < args.nInputFiles:
            
            print("[ ERROR ]: Too few resampling options set.")
            print("[ ERROR ]: Using default: %s" % args.defaultResampling)
            
            while args.nResamplingOptions < args.nInputFiles:
                args.resampling.append(args.defaultResampling)
                args.nResamplingOptions += 1
                
        # if too many options set, use the first [n input files] specified
        else:
            args.resampling = args.resampling[:args.nInputFiles-1]
            print("[ ERROR ]: Too many resampling options set.")
            print("[ ERROR ]: Using the following: %s" % aei.fn.strJoin(args.resampling, ", "))
    
    # check number of output files
    args.nOutputFiles = len(args.outputFiles)
    if args.nOutputFiles == 0:
        usage()
        print("[ ERROR ]: No output file specified")
        sys.exit(1)
    
    # if -stack is set, only one output should be specified. if -stack is not set, multiple output files must be specified. 
    if args.stack:
        if args.nOutputFiles != 1:
            args.outputFiles = [args.outputFiles[0]]
            args.nOutputFiles = 1
            print("[ ERROR ]: Cannot specify more than one output file with -stack")
            print("[ ERROR ]: Using the first specified output file: %s" % args.outputFiles[0])
        
        # create a series of temp files to hold the resampled outputs prior to stacking
        for i in range(args.nInputFiles):
            args.tempResampleFiles.append(aei.params.scratchdir + aei.params.pathsep + \
                args.tempResampleBase + "%03d.tif" % (i + 1))
            
    else:
        if args.nOutputFiles != args.nInputFiles:
            usage()
            print("[ ERROR ]: Number of input files does not match number of output files")
            sys.exit(1)
    
    # return fixed arguments
    return(args)
            
def usage(exit=False):
    """
    describes the aei-align.py procedure in case of incorrect parameter calls
    
    syntax: usage()
    """
    print(
        """
$ aei-align.py -ref baseReferenceFile -i inputFileList -o outputFileList 
      [-overwrite] [-stack] [-r resamplingOptionList]
        """
        )
    if exit:
        sys.exit(1)

def main():
    """
    the main program for aei-align.py
    
    the order of operations for this procedure is to
    1) read georeferencing data for a base reference file
    2) use gdalwarp to resample the input files to the same extent
       and resolution as the reference file
    3) optionally stack the reference and input files into a single output raster
    
    syntax: main()
    """
    # parse the argument list
    args = parse_args(sys.argv)
    
    # check the argument list to ensure consistent arguments are set
    arge = check_args(args)
    
    # report starting
    print("[ STATUS ]: Starting aei-align.py")
    print("[ STATUS ]: Reference file : %s" % args.refFile)
    print("[ STATUS ]: Input files    : %s" % aei.fn.strJoin(args.inputFiles))
    print("[ STATUS ]: Output files   : %s" % aei.fn.strJoin(args.outputFiles))
    print("[ STATUS ]: Resampling     : %s" % aei.fn.strJoin(args.resampling))
    print("[ STATUS ]: ----------")
    print("[ STATUS ]: Reading georeferencing data")
    
    # open the reference file
    refFile = gdal.Open(args.refFile)
    
    # get the georeferencing information
    xmin, xps, xoff, ymax, yoff, yps = refFile.GetGeoTransform()
    nx = refFile.RasterXSize
    ny = refFile.RasterYSize
    
    # calculate bounding box
    xmax = xmin + xoff + (nx * xps)
    ymin = ymax + yoff + (ny * yps)
    
    # get projection in EPSG
    refProj = refFile.GetProjection()
    srs.ImportFromWkt(refProj)
    refPROJ = '"' + srs.ExportToProj4() + '"'
    
    # report info
    print("[ STATUS ]: Projection     : %s" % refPROJ)
    print("[ STATUS ]: Bounding box   : [%0.2f, %0.2f, %0.2f, %0.2f]" % (xmin, ymin, xmax, ymax))
    print("[ STATUS ]: Pixel size     : %0.2f, %0.2f" % (abs(xps), abs(yps)))
    print("[ STATUS ]: ----------")
    print("[ STATUS ]: Beginning resampling")
    
    # set the reference for output file in resampling. This will be a list of temp files
    #  if -stack is set
    if args.stack:
        outputFiles = args.tempResampleFiles
    else:
        outputFiles = args.outputFiles
        
    # loop through each input and resample
    for i in range(args.nInputFiles):
        
        # report iteration
        print("[ STATUS ]: Resampling image %03d: %s" % (i+1, args.inputFiles[i]))
        
        # call gdalwarp to resample
        aei.cmd.gdalwarp(args.inputFiles[i], outputFiles[i], multi=True, 
            tr=[abs(xps), abs(yps)], te=[xmin, ymin, xmax, ymax],
            t_srs=refPROJ, r=args.resampling[i], overwrite=args.overwrite)
            
    # stack 'em up if set
    if args.stack:
        print("[ STATUS ]: ----------")
        print("[ STATUS ]: Stacking aligned images")
        aei.cmd.otb.ConcatenateImages(args.tempResampleFiles, args.outputFiles[0])
        
        # delete the temporary files
        print("[ STATUS ]: Deleting temporary files")
        for i in range(args.nInputFiles):
            os.remove(args.tempResampleFiles[i])
            
    # report finished
    print("[ STATUS ]: ----------")
    print("[ STATUS ]: Finished aligning raster layers!")
    print("[ STATUS ]: Please see output file(s) : %s" % aei.fn.strJoin(args.outputFiles))
    
# call the aain routine when run from command lne
if __name__ == "__main__":
    main()