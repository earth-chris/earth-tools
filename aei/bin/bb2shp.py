#!/usr/bin/python
##
## converts the bounding box of a raster file to a shape file
##
####################

import aei as aei

# define global parameters
class default_params:
    def __init__(self):
        
        import datetime
        now = datetime.datetime.now()
        
        ### OUTPUT / PROCESSING PARAMETERS
        
        # set the input file parameter
        self.inputFiles = None
        
        # set the output file parameters (if single file)
        # logic for output files is:
        #  if single input file: single output file
        #    if -o set: use -o parameter
        #    if not set: use base name + dir of input file
        #  if multi input file: 
        #    if -o set: single output file
        #    if -odir set: multi output file
        #    if not specified: multi output file in directory of input files
        self.outputFiles = None
        self.outputDir = None
        self.singleOutput = True
        
        # set this to calculate raster statistics and include in attribute table.
        self.getStats = False

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
                
                # if glob finds files, use what was set in command line
                if len(glob.glob(arg)) > 0:
                    params.inputFiles = arg
                
                # if glob fails to find files, try splitting up if there
                #  are spaces between files
                else:
                    split = arg.split()
                    for files in split:
                        if not aei.checkFile(files):
                            aei.params.sys.exit(1)
                            
                    params.inputFiles = arg
                    
            # check output dir flag
            elif arg.lower() == '-odir':
                i += 1
                arg = arglist[i]
                
                params.outputDir = arg
                params.singleOutput = False
                
            # check single output file flag
            elif arg.lower() == '-o':
                i += 1
                arg = arglist[i]
                
                params.outputFile = arg
                params.singleOutput = True
                
            elif arg.lower() == '-get_stats':
                params.getStats = True
                
            else:
                usage()
                print("[ ERROR ]: Unrecognized argument: %s" % arg)
                aei.params.sys.exit(1)
                
            i += 1
                    
# updates params variable to reflect updates in processing options
def update_params(params):
    """
    this procedure takes the input parameters as set at 
      command line and updates within script
    """
    import gdal as gdal
    
    # go through logic of input and output files and update parameters
    #  check how many input files are specified. If more than one...
    if len(params.inputFiles) > 1:
        # check if single output
        if params.singleOutput:
            
        # or multi output
        else:
            # check if output dir specified
            if params.outputDir:
                
            # otherwise use from each input file
            else:
                
    # if only one input file    
    else:
        # if -o flag set, nothing needs to be updated
        if params.singleOutput:
            pass
        
        # otherwise, set output file to match the input file
        else:
            # use gdal to get the file name
            gdalRef = gdal.Open(params.inputFiles)
            params.outputFiles = gdalRef.getName()
            
            # release gdal reference
            gdalRef = None

# set usage string    
def usage(exit=False):
    """
    describes the bb2shp.py procedure in case of incorrect parameter calls
    
    syntax: usage(exit=False)
    """

    print(
        """
$ bb2shp.py -i input_files [-o output file] [-odir output_directory]
    [-of output_format] [-get_stats]
        """
        )
    
    if exit:    
        aei.params.sys.exit(1)
        
def main ():
    """
    the main program for bb2shp.py
    
    syntax: main()
    """
    
    # first, read the default parameters to get the processing object
    params = default_params()
    
    # then parse the arguments passed via command line
    args = aei.params.sys.argv
    parse_args(args, params)
    
    # check that input files were specified
    if not params.inputFiles:
        print("[ ERROR ]: No input file(s) specified")
        aei.params.sys.exit(1)
    
    # update the params class for any dependencies set at runtime
    update_params(params)
    
if __name__ == "__main__":
    main()