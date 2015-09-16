#!/bin/env/python2.7
##
## full processing stream for las/laz files to raster data products
##
####################
##
import aei as aei
import liblas as ll
import gdal as gdal
import mutliprocessing

# define global parameters
params = {
    'buffer': 20, # buffer distance in units of output projection
    'cores' : multiprocessing.cpu_count()/-1, # default to all cores but 1
    'lt' : 'laz', # las format for output
    'ot' : 'Float32', # raster output format
    'max_tch' : 65, # max height above which to clip spurious (noise) points
    'outdir' : aei.get_outdir(), # the output directory for final products
    'outfile' : 'default', # the default output basename
    'overwrite' : False, # default to not overwrite any existing files
    'procdir' : aei.get_procdir(), # the output directory for scratch files
    'res' : 1, # default raster resolution
    'slicer_res' : 5, # slicer resolution (usually 5x raster res)
    'tile_size' : 500 # tile size in units of output projection
    }

def usage():
    """
    describes the aeilas.py procedure in case of incorrect parameter calls
    
    syntax: usage()
    """

    print(
        """
        $ aeilas.py [-cores cores] [-lt las_type][-max_tch max_tch] 
        [-o output_file] [-overwrite] [-proj proj4string]
        [-res resolution] [-sres slicer_resolution] [-tile_size tile_size] 
        input_files
        
        
        """
        )
        
    sys.exit(1)

step = 0
def resume():
    """
    resumes the aeilas.py procedure when flagged at runtime
    
    syntax: resume()
    """
    print(steps)
    while True:
        try:
            step=int(raw_input('Select the step to resume: '))
            break
        except ValueError:
            print('Invalid entry. Please enter number from 1 to %d' % (len(steps)))
        
def main ():
    """
    the main program for aeilas.py
    
    syntax: main()
    """
    
    # read arguments from command line
    i = 1
    while i < len(sys.argv):
        arg = sys.argv[i]
        
    if arg.lower() == '-cores':
        i += 1
        arg = sys.argv[i]
        try:
            if arg < multiprocessing.cpu_count():
                params['cores'] = arg
        except ValueError:
            print("Could not parse -cores argument: %s" % (arg))
            usage()
        
    elif arg.lower() == '-lt':
        i += 1
        arg = sys.argv[i]
        try:
            if (arg.lower() != 'las') or (arg.lower() != 'laz'):
                params['lt'] = arg
        except ValueError:
            print("Could not parse -lt argument: %s" % (arg))
            usage()
            
    elif arg.lower() == '-max_tch':
        i += 1
        arg = sys.argv[i]
        try:
            if arg > 0:
                params['max_tch'] = arg
        except ValueError:
            print("Could not parse -max_tch argument: %s" % (arg))
            usage()
            
    elif arg.lower() == '-o':
        i += 1
        arg = sys.argv[i]
        
    elif arg.lower() == '-overwrite':
        params['overwrite'] = True
        
    elif arg.lower() == '-proj':
        i += 1
        arg = sys.argv[i]
        
    elif arg.lower() == '-res':
        i += 1
        arg = sys.argv[i]
        try:
            if arg > 0:
                params['res'] = arg
        except ValueError:
            usage()
            print('Could not parse -res argument: %s' % ( arg))
        
    elif arg.lower() == '-sres':
        i += 1
        arg = sys.argv[i]
        try:
            if arg > 0:
                params['sres'] = arg
        except ValueError:
            usage()
            print('Could not parse -sres argument: %s' % ( arg))
        
    elif arg.lower() == '-tile_size':
        i += 1
        arg = sys.argv[i]
        try:
            if arg > 0:
                params['tile_size'] = arg
        except ValueError:
            usage()
            print('Could not parse -tile_size argument: %s' % ( arg))
        
if __name__ == "__main__":
    main()