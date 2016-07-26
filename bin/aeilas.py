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
    'cores' : multiprocessing.cpu_count()-1, # default to all cores but 1
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

class parse_args:
    def __init__(self, arglist):
        
        # set up main variables and defaults to parse
        self.infile = ''
        self.outfile = ''
        self.spectral_libs = []
        self.n = 20
        self.bands = []
        self.of = 'GTiff'
        self.normalize = False
        
        # exit if no arguments passed
        if len(arglist) == 1:
            usage(exit=True)

        # read arguments from command line
        i = 1
        while i < len(arglist):
            arg = arglist[i]
        
            # check input flag    
            if arg.lower() == '-i':
                i += 1
                arg = arglist[i]
                
                if type(arg) is str:
                    self.infile = arg
                    if not aei.checkFile(self.infile, quiet = True):
                        usage()
                        aei.checkFile(self.infile)
                        sys.exit(1)
            
            # check output flag
            if arg.lower() == '-o':
                i += 1
                arg = arglist[i]
                
                if type(arg) is str:
                    self.outfile = arg
                    outpath = os.path.dirname(self.outfile)
                    if outpath == '':
                        outpath = '.'
                    if not os.access(outpath, os.W_OK):
                        usage()
                        print("[ ERROR ]: unable to write to output path: %s" % outpath)
                        sys.exit(1)
    
def usage(exit=False):
    """
    describes the aeilas.py procedure in case of incorrect parameter calls
    
    syntax: usage(exit=False)
    """

    print(
        """
        $ aeilas.py [-cores cores] [-lt las_type][-max_tch max_tch] 
        [-o output_file] [-overwrite] [-proj proj4string]
        [-res resolution] [-sres slicer_resolution] [-tile_size tile_size] 
        -i input_files -o output_file
        
        
        """
        )
    
    if exit:    
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
            usage()
            print("Could not parse -cores argument: %s" % (arg))
        
    elif arg.lower() == '-lt':
        i += 1
        arg = sys.argv[i]
        try:
            if (arg.lower() != 'las') or (arg.lower() != 'laz'):
                params['lt'] = arg
        except ValueError:
            usage()
            print("Could not parse -lt argument: %s" % (arg))
            
    elif arg.lower() == '-max_tch':
        i += 1
        arg = sys.argv[i]
        try:
            if arg > 0:
                params['max_tch'] = arg
        except ValueError:
            usage()
            print("Could not parse -max_tch argument: %s" % (arg))
            
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

#################################
::
:: aei_merge_las
::
:: windows batch script for converting to UTM coordinates, setting xyz precision, 
::  then merging all las files into a single laz file
::

:: if you are resuming this batch from a failed point, uncomment the goto option you would like
::
::goto change_projection
::goto merge_to_unclassified
::goto remove_noise
::goto tile_unclassified
::goto classify_ground
::goto classify_tiles
::goto normalize_height
::goto remove_buffers
::goto create_ground_model
::goto create_surface_model
::goto create_tch
::goto create_slicer

:: set your working path, number of cores
set PATH=%PATH%;C:\software\lastools\bin;
set N_CORES=1

:: specify UTM zone (eg: 10S for CA, 18M for northern Peru, auto for auto-determination), 
::  resolution (m), buffer distance (m), tile size (m), max height (m)
set ZONE=auto
set RES=1
set SLICER_RES=10
set BUFFER=20
set TILE_SIZE=500
set MAX_HEIGHT=60

:: specify your processing, temp, and output directories (can be relative or absolute). Site should
::  change for each batch to ensure no overlap of files
set SITE=big_basin
set FORMAT=laz
set INPUT_DIR=G:\AEI\Data\raw\lidar\california\big_basin
set OUTPUT_DIR=G:\AEI\Data\processed\%SITE%
set PROCESSING_DIR=G:\AEI\Data\scratch\%SITE%

:: make your directories
mkdir %PROCESSING_DIR%
mkdir %PROCESSING_DIR%\tmp_reprojected
mkdir %PROCESSING_DIR%\tiles_raw
mkdir %PROCESSING_DIR%\tiles_ground
mkdir %PROCESSING_DIR%\tiles_classified
mkdir %PROCESSING_DIR%\tiles_normalized 
mkdir %PROCESSING_DIR%\tiles_final

:change_projection
:: change the projection of your las files
las2las -i %INPUT_DIR%\*.%FORMAT% -olaz -odir %PROCESSING_DIR%\tmp_reprojected -target_utm %ZONE% -set_classification 0

:merge_to_unclassified
:: merge to single unclassified laz file
lasmerge -i %PROCESSING_DIR%\tmp_reprojected\*.laz -o %PROCESSING_DIR%\%SITE%_unclassified.laz -rescale 0.01 0.01 0.001

:: delete your temp reclassification products
::rmdir %PROCESSING_DIR%\%SITE%\tmp_reprojected /s /q

:remove_noise
:: automatically remove noise from unclassified laz file
lasnoise -i %PROCESSING_DIR%\%SITE%_unclassified.laz -o %PROCESSING_DIR%\%SITE%_unclassified_denoised.laz -step 3 -isolated 2

:tile_unclassified
:: tile out your unclassified laz file
lastile -i %PROCESSING_DIR%\%SITE%_unclassified_denoised.laz -olaz -odir %PROCESSING_DIR%\tiles_raw -buffer %BUFFER% -tile_size %TILE_SIZE% -reversible -extra_pass

:classify_ground
:: classify the ground points in each of the tiles. compute height AGL on the fly so we don't need lasheight
lasground -i %PROCESSING_DIR%\tiles_raw\*.laz -olaz -odir %PROCESSING_DIR%\tiles_ground -town -fine -cores %N_CORES% -compute_height

:classify_tiles
:: classify your height normalized tiles
lasclassify -i %PROCESSING_DIR%\tiles_ground\*.laz -olaz -odir %PROCESSING_DIR%\tiles_classified -cores %N_CORES%

:normalize_height
:: normalize height to m AGL for lascanopy inputs
lasheight -i %PROCESSING_DIR%\tiles_classified\*.laz -olaz -odir %PROCESSING_DIR%\tiles_normalized -replace_z -cores %N_CORES%

:remove_buffers
:: remove buffers from your tiles
lastile -i %PROCESSING_DIR%\tiles_classsified\*.laz -olaz -odir %PROCESSING_DIR%\tiles_final -remove_buffer

::
:: begin generating raster products
::

:create_ground_model
:: create ground models
las2dem -i %PROCESSING_DIR%\tiles_final\*.laz -obil -odir %PROCESSING_DIR%\tiles_final -odix _ground -keep_class 2 -step %RES% -use_tile_bb -extra_pass
gdalwarp -of GTiff %PROCESSING_DIR%\tiles_final\*_ground.bil %OUTPUT_DIR%\%SITE%_ground.tif -co "COMPRESS=LZW" -multi -dstnodata -9999 -s_srs "+proj=utm +zone=10 +north +datum=nad83" -t_srs "+proj=utm +zone=10 +north +datum=wgs84"

:create_surface_model
:: create surface model
las2dem -i %PROCESSING_DIR%\tiles_final\*.laz -obil -odir %PROCESSING_DIR%\tiles_final -odix _surface -first_only -step %RES% -use_tile_bb -extra_pass
gdalwarp -of GTiff %PROCESSING_DIR%\tiles_final\*_surface.bil %OUTPUT_DIR%\%SITE%_surface.tif -co "COMPRESS=LZW" -multi -dstnodata -9999 -s_srs "+proj=utm +zone=10 +north +datum=nad83" -t_srs "+proj=utm +zone=10 +north +datum=wgs84"

:create_tch
:: create tch model
lascanopy -i %PROCESSING_DIR%\tiles_normalized\*.laz -obil -odir %PROCESSING_DIR%\tiles_final -odix _tch -step %RES% -max
gdalwarp -of GTiff %PROCESSING_DIR%\tiles_final\*_tch.bil %OUTPUT_DIR%\%SITE%_tch.tif -co "COMPRESS=LZW" -multi -dstnodata -9999 -s_srs "+proj=utm +zone=10 +north +datum=nad83" -t_srs "+proj=utm +zone=10 +north +datum=wgs84"

:create_slicer
:: create the slice layers
lascanopy -i %PROCESSING_DIR%\tiles_normalized\*.laz -obil -odir %PROCESSING_DIR%\tiles_final -odix _slice -step %SLICER_RES% -d 0.5 2 4 6 8 10 12 14 16 18 20 22 24 26 28 30 32 34 36 38 40

#################################
        
if __name__ == "__main__":
    main()