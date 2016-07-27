#!/bin/env/python2.7
##
## full processing stream for las/laz files to raster data products
##
####################
##
import aei as aei

# define global parameters
class default_params:
    def __init__(self):
        
        import datetime
        now = datetime.datetime.now()
        
        ### OUTPUT / PROCESSING PARAMETERS
        
        # default raster resolution
        self.res = 1 
        
        # slicer resolution (usually 5x raster res)
        self.slicer_res = 5 
        
        # tile size in units of output projection
        self.tile_size = 500 
        
        # buffer distance in units of output projection
        self.buffer = 50 
        
        # max height above which to clip spurious (noise) points
        self.max_tch = 50
        
        # default to all cores but 1
        self.cores = aei.params.cores - 1
        
        ### FILE FORMAT PARAMETERS
        
        # default format for output
        self.lt = 'laz' 
        
        # raster output format
        self.ot = 'Float32' 
        
        ### DEFAULT DATA PRODUCTS
        
        # output ground model
        self.ground = True
        
        # output surface model
        self.surface = True
        
        # output slicer
        self.slicer = True
        
        # output shape file
        self.shape = True
        
        # output merged las file
        self.merged = True
        
        ### OUTPUT DIRECTORIES AN DEFAULT FILE NAMES
        
        # the input file list
        self.input_files = []
        
        # the site name that will have everything appended to it e.g. site_tch.tif, site_merged.laz
        self.sitename = 'aeilas'
        
        # the output directory for final products
        #  default is to use outdir/YYYY-MM-DD_HH-MM/ if not set
        datestring = ("%04d-%02d-%02d_%02d%02d" % (now.year, now.month, now.day, now.hour, now.minute))
        self.outdir = aei.params.outdir + aei.params.pathsep + \
          datestring + aei.params.pathsep
        
        # the output directory for intermediate products
        self.scratchdir = aei.params.scratchdir
        
        # set directory names for temp outputs
        self.dir_unclassified = 'unclassified'
        self.dir_ground = 'ground'
        self.dir_classified = 'classified'
        self.dir_height = 'height_normalized'
        
        # the scratch directory for unclassified tiles
        self.tiles_unclassified = self.scratchdir + aei.params.pathsep + self.dir_unclassified
        
        # the scratch directory for ground classification
        self.tiles_ground = self.scratchdir + aei.params.pathsep + self.dir_ground
        
        # the scratch directory for classified tiles
        self.tiles_classified = self.scratchdir + aei.params.pathsep + self.dir_classified
        
        # scratch directory for height-normalized tiles
        self.tiles_neight = self.scratchdir + aei.params.pathsep + self.dir_height
        
        # default to overwrite any existing files
        self.overwrite = True 
        
        # default to clean up temp files
        self.cleanup = True

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
        while i < len(arglist):
            arg = arglist[i]
        
            # check input flag    
            if arg.lower() == '-i':
                i += 1
                arg = arglist[i]
                
                # if glob finds files, use what was set in command line
                if len(glob.glob(arg)) > 0:
                    params.input_files = arg
                    
                # if glob fails to find files, try splitting up if there
                #  are spaces between files
                else:
                    split = arg.split()
                    for files in split:
                        if not aei.checkFile(files):
                            sys.exit(1)
                            
                    params.input_files = arg
                
                if type(arg) is str:
                    self.infile = arg
                    if not aei.checkFile(self.infile, quiet = True):
                        usage()
                        aei.checkFile(self.infile)
                        sys.exit(1)
            
            # check output flag
            if arg.lower() == '-odir':
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
            
            if arg.lower() == '-scratchdir':
            
            if arg.lower() == '-name':
            
            if arg.lower() == '-keep_temp_files':
                
            if arg.lower() == '-res':
                
            if arg.lower() == '-sres':
                
            if arg.lower() == '-max_tch':
                
            if arg.lower() == '-tile_size':
                
            if arg.lower() == '-buffer_size':
                
            if arg.lower() == '-start_step':
                
            if arg.lower() == '-no_ground':
                
            if arg.lower() == '-no_surface':
                
            if arg.lower() == '-no_tch':
                
            if arg.lower() == '-no_merged':
                
            if arg.lower() == '-no_shape':
                
            if arg.lower() == '-no_slicer':
                
            if arg.lower() == '-no_density':
                
            
            i += 1
            
        return params
    
def usage(exit=False):
    """
    describes the aeilas.py procedure in case of incorrect parameter calls
    
    syntax: usage(exit=False)
    """

    print(
        """
        $ aeilas.py -i input_files [-odir output_directory]
          [-scratchdir scratch_directory] [-name output_basename]
          [-cores cores] [-max_tch max_tch] 
          [-res resolution] [-sres slicer_resolution] [-tile_size tile_size] 
          []
        """
        )
    
    if exit:    
        sys.exit(1)
        
def main ():
    """
    the main program for aeilas.py
    the order of operations for las processing is as follows:
    
    01. merge the input data files (assumes multiple *.las or *.laz files set)
    02. tile the data
    03. classify ground points
    04. calculate height agl (in place)
    05. classify the tiles
    06. create height-normalized tiles
    07. create ground models and mosaic
    08. create surface models and mosaic
    09. create tree-height models and mosaic
    10. create canopy density models
    11. create slicer models
    
    syntax: main()
    """
    
    # first, read the default parameters to get the processing object
    default_params = default_params()
    
    # then parse the arguments passed via command line
    args = sys.argv
    params = parse_args(args, default_params)
    
    # check that input files were specified
    if not params.input_files:
        print("[ ERROR ]: No input files specified")
        sys.exit(1)
    
    # report that we're starting, and the parameters used
    print("[ STATUS ]: Beginning aeilas.py")
    print("[ STATUS ]: Parameters set:")
    print("[ STATUS ]: Site name        : %s" % params.name)
    print("[ STATUS ]: Output directory : %s" % params.odir)
    
    # based on these params, proces through each step
    if params.step_1:
        step_1(params)
        
    if params.step_2:
        step_2(params)
    
    if params.step_3:
        step_3(params)
        
    if params.step_4:
        step_4(params)
        
    if params.step_5:
        step_5(params)
        
    if params.step_6:
        step_6(params)
    
    if params.step_7:
        step_7(params)
        
    if params.step_8:
        step_8(params)
        
    if params.step_9:
        step_9(params)
        
    if params.step_10:
        step_10(params)
    
    if params.step_11:
        step_11(params)

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