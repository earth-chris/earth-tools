#!/usr/bin/python
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
        self.sres = 5 
        
        # tile size in units of output projection
        self.tile_size = 500 
        
        # buffer distance in units of output projection
        self.buffer_size = 50 
        
        # max height above which to clip spurious (noise) points
        self.max_tch = 50
        
        # set canopy cover height threshold
        self.cover_cutoff = 4.0
        
        # default to all cores but 1
        self.cores = aei.params.cores - 1
        
        # default no-data value
        self.nodata = -9999
        
        # set up some processing keywords
        self.city = False
        self.town = False
        self.fine = False
        self.extra_fine = False
        self.coarse = False
        self.extra_coarse = False
        
        ### FILE FORMAT PARAMETERS
        
        # default format for output
        self.lt = '.laz' 
        self.olaz = True
        
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
        self.input_files = ''
        
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
        self.tiles_height = self.scratchdir + aei.params.pathsep + self.dir_height
        
        # default to overwrite any existing files
        self.overwrite = True 
        
        # default to clean up temp files
        self.cleanup = True
        
        ### ASSUME ALL DATA PRODUCTS ON
        
        # set the starting step
        self.start_step = 1
        
        # and all steps are go
        self.step_1 = True
        self.step_2 = True
        self.step_3 = True
        self.step_4 = True
        self.step_5 = True
        self.step_6 = True
        self.step_7 = True
        self.step_8 = True
        self.step_9 = True
        self.step_10 = True
        self.step_11 = True
        self.step_12 = True
        self.step_13 = True
        
        # and create parameters for each step
        self.input_step_1 = self.input_files
        self.input_step_2 = self.input_files
        self.input_step_3 = self.input_files
        self.input_step_4 = self.input_files
        self.input_step_5 = self.input_files
        self.input_step_6 = self.input_files
        self.input_step_7 = self.input_files
        self.input_step_8 = self.input_files
        self.input_step_9 = self.input_files
        self.input_step_10 = self.input_files
        self.input_step_11 = self.input_files
        self.input_step_12 = self.input_files
        self.input_step_13 = self.input_files

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
                    params.input_files = arg
                    
                # if glob fails to find files, try splitting up if there
                #  are spaces between files
                else:
                    split = arg.split()
                    for files in split:
                        if not aei.checkFile(files):
                            aei.params.sys.exit(1)
                            
                    params.input_files = arg
            
            # check output flag
            elif arg.lower() == '-odir':
                i += 1
                arg = arglist[i]
                
                #if not aei.params.os.access(arg, aei.params.os.W_OK):
                #    usage()
                #    print("[ ERROR ]: unable to write to output path: %s" % arg)
                #    aei.params.sys.exit(1)
                #    
                #else: 
                params.odir = arg
            
            # check scratchdir flag
            elif arg.lower() == '-scratchdir':
                i += 1
                arg = arglist[i]
                
                #if not aei.params.os.access(arg, aei.params.os.W_OK):
                #    usage()
                #    print("[ ERROR ]: unable to write to output path: %s" % arg)
                #    aei.params.sys.exit(1)
                #    
                #else: 
                params.scratchdir = arg
            
            # check -name flag
            elif arg.lower() == '-name':
                i += 1
                arg = arglist[i]
                
                params.name = str(arg)
            
            # check if we keep temp files    
            elif arg.lower() == '-keep_temp_files':
                params.cleanup = False
            
            # check the resolution set    
            elif arg.lower() == '-res':
                i += 1
                arg = arglist[i]
                
                try:
                    params.res = float(arg)
                except ValueError:
                    print("[ ERROR ]: unable to set res using arg: %s" % arg)
                    aei.params.sys.exit(1)
            
            # check the slicer resolution set    
            elif arg.lower() == '-sres':
                i += 1
                arg = arglist[i]
                
                try:
                    params.sres = float(arg)
                except ValueError:
                    print("[ ERROR ]: unable to set sres using arg: %s" % arg)
                    aei.params.sys.exit(1)
                
            # check the max tch value set
            elif arg.lower() == '-max_tch':
                i += 1
                arg = arglist[i]
                
                try:
                    params.max_tch = float(arg)
                except ValueError:
                    print("[ ERROR ]: unable to set max_tch using arg: %s" % arg)
                    aei.params.sys.exit(1)
                
            # check the tile size set
            elif arg.lower() == '-tile_size':
                i += 1
                arg = arglist[i]
                
                try:
                    params.tile_size = float(arg)
                except ValueError:
                    print("[ ERROR ]: unable to set tile_size using arg: %s" % arg)
                    aei.params.sys.exit(1)
                
            # check the buffer size set
            elif arg.lower() == '-buffer_size':
                i += 1
                arg = arglist[i]
                
                try:
                    params.buffer_size = float(arg)
                except ValueError:
                    print("[ ERROR ]: unable to set buffer_size using arg: %s" % arg)
                    aei.params.sys.exit(1)
                
            # check the starting step
            elif arg.lower() == '-start_step':
                i += 1
                arg = arglist[i]
                
                try:
                    params.start_step = int(arg)
                except ValueError:
                    print("[ ERROR ]: unable to set start_step using arg: %s" % arg)
                    aei.params.sys.exit(1)
                    
                # turn off all steps before the start step
                for j in range(1, params.start_step):
                    params.__dict__["step_" + str(j)] = False
                
            # check if any outputs are turned off
            elif arg.lower() == '-no_ground':
                params.step_7 = False
                
            elif arg.lower() == '-no_surface':
                params.step_8 = False
                
            elif arg.lower() == '-no_tch':
                params.step_9 = False
                
            elif arg.lower() == '-no_density':
                params.step_10 = False
                
            elif arg.lower() == '-no_slicer':
                params.step_11 = False
                
            elif arg.lower() == '-no_merged':
                params.step_12 = False
                
            elif arg.lower() == '-no_shape':
                params.step_13 = False
                
            elif arg.lower() == '-city':
                params.city = False
                
            elif arg.lower() == '-town':
                params.town = False
                
            elif arg.lower() == '-fine':
                params.fine = False
                
            elif arg.lower() == '-extra_fine':
                params.extra_fine = False
                
            elif arg.lower() == '-coarse':
                params.coarse = True
                
            elif arg.lower() == '-extra_coarse':
                params.extra_coarse = False
                
            elif arg.lower() == '-cores':
                i += 1
                arg = arglist[i]
                params.cores = int(arg)
                
            else:
                usage()
                print("[ ERROR ]: Unrecognized argument: %s" % arg)
                print("%s" % i)
                aei.params.sys.exit(1)
                
            i += 1
        
def update_params(params):
    """
    this procedure takes the input parameters as set and updates
      which files should be processed in which step. these will be
      updated later by various steps
    """
    
    # assume input for each step will be the input files listed. 
    #  these will be updated by the various steps below
    params.input_step_1 = params.input_files
    params.input_step_2 = params.input_files
    params.input_step_3 = params.input_files
    params.input_step_4 = params.input_files
    params.input_step_5 = params.input_files
    params.input_step_6 = params.input_files
    params.input_step_7 = params.input_files
    params.input_step_8 = params.input_files
    params.input_step_9 = params.input_files
    params.input_step_10 = params.input_files
    params.input_step_11 = params.input_files
    params.input_step_12 = params.input_files
    params.input_step_13 = params.input_files
    
    # strip the path separater if there
    if params.outdir[-1] == aei.params.pathsep:
        params.outdir = params.outdir[:-1]
        
    if params.scratchdir[-1] == aei.params.pathsep:
        params.scratchdir = params.scratchdir[:-1]
    
    # update the scratch directory for unclassified tiles
    params.tiles_unclassified = params.scratchdir + aei.params.pathsep + params.dir_unclassified
    
    # update the scratch directory for ground classification
    params.tiles_ground = params.scratchdir + aei.params.pathsep + params.dir_ground
    
    # update the scratch directory for classified tiles
    params.tiles_classified = params.scratchdir + aei.params.pathsep + params.dir_classified
    
    # update the scratch directory for height-normalized tiles
    params.tiles_height = params.scratchdir + aei.params.pathsep + params.dir_height
    
def create_directories(params):
    """
    creates the directories necessary for processing if not already there
    """
    
    # check output directory
    if not aei.params.os.path.exists(params.odir):
        try:
            aei.params.os.makedirs(params.odir)
        except:
            print("[ ERROR ]: Unable to create output directory %s" % params.odir)
            aei.params.sys.exit(1)
            
    # check scratch directory
    if not aei.params.os.path.exists(params.scratchdir):
        try:
            aei.params.os.makedirs(params.scratchdir)
        except:
            print("[ ERROR ]: Unable to create scratch directory %s" % params.scratchdirdir)
            aei.params.sys.exit(1)
            
    # check unclassified
    if not aei.params.os.path.exists(params.tiles_unclassified):
        try:
            aei.params.os.makedirs(params.tiles_unclassified)
        except:
            print("[ ERROR ]: Unable to create unclassified directory %s" % params.tiles_unclassified)
            aei.params.sys.exit(1)
    
    # check ground
    if not aei.params.os.path.exists(params.tiles_ground):
        try:
            aei.params.os.makedirs(params.tiles_ground)
        except:
            print("[ ERROR ]: Unable to create output directory %s" % params.tiles_ground)
            aei.params.sys.exit(1)
    
    # check classified
    if not aei.params.os.path.exists(params.tiles_classified):
        try:
            aei.params.os.makedirs(params.tiles_classified)
        except:
            print("[ ERROR ]: Unable to create output directory %s" % params.tiles_classified)
            aei.params.sys.exit(1)
        
    # check height
    if not aei.params.os.path.exists(params.tiles_height):
        try:
            aei.params.os.makedirs(params.tiles_height)
        except:
            print("[ ERROR ]: Unable to create output directory %s" % params.tiles_height)
            aei.params.sys.exit(1)

def step_1(params):
    """
    merges the input data files
    """
    # report starting
    print("[ STATUS ]: Starting step 1")
    
    # set up the output file
    output_file = params.scratchdir + aei.params.pathsep + \
      params.name + "_merged.laz"
      
    # run the command  
    aei.cmd.lasmerge(params.input_step_1, output_file)
    
    # update the inputs for next steps
    params.input_step_2 = output_file
    
def step_2(params):
    """
    tiles the data
    """
    # report starting
    print("[ STATUS ]: Starting step 2")
    
    # set up the output directory
    output_directory = params.tiles_unclassified
    
    # run the command
    aei.cmd.lastile(params.input_step_2, odir = output_directory,
      tile_size = params.tile_size, buffer = params.buffer_size,
      olaz = params.olaz, cores = params.cores)
      
    # update the inputs for next steps
    params.input_step_3 = output_directory + aei.params.pathsep + "*" + params.lt
    
def step_3(params):
    """
    classifies ground points
    """
    # report starting
    print("[ STATUS ]: Starting step 3")
    
    # set up output directory
    output_directory = params.tiles_ground
    
    # check if step 4 is to be run - if so, we can perform both
    #  in a single command and skip step 4 below
    compute_height = False
    
    if params.step_4:
        compute_height = True
        params.step_4 = False
        print("[ STATUS ]: Starting step 4")
        
    # run the command
    aei.cmd.lasground(params.input_step_3, odir = output_directory,
      olaz = params.olaz, city = params.city, town = params.town,
      compute_height = compute_height, fine = params.fine,
      extra_fine = params.extra_fine, coarse = params.coarse,
      extra_coarse = params.extra_coarse, cores = params.cores)
      
    # update inputs for next steps
    params.input_step_4 = output_directory + aei.params.pathsep + "*" + params.lt
    params.input_step_5 = output_directory + aei.params.pathsep + "*" + params.lt
    
def step_4(params):
    """
    calculates height agl
    """
    # report starting
    print("[ STATUS ]: Starting step 4")
    
    # no need to set up output directory, calculates height in place
    
    # run the command
    aei.cmd.lasheight(params.input_step_4, cores=params.cores)
    
    # update inputs for next steps
    params.input_step_5 = params.input_step_4
    
def step_5(params):
    """
    classifies veg & buildings
    """
    # report starting
    print("[ STATUS ]: Starting step 5")
    
    # set up output directory
    output_directory = params.tiles_classified
    
    # run the command
    aei.cmd.lasclassify(params.input_step_5, odir = output_directory,
      olaz = params.olaz, cores = params.cores)
      
    # update inputs for next steps
    params.input_step_6 = output_directory + aei.params.pathsep + "*" + params.lt
    params.input_step_7 = output_directory + aei.params.pathsep + "*" + params.lt
    params.input_step_8 = output_directory + aei.params.pathsep + "*" + params.lt
    params.input_step_12 = output_directory + aei.params.pathsep + "*" + params.lt
    
def step_6(params):
    """
    height-normalizes the data
    """
    # report starting
    print("[ STATUS ]: Starting step 6")
    
    # set up output directory
    output_directory = params.tiles_height
    
    # run the command
    aei.cmd.lasheight(params.input_step_6, odir = output_directory,
      replace_z = True, olaz = params.olaz, cores=params.cores)
      
    # update inputs for next steps
    params.input_step_9 = output_directory + aei.params.pathsep + "*" + params.lt
    params.input_step_10 = output_directory + aei.params.pathsep + "*" + params.lt
    params.input_step_11 = output_directory + aei.params.pathsep + "*" + params.lt
    
def step_7(params):
    """
    creates ground models and mosaics
    """
    # report starting
    print("[ STATUS ]: Starting step 7")
    
    # set up output directory
    output_directory = params.tiles_classified
    
    # run the command
    aei.cmd.las2dem(params.input_step_7, odir = output_directory,
      step = params.res, nodata = params.nodata, ground = True,
      otif = True, use_tile_bb = True, odix = "_ground",
      cores = params.cores)
      
    # mosaic the individual tiles into an output file
    tiles = output_directory + aei.params.pathsep + "*_ground*.tif"
    mosaic_file = params.odir + aei.params.pathsep + \
      params.name + "_ground.tif"
    
    # run the command  
    aei.cmd.gdalwarp(tiles, mosaic_file, dstnodata=params.nodata,
      n_threads = params.cores, overwrite = params.overwrite)
      
    # no need to update inputs
    
def step_8(params):
    """
    creates surface models and mosaics
    """
    # report starting
    print("[ STATUS ]: Starting step 8")
    
    # set up the output directory
    output_directory = params.tiles_classified
    
    # run the command
    aei.cmd.las2dem(params.input_step_8, odir = output_directory,
      odix = "_surface", step = params.res, nodata = params.nodata,
      ground = False, otif = True, use_tile_bb = True,
      cores = params.cores)
      
    # mosaic the tiles to an output file
    tiles = output_directory + aei.params.pathsep + "*_surface*.tif"
    mosaic_file = params.odir + aei.params.pathsep + \
      params.name + "_surface.tif"
    
    # run the command  
    aei.cmd.gdalwarp(tiles, mosaic_file, dstnodata = params.nodata,
      n_threads = params.cores, overwrite = params.overwrite,
      srcnodata = params.nodata)
      
    # no need to update inputs
    
def step_9(params):
    """
    creates tree-height models and mosaics
    """
    # report starting
    print("[ STATUS ]: Starting step 9")
    
    # set up the output directory
    output_directory = params.tiles_height
    
    # run the command 
    aei.cmd.lascanopy(params.input_step_9, odir = output_directory,
      height_cutoff = 0, step = params.res, max = True, otif = True,
      use_tile_bb = True, cores = params.cores)
      
    # mosaic the tiles to an output file
    tiles = output_directory + aei.params.pathsep + "*_max.tif"
    mosaic_file = params.odir + aei.params.pathsep + \
      params.name + "_tch.tif"
    
    # run the command  
    aei.cmd.gdalwarp(tiles, mosaic_file, dstnodata=params.nodata,
      n_threads = params.cores, overwrite = params.overwrite)
      
    # no need to update inputs
    
def step_10(params):
    """
    creates canopy density models and mosaics
    """
    # report starting
    print("[ STATUS ]: Starting step 10")
    
    # set up the output directory
    output_directory = params.tiles_height
    
    # run the command 
    aei.cmd.lascanopy(params.input_step_9, odir = output_directory,
      height_cutoff = 0, step = params.res, dns = True, 
      otif = True, cover_cutoff = params.cover_cutoff,
      use_tile_bb = True, cores = params.cores)
      
    # mosaic the tiles to an output file
    tiles = output_directory + aei.params.pathsep + "*_dns.tif"
    mosaic_file = params.odir + aei.params.pathsep + \
      params.name + "_dns.tif"
    
    # run the command  
    aei.cmd.gdalwarp(tiles, mosaic_file, dstnodata=params.nodata,
      n_threads = params.cores, overwrite = params.overwrite)
      
    # no need to update inputs
    
def step_11(params):
    """
    creates slicer models and mosaics
    """
    import glob
    
    # report starting
    print("[ STATUS ]: Starting step 11")
    
    # set up the output directory
    output_directory = params.tiles_height
    
    # run the command 
    aei.cmd.lascanopy(params.input_step_9, odir = output_directory,
      height_cutoff = 0, step = params.res, 
      density = range(0,int(params.max_tch+1)), 
      use_tile_bb = True, otif = True, cores = params.cores)
      
    # mosaic the tiles to output files
    for i in range(0, int(params.max_tch+1)):    
        tiles = str(output_directory + aei.params.pathsep + \
          "*_d%02d.tif" % i)
        mosaic_file = str(params.odir + aei.params.pathsep + \
          params.name + "_slicer_%02d.tif" % i)
    
        # run the command  
        aei.cmd.gdalwarp(tiles, mosaic_file, dstnodata=params.nodata,
          n_threads = params.cores, overwrite = params.overwrite)
    
    # stack the slicer bands into a single output file
    individual_files = params.odir + aei.params.pathsep + \
      params.name + "_slicer*.tif"
    file_list = glob.glob(individual_files)
    
    vrt_file = params.odir + aei.params.pathsep + \
      params.name + "_slicer.vrt"
    output_stack = params.odir + aei.params.pathsep + \
      params.name + "_slicer.tif"
    
    # build a vrt
    aei.cmd.gdalbuildvrt(individual_files, vrt_file, separate=True)
    
    # translate to a stacked tif
    aei.cmd.gdal_translate(vrt_file, output_stack)
    
    # clean up individual files
    aei.params.os.remove(vrt_file)
    for file in file_list:
        aei.params.os.remove(file)
      
    # no need to update inputs
    
def step_12(params):
    """
    creats a final classified las file
    """
    # report starting
    print("[ STATUS ]: Starting step 12")
    
    # set up output file
    output_file = params.odir + aei.params.pathsep + \
      params.name + "_merged" + params.lt
      
    # run the command
    aei.cmd.lasmerge(params.input_step_12, output_file,
      rescale = False, olaz = params.olaz)
      
    # index the files
    aei.cmd.lasindex(output_file)
      
    # update inputs
    params.input_step_13 = output_file
    
def step_13(params):
    """
    creates a bounding shape file
    """
    # report starting
    print("[ STATUS ]: Starting step 13")
    
    # set up output file
    output_file = params.odir + aei.params.pathsep + \
    params.name + "_boundary.shp"
    
    # run the command
    aei.cmd.lasboundary(params.input_step_13, output = output_file,
      disjoint = True, cores=params.cores)
      
    # no need to update inputs
    
def cleanup_files(params):
    """
    cleans up temp files and directories
    """
    import shutil
    
    # delete unclassified
    shutil.rmtree(params.tiles_unclassified)
    
    # update the scratch directory for ground classification
    shutil.rmtree(params.tiles_ground)
    
    # update the scratch directory for classified tiles
    shutil.rmtree(params.tiles_classified)
    
    # update the scratch directory for height-normalized tiles
    shutil.rmtree(params.tiles_height)
    
def usage(exit=False):
    """
    describes the aeilas.py procedure in case of incorrect parameter calls
    
    syntax: usage(exit=False)
    """

    print(
        """
$ aeilas.py -i input_files [-odir output_directory]
  [-scratchdir scratch_directory] [-name output_basename]
  [-keep_temp_files] [-res resolution] [-sres slicer_resolution] 
  [-max_tch max_tch] [-tile_size tile_size] [-buffer_size buffer_size]
  [-start_step step] [-no_ground] [-no_surface] [-no_tch]
  [-no_merged] [-no_shape] [-no_slicer] [-no_density] [-cores cores]
  [-city] [-town] [-fine] [-extra_fine] [-coarse] [-extra_coarse]
        """
        )
    
    if exit:    
        aei.params.sys.exit(1)
        
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
    12. create a merged laz file
    13. create shape file
    
    syntax: main()
    """
    
    # first, read the default parameters to get the processing object
    params = default_params()
    
    # then parse the arguments passed via command line
    args = aei.params.sys.argv
    parse_args(args, params)
    
    # check that input files were specified
    if not params.input_files:
        print("[ ERROR ]: No input files specified")
        aei.params.sys.exit(1)
    
    # update the params class for any dependencies set at runtime
    update_params(params)
    
    # create the output directories necessary 
    create_directories(params)
    
    # report that we're starting, and the parameters used
    print("[ STATUS ]: Beginning aeilas.py")
    print("[ STATUS ]: Parameters set:")
    print("[ STATUS ]: Site name        : %s" % params.name)
    print("[ STATUS ]: Output directory : %s" % params.odir)
    
    # based on these params, process each step
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
    
    if params.step_12:
        step_12(params)
        
    if params.step_13:
        step_13(params)
        
    # clean up, if set
    if params.cleanup:
        cleanup_files(params)
        
    # report finished up    
    print("[ STATUS ]: aeilas.py processing complete")
    print("[ STATUS ]: Output project  : %s" % params.name)
    print("[ STATUS ]: Output directory: %s" % params.odir)

if __name__ == "__main__":
    main()