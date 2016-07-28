#####
# sets up paths and calls to external executables.
# c. 2016 Christopher Anderson
#####

# set up commands to call binaries
class getPaths:
    def __init__(self):
        from aei import params
        
        # lastools commands
        self.las2dem = ' '.join({
            'Windows' : [params.pathcat(params.ltbase, 'las2dem.exe')], 
            'Linux'   : ['wine', params.pathcat(params.ltbase, 'las2dem.exe')]
            }[params.system])
        self.las2las = ' '.join({
            'Windows' : [params.pathcat(params.ltbase, 'las2las.exe')], 
            'Linux'   : ['wine', params.pathcat(params.ltbase, 'lasmerge.exe')]
            }[params.system])
        self.lasboundary = ' '.join({
            'Windows' : [params.pathcat(params.ltbase, 'lasboundary.exe')], 
            'Linux'   : ['wine', params.pathcat(params.ltbase, 'lasboundary.exe')]
            }[params.system])
        self.lascanopy = ' '.join({
            'Windows' : [params.pathcat(params.ltbase, 'lascanopy.exe')], 
            'Linux'   : ['wine', params.pathcat(params.ltbase, 'lascanopy.exe')]
            }[params.system])
        self.lasclassify = ' '.join({
            'Windows' : [params.pathcat(params.ltbase, 'lasclassify.exe')],
            'Linux'   : ['wine', params.pathcat(params.ltbase, 'lasclassify.exe')]
            }[params.system])
        self.lasclip = ' '.join({
            'Windows' : [params.pathcat(params.ltbase, 'lasclip.exe')], 
            'Linux'   : ['wine', params.pathcat(params.ltbase, 'lasclip.exe')]
            }[params.system])
        self.lasgrid = ' '.join({
            'Windows' : [params.pathcat(params.ltbase, 'lasgrid.exe')], 
            'Linux'   : ['wine', params.pathcat(params.ltbase, 'lasgrid.exe')]
            }[params.system])
        self.lasground = ' '.join({
            'Windows' : [params.pathcat(params.ltbase, 'lasground.exe')], 
            'Linux'   : ['wine', params.pathcat(params.ltbase, 'lasground.exe')]
            }[params.system])
        self.lasheight = ' '.join({
            'Windows' : [params.pathcat(params.ltbase, 'lasheight.exe')], 
            'Linux'   : ['wine', params.pathcat(params.ltbase, 'lasheight.exe')]
            }[params.system])
        self.lasinfo = ' '.join({
            'Windows' : [params.pathcat(params.ltbase, 'lasinfo.exe')], 
            'Linux'   : ['wine', params.pathcat(params.ltbase, 'lasinfo.exe')]
            }[params.system])
        self.lasindex = ' '.join({
            'Windows' : [params.pathcat(params.ltbase, 'lasindex.exe')], 
            'Linux'   : ['wine', params.pathcat(params.ltbase, 'lasindex.exe')]
            }[params.system])
        self.lasmerge = ' '.join({
            'Windows' : [params.pathcat(params.ltbase, 'lasmerge.exe')], 
            'Linux'   : ['wine', params.pathcat(params.ltbase, 'lasmerge.exe')]
            }[params.system])
        self.lasnoise = ' '.join({
            'Windows' : [params.pathcat(params.ltbase, 'lasnoise.exe')],
            'Linux'   : ['wine', params.pathcat(params.ltbase, 'lasnoise.exe')]
            }[params.system])
        self.lastile = ' '.join({
            'Windows' : [params.pathcat(params.ltbase, 'lastile.exe')],
            'Linux'   : ['wine', params.pathcat(params.ltbase, 'lastile.exe')]
            }[params.system])
        self.laszip = ' '.join({
            'Windows' : [params.pathcat(params.ltbase, 'laszip.exe')], 
            'Linux'   : ['wine', params.pathcat(params.ltbase, 'laszip.exe')]
            }[params.system])
        
        # ogr/gdal commands
        self.ogr2ogr = ' '.join({
            'Windows' : [params.pathcat(params.gdalbase, 'ogr2ogr.exe')],
            'Linux'   : [params.pathcat(params.gdalbase, 'ogr2ogr')]
            }[params.system])
        self.gdalbuildvrt = ' '.join({
            'Windows' : [params.pathcat(params.gdalbase, 'gdalbuildvrt.exe')],
            'Linux'   : [params.pathcat(params.gdalbase, 'gdalbuildvrt')]
            }[params.system])
        self.gdalwarp = ' '.join({
            'Windows' : [params.pathcat(params.gdalbase, 'gdalwarp.exe')],
            'Linux'   : [params.pathcat(params.gdalbase, 'gdalwarp')]
            }[params.system])
        self.gdal_translate = ' '.join({
            'Windows' : [params.pathcat(params.gdalbase, 'gdal_translate.exe')],
            'Linux'   : [params.pathcat(params.gdalbase, 'gdal_translate')]
            }[params.system])
        self.gdal_rasterize = ' '.join({
            'Windows' : [params.pathcat(params.gdalbase, 'gdal_rasterize.exe')],
            'Linux'   : [params.pathcat(params.gdalbase, 'gdal_rasterize')]
            }[params.system])

raw = getPaths()

###
# GDAL commands
###

# gdalwarp
def gdalwarp(inputs='', output='', etc='', dstnodata=False, 
        multi=False, n_threads=False, overwrite=False,
        of=False, ot=False, r='', srcnodata=False,
        s_srs='', t_srs='', te=[], tr=[]):
    """
    mosaics, reprojects or warps raster imagery
    
    syntax: gdalwarp(inputs, output, etc=etc, dstnodata=dstnodata, multi=multi, n_threads=n_threads, overwrite=overwrite, of=of, ot=ot, r=r, srcnodata=srcnodata, s_srs=s_srs, t_srs=t_srs, te=te, tr=tr)
    
    inputs [string] - the input raster file(s). accepts wild
                      cards. enter multiple files as one string.
    output [string] - the output raster file
    etc    [string] - additional command line params. enter all
                      as a scalar string. 
    dstnodata [num] - output no data value
    multi     [T/F] - set multithreading. True/False
    n_threads [num] - set the number of threads. use in place of
                      the multi=True command. 
    overwrite [T/F] - set to overwrite output file. True/False
    of     [string] - output format (e.g. 'ENVI')
    ot      [multi] - output data type (idl or python style)
    r      [string] - resampling method (e.g. 'average')
    srcnodata [num] - input no data value
    s_srs  [string] - the input projection
    t_srs  [string] - the output projection 
    te       [list] - output extent [xmin, ymin, xmax, ymax]
    tr       [list] - output resolution [xres, yres]
    """
    import os
    import aei as aei
    
    # return the docstring if user does not set inputs correctly
    if not (inputs or output):
        print(gdalwarp.__doc__)
        return -1
        
    # parse input params
    fparams = []
    
    if dstnodata:
        if not isinstance(dstnodata, (int, long, float, complex,)):
            print('[ ERROR ]: Unable to parse -dstnodata option %s' % (dstnodata))
            print('[ ERROR ]: Must be scalar int, long or float')
            return -1
        else:
            fparams.append(aei.strJoin(['-dstnodata', dstnodata]))
    
    if srcnodata:
        if not isinstance(srcnodata, (int, long, float, complex)):
            print('[ ERROR ]: Unable to parse -srcnodata option %s' % (srcnodata))
            print('[ ERROR ]: Must be scalar int, long or float')
            return -1
        else:
            fparams.append(aei.strJoin(['-srcnodata', srcnodata]))
    
    if n_threads:
        if not isinstance(n_threads, (int, long, float)):
            print('[ ERROR ]: Unable to parse -n_threads option %s' % (n_threads))
            print('[ ERROR ]: Must be scalar int, long or float')
            return -1
        else:
            fparams.append(aei.strJoin(['-multi -wo NUM_THREADS',n_threads], join='='))
        multi = False
    
    if multi:
        fparams.append('-multi')
    
    if overwrite:
        fparams.append('-overwrite')
    
    if of:
        if type(of) is not str:
            print('[ ERROR ]: Unable to parse -of option %s' % (of))
            print('[ ERROR ]: Must be scalar string')
            return -1
        else:
            fparams.append(aei.strJoin(['-of', of]))
    if ot:
        ot=read.ot(ot)
        fparams.append(aei.strJoin(['-ot', ot]))
    
    if s_srs:
        if type(s_srs) is not str:
            print('[ ERROR ]: Unable to parse -s_srs option %s' % (s_srs))
            print('[ ERROR ]: Must be scalar string')
            return -1
        else:
            fparams.append(aei.strJoin(['-s_srs', s_srs]))
    
    if t_srs:
        if type(t_srs) is not str:
            print('[ ERROR ]: Unable to parse -t_srs option %s' % (t_srs))
            print('[ ERROR ]: Must be a scalar string')
            return -1
        else:
            fparams.append(aei.strJoin(['-t_srs', t_srs]))
    
    if te:
        if type(te) is not list and len(te) is not 4:
            print('[ ERROR ]: Unable to parse -te option %s' % (te))
            print('[ ERROR ]: Must be a 4-element list')
            return -1
        else:
            fparams.append(aei.strJoin(['-te', aei.strJoin(te)]))
    
    if tr:
        if type(te) is not list and len (tr) is not 2:
            print('[ ERROR ]: Unable to parse -tr option %s' % (tr))
            print('[ ERROR ]: Must be a 2-element list')
            return -1
        else:
            fparams.append(aei.strJoin(['-tr', aei.strJoin(tr)]))

    fparams.append(aei.read.etc(etc))
    fparams=aei.strJoin(fparams)
    
    # join and run the command
    ocmd = aei.strJoin([raw.gdalwarp, fparams, inputs, output])
    print("[ STATUS ]: Running gdalwarp")
    print("[ STATUS ]: %s" % ocmd)
    os.system(ocmd)
    
    return(ocmd)
    
# gdal_translate
def gdal_translate(inputs='', output='', etc='', dstnodata=False, 
        of=False, ot=False, srcnodata=False):
    """
    mosaics, reprojects or warps raster imagery
    
    syntax: gdal_translate(inputs, output, etc=etc, dstnodata=dstnodata, of=of, ot=ot, r=r, srcnodata=srcnodata)
    
    inputs [string] - the input raster file(s). accepts wild
                      cards. enter multiple files as one string.
    output [string] - the output raster file
    etc    [string] - additional command line params. enter all
                      as a scalar string. 
    dstnodata [num] - output no data value
    of     [string] - output format (e.g. 'ENVI')
    ot      [multi] - output data type (idl or python style)
    srcnodata [num] - input no data value
    """
    import os
    import aei as aei
    
    # return the docstring if user does not set inputs correctly
    if not (inputs or output):
        print(gdal_translate.__doc__)
        return -1
        
    # parse input params
    fparams = []
    
    if dstnodata:
        if not isinstance(dstnodata, (int, long, float, complex,)):
            print('[ ERROR ]: Unable to parse -dstnodata option %s' % (dstnodata))
            print('[ ERROR ]: Must be scalar int, long or float')
            return -1
        else:
            fparams.append(aei.strJoin(['-dstnodata', dstnodata]))
    
    if srcnodata:
        if not isinstance(srcnodata, (int, long, float, complex)):
            print('[ ERROR ]: Unable to parse -srcnodata option %s' % (srcnodata))
            print('[ ERROR ]: Must be scalar int, long or float')
            return -1
        else:
            fparams.append(aei.strJoin(['-srcnodata', srcnodata]))
    
    if of:
        if type(of) is not str:
            print('[ ERROR ]: Unable to parse -of option %s' % (of))
            print('[ ERROR ]: Must be scalar string')
            return -1
        else:
            fparams.append(aei.strJoin(['-of', of]))
    if ot:
        ot=read.ot(ot)
        fparams.append(aei.strJoin(['-ot', ot]))

    if etc:
        fparams.append(aei.read.etc(etc))
    
    fparams=aei.strJoin(fparams)
    
    # join and run the command
    ocmd = aei.strJoin([raw.gdal_translate, fparams, inputs, output])
    print("[ STATUS ]: Running gdal_translate")
    print("[ STATUS ]: %s" % ocmd)
    os.system(ocmd)
    
    return(ocmd)
    
# gdalbuildvrt
def gdalbuildvrt(inputs='', output='', etc='', vrtnodata=False, 
        of=False, ot=False, srcnodata=False, separate=False):
    """
    mosaics, reprojects or warps raster imagery
    
    syntax: gdaluildvrt(inputs, output, etc=etc, dstnodata=dstnodata, multi=multi, n_threads=n_threads, overwrite=overwrite, of=of, ot=ot, r=r, srcnodata=srcnodata, s_srs=s_srs, t_srs=t_srs, te=te, tr=tr)
    
    inputs [string] - the input raster file(s). accepts wild
                      cards. enter multiple files as one string.
    output [string] - the output raster file
    etc    [string] - additional command line params. enter all
                      as a scalar string. 
    vrtnodata [num] - output no data value
    of     [string] - output format (e.g. 'ENVI')
    ot      [multi] - output data type (idl or python style)
    srcnodata [num] - input no data value
    separate [bool] - set inputs as separate bands
    """
    import os
    import aei as aei
    
    # return the docstring if user does not set inputs correctly
    if not (inputs or output):
        print(gdal_translate.__doc__)
        return -1
        
    # parse input params
    fparams = []
    
    if vrtnodata:
        if not isinstance(vrtnodata, (int, long, float, complex,)):
            print('[ ERROR ]: Unable to parse -vrtnodata option %s' % (dstnodata))
            print('[ ERROR ]: Must be scalar int, long or float')
            return -1
        else:
            fparams.append(aei.strJoin(['-vrtnodata', dstnodata]))
    
    if srcnodata:
        if not isinstance(srcnodata, (int, long, float, complex)):
            print('[ ERROR ]: Unable to parse -srcnodata option %s' % (srcnodata))
            print('[ ERROR ]: Must be scalar int, long or float')
            return -1
        else:
            fparams.append(aei.strJoin(['-srcnodata', srcnodata]))
    
    if of:
        if type(of) is not str:
            print('[ ERROR ]: Unable to parse -of option %s' % (of))
            print('[ ERROR ]: Must be scalar string')
            return -1
        else:
            fparams.append(aei.strJoin(['-of', of]))
    if ot:
        ot=read.ot(ot)
        fparams.append(aei.strJoin(['-ot', ot]))
        
    if separate:
        fparams.append('-separate')

    if etc:
        fparams.append(aei.read.etc(etc))
    
    fparams=aei.strJoin(fparams)
    
    # join and run the command
    ocmd = aei.strJoin([raw.gdalbuildvrt, fparams, output, inputs])
    print("[ STATUS ]: Running gdalbuildvrt")
    print("[ STATUS ]: %s" % ocmd)
    os.system(ocmd)
    
    return(ocmd)

###
# LAStools commands
###

# las2dem
def las2dem(inputs, output='', etc='', odir='', odix='',
    step=2.0, nodata=-9999, ground=False, otif=True, 
    use_tile_bb=False, cores=1):
    """
    interpolates, grids and rasterizes las files
    
    syntax: las2dem(inputs, output='', etc='', odir='', odix='',
      step=2.0, fill=5, subcircle=0.4, nodata=-9999, elevation=True,
      ground=False, lowest=False, highest=False, wgs84=True, 
      utm='', otif=True, use_tile_bb=False, cores=1)
    
    inputs   [string] - the input las/laz file(s). accepts wild
                        cards. enter multiple files as one string.
    output   [string] - the output merged las/laz file. 
    etc      [string] - additional command line params. enter all
                        as a scalar string. 
    odir     [string] - the output directory. superceded by using 'output' variable
    odix     [string] - a string to append to each output. superceded by using 'output' variable
    step      [float] - the raster grid size.
    nodata    [float] - the value where no points are found
    ground     [bool] - create ground model
    otif       [bool] - output as tif format
    use_tile_bb[bool] - set to output a grid only to the extent of the tile
    cores     [float] - number of processors to use
    """
    import aei as aei
    
    # parse input params
    fparams = []
    
    # determine if -o, -odir, or -odix should be set
    #  if -o
    if output:
        fparams.append(aei.strJoin(['-o', output]))
    
    else:
        # if -odir
        if odir:
            fparams.append(aei.strJoin(['-odir', odir]))
        else:
            fparams.append(aei.strJoin(['-odir', aei.params.outdir]))
        
        # if -odix
        if odix:
            fparams.append(aei.strJoin(['-odix', odix]))
            
    # set the step size
    fparams.append(aei.strJoin(['-step', step]))
    
    # set the nodata value
    fparams.append(aei.strJoin(['-nodata', nodata]))
    
    # set ground
    if ground:
        fparams.append(aei.strJoin(['-keep_class', 2]))
        lowest = True
        
    # set tif as output format
    if otif:
        fparams.append('-otif')
    
    # set bounding box
    if use_tile_bb:
        fparams.append('-use_tile_bb')
    
    # set cores
    fparams.append(aei.strJoin(['-cores', cores]))
    
    # add additional parameters passed through etc keyword
    if etc:
        fparams.append(aei.read.etc(etc))
    
    # concatenate params list to string
    fparams = aei.strJoin(fparams)
    
    # join and run the command
    ocmd = aei.strJoin([raw.las2dem, '-i', inputs, fparams])
    
    # report status
    print("[ STATUS ]: Running las2dem")
    print("[ STATUS ]: %s" % ocmd)
    
    aei.params.os.system(ocmd)

# las2las

# lasboundary
def lasboundary(inputs, output='', etc='', odir='', odix='',
    use_bb=False, disjoint=False, concavity=50, oshp=True,
    cores=1):
    """
    creates a boundary polygon for las files
    
    syntax: lasboundary(inputs, output='', etc='', odir='', odix='',
      use_bb=False, disjoint=False, concavity=50, oshp=True
      cores=1)
    
    inputs   [string] - the input las/laz file(s). accepts wild
                        cards. enter multiple files as one string.
    output   [string] - the output merged las/laz file. 
    etc      [string] - additional command line params. enter all
                        as a scalar string. 
    odir     [string] - the output directory. superceded by using 'output' variable
    odix     [string] - a string to append to each output. superceded by using 'output' variable
    use_bb   [string] - set to output the bounding box- not bounds of individual pts
    disjoint   [bool] - set to allow multiple hulls (i.e. islands of pts)
    concavity [float] - distance in meters of largest internal hole
    oshp       [bool] - ensures output format is shp
    cores     [float] - number of processors to use
    """
    import aei as aei
    
    # parse input params
    fparams = []
    
    # determine if -o, -odir, or -odix should be set
    #  if -o
    if output:
        fparams.append(aei.strJoin(['-o', output]))
    
    else:
        # if -odir
        if odir:
            fparams.append(aei.strJoin(['-odir', odir]))
        else:
            fparams.append(aei.strJoin(['-odir', aei.params.outdir]))
        
        # if -odix
        if odix:
            fparams.append(aei.strJoin(['-odix', odix]))
            
    # check if only bounding box is to be used
    if use_bb:
        fparams.append('-use_bb')
        
    # check if allowing disjointed hulls
    if disjoint:
        fparams.append('-disjoint')
        
    # set concavity parameter
    fparams.append(aei.strJoin(['-concavity', concavity]))
    
    # check output format
    if oshp:
        fparams.append('-oshp')
    
    # set number of cores
    fparams.append(aei.strJoin(['-cores', cores]))
    
    # add additional parameters passed through etc keyword
    if etc:
        fparams.append(aei.read.etc(etc))
    
    # concatenate params list to string
    fparams = aei.strJoin(fparams)
    
    # join and run the command
    ocmd = aei.strJoin([raw.lasboundary, '-i', inputs, fparams])
    
    # report status
    print("[ STATUS ]: Running lasboundary")
    print("[ STATUS ]: %s" % ocmd)
    
    aei.params.os.system(ocmd)

# lascanopy
def lascanopy(inputs, output='', etc='', odir='', odix='', merged=False,
    height_cutoff=1.37, step=20, file_are_plots=False,
    loc='', min=False, max=False, avg=False, std=False,
    ske=False, kur=False, cov=False, dns=False, gap=False,
    cover_cutoff=5.0, count=[], density=[], fractions=False, 
    use_tile_bb = False, otif=True, cores=1):
    """
    calculates and rasterizes forestry/canopy metrics using las files. 
      must be run using a height-normalized input file.
    
    syntax: lascanopy(inputs, output='', etc='', odir='', odix='', merged=False,
      height_cutoff=1.37, step=20, file_are_plots=False,
      loc='', min=False, max=False, avg=False, std=False,
      ske=False, kur=False, cov=False, dns=False, gap=False,
      cover_cutoff=5.0, count=[], density=[], fractions=False, 
      use_tile_bb=False, otif=True, cores=1)
    
    inputs        [string] - the input las/laz file(s). accepts wild
                             cards. enter multiple files as one string.
    output        [string] - the output merged las/laz file. 
    etc           [string] - additional command line params. enter all
                             as a scalar string. 
    odir          [string] - the output directory. superceded by using 'output' variable
    odix          [string] - a string to append to each output. superceded by using 'output' variable
    merged          [bool] - set if specifying multiple input files that should be merged
    height_cutoff  [float] - the minimum z value to include (i.e. above dbh)
    step           [float] - the raster grid size.
    files_are_plots [bool] - returns results for the full las file, instead
                             of gridding output using the -step parameter
    loc           [string] - use an ascii file with each line formatted as 
                             [center_x center_y radius] instead of gridding.
    min/max/avg/std [bool] - calculate the min/max/average/standard dev for each grid
    ske/kur         [bool] - calculate the skewness/kurtosis for each grid
    cov             [bool] - calculate canopy cover - 
                             100 * (n 1st returns above canopy_cutoff) / (n 1st returns)
    dns             [bool] - calculate canopy density - 
                             100 * (n points above canopy cutoff) / (n points)
    gap             [bool] - calculate canopy gaps - the inverse of cover or density
    count           [list] - counts number of points within intervals of the
                             list provided (e.g. count = [2, 4, 6, 8] will produce
                             three rasters in intervals [2,4), [4, 6), [6, 8)
    density         [list] - similar to 'count' but reports as normalized density
    fractions       [bool] - report output in fractions (0-1) instead of percent (0-100)
    use_tile_bb     [bool] - set to use the tile bounding box, not full file
    otif            [bool] - ensures output format is tif.
    cores          [float] - set number of processors
    """
    import aei as aei
    
    # parse input params
    fparams = []
    
    # determine if -o, -odir, or -odix should be set
    #  if -o
    if output:
        fparams.append(aei.strJoin(['-o', output]))
    
    else:
        # if -odir
        if odir:
            fparams.append(aei.strJoin(['-odir', odir]))
        else:
            fparams.append(aei.strJoin(['-odir', aei.params.outdir]))
        
        # if -odix
        if odix:
            fparams.append(aei.strJoin(['-odix', odix]))
        else:
            fparams.append(aei.strJoin(['-odix', '_lascanopy']))
    
    # set  merging option
    if merged:
        fparams.append('-merged')
        
    # set height cutoff
    fparams.append(aei.strJoin(['-height_cutoff', height_cutoff]))
    
    # set cover cutoff
    fparams.append(aei.strJoin(['-cover_cutoff', cover_cutoff]))
    
    # determine if gridding, using full file, or using circles
    #  if using full file
    if file_are_plots:
        fparams.append('-files_are_plots')
    
    else:
        # if using cirles with radii
        if loc:
            fparams.append(aei.strJoin(['-loc', loc]))
        
        # if gridding
        else:
            fparams.append(aei.strJoin(['-step', step]))
            
    # set metrics min/max/avg/std/ske/kur/cov/dns/gap
    if min:
        fparams.append('-min')
    if max:
        fparams.append('-max')
    if avg:
        fparams.append('-avg')
    if std:
        fparams.append('-std')
    if ske:
        fparams.append('-ske')
    if kur:
        fparams.append('-kur')
    if cov:
        fparams.append('-cov')
    if dns:
        fparams.append('-dns')
        
    if gap:
        fparams.append('-gap')
        
        # -gap depends on cov or dns being set, so add dns if necessary
        if not cov or dns:
            fparams.append('-dns')
            
    # set count parameters
    if count:
        fparams.append(aei.strJoin(['-c', aei.strJoin(count)]))
    
    # set density parameters    
    if density:
        fparams.append(aei.strJoin(['-d', aei.strJoin(density)]))
    
    # set if returning as fraction
    if fractions:
        fparams.append('-fractions')
        
    # set if output raster format is tif
    if otif:
        fparams.append('-otif')
    
    # set number of cores
    fparams.append(aei.strJoin(['-cores', cores]))
    
    # add additional parameters passed through etc keyword
    if etc:
        fparams.append(aei.read.etc(etc))
    
    # concatenate params list to string
    fparams = aei.strJoin(fparams)
    
    # join and run the command
    ocmd = aei.strJoin([raw.lascanopy, '-i', inputs, fparams])
    
    # report status
    print("[ STATUS ]: Running lascanopy")
    print("[ STATUS ]: %s" % ocmd)
    
    aei.params.os.system(ocmd)

# lasclassify
def lasclassify(inputs, output='', etc='', odir='', odix='',
    planar=0.1, ground_offset=2.0, olaz=True, cores=1):
    """
    classifies building and vegetation in las files. requires ground
      points and heights calculated in input file.
    
    syntax: lasclassify(input[s], output='', etc='', odir='', odix='',
      planar=0.1, ground_offset=2.0, olaz=True, cores=1)
    
    inputs       [string] - the input las/laz file(s). accepts wild
                            cards. enter multiple files as one string.
    output       [string] - the output merged las/laz file. 
    etc          [string] - additional command line params. enter all
                            as a scalar string. 
    odir         [string] - the output directory. superceded by using 'output' variable
    odix         [string] - a string to append to each output. superceded by using 'output' variable
    planar        [float] - threshold for determining planar surfaces.
                            increasing the value allows for fitting planes
                            to noisier (e.g. less-planar) data.
    ground_offset [float] - sets the minimum height above which veg is
                            classified. default = 2.0
    olaz           [bool] - ensures output format is laz.
    cores         [float] - number of processors to use
    """
    import aei as aei
    
    # parse input params
    fparams = []
    
    # determine if -o, -odir, or -odix should be set
    #  if -o
    if output:
        fparams.append(aei.strJoin(['-o', output]))
    
    else:
        # if -odir
        if odir:
            fparams.append(aei.strJoin(['-odir', odir]))
        else:
            fparams.append(aei.strJoin(['-odir', aei.params.outdir]))
        
        # if -odix
        if odix:
            fparams.append(aei.strJoin(['-odix', odix]))
    
    # add planar parameter
    fparams.append(aei.strJoin(['-planar', planar]))
    
    # add ground_offset parameter
    fparams.append(aei.strJoin(['-ground_offset', ground_offset]))
    
    # add olaz parameter if set
    if olaz:
        fparams.append('-olaz')
    
    # set number of cores
    fparams.append(aei.strJoin(['-cores', cores]))
    
    # add additional parameters passed through etc keyword
    if etc:
        fparams.append(aei.read.etc(etc))
    
    # concatenate params list to string
    fparams = aei.strJoin(fparams)
    
    # join and run the command
    ocmd = aei.strJoin([raw.lasclassify, '-i', inputs, fparams])
    
    # report status
    print("[ STATUS ]: Running lasclassify")
    print("[ STATUS ]: %s" % ocmd)
    
    aei.params.os.system(ocmd)

# lasclip

# lasgrid
def lasgrid(inputs, output='', etc='', odir='', odix='',
    step=2.0, fill=5, subcircle=0.4, nodata=-9999, elevation=True,
    ground=False, lowest=False, highest=False, wgs84=True, 
    utm='', otif=True, use_tile_bb=False, cores=1):
    """
    grids and rasterizes las files
    
    syntax: lasgrid(inputs, output='', etc='', odir='', odix='',
      step=2.0, fill=5, subcircle=0.4, nodata=-9999, elevation=True,
      ground=False, lowest=False, highest=False, wgs84=True, 
      utm='', otif=True, use_tile_bb=False, cores=1)
    
    inputs   [string] - the input las/laz file(s). accepts wild
                        cards. enter multiple files as one string.
    output   [string] - the output merged las/laz file. 
    etc      [string] - additional command line params. enter all
                        as a scalar string. 
    odir     [string] - the output directory. superceded by using 'output' variable
    odix     [string] - a string to append to each output. superceded by using 'output' variable
    step      [float] - the raster grid size.
    fill      [float] - fills voids in the grid within search radius
    subcircle [float] - 'splats' the points to larger radius
    nodata    [float] - the value where no points are found
    elevation  [bool] - use elevation values
    ground     [bool] - create ground model
    lowest     [bool] - grid the lowest points
    highest    [bool] - grid the highest points
    wgs84      [bool] - use the wgs84 ellipsoid
    otif       [bool] - output as tif format
    use_tile_bb[bool] - set to output a grid only to the extent of the tile
    cores     [float] - set number of processors to use
    """
    import aei as aei
    
    # parse input params
    fparams = []
    
    # determine if -o, -odir, or -odix should be set
    #  if -o
    if output:
        fparams.append(aei.strJoin(['-o', output]))
    
    else:
        # if -odir
        if odir:
            fparams.append(aei.strJoin(['-odir', odir]))
        else:
            fparams.append(aei.strJoin(['-odir', aei.params.outdir]))
        
        # if -odix
        if odix:
            fparams.append(aei.strJoin(['-odix', odix]))
            
    # set the step size
    fparams.append(aei.strJoin(['-step', step]))
    
    # set the fill size
    fparams.append(aei.strJoin(['-fill', fill]))
    
    # set the subcircle size
    fparams.append(aei.strJoin(['-subcircle', subcircle]))
    
    # set the nodata value
    fparams.append(aei.strJoin(['-nodata', nodata]))
    
    # set elevation as target value
    if elevation:
        fparams.append('-elevation')
    
    # set ground
    if ground:
        fparams.append(aei.strJoin(['-keep_class', 2]))
        lowest = True
    
    # set if highest or lowest values are returned
    if highest:
        fparams.append('-highest')
    else:
        if lowest:
            fparams.append('-lowest')
        
    # set wgs84 as default ellipsoid
    if wgs84:
        fparams.append('-wgs84')
        
    # set tif as output format
    if otif:
        fparams.append('-otif')
    
    # set bounding box
    if use_tile_bb:
        fparams.append('-use_tile_bb')
    
    # set number of cores to use
    fparams.append(aei.strJoin(['-cores', cores]))
    
    # add additional parameters passed through etc keyword
    if etc:
        fparams.append(aei.read.etc(etc))
    
    # concatenate params list to string
    fparams = aei.strJoin(fparams)
    
    # join and run the command
    ocmd = aei.strJoin([raw.lasgrid, '-i', inputs, fparams])
    
    # report status
    print("[ STATUS ]: Running lasgrid")
    print("[ STATUS ]: %s" % ocmd)
    
    aei.params.os.system(ocmd)

# lasground
def lasground(inputs, output='', etc='', odir='', odix='',
    olaz=True, not_airborne=False, city=False, town=False,
    compute_height=True, replace_z=False,
    fine=False, extra_fine=False,
    coarse=False, extra_coarse=False, cores=1,
    non_ground_unchanged=True):
    """
    calculates ground points for las files
    
    syntax: lasground(input[s], output='', etc='',
      odir='', odix='', replace_z=False, olaz=True,
      city=False, town=False, compute_height=True,
      fine=False, extra_fine=False, coarse=false, 
      extra_coarse=False, non_ground_unchanged=True)
    
    inputs       [string] - the input las/laz file(s). accepts wild
                            cards. enter multiple files as one string.
    output       [string] - the output merged las/laz file. 
    etc          [string] - additional command line params. enter all
                            as a scalar string. 
    odir         [string] - the output directory. superceded by using 'output' variable
    odix         [string] - a string to append to each output. 
                            superceded by using 'output' variable
    olaz           [bool] - ensures output format is laz. default = True
    not_airborne   [bool] - set this flag if using ground LiDAR data.
    city           [bool] - uses widest radius to find ground
    town           [bool] - uses wide radius to find ground
    fine           [bool] - use in steep terrain
    extra_fine     [bool] - use in extra steep terrain
    coarse         [bool] - use in flat terrain
    extra_coarse   [bool] - use in extra flat terrain
    compute_height [bool] - use to compute height AGL, so you don't 
                            need lasheight
    replace_z      [bool] - replaces the z value with the height above
                            ground. use in conjunction with compute_height
    cores         [float] - set the number of processors to use
    non_ground_unchanged [bool] - set this to keep original classification
                                  in place and only re-classify ground.
    """
    import aei as aei
    
    # parse input params
    fparams = []
    
    # determine if -o, -odir, or -odix should be set
    #  if -o
    if output:
        fparams.append(aei.strJoin(['-o', output]))
    
    else:
        # if -odir
        if odir:
            fparams.append(aei.strJoin(['-odir', odir]))
        else:
            fparams.append(aei.strJoin(['-odir', aei.params.outdir]))
        
        # if -odix
        if odix:
            fparams.append(aei.strJoin(['-odix', odix]))
    
    # add olaz parameter if set
    if olaz:
        fparams.append('-olaz')
    
    # allow non-airborne lidar
    if not_airborne:
        fparams.append('-not_airborne')
        
    # use largest search radius    
    if city:
        fparams.append('-city')
    
    # use large search radius    
    if town:
        fparams.append('-town')
        
    # compute height AGL    
    if compute_height:
        fparams.append('-compute_height')
        
    # add replace_z to param list if set
    if replace_z:
        fparams.append('-replace_z')
        
        # double check that compute height is also set
        if not compute_height:
            fparams.append('-compute_height')
    
    # use small search radius
    if fine:
        fparams.append('-fine')
    
    # use extra small search radius    
    if extra_fine:
        fparams.append('-extra_fine')
    
    # for flat areas    
    if coarse:
        fparams.append('-coarse')
    
    # for extra flat areas    
    if extra_coarse:
        fparams.append('-extra_coarse')
        
    # set number of processors to use
    fparams.append(aei.strJoin(['-cores', cores]))
    
    # keep default classification info
    if non_ground_unchanged:
        fparams.append('-non_ground_unchanged')
    
    # add additional parameters passed through etc keyword
    if etc:
        fparams.append(aei.read.etc(etc))
    
    # concatenate params list to string
    fparams = aei.strJoin(fparams)
    
    # join and run the command
    ocmd = aei.strJoin([raw.lasground, '-i', inputs, fparams])
    
    # report status
    print("[ STATUS ]: Running lasground")
    print("[ STATUS ]: %s" % ocmd)
    
    aei.params.os.system(ocmd)

# lasheight
def lasheight(inputs, output='', etc='', odir='', odix='',
    replace_z=False, olaz=True, cores=1):
    """
    calculates height of las files
    
    syntax: lasheight(input[s], output='', etc='', odir='',
      odix='', replace_z=False, olaz=True, cores=1)
    
    inputs  [string] - the input las/laz file(s). accepts wild
                       cards. enter multiple files as one string.
    output  [string] - the output merged las/laz file. 
    etc     [string] - additional command line params. enter all
                       as a scalar string. 
    odir    [string] - the output directory. superceded by using 'output' variable
    odix    [string] - a string to append to each output. superceded by using 'output' variable
    replace_z [bool] - replaces the z value with the height above
                       ground. default = False
    olaz      [bool] - ensures output format is laz. default = True
    cores    [float] - sets the number of processors to use
    """
    import aei as aei
    
    # parse input params
    fparams = []
    
    # determine if -o, -odir, or -odix should be set
    #  if -o
    if output:
        fparams.append(aei.strJoin(['-o', output]))
    
    else:
        # if -odir
        if odir:
            fparams.append(aei.strJoin(['-odir', odir]))
        else:
            fparams.append(aei.strJoin(['-odir', aei.params.outdir]))
        
        # if -odix
        if odix:
            fparams.append(aei.strJoin(['-odix', odix]))
    
    # add replace_z to param list if set
    if replace_z:
        fparams.append('-replace_z')
    
    # add olaz parameter if set
    if olaz:
        fparams.append('-olaz')
    
    # set number of cores to use
    fparams.append(aei.strJoin(['-cores', cores]))
    
    # add additional parameters passed through etc keyword
    if etc:
        fparams.append(aei.read.etc(etc))
    
    # concatenate params list to string
    fparams = aei.strJoin(fparams)
    
    # join and run the command
    ocmd = aei.strJoin([raw.lasheight, '-i', inputs, fparams])
    
    # report status
    print("[ STATUS ]: Running lasheight")
    print("[ STATUS ]: %s" % ocmd)
    
    aei.params.os.system(ocmd)
    
# lasindex
def lasindex(inputs, etc=''):
    """
    creates an index (*.lax) for las files
    
    syntax: lasindex(input[s], etc=etc)
    
    inputs  [string] - the input las/laz file(s). accepts wild
                       cards. enter multiple files as one string.
    etc     [string] - additional command line params. enter all
                       as a scalar string. 
    """
    import aei as aei
    
    # parse input params
    fparams = []
    
    # add additional parameters passed through etc keyword
    if etc:
        fparams.append(aei.read.etc(etc))
    
    # concatenate params list to string
    fparams = aei.strJoin(fparams)
    
    # join and run the command
    ocmd = aei.strJoin([raw.lasindex, '-i', inputs, fparams])
    
    # report status
    print("[ STATUS ]: Running lasindex")
    print("[ STATUS ]: %s" % ocmd)
    
    aei.params.os.system(ocmd)

# lasinfo

# lasmerge
def lasmerge(inputs, output, etc='',
    rescale=[0.01, 0.01, 0.001], olaz=True):
    """
    merges and rescales las files
    
    syntax: lasmerge(inputs, output, etc=etc, rescale=rescale)
    
    inputs [string] - the input las/laz file(s). accepts wild
                      cards. enter multiple files as one string.
    output [string] - the output merged las/laz file. 
    etc    [string] - additional command line params. enter all
                      as a scalar string. 
    rescale  [list] - the rescale/precision parameters for the
                      output merged output file. default: 
                      [0.01, 0.01, 0.001]. set to False to not
                      rescale
    """
    import aei as aei
    
    # parse input params
    fparams = []
    
    # add rescaling to param list if set
    if rescale:
        fparams.append(aei.strJoin(['-rescale', aei.strJoin(rescale)]))
    
    # add olaz parameter if set
    if olaz:
        fparams.append('-olaz')
    
    # add additional parameters passed through etc keyword
    if etc:
        fparams.append(aei.read.etc(etc))
    
    # concatenate params list to string
    fparams = aei.strJoin(fparams)
    
    # join and run the command
    ocmd = aei.strJoin([raw.lasmerge, '-i', inputs, '-o', output,
            fparams])
    
    # report status
    print("[ STATUS ]: Running lasmerge")
    print("[ STATUS ]: %s" % ocmd)
    
    aei.params.os.system(ocmd)
    
# lastile
def lastile(inputs, odir, etc='', odix='', tile_size=1000, buffer=100,
    olaz=True, cores=1):
    """
    tiles las files
    
    syntax: lastile(inputs, odir, etc=etc, odix=odix, tile_size=tile_size)
    
    inputs   [string] - the input las/laz file(s). accepts wild
                        cards. enter multiple files as one string.
    odir     [string] - the output directory for the tiles. if not set,
                        outputs tiles to the default output directory
                        (set under aei.params)
    odix     [string] - a string to append to each output. superceded by using 'output' variable
    etc      [string] - additional command line params. enter all
                        as a scalar string. 
    tile_size [float] - the output tile size in units of the horizontal
                        projection. default = 1000
    buffer    [float] - the buffer surrounding each tile. default = 100
    cores     [float] - the number of processors to use
    """
    import aei as aei
    
    # parse input params
    fparams = []
    
    # if -odix
    if odix:
        fparams.append(aei.strJoin(['-odix', odix]))

    # add tile_size to parameters list
    fparams.append(aei.strJoin(['-tile_size', tile_size]))
    
    # add buffer distance to parameters list
    fparams.append(aei.strJoin(['-buffer', buffer]))
    
    # add olaz parameter if set
    if olaz:
        fparams.append('-olaz')
        
    # add n cores
    fparams.append(aei.strJoin(['-cores', cores]))
    
    # add additional parameters passed through etc keyword
    if etc:
        fparams.append(aei.read.etc(etc))
    
    # concatenate params list to string
    fparams = aei.strJoin(fparams)
    
    # join and run the command
    ocmd = aei.strJoin([raw.lastile, '-i', inputs, '-odir', odir,
            fparams])
    
    # report status
    print("[ STATUS ]: Running lastile")
    print("[ STATUS ]: %s" % ocmd)
    
    aei.params.os.system(ocmd)

# laszip