#####
# sets up paths and calls to external executables.
#
# c. 2016-2017 Christopher Anderson
#####
import subprocess

# set up a function to run commands and return stdout/stderr
def run(cmd, stderr=True):
    """
    """
    # set whether or not stderr is included in return or just stdout
    if stderr:
        se = subprocess.STDOUT
    else:
        se = None

    # run the command, and return stdout as a list
    try:
        subprocess.check_output(cmd, shell=True, stderr=se)

    # raise an exception and print the error if the command fails
    except subprocess.CalledProcessError as e:
        output = e.output.strip()
        sp = output.find(b":") + 2
        print(output[sp:])


# set up commands to call binaries
class getPaths:
    def __init__(self):
        from earthtools import params

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
        self.gdal_proximity = ' '.join({
            'Windows' : [params.pathcat(params.gdalbase, 'gdal_proximity.py')],
            'Linux'   : [params.pathcat(params.gdalbase, 'gdal_proximity.py')]
            }[params.system])

        # orfeo toolbox commands
        self.otbBandMath = ''.join({
            'Windows' : ['otbcli_BandMath'],
            'Linux'   : ['otbcli_BandMath']
            }[params.system])
        self.otbConcatenateImages = ''.join({
            'Windows' : ['otbcli_ConcatenateImages'],
            'Linux'   : ['otbcli_ConcatenateImages']
            }[params.system])
        self.otbHyperspectralUnmixing = ''.join({
            'Windows' : ['otbcli_HyperspectralUnmixing'],
            'Linux'   : ['otbcli_HyperspectralUnmixing']
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
    import earthtools as et

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
            fparams.append(et.fn.strJoin(['-dstnodata', dstnodata]))

    if srcnodata:
        if not isinstance(srcnodata, (int, long, float, complex)):
            print('[ ERROR ]: Unable to parse -srcnodata option %s' % (srcnodata))
            print('[ ERROR ]: Must be scalar int, long or float')
            return -1
        else:
            fparams.append(et.fn.strJoin(['-srcnodata', srcnodata]))

    if n_threads:
        if not isinstance(n_threads, (int, long, float)):
            print('[ ERROR ]: Unable to parse -n_threads option %s' % (n_threads))
            print('[ ERROR ]: Must be scalar int, long or float')
            return -1
        else:
            fparams.append(et.fn.strJoin(['-multi -wo NUM_THREADS',n_threads], join='='))
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
            fparams.append(et.fn.strJoin(['-of', of]))
    if ot:
        ot=read.ot(ot)
        fparams.append(et.fn.strJoin(['-ot', ot]))

    if s_srs:
        if type(s_srs) is not str:
            print('[ ERROR ]: Unable to parse -s_srs option %s' % (s_srs))
            print('[ ERROR ]: Must be scalar string')
            return -1
        else:
            fparams.append(et)

    if t_srs:
        if type(t_srs) is not str:
            print('[ ERROR ]: Unable to parse -t_srs option %s' % (t_srs))
            print('[ ERROR ]: Must be a scalar string')
            return -1
        else:
            fparams.append(et.fn.strJoin(['-t_srs', t_srs]))

    if te:
        if type(te) is not list and len(te) is not 4:
            print('[ ERROR ]: Unable to parse -te option %s' % (te))
            print('[ ERROR ]: Must be a 4-element list')
            return -1
        else:
            fparams.append(et.fn.strJoin(['-te', et.fn.strJoin(te)]))

    if tr:
        if type(te) is not list and len (tr) is not 2:
            print('[ ERROR ]: Unable to parse -tr option %s' % (tr))
            print('[ ERROR ]: Must be a 2-element list')
            return -1
        else:
            fparams.append(et.fn.strJoin(['-tr', et.fn.strJoin(tr)]))

    if etc:
        fparams.append(et.read.etc(etc))

    fparams=et.fn.strJoin(fparams)

    # join and run the command
    ocmd = et.fn.strJoin([raw.gdalwarp, fparams, inputs, output])
    print("[ STATUS ]: Running gdalwarp")
    print("[ STATUS ]: %s" % ocmd)

    # run the command
    run(ocmd)

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
    import earthtools as et

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
            fparams.append(et.fn.strJoin(['-dstnodata', dstnodata]))

    if srcnodata:
        if not isinstance(srcnodata, (int, long, float, complex)):
            print('[ ERROR ]: Unable to parse -srcnodata option %s' % (srcnodata))
            print('[ ERROR ]: Must be scalar int, long or float')
            return -1
        else:
            fparams.append(et.fn.strJoin(['-srcnodata', srcnodata]))

    if of:
        if type(of) is not str:
            print('[ ERROR ]: Unable to parse -of option %s' % (of))
            print('[ ERROR ]: Must be scalar string')
            return -1
        else:
            fparams.append(et.fn.strJoin(['-of', of]))
    if ot:
        ot=et.read.ot(ot)
        fparams.append(et.fn.strJoin(['-ot', ot]))

    if etc:
        fparams.append(et.read.etc(etc))

    fparams=et.fn.strJoin(fparams)

    # join the command
    ocmd = et.fn.strJoin([raw.gdal_translate, fparams, inputs, output])
    print("[ STATUS ]: Running gdal_translate")
    print("[ STATUS ]: %s" % ocmd)

    # run the command
    run(ocmd)

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
    import earthtools as et

    # return the docstring if user does not set inputs correctly
    if not (inputs or output):
        print(gdalbuildvrt.__doc__)
        return -1

    # parse input params
    fparams = []

    if vrtnodata:
        if not isinstance(vrtnodata, (int, long, float, complex,)):
            print('[ ERROR ]: Unable to parse -vrtnodata option %s' % (dstnodata))
            print('[ ERROR ]: Must be scalar int, long or float')
            return -1
        else:
            fparams.append(et.fn.strJoin(['-vrtnodata', dstnodata]))

    if srcnodata:
        if not isinstance(srcnodata, (int, long, float, complex)):
            print('[ ERROR ]: Unable to parse -srcnodata option %s' % (srcnodata))
            print('[ ERROR ]: Must be scalar int, long or float')
            return -1
        else:
            fparams.append(et.fn.strJoin(['-srcnodata', srcnodata]))

    if of:
        if type(of) is not str:
            print('[ ERROR ]: Unable to parse -of option %s' % (of))
            print('[ ERROR ]: Must be scalar string')
            return -1
        else:
            fparams.append(et.fn.strJoin(['-of', of]))
    if ot:
        ot=read.ot(ot)
        fparams.append(et.fn.strJoin(['-ot', ot]))

    if separate:
        fparams.append('-separate')

    if etc:
        fparams.append(et.read.etc(etc))

    fparams=et.fn.strJoin(fparams)

    # join the command
    ocmd = et.fn.strJoin([raw.gdalbuildvrt, fparams, output, inputs])
    print("[ STATUS ]: Running gdalbuildvrt")
    print("[ STATUS ]: %s" % ocmd)

    # run the command
    run(ocmd)

# gdal_proximity
def gdal_proximity(inputs='', output='', etc='', srcband=None,
        dstband=None, of="GTiff", ot=False, co=None, distunits="GEO",
        maxdist=None, nodata=None, use_input_nodata=True):
    """
    generates a raster proximity map

    syntax: gdal_proximity(inputs, output, etc=etc, dstnodata=dstnodata, multi=multi, n_threads=n_threads, overwrite=overwrite, of=of, ot=ot, r=r, srcnodata=srcnodata, s_srs=s_srs, t_srs=t_srs, te=te, tr=tr)

    inputs [string] - the input raster file(s). accepts wild
                      cards. enter multiple files as one string.
    output [string] - the output raster file
    etc    [string] - additional command line params. enter all
                      as a scalar string.
    of     [string] - output format (e.g. 'ENVI')
    ot      [multi] - output data type (idl or python style)
    nodata    [num] - output no data value
    """
    import os
    import earthtools as et

    # return the docstring if user does not set inputs correctly
    if not (inputs or output):
        print(gdal_proximity.__doc__)
        return -1

    # parse input params
    fparams = []

    if srcband:
        fparams.append(et.fn.strJoin(['-srcband', srcband]))

    if dstband:
        fparams.append(et.fn.strJoin(['-dstband', dstband]))

    if of:
        if type(of) is not str:
            print('[ ERROR ]: Unable to parse -of option %s' % (of))
            print('[ ERROR ]: Must be scalar string')
            return -1
        else:

            # set default compression if using geotiff format
            if of == "GTiff":
                co = "COMPRESS=LZW"

            fparams.append(et.fn.strJoin(['-of', of]))
    if ot:
        ot=read.ot(ot)
        fparams.append(et.fn.strJoin(['-ot', ot]))

    if co:
        fparams.append(et.fn.strJoin(['-co', co]))

    if distunits:
        if not distunits.upper() in ["PIXEL", "GEO"]:
            print('[ ERROR ]: Unrecognized distance units: %s' % dstunits)
            print('[ ERROR ]: Using GEO as default')
            distunits = "GEO"
        fparams.append(et.fn.strJoin(['-distunits', distunits]))

    if maxdist:
        fparams.append(et.fn.strJoin(['-maxdist', maxdist]))

    if use_input_nodata:
        fparams.append(et.fn.strJoin(['-use_input_nodata', 'YES']))

    if etc:
        fparams.append(et.read.etc(etc))

    fparams=et.fn.strJoin(fparams)

    # join the command
    ocmd = et.fn.strJoin([raw.gdal_proximity, fparams, inputs, output])
    print("[ STATUS ]: Running gdal_proximity")
    print("[ STATUS ]: %s" % ocmd)

    # run the command
    run(ocmd)

###
# LAStools commands
###

# las2dem
def las2dem(inputs, output='', etc='', odir='', odix='',
    step=2.0, nodata=-9999, ground=False, otif=True,
    use_tile_bb=False, cores=1, drop_above=None):
    """
    interpolates, grids and rasterizes las files

    syntax: las2dem(inputs, output='', etc='', odir='', odix='',
      step=2.0, fill=5, subcircle=0.4, nodata=-9999, elevation=True,
      ground=False, lowest=False, highest=False, wgs84=True,
      utm='', otif=True, use_tile_bb=False, cores=1, drop_above=None)

    inputs    [string] - the input las/laz file(s). accepts wild
                         cards. enter multiple files as one string.
    output    [string] - the output merged las/laz file.
    etc       [string] - additional command line params. enter all
                         as a scalar string.
    odir      [string] - the output directory. superceded by using 'output' variable
    odix      [string] - a string to append to each output. superceded by using 'output' variable
    step       [float] - the raster grid size.
    nodata     [float] - the value where no points are found
    ground      [bool] - create ground model
    otif        [bool] - output as tif format
    use_tile_bb [bool] - set to output a grid only to the extent of the tile
    drop_above [float] - set this to drop points above a certain height
    cores      [float] - number of processors to use
    """
    import earthtools as et

    # parse input params
    fparams = []

    # determine if -o, -odir, or -odix should be set
    #  if -o
    if output:
        fparams.append(et.fn.strJoin(['-o', output]))

    else:
        # if -odir
        if odir:
            fparams.append(et.fn.strJoin(['-odir', odir]))
        else:
            fparams.append(et.fn.strJoin(['-odir', et.params.outdir]))

        # if -odix
        if odix:
            fparams.append(et.fn.strJoin(['-odix', odix]))

    # set the step size
    fparams.append(et.fn.strJoin(['-step', step]))

    # set the nodata value
    fparams.append(et.fn.strJoin(['-nodata', nodata]))

    # set ground
    if ground:
        fparams.append(et.fn.strJoin(['-keep_class', 2]))
        lowest = True

    # set tif as output format
    if otif:
        fparams.append('-otif')

    # set bounding box
    if use_tile_bb:
        fparams.append('-use_tile_bb')

    # set cores
    fparams.append(et.fn.strJoin(['-cores', cores]))

    # set drop_above
    if drop_above:
        fparams.append(et.fn.strJoin(['-drop_above', drop_above]))

    # add additional parameters passed through etc keyword
    if etc:
        fparams.append(et.read.etc(etc))

    # concatenate params list to string
    fparams = et.fn.strJoin(fparams)

    # join and run the command
    ocmd = et.fn.strJoin([raw.las2dem, '-i', inputs, fparams])

    # report status
    print("[ STATUS ]: Running las2dem")
    print("[ STATUS ]: %s" % ocmd)

    # run the command
    run(ocmd)

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
    import earthtools as et

    # parse input params
    fparams = []

    # determine if -o, -odir, or -odix should be set
    #  if -o
    if output:
        fparams.append(et.fn.strJoin(['-o', output]))

    else:
        # if -odir
        if odir:
            fparams.append(et.fn.strJoin(['-odir', odir]))
        else:
            fparams.append(et.fn.strJoin(['-odir', et.params.outdir]))

        # if -odix
        if odix:
            fparams.append(et.fn.strJoin(['-odix', odix]))

    # check if only bounding box is to be used
    if use_bb:
        fparams.append('-use_bb')

    # check if allowing disjointed hulls
    if disjoint:
        fparams.append('-disjoint')

    # set concavity parameter
    fparams.append(et.fn.strJoin(['-concavity', concavity]))

    # check output format
    if oshp:
        fparams.append('-oshp')

    # set number of cores
    fparams.append(et.fn.strJoin(['-cores', cores]))

    # add additional parameters passed through etc keyword
    if etc:
        fparams.append(et.read.etc(etc))

    # concatenate params list to string
    fparams = et.fn.strJoin(fparams)

    # join and run the command
    ocmd = et.fn.strJoin([raw.lasboundary, '-i', inputs, fparams])

    # report status
    print("[ STATUS ]: Running lasboundary")
    print("[ STATUS ]: %s" % ocmd)

    # run the command
    run(ocmd)

# lascanopy
def lascanopy(inputs, output='', etc='', odir='', odix='', merged=False,
    height_cutoff=1.37, step=20, file_are_plots=False,
    loc='', min=False, max=False, avg=False, std=False,
    ske=False, kur=False, cov=False, dns=False, gap=False,
    cover_cutoff=5.0, count=[], density=[], fractions=False,
    use_tile_bb = False, otif=True, cores=1, veg=False):
    """
    calculates and rasterizes forestry/canopy metrics using las files.
      must be run using a height-normalized input file.

    syntax: lascanopy(inputs, output='', etc='', odir='', odix='', merged=False,
      height_cutoff=1.37, step=20, file_are_plots=False,
      loc='', min=False, max=False, avg=False, std=False,
      ske=False, kur=False, cov=False, dns=False, gap=False,
      cover_cutoff=5.0, count=[], density=[], fractions=False,
      use_tile_bb=False, otif=True, cores=1, veg=True)

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
    veg             [bool] - masks non-veg data.
    """
    import earthtools as et

    # parse input params
    fparams = []

    # determine if -o, -odir, or -odix should be set
    #  if -o
    if output:
        fparams.append(et.fn.strJoin(['-o', output]))

    else:
        # if -odir
        if odir:
            fparams.append(et.fn.strJoin(['-odir', odir]))
        else:
            fparams.append(et.fn.strJoin(['-odir', et.params.scratchdir]))

        # if -odix
        if odix:
            fparams.append(et.fn.strJoin(['-odix', odix]))

    # set  merging option
    if merged:
        fparams.append('-merged')

    # set height cutoff
    fparams.append(et.fn.strJoin(['-height_cutoff', height_cutoff]))

    # set cover cutoff
    fparams.append(et.fn.strJoin(['-cover_cutoff', cover_cutoff]))

    # determine if gridding, using full file, or using circles
    #  if using full file
    if file_are_plots:
        fparams.append('-files_are_plots')

    else:
        # if using cirles with radii
        if loc:
            fparams.append(et.fn.strJoin(['-loc', loc]))

        # if gridding
        else:
            fparams.append(et.fn.strJoin(['-step', step]))

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
        fparams.append(et.fn.strJoin(['-c', et.fn.strJoin(count)]))

    # set density parameters
    if density:
        fparams.append(et.fn.strJoin(['-d', et.fn.strJoin(density)]))

    # set if returning as fraction
    if fractions:
        fparams.append('-fractions')

    # set if output raster format is tif
    if otif:
        fparams.append('-otif')

    # set number of cores
    fparams.append(et.fn.strJoin(['-cores', cores]))

    # set params to run on veg only
    if veg:
        fparams.append(et.fn.strJoin(['-keep_class', '2 3 4 5']))

    # add additional parameters passed through etc keyword
    if etc:
        fparams.append(et.read.etc(etc))

    # concatenate params list to string
    fparams = et.fn.strJoin(fparams)

    # join and run the command
    ocmd = et.fn.strJoin([raw.lascanopy, '-i', inputs, fparams])

    # report status
    print("[ STATUS ]: Running lascanopy")
    print("[ STATUS ]: %s" % ocmd)

    # run the command
    run(ocmd)

# lasclassify
def lasclassify(inputs, output='', etc='', odir='', odix='',
    planar=0.1, ground_offset=2.0, olaz=True, cores=1, step = 4):
    """
    classifies building and vegetation in las files. requires ground
      points and heights calculated in input file.

    syntax: lasclassify(input[s], output='', etc='', odir='', odix='',
      planar=0.1, ground_offset=2.0, olaz=True, cores=1, res = 4)

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
    res           [float] - the search radius to find similar points
    """
    import earthtools as et

    # parse input params
    fparams = []

    # determine if -o, -odir, or -odix should be set
    #  if -o
    if output:
        fparams.append(et.fn.strJoin(['-o', output]))

    else:
        # if -odir
        if odir:
            fparams.append(et.fn.strJoin(['-odir', odir]))
        else:
            fparams.append(et.fn.strJoin(['-odir', et.params.outdir]))

        # if -odix
        if odix:
            fparams.append(et.fn.strJoin(['-odix', odix]))

    # add planar parameter
    fparams.append(et.fn.strJoin(['-planar', planar]))

    # add ground_offset parameter
    fparams.append(et.fn.strJoin(['-ground_offset', ground_offset]))

    # add olaz parameter if set
    if olaz:
        fparams.append('-olaz')

    # set number of cores
    fparams.append(et.fn.strJoin(['-cores', cores]))

    # add additional parameters passed through etc keyword
    if etc:
        fparams.append(et.read.etc(etc))

    # concatenate params list to string
    fparams = et.fn.strJoin(fparams)

    # join and run the command
    ocmd = et.fn.strJoin([raw.lasclassify, '-i', inputs, fparams])

    # report status
    print("[ STATUS ]: Running lasclassify")
    print("[ STATUS ]: %s" % ocmd)

    # run the command
    run(ocmd)

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
    import earthtools as et

    # parse input params
    fparams = []

    # determine if -o, -odir, or -odix should be set
    #  if -o
    if output:
        fparams.append(et.fn.strJoin(['-o', output]))

    else:
        # if -odir
        if odir:
            fparams.append(et.fn.strJoin(['-odir', odir]))
        else:
            fparams.append(et.fn.strJoin(['-odir', et.params.outdir]))

        # if -odix
        if odix:
            fparams.append(et.fn.strJoin(['-odix', odix]))

    # set the step size
    fparams.append(et.fn.strJoin(['-step', step]))

    # set the fill size
    fparams.append(et.fn.strJoin(['-fill', fill]))

    # set the subcircle size
    fparams.append(et.fn.strJoin(['-subcircle', subcircle]))

    # set the nodata value
    fparams.append(et.fn.strJoin(['-nodata', nodata]))

    # set elevation as target value
    if elevation:
        fparams.append('-elevation')

    # set ground
    if ground:
        fparams.append(et.fn.strJoin(['-keep_class', 2]))
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
    fparams.append(et.fn.strJoin(['-cores', cores]))

    # add additional parameters passed through etc keyword
    if etc:
        fparams.append(et.read.etc(etc))

    # concatenate params list to string
    fparams = et.fn.strJoin(fparams)

    # join and run the command
    ocmd = et.fn.strJoin([raw.lasgrid, '-i', inputs, fparams])

    # report status
    print("[ STATUS ]: Running lasgrid")
    print("[ STATUS ]: %s" % ocmd)

    # run the command
    run(ocmd)

# lasground
def lasground(inputs, output='', etc='', odir='', odix='',
    olaz=True, not_airborne=False, city=False, town=False,
    metro=False, compute_height=True, replace_z=False,
    fine=False, extra_fine=False,
    coarse=False, extra_coarse=False, cores=1,
    non_ground_unchanged=False):
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
    city           [bool] - uses very wide radius to find ground
    town           [bool] - uses wide radius to find ground
    metro          [bool] - uses widest radius to find ground
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
    import earthtools as et

    # parse input params
    fparams = []

    # determine if -o, -odir, or -odix should be set
    #  if -o
    if output:
        fparams.append(et.fn.strJoin(['-o', output]))

    else:
        # if -odir
        if odir:
            fparams.append(et.fn.strJoin(['-odir', odir]))
        else:
            fparams.append(et.fn.strJoin(['-odir', et.params.outdir]))

        # if -odix
        if odix:
            fparams.append(et.fn.strJoin(['-odix', odix]))

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
    fparams.append(et.fn.strJoin(['-cores', cores]))

    # keep default classification info
    if non_ground_unchanged:
        fparams.append('-non_ground_unchanged')

    # add additional parameters passed through etc keyword
    if etc:
        fparams.append(et.read.etc(etc))

    # concatenate params list to string
    fparams = et.fn.strJoin(fparams)

    # join and run the command
    ocmd = et.fn.strJoin([raw.lasground, '-i', inputs, fparams])

    # report status
    print("[ STATUS ]: Running lasground")
    print("[ STATUS ]: %s" % ocmd)

    # run the command
    run(ocmd)

# lasheight
def lasheight(inputs, output='', etc='', odir='', odix='',
    replace_z=False, olaz=True, cores=1, drop_above=None):
    """
    calculates height of las files

    syntax: lasheight(input[s], output='', etc='', odir='',
      odix='', replace_z=False, olaz=True, cores=1)

    inputs    [string] - the input las/laz file(s). accepts wild
                         cards. enter multiple files as one string.
    output    [string] - the output merged las/laz file.
    etc       [string] - additional command line params. enter all
                         as a scalar string.
    odir      [string] - the output directory. superceded by using 'output' variable
    odix      [string] - a string to append to each output. superceded by using 'output' variable
    replace_z   [bool] - replaces the z value with the height above
                         ground. default = False
    olaz        [bool] - ensures output format is laz. default = True
    cores      [float] - sets the number of processors to use
    drop_above [float] - set this to drop points above a certain height
    """
    import earthtools as et

    # parse input params
    fparams = []

    # determine if -o, -odir, or -odix should be set
    #  if -o
    if output:
        fparams.append(et.fn.strJoin(['-o', output]))

    else:
        # if -odir
        if odir:
            fparams.append(et.fn.strJoin(['-odir', odir]))
        else:
            fparams.append(et.fn.strJoin(['-odir', et.params.outdir]))

        # if -odix
        if odix:
            fparams.append(et.fn.strJoin(['-odix', odix]))

    # add replace_z to param list if set
    if replace_z:
        fparams.append('-replace_z')

    # add olaz parameter if set
    if olaz:
        fparams.append('-olaz')

    # set to drop certain points
    if drop_above:
        fparams.append(et.fn.strJoin(['-drop_above', drop_above]))

    # set number of cores to use
    fparams.append(et.fn.strJoin(['-cores', cores]))

    # add additional parameters passed through etc keyword
    if etc:
        fparams.append(et.read.etc(etc))

    # concatenate params list to string
    fparams = et.fn.strJoin(fparams)

    # join and run the command
    ocmd = et.fn.strJoin([raw.lasheight, '-i', inputs, fparams])

    # report status
    print("[ STATUS ]: Running lasheight")
    print("[ STATUS ]: %s" % ocmd)

    # run the command
    run(ocmd)

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
    import earthtools as et

    # parse input params
    fparams = []

    # add additional parameters passed through etc keyword
    if etc:
        fparams.append(et.read.etc(etc))

    # concatenate params list to string
    fparams = et.fn.strJoin(fparams)

    # join and run the command
    ocmd = et.fn.strJoin([raw.lasindex, '-i', inputs, fparams])

    # report status
    print("[ STATUS ]: Running lasindex")
    print("[ STATUS ]: %s" % ocmd)

    # run the command
    run(ocmd)

# lasinfo

# lasmerge
def lasmerge(inputs, output, etc='',
    rescale=[0.01, 0.01, 0.01], olaz=True):
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
    import earthtools as et

    # parse input params
    fparams = []

    # add rescaling to param list if set
    if rescale:
        fparams.append(et.fn.strJoin(['-rescale', et.fn.strJoin(rescale)]))

    # add olaz parameter if set
    if olaz:
        fparams.append('-olaz')

    # add additional parameters passed through etc keyword
    if etc:
        fparams.append(et.read.etc(etc))

    # concatenate params list to string
    fparams = et.fn.strJoin(fparams)

    # join and run the command
    ocmd = et.fn.strJoin([raw.lasmerge, '-i', inputs, '-o', output,
            fparams])

    # report status
    print("[ STATUS ]: Running lasmerge")
    print("[ STATUS ]: %s" % ocmd)

    # run the command
    run(ocmd)

# lastile
def lastile(inputs, output='', odir='', etc='', odix='', tile_size=1000,
    buffer=100, olaz=True, cores=1, reversible=True, reverse_tiling=False):
    """
    tiles las files

    syntax: lastile(inputs, odir, etc=etc, odix=odix, tile_size=1000
      buffer=100, olaz=True, cores=1, reversible=True, reverse_tiling=False)

    inputs   [string] - the input las/laz file(s). accepts wild
                        cards. enter multiple files as one string.
    odir     [string] - the output directory for the tiles. if not set,
                        outputs tiles to the default output directory
                        (set under et.params)
    odix     [string] - a string to append to each output. superceded by using 'output' variable
    etc      [string] - additional command line params. enter all
                        as a scalar string.
    tile_size [float] - the output tile size in units of the horizontal
                        projection. default = 1000
    buffer    [float] - the buffer surrounding each tile. default = 100
    cores     [float] - the number of processors to use
    reversible [bool] - set this if you want to be able to re-assemble the
                        original file after tiling
    reverse_tiling [bool] - use this to re-assemble the original file from
                        tiles using the -reversible option
    """
    import earthtools as et

    # parse input params
    fparams = []

    # check if output, odir or odix are set
    if output:
        fparams.append(et.fn.strJoin(['-o', output]))

    else:
        # if -odir
        if odir:
            fparams.append(et.fn.strJoin(['-odir', odir]))
        else:
            fparams.append(et.fn.strJoin(['-odir', et.params.outdir]))

        # if -odix
        if odix:
            fparams.append(et.fn.strJoin(['-odix', odix]))

    # add tile_size to parameters list
    if tile_size:
        fparams.append(et.fn.strJoin(['-tile_size', tile_size]))

    # add buffer distance to parameters list
    if buffer:
        fparams.append(et.fn.strJoin(['-buffer', buffer]))

    # add olaz parameter if set
    if olaz:
        fparams.append('-olaz')

    # add n cores
    fparams.append(et.fn.strJoin(['-cores', cores]))

    # check if reversible
    if reversible:
        fparams.append('-reversible')

    # check if re-assembling tiles
    if reverse_tiling:
        fparams.append('-reverse_tiling')

    # add additional parameters passed through etc keyword
    if etc:
        fparams.append(et.read.etc(etc))

    # concatenate params list to string
    fparams = et.fn.strJoin(fparams)

    # join and run the command
    ocmd = et.fn.strJoin([raw.lastile, '-i', inputs, fparams])

    # report status
    print("[ STATUS ]: Running lastile")
    print("[ STATUS ]: %s" % ocmd)

    # run the command
    run(ocmd)

# laszip

###
# Orfeo Toolbox (OTB) commands
###

class otb:
    from earthtools import params

    def __init__(self):
        pass

    # BandMath
    @staticmethod
    def BandMath(inputs, output, exp, etc = '',
        ram = 2048, inxml = None):
        """
        performs band math on an image file

        syntax: otb.BandMath(inputs, output, exp,
          ram = 2048, inxml = None, etc=etc)

        inputs [list/string] - the input raster file
        output      [string] - the output unmixed file
        exp         [string] - the expression to use
        ram         [number] - the amount of ram to use. default is half system ram
        inxml       [string] - an optional xml application file
        etc         [string] - additional command line params. enter all
                               as a scalar string.

        full doc: http://www.orfeo-toolbox.org/Applications/BandMath.html
        """
        import earthtools as et

        # parse input params
        fparams = []

        # concatenate input list to string if not already done
        if type(inputs) is list:
            inputs = et.fn.strJoin(inputs)

        # specify how much memory to use
        if ram:
            fparams.append(et.fn.strJoin(['-ram', '%s' % ram]))

        # add xml app if set
        if inxml:
            fparams.append(et.fn.strJoin(['-inxml', inxml]))

        # add additional parameters passed through etc keyword
        if etc:
            fparams.append(et.read.etc(etc))

        # concatenate params list to string
        fparams = et.fn.strJoin(fparams)

        # join and run the command
        ocmd = et.fn.strJoin([raw.otbBandMath, '-il', inputs, '-out',
                output, '-exp', exp, fparams])

        # report status
        print("[ STATUS ]: Running OTB Band Math")
        print("[ STATUS ]: %s" % ocmd)

        # run the command
        run(ocmd)

    # ConcatenateImages
    @staticmethod
    def ConcatenateImages(inputs, output, etc = '',
        ram = 2058, inxml = None):
        """
        stacks raster bands

        syntax: otb.ConcatenateImages(inputs, output,
          ram = 2048, inxml = None, etc=etc)

        inputs [list/string] - the input raster files
        output      [string] - the output unmixed file
        ram         [number] - the amount of ram to use. default is half system ram
        inxml       [string] - an optional xml application file
        etc         [string] - additional command line params. enter all
                               as a scalar string.

        full doc: http://www.orfeo-toolbox.org/Applications/ConcatenateImages.html
        """
        import earthtools as et

        # parse input params
        fparams = []

        # concatenate input list to string if not already done
        if type(inputs) is list:
            inputs = et.fn.strJoin(inputs)

        # specify how much memory to use
        if ram:
            fparams.append(et.fn.strJoin(['-ram', '%s' % ram]))

        # add xml app if set
        if inxml:
            fparams.append(et.fn.strJoin(['-inxml', inxml]))

        # add additional parameters passed through etc keyword
        if etc:
            fparams.append(et.read.etc(etc))

        # concatenate params list to string
        fparams = et.fn.strJoin(fparams)

        # join and run the command
        ocmd = et.fn.strJoin([raw.otbConcatenateImages, '-il', inputs, '-out',
                output, fparams])

        # report status
        print("[ STATUS ]: Running OTB Concatenate Images")
        print("[ STATUS ]: %s" % ocmd)

        # run the command
        run(ocmd)

    # HyperspectralUnmixing
    @staticmethod
    def HyperspectralUnmixing(inputs, output, endmembers, etc='',
        ua='ucls', inxml = None):
        """
        performs spectral unmixing on an image file

        syntax: otb.HyperspectralUnmixing(inputs, output,
          endmembers, ua = 'ucls', inxml = None, etc=etc)

        inputs     [string] - the input raster file
        output     [string] - the output unmixed file
        endmembers [string] - the endmember image file
        ua         [string] - the unmixing algorithm (ucls/nnls/isra/mdmdnmf)
        inxml      [string] - an optional xml application file
        etc        [string] - additional command line params. enter all
                              as a scalar string.

        full doc: http://www.orfeo-toolbox.org/Applications/HyperspectralUnmixing.html
        """
        import earthtools as et

        # parse input params
        fparams = []

        # add xml app parameter if set
        if inxml:
            fparams.append(et.fn.strJoin(['-inxml', inxml]))

        # add additional parameters passed through etc keyword
        if etc:
            fparams.append(et.read.etc(etc))

        # concatenate params list to string
        fparams = et.fn.strJoin(fparams)

        # join and run the command
        ocmd = et.fn.strJoin([raw.otbHyperspectralUnmixing, '-in', inputs, '-out',
                output, '-ie', endmembers, '-ua', ua, fparams])

        # report status
        print("[ STATUS ]: Running OTB Spectral Unmixing")
        print("[ STATUS ]: %s" % ocmd)

        # run the command
        run(ocmd)

###
# SAGA commands
###

class saga:
    def __init__(self):
        pass

    @staticmethod
    def grid_calculs(module, parameters, etc=''):
        """
        """

        # add additional parameters passed through etc keyword
        if etc:
            fparams.append(et.read.etc(etc))

    @staticmethod
    def grid_filter(module, parameters, etc=''):
        """
        """

        # add additional parameters passed through etc keyword
        if etc:
            fparams.append(et.read.etc(etc))

    @staticmethod
    def imagery_segmentation(module, parameters, etc=''):
        """
        """

        # add additional parameters passed through etc keyword
        if etc:
            fparams.append(et.read.etc(etc))

    @staticmethod
    def io_gdal(module, parameters, etc=''):
        """
        """

        # add additional parameters passed through etc keyword
        if etc:
            fparams.append(et.read.etc(etc))

    @staticmethod
    def shapes_grid(module, parameters, etc=''):
        """
        """

        # add additional parameters passed through etc keyword
        if etc:
            fparams.append(et.read.etc(etc))

    @staticmethod
    def shapes_polygons(module, parameters, etc=''):
        """
        """

        # add additional parameters passed through etc keyword
        if etc:
            fparams.append(et.read.etc(etc))

    @staticmethod
    def shapes_tools(module, parameters, etc=''):
        """
        """

        # add additional parameters passed through etc keyword
        if etc:
            fparams.append(et.read.etc(etc))

    @staticmethod
    def ta_lighting(module, parameters, etc=''):
        """
        """

        # add additional parameters passed through etc keyword
        if etc:
            fparams.append(et.read.etc(etc))
