#####
# sets up paths and calls to external executables.
# c. 2016 Christopher Anderson
#####

# set up commands to call binaries
class getPaths:
    def __init__(self):
        
        from aei import params
        
        self.lasmerge = {
            'Windows' : [params.pathcat(params.ltbase, 'lasmerge.exe')], 
            'Linux'   : ['wine', params.pathcat(params.ltbase, 'lasmerge.exe')]
            }[params.system]
        self.las2las = {
            'Windows' : [params.pathcat(params.ltbase, 'lasmerge.exe')], 
            'Linux'   : ['wine', params.pathcat(params.ltbase, 'lasmerge.exe')]
            }[params.system]
        self.lasnoise = {
            'Windows' : [params.pathcat(params.ltbase, 'lasnoise.exe')],
            'Linux'   : ['wine', params.pathcat(params.ltbase, 'lasnoise.exe')]
            }[params.system]
        self.lastile = {
            'Windows' : [params.pathcat(params.ltbase, 'lastile.exe')],
            'Linux'   : ['wine', params.pathcat(params.ltbase, 'lastile.exe')]
            }[params.system]
        self.lasclassify = {
            'Windows' : [params.pathcat(params.ltbase, 'lasclassify.exe')],
            'Linux'   : ['wine', params.pathcat(params.ltbase, 'lasclassify.exe')]
            }[params.system]
        self.lasboundary = {
            'Windows' : [params.pathcat(params.ltbase, 'lasboundary.exe')], 
            'Linux'   : ['wine', params.pathcat(params.ltbase, 'lasboundary.exe')]
            }[params.system]
        self.lasclip = {
            'Windows' : [params.pathcat(params.ltbase, 'lasclip.exe')], 
            'Linux'   : ['wine', params.pathcat(params.ltbase, 'lasclip.exe')]
            }[params.system]
        self.lasheight = {
            'Windows' : [params.pathcat(params.ltbase, 'lasheight.exe')], 
            'Linux'   : ['wine', params.pathcat(params.ltbase, 'lasheight.exe')]
            }[params.system]
        self.lasground = {
            'Windows' : [params.pathcat(params.ltbase, 'lasground.exe')], 
            'Linux'   : ['wine', params.pathcat(params.ltbase, 'lasground.exe')]
            }[params.system]
        self.lascanopy = {
            'Windows' : [params.pathcat(params.ltbase, 'lascanopy.exe')], 
            'Linux'   : ['wine', params.pathcat(params.ltbase, 'lascanopy.exe')]
            }[params.system]
        self.las2dem = {
            'Windows' : [params.pathcat(params.ltbase, 'las2dem.exe')], 
            'Linux'   : ['wine', params.pathcat(params.ltbase, 'las2dem.exe')]
            }[params.system]
        self.ogr2ogr = {
            'Windows' : [params.pathcat(params.gdalbase, 'ogr2ogr.exe')],
            'Linux'   : [params.pathcat(params.gdalbase, 'ogr2ogr')]
            }[params.system]
        self.gdalbuildvrt = {
            'Windows' : [params.pathcat(params.gdalbase, 'gdalbuildvrt.exe')],
            'Linux'   : [params.pathcat(params.gdalbase, 'gdalbuildvrt')]
            }[params.system]
        self.gdalwarp = {
            'Windows' : [params.pathcat(params.gdalbase, 'gdalwarp.exe')],
            'Linux'   : [params.pathcat(params.gdalbase, 'gdalwarp')]
            }[params.system]
        self.gdal_translate = {
            'Windows' : [params.pathcat(params.gdalbase, 'gdal_translate.exe')],
            'Linux'   : [params.pathcat(params.gdalbase, 'gdal_translate')]
            }[params.system]
        self.gdal_rasterize = {
            'Windows' : [params.pathcat(params.gdalbase, 'gdal_rasterize.exe')],
            'Linux'   : [params.pathcat(params.gdalbase, 'gdal_rasterize')]
            }[params.system]

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
    import aei
    
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
            fparams.append(strJoin(['-dstnodata', dstnodata]))
    
    if srcnodata:
        if not isinstance(srcnodata, (int, long, float, complex)):
            print('[ ERROR ]: Unable to parse -srcnodata option %s' % (srcnodata))
            print('[ ERROR ]: Must be scalar int, long or float')
            return -1
        else:
            fparams.append(strJoin(['-srcnodata', srcnodata]))
    
    if n_threads:
        if not isinstance(n_threads, (int, long, float)):
            print('[ ERROR ]: Unable to parse -n_threads option %s' % (n_threads))
            print('[ ERROR ]: Must be scalar int, long or float')
            return -1
        else:
            fparams.append(strJoin(['-multi -wo NUM_THREADS',n_threads], join='='))
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
            fparams.append(strJoin(['-of', of]))
    if ot:
        ot=read.ot(ot)
        fparams.append(strJoin(['-ot', ot]))
    
    if s_srs:
        if type(s_srs) is not str:
            print('[ ERROR ]: Unable to parse -s_srs option %s' % (s_srs))
            print('[ ERROR ]: Must be scalar string')
            return -1
        else:
            fparams.append(strJoin(['-s_srs', s_srs]))
    
    if t_srs:
        if type(t_srs) is not str:
            print('[ ERROR ]: Unable to parse -t_srs option %s' % (t_srs))
            print('[ ERROR ]: Must be a scalar string')
            return -1
        else:
            fparams.append(strJoin(['-t_srs', t_srs]))
    
    if te:
        if type(te) is not list and len(te) is not 4:
            print('[ ERROR ]: Unable to parse -te option %s' % (te))
            print('[ ERROR ]: Must be a 4-element list')
            return -1
        else:
            fparams.append(strJoin(['-te', strJoin(te)]))
    
    if tr:
        if type(te) is not list and len (tr) is not 2:
            print('[ ERROR ]: Unable to parse -tr option %s' % (tr))
            print('[ ERROR ]: Must be a 2-element list')
            return -1
        else:
            fparams.append(strjoin(['-tr', strJoin(tr)]))

    fparams.append(aei.read.etc(etc))
    fparams=aei.strJoin(fparams)
    
    # join and run the command
    ocmd = aei.strJoin([raw.gdalwarp[0], fparams, inputs, output])
    print("[ STATUS ]: Running gdalwarp")
    print("[ STATUS ]: %s" % ocmd)
    os.system(ocmd)
    
    return(ocmd)

###
# LAStools commands
###
las2dem
las2las
lasclassify
lasgrid
lasheight

def lasmerge(inputs, output, etc='',
    rescale=[0.01, 0.01, 0.01]
    ):
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
                      [0.01, 0.01, 0.01]. set to False to not
                      rescale
    """
    import aei as aei
    
    # parse input params
    fparams = []
    
    # add rescaling to param list if set
    if rescale:
        fparams.append(aei.strJoin(['-rescale', aei.strJoin(rescale)]))
    
    # add additional parameters passed through etc keyword
    fparams.append(aei.read.etc(etc))
    
    # concatenate params list to string
    fparams = aei.strJoin(fparams)
    
    # join and run the command
    ocmd = aei.strjoin(raw.lasmerge,['-i', inputs, '-o', output,
            fparams])
    
    # report status
    print("[ STATUS ]: Running lasmerge")
    print("[ STATUS ]: %s" % ocmd)
    
    os.system(ocmd)
    
    return ocmd

def lastile(inputs, outdir='', etc='', tile_size=1000, buffer=100
    ):
    """
    tiles las files
    
    syntax: lastile(inputs, outdir=outdir, etc=etc, tile_size=tile_size)
    
    inputs   [string] - the input las/laz file(s). accepts wild
                        cards. enter multiple files as one string.
    outdir   [string] - the output directory for the tiles. if not set,
                        outputs tiles to the default output directory
                        (set under aei.params)
    etc      [string] - additional command line params. enter all
                        as a scalar string. 
    tile_size [float] - the output tile size in units of the horizontal
                        projection. default = 1000
    buffer    [float] - the buffer surrounding each tile. default = 100
    """
    import aei as aei
    
    # parse input params
    fparams = []
                
    # add tile_size to parameters list
    fparams.append(aei.strJoin(['-tile_size', tile_size]))
    
    # add buffer distance to parameters list
    fparams.append(aei.strJoin(['-buffer', buffer]))
    
    # add additional parameters passed through etc keyword
    fparams.append(aei.read.etc(etc))
    
    # concatenate params list to string
    fparams = aei.strJoin(fparams)
    
    # join and run the command
    ocmd = aei.strjoin(raw.lasmerge,['-i', inputs, '-o', output,
            fparams])
    
    # report status
    print("[ STATUS ]: Running lasmerge")
    print("[ STATUS ]: %s" % ocmd)
    
    os.system(ocmd)
    
    return ocmd