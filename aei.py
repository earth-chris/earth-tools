#!/bin/env/python2.7
##
## The wonderous objects and functions for [AEI]
##
###################
##
import sys
import platform
import threading
import multiprocessing
from subprocess import call, STDOUT
import gdal as gdal
import numpy as np

## allow python exceptions
gdal.UseExceptions()

## get environment variables
environ = platform.os.environ
pathsep = platform.os.pathsep
pathcat = platform.os.path.join
system  = platform.system()

try:
    ltbase = environ['LTBASE']
except KeyError:
    ltbase = ''
    
try:
    gdalbase = environ['GDALBASE']
except KeyError:
    gdalbase = ''
    
## set up commands to call binaries
cmd_lasmerge = {
    'Windows' : [pathcat(ltbase, 'lasmerge.exe')], 
    'Linux'   : ['wine', pathcat(ltbase, 'lasmerge.exe')]
    }[system]
cmd_las2las = {
    'Windows' : [pathcat(ltbase, 'las2las.exe')],
    'Linux'   : ['wine', pathcat(ltbase, 'las2las.exe')]
    }[system]
cmd_lasnoise = {
    'Windows' : [pathcat(ltbase, 'lasnoise.exe')],
    'Linux'   : ['wine', pathcat(ltbase, 'lasnoise.exe')]
    }[system]
cmd_lastile = {
    'Windows' : [pathcat(ltbase, 'lastile.exe')],
    'Linux'   : ['wine', pathcat(ltbase, 'lastile.exe')]
    }[system]
cmd_lasclassify = {
    'Windows' : [pathcat(ltbase, 'lasclassify.exe')],
    'Linux'   : ['wine', pathcat(ltbase, 'lasclassify.exe')]
    }[system]
cmd_lasboundary = {
    'Windows' : [pathcat(ltbase, 'lasboundary.exe')], 
    'Linux'   : ['wine', pathcat(ltbase, 'lasboundary.exe')]
    }[system]
cmd_lasclip = {
    'Windows' : [pathcat(ltbase, 'lasclip.exe')], 
    'Linux'   : ['wine', pathcat(ltbase, 'lasclip.exe')]
    }[system]
cmd_lasheight = {
    'Windows' : [pathcat(ltbase, 'lasheight.exe')], 
    'Linux'   : ['wine', pathcat(ltbase, 'lasheight.exe')]
    }[system]
cmd_lasground = {
    'Windows' : [pathcat(ltbase, 'lasground.exe')], 
    'Linux'   : ['wine', pathcat(ltbase, 'lasground.exe')]
    }[system]
cmd_lascanopy = {
    'Windows' : [pathcat(ltbase, 'lascanopy.exe')], 
    'Linux'   : ['wine', pathcat(ltbase, 'lascanopy.exe')]
    }[system]
cmd_las2dem = {
    'Windows' : [pathcat(ltbase, 'las2dem.exe')], 
    'Linux'   : ['wine', pathcat(ltbase, 'las2dem.exe')]
    }[system]
cmd_ogr2ogr = {
    'Windows' : [pathcat(gdalbase, 'ogr2ogr.exe')],
    'Linux'   : [pathcat(gdalbase, 'ogr2ogr')]
    }[system]
cmd_gdalbuildvrt = {
    'Windows' : [pathcat(gdalbase, 'gdalbuildvrt.exe')],
    'Linux'   : [pathcat(gdalbase, 'gdalbuildvrt')]
    }[system]
cmd_gdalwarp = {
    'Windows' : [pathcat(gdalbase, 'gdalwarp.exe')],
    'Linux'   : [pathcat(gdalbase, 'gdalwarp')]
    }[system]
cmd_gdal_translate = {
    'Windows' : [pathcat(gdalbase, 'gdal_translate.exe')],
    'Linux'   : [pathcat(gdalbase, 'gdal_translate')]
    }[system]
cmd_gdal_rasterize = {
    'Windows' : [pathcat(gdalbase, 'gdal_rasterize.exe')],
    'Linux'   : [pathcat(gdalbase, 'gdal_rasterize')]
    }[system]

###################
## 
## user editable configuration params and defaults
##
###################
##
params = {
    'cores' : multiprocessing.cpu_count()/-1, # num threads for processing
    'gdalbase' : gdalbase, # gdal path
    'lt' : 'laz', # las/laz type
    'ltbase' : ltbase, # lastools path
    'outdir' : 'G:\\AEI\\Data\\processed\\default\\', # default final product directory
    'ot' : 'Float32', # data type
    'procdir' : 'G:\\AEI\\Data\\scratch\\default\\' # default scratch product directory
    }
    
## set command class
class command(list):
    def __init__(self, *args, **kwargs):
        super(command, self).__init__(args[0])
        return
    def __str__(self):
        return "%s\n"%(" ".join(self))
    def write(self,outfile=sys.stdout):
        outfile.write(str(self))
        return
    def run(self,outfile=sys.stdout):
        try:
            ret = call(strjoin(self),stdout=outfile,stderr=STDOUT)
        except:
            print('Unable to call command %s'%(self))
        outfile.flush()
        if ret > 0:
            raise Exception("Failed at command (return code %d): %s"%(ret,str(self)))
        return ret
        
# set command merger
def cat_cmd(cmd, args):
    parts = command(cmd)
    parts.extend(args)
    return parts
    
# set basic tools
def strjoin(list, join=' '):
    outstr = join.join(map(str,list))
    return outstr

def parse_etc(etc): 
    if type(etc) is list:
        etc=strjoin(etc)
    if type(etc) is not str:
        print('Unable to parse the etc argument.')
        type(etc)
        etc=''
    return etc

#########################################
## build executable commands
#########################################

# LASTOOLS stuff
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
    # parse input params
    fparams = []
    if rescale:
        fparams.append(strjoin(['-rescale', strjoin(rescale)]))
    fparams.append(parse_etc(etc))
    fparams = strjoin(fparams)
    cmd = cat_cmd(cmd_lasmerge,['-i', inputs, '-o', output,
            fparams])
    ret=command.run(cmd)
    return ret

# GDAL stuff    
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
    # set up return string if user does not set inputs correctly
    if not (inputs or output):
        print("""
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
              )
        return -1
        
    # add quotes to input and output strings to handle wildcards/escape chars
    inputs = strjoin(['"', inputs, '"'], join='')
    output = strjoin(['"', output, '"'], join='')
    
    # parse input params
    fparams = []
    if dstnodata:
        if not isinstance(dstnodata, (int, long, float, complex)):
            print('Unable to parse -dstnodata option %s' % (dstnodata))
            print('Must be scalar int, long or float')
            return 01
        else:
            fparams.append(strjoin(['-dstnodata', dstnodata]))
    if srcnodata:
        if not isinstance(srcnodata, (int, long, float, complex)):
            print('Unable to parse -srcnodata option %s' % (srcnodata))
            print('Must be scalar int, long or float')
            return -1
        else:
            fparams.append(strjoin(['-srcnodata', srcnodata]))
    if n_threads:
        if not isinstance(n_threads, (int, long, float, complex)):
            print('Unable to parse -n_threads option %s' % (n_threads))
            print('Must be scalar int, long or float')
            return -1
        else:
            fparams.append(strjoin(['-multi -wo NUM_THREADS',n_threads], join='='))
        multi = False
    if multi:
        fparams.append('-multi')
    if overwrite:
        fparams.append('-overwrite')
    if of:
        if type(of) is not str:
            print('Unable to parse -of option %s' % (of))
            print('Must be scalar string')
            return -1
        else:
            fparams.append(strjoin(['-of', of]))
    if ot:
        ot=get_ot(ot)
        fparams.append(strjoin(['-ot', ot]))
    if s_srs:
        if type(s_srs) is not str:
            print('Unable to parse -s_srs option %s' % (s_srs))
            print('Must be scalar string')
            return -1
        else:
            fparams.append(strjoin(['-s_srs', s_srs]))
    if t_srs:
        if type(t_srs) is not str:
            print('Unable to parse -t_srs option %s' % (t_srs))
            print('Must be a scalar string')
            return -1
        else:
            fparams.append(strjoin(['-t_srs', t_srs]))
    if te:
        if type(te) is not list and len(te) is not 4:
            print('Unable to parse -te option %s' % (te))
            print('Must be a 4-element list')
            return -1
        else:
            fparams.append(strjoin(['-te', strjoin(te)]))
    if tr:
        if type(te) is not list and len (tr) is not 2:
            print('Unable to parse -tr option %s' % (tr))
            print('Must be a 2-element list')
            return -1
        else:
            fparams.append(strjoin(['-tr', strjoin(tr)]))

    fparams.append(parse_etc(etc))
    fparams=strjoin(fparams)
    
    cmd = [cat_cmd(cmd_gdalwarp,[fparams, inputs, output])]
    ret=command.run(cmd)
    return ret

# set up the return functions from params. can be done by param
#  or with just returning all params
def get_params():
    return params

def get_ot(arg):
    """
    reads the argument passed by -ot and returns the output
data type for the params dictionary.

    accepts both python style data type codes (e.g. Float32, Int16) and IDL
    style data type codes (e.g. 4 for Float32, 1 for Byte)

    syntax: read_ot(arg)
    """
    # Python style data type flags
    try:
        if (str(arg).lower() == 'byte'):
            params['ot'] = 'Byte'
        elif (str(arg).lower() == 'int16'):
            params['ot'] = 'Int16'
        elif (str(arg).lower() == 'int32'):
            params['ot'] = 'Int32'
        elif (str(arg).lower() == 'uint16'):
            params['ot'] = 'UInt16'
        elif (str(arg).lower() == 'uint32'):
            params['ot'] = 'UInt32'
        elif (str(arg).lower() == 'float32'):
            params['ot'] = 'Float32'
        elif (str(arg).lower() == 'float64'):
            params['ot'] = 'Float64'
        elif (str(arg).lower() == 'cint16'):
            params['ot'] = 'CInt16'
        elif (str(arg).lower() == 'cint32'):
            params['ot'] = 'CInt32'
        elif (str(arg).lower() == 'cfloat32'):
            params['ot'] = 'CFloat32'
        elif (str(arg).lower() == 'cfloat64'):
            params['ot'] = 'CFloat64'
        # IDL style data type flags
        elif (str(arg).lower() == '1'):
            params['ot'] = 'Byte'
        elif (str(arg).lower() == '2'):
            params['ot'] = 'Int16'
        elif (str(arg).lower() == '3'):
            params['ot'] = 'Int32'
        elif (str(arg).lower() == '4'):
            params['ot'] = 'Float32'
        elif (str(arg).lower() == '5'):
            params['ot'] = 'Float64'
        elif (str(arg).lower() == '6'):
            params['ot'] = 'CFloat32'
        elif (str(arg).lower() == '7'):
            print('Invalid -ot argument: %s, unable to output using String format' % (arg))
        elif (str(arg).lower() == '8'):
            print('Invalid -ot argument: %s, unable to output using Structure format' % (arg))
        elif (str(arg).lower() == '9'):
            params['ot'] = 'CFloat64'
        elif (str(arg).lower() == '10'):
            print('Invalid -ot argument: %s, unable to output using Pointer format' % (arg))
        elif (str(arg).lower() == '11'):
            print('Invalid -ot argument: %s unable to output using Object format' % (arg))
        elif (str(arg).lower() == '12'):
            params['ot'] = 'UInt16'
        elif (str(arg).lower() == '13'):
            params['ot'] = 'UInt32'
        elif (str(arg).lower() == '14'):
            params['ot'] = 'Int64'
        elif (str(arg).lower() == '15'):
            params['ot'] = 'UInt64'
        # LAS data type flags
        elif (str(arg).lower() == 'las'):
            params['lt'] = 'las'
        elif (str(arg).lower() == 'laz'):
            params['lt'] = 'laz'
    except ValueError:
        print('Unable to parse -ot argument: %s' % (arg))
    return params['ot']
   
def get_outdir():
    return params['outdir']
    
def get_procdir():
    return params['procdir']
    
if __name__ == "__main__":
    main()