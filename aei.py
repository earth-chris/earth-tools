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
        ret = call(strjoin(self),stdout=outfile,stderr=STDOUT)
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

# build executable commands
def lasmerge(inputs, outputs,
    rescale=[0.01, 0.01, 0.01],
    etc=''
    ):
    """
    merges and rescales las files

    syntax: lasmerge(inputs, output, etc, rescale=rescale)
    """
    if rescale == False:
        rs = ''
    if rescale:
        rs=strjoin(['-rescale',strjoin(rescale)])
    if type(etc) is not str:
        print('Unable to parse the etc argument.')
        type(etc)
        etc=''
    cmd = cat_cmd(cmd_lasmerge,['-i', inputs, '-o', outputs,
            rs, etc])
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