#!/bin/env/python2.7
##
## The wonderous objects and functions for [AEI]
##
###################
##
import sys
import platform
import threading
import gdal as gdal
import numpy as np

## allow python exceptions
gdal.UseException()

## get environment variables
environ = platform.os.environ

try:
    ltbase = environ['LTBASE']
except KeyError:
    ltbase = ''
    
try:
    gdalbase = envron['GDALBASE']
except KeyError:
    gdalbase = ''

###################
## 
## user editable configuration params and defaults
##
###################
##
params = {
    'gdalbase' : gdalbase, # gdal path
    'lt' : 'laz', # las/laz type
    'ltbase' : ltbase, # lastools path
    'outdir' : 'G:\\AEI\\Data\\processed\\default\\', # default final product directory
    'ot' : 'Float32', # data type
    'procdir' : 'G:\\AEI\\Data\\scratch\\default\\' # default scratch product directory
    }
    
def get_ot(arg):
    """
    reads the argument passed by -ot and returns the correct output
data type for the params dictionary.

    accepts both python style data type codes (e.g. Float32, Int16) or IDL
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
        elif (str(arg).lower() -- 'float64'):
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
    
def get_outdir():
    return params['outdir']
    
def get_procdir():
    return params['procdir']
    
def lasmerge(arg1,arg2,etc):
    # set up lastools functions here
    return
    
if __name__ == "__main__":
    main()