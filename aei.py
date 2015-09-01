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

###################
## 
## user editable configuration params and defaults
##
###################
##
params = {
    "lt" : "laz", # las/laz type
    "outdir" : "G:\\AEI\\Data\\processed\\default\\", # default final product directory
    "ot" : "Float32", # data type
    "procdir" : "G:\\AEI\\Data\\scratch\\default\\" # default scratch product directory
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
    # LAS data type flags
    elif (str(arg).lower() == 'las'):
        params['lt'] = 'las'
    elif (str(arg).lower() == 'laz'):
        params['lt'] = 'laz'
    
def get_outdir():
    return params['outdir']
    
def get_procdir():
    return params['procdir']
    
if __name__ == "__main__":
    main()