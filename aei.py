#!/bin/env/python2.7
##
## The wonderous objects and functions for [AEI]
##
###################
##
import os
import sys
import platform
import threading
import multiprocessing
from psutil import virtual_memory
from subprocess import call, STDOUT
import gdal as gdal
import numpy as np

## allow python exceptions
gdal.UseExceptions()

## get environment variables
environ = platform.os.environ
pathsep = platform.os.path.sep
pathcat = platform.os.path.join
system  = platform.system()

###################
## set up system params and defaults
###################

# set up class to return os/fs/computing defaults
class getParams:
    def __init__(self):
        try:
            ltbase = environ['LTBASE']
        except KeyError:
          ltbase = ''

        try:
            gdalbase = environ['GDALBASE']
        except KeyError:
            gdalbase = ''
        
        try:
            outdir = environ['OUTPUT_DIR']
        except KeyError:
            outdir = './'
        
        try:
            scratchdir = environ['SCRATCH_DIR']
        except KeyError:
            scratchdir = './'

        self.cores = multiprocessing.cpu_count()
        self.gdalbase = gdalbase
        self.lt = 'laz'
        self.ltbase = ltbase
        self.mem = virtual_memory().total
        self.outdir = outdir
        self.ot = 'Float32'
        self.scratchdir = scratchdir

params = getParams()

## set up commands to call binaries
class getCommands:
    def __init__(self):
        self.lasmerge = {
            'Windows' : [pathcat(params.ltbase, 'lasmerge.exe')], 
            'Linux'   : ['wine', pathcat(params.ltbase, 'lasmerge.exe')]
            }[system]
        self.las2las = {
            'Windows' : [pathcat(params.ltbase, 'lasmerge.exe')], 
            'Linux'   : ['wine', pathcat(params.ltbase, 'lasmerge.exe')]
            }[system]
        self.lasnoise = {
            'Windows' : [pathcat(params.ltbase, 'lasnoise.exe')],
            'Linux'   : ['wine', pathcat(params.ltbase, 'lasnoise.exe')]
            }[system]
        self.lastile = {
            'Windows' : [pathcat(params.ltbase, 'lastile.exe')],
            'Linux'   : ['wine', pathcat(params.ltbase, 'lastile.exe')]
            }[system]
        self.lasclassify = {
            'Windows' : [pathcat(params.ltbase, 'lasclassify.exe')],
            'Linux'   : ['wine', pathcat(params.ltbase, 'lasclassify.exe')]
            }[system]
        self.lasboundary = {
            'Windows' : [pathcat(params.ltbase, 'lasboundary.exe')], 
            'Linux'   : ['wine', pathcat(params.ltbase, 'lasboundary.exe')]
            }[system]
        self.lasclip = {
            'Windows' : [pathcat(params.ltbase, 'lasclip.exe')], 
            'Linux'   : ['wine', pathcat(params.ltbase, 'lasclip.exe')]
            }[system]
        self.lasheight = {
            'Windows' : [pathcat(params.ltbase, 'lasheight.exe')], 
            'Linux'   : ['wine', pathcat(params.ltbase, 'lasheight.exe')]
            }[system]
        self.lasground = {
            'Windows' : [pathcat(params.ltbase, 'lasground.exe')], 
            'Linux'   : ['wine', pathcat(params.ltbase, 'lasground.exe')]
            }[system]
        self.lascanopy = {
            'Windows' : [pathcat(params.ltbase, 'lascanopy.exe')], 
            'Linux'   : ['wine', pathcat(params.ltbase, 'lascanopy.exe')]
            }[system]
        self.las2dem = {
            'Windows' : [pathcat(params.ltbase, 'las2dem.exe')], 
            'Linux'   : ['wine', pathcat(params.ltbase, 'las2dem.exe')]
            }[system]
        self.ogr2ogr = {
            'Windows' : [pathcat(params.gdalbase, 'ogr2ogr.exe')],
            'Linux'   : [pathcat(params.gdalbase, 'ogr2ogr')]
            }[system]
        self.gdalbuildvrt = {
            'Windows' : [pathcat(params.gdalbase, 'gdalbuildvrt.exe')],
            'Linux'   : [pathcat(params.gdalbase, 'gdalbuildvrt')]
            }[system]
        self.gdalwarp = {
            'Windows' : [pathcat(params.gdalbase, 'gdalwarp.exe')],
            'Linux'   : [pathcat(params.gdalbase, 'gdalwarp')]
            }[system]
        self.gdal_translate = {
            'Windows' : [pathcat(params.gdalbase, 'gdal_translate.exe')],
            'Linux'   : [pathcat(params.gdalbase, 'gdal_translate')]
            }[system]
        self.gdal_rasterize = {
            'Windows' : [pathcat(params.gdalbase, 'gdal_rasterize.exe')],
            'Linux'   : [pathcat(params.gdalbase, 'gdal_rasterize')]
            }[system]

cmd = getCommands()

# set function to return a standardized data type used by gdal, based on idl/python inputs
def getOt(arg):
    """
    reads the argument passed by -ot and returns the output data type for the params dictionary.

    accepts both python style data type codes (e.g. Float32, Int16) and IDL
    style data type codes (e.g. 4 for Float32, 1 for Byte)

    syntax: getOt(arg)
    """
    # Python style data type flags
    try:
        if (str(arg).lower() == 'byte'):
            params.ot = 'Byte'
        elif (str(arg).lower() == 'int16'):
            params.ot = 'Int16'
        elif (str(arg).lower() == 'int32'):
            params.ot = 'Int32'
        elif (str(arg).lower() == 'uint16'):
            params.ot = 'UInt16'
        elif (str(arg).lower() == 'uint32'):
            params.ot = 'UInt32'
        elif (str(arg).lower() == 'float32'):
            params.ot = 'Float32'
        elif (str(arg).lower() == 'float64'):
            params.ot = 'Float64'
        elif (str(arg).lower() == 'cint16'):
            params.ot = 'CInt16'
        elif (str(arg).lower() == 'cint32'):
            params.ot = 'CInt32'
        elif (str(arg).lower() == 'cfloat32'):
            params.ot = 'CFloat32'
        elif (str(arg).lower() == 'cfloat64'):
            params.ot = 'CFloat64'
        # IDL style data type flags
        elif (str(arg).lower() == '1'):
            params.ot = 'Byte'
        elif (str(arg).lower() == '2'):
            params.ot = 'Int16'
        elif (str(arg).lower() == '3'):
            params.ot = 'Int32'
        elif (str(arg).lower() == '4'):
            params.ot = 'Float32'
        elif (str(arg).lower() == '5'):
            params.ot = 'Float64'
        elif (str(arg).lower() == '6'):
            params.ot = 'CFloat32'
        elif (str(arg).lower() == '7'):
            print('Invalid -ot argument: %s, unable to output using String format' % (arg))
        elif (str(arg).lower() == '8'):
            print('Invalid -ot argument: %s, unable to output using Structure format' % (arg))
        elif (str(arg).lower() == '9'):
            params.ot = 'CFloat64'
        elif (str(arg).lower() == '10'):
            print('Invalid -ot argument: %s, unable to output using Pointer format' % (arg))
        elif (str(arg).lower() == '11'):
            print('Invalid -ot argument: %s unable to output using Object format' % (arg))
        elif (str(arg).lower() == '12'):
            params.ot = 'UInt16'
        elif (str(arg).lower() == '13'):
            params.ot = 'UInt32'
        elif (str(arg).lower() == '14'):
            params.ot = 'Int64'
        elif (str(arg).lower() == '15'):
            params.ot = 'UInt64'
        # LAS data type flags
        elif (str(arg).lower() == 'las'):
            params.lt = 'las'
        elif (str(arg).lower() == 'laz'):
            params.lt = 'laz'
    except ValueError:
        print('Unable to parse -ot argument: %s' % (arg))
    return params.ot
   
# set command class
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

# set up class for spectral objects
class spectralObject:
    def __init__(self, n_spectra=1, n_wl=2151, type='', band_unit = '',
            band_quantity = '', band_centers = []):
        
        # set to asd type if no params set to change n_wl
        if n_wl == 2151:
            type = 'asd'
        
        # set up pre-defined types
        if (str(type).lower()) == 'asd':
            n_wl = 2151
            band_unit = 'Nanometers'
            band_quantity = 'Wavelength'
            band_centers = np.arange(350.,2501.)
        
        # return a list same size as number of spectra
        names = []
        for i in range(n_spectra):
            names.append('Spectrum ' + str(i))
        self.names = names
        
        # set up the band definitions
        try:
            self.band_unit = band_unit
        except NameError:
            pass
        try:
            self.band_quantity = band_quantity
        except NameError:
            pass
        try:
            self.band_centers = band_centers
        except NameError:
            pass
        #if band_unit:
        #    self.band_unit = band_unit
        #if band_quantity:
        #    self.band_quantity = band_quantity
        #if 'band_centers' in local():
        #    self.band_centers = band_centers
        #    n_wl = len(band_centers)
        
        # return an np array size of n spectra x n wavelengths
        self.spectra = np.zeros([n_spectra, n_wl])
        
    def remove_water_bands(self):
        """ 
        a function to set data from water absorption bands,
        i.e. 1350 - 1460 nm and 1790 - 1960 nm, to zero
        
        usage: self.remove_water_bands()
        """
        water_bands = [[1350.0, 1460.0], [1790.0, 1960.0]]
        
        # start with nir-swir1 transition
        gt = np.where(self.band_centers > water_bands[0][0])
        lt = np.where(self.band_centers < water_bands[0][1])
        nd = np.intersect1d(gt[0], lt[0])
        self.spectra[:,nd] = 0.0
        
        # then swir1-swir2 transition
        gt = np.where(self.band_centers > water_bands[1][0])
        lt = np.where(self.band_centers < water_bands[1][1])
        nd = np.intersect1d(gt[0], lt[0])
        self.spectra[:,nd] = 0.0
    
    def get_shortwave_bands(self, bands=[]):
        """
        a function that returns an index of the bands that encompass
        the shortwave range (350 - 2500 nm)
        
        usage: self.get_shortwave_bands(bands=[])
        
        returns: an index of bands to subset to the shortwave range
        """
        # set range to return in nanometers
        shortwave_range = [350., 2500.]
        
        # normalize if wavelength units are different
        if self.band_unit == 'Micrometers':
            shortwave_range /= 1000.
            
        # find overlapping range
        gt = np.where(self.band_centers > shortwave_range[0])
        lt = np.where(self.band_centers < shortwave_range[1])
        overlap = np.intersect1d(gt[0], lt[0])
        
        # return output
        return overlap
        
    def plot(self, inds = []):
        """
        plots the spectra using a standard plot format
        can be set to only plot selected spectra
        
        usage: self.plot(inds = [])
          where inds = the 0-based indices for spectra to plot
        """
        # import pyplot
        import matplotlib.pyplot as plt
        
        # set basic parameters
        plt.xlim((self.band_centers.min(), self.band_centers.max()))
        plt.xlabel('Wavelength (' + self.band_unit + ')')
        plt.ylabel('Reflectance (%)')
        plt.title('spectralObjet plot')
        
        # check if indices were set and valid. if not, plot all items
        if inds:
            if max(inds) > len(self.names):
                inds = range(0, len(self.names))
                print("invalid range set. using all spectra")
            if min(inds) < 0:
                inds = range(0, len(self.names))
                print("invalid range set. using all spectra")
        else:
            inds = range(0, len(self.names))
            
        # loop through each item to plot
        for i in inds:
            plt.plot(self.band_centers, self.spectra[i,:], 
                label = self.names[i])
            
        # add the legend with each spectrum's name
        plt.legend()
        
        # display the plot
        plt.show()
        
# set command merger
def catCmd(executable, args):
    parts = command(executable)
    parts.extend(args)
    return parts
    
# set function to join strings
def strJoin(list, join=' '):
    outstr = join.join(map(str,list))
    return outstr

# set function to join miscellaneous arguments
def parseEtc(etc): 
    if type(etc) is list:
        etc=strJoin(etc)
    if type(etc) is not str:
        print('Unable to parse the etc argument.')
        type(etc)
        etc=''
    return etc
    
# set function to read the ascii spectra from the joint fire science
#  program (https://www.frames.gov/partner-sites/assessing-burn-severity/spectral/spectral-library-southern-california/)
def readJFSC(infile):
    """
    reads the ascii format spectral data from the joint-fire-science-program
    and returns an object with the mean and +/- standard deviation
    reflectance data
    """
    if os.path.isfile(infile) and os.access(infile, os.R_OK):
        
        # create the spectral object
        s = spectralObject(1,type='asd')
        s.spectra_stdevm = np.zeros(s.spectra.shape)
        s.spectra_stdevp = np.zeros(s.spectra.shape)
        
        # open the file and read the data out
        f = open(infile, 'r')
        header = f.readline()
        i = 0
        for line in f:
            line = line.strip().split()
            s.spectra[0,i] = line[1]
            s.spectra_stdevp[0,i] = line[2]
            s.spectra_stdevm[0,i] = line[3]
            i += 1
        
        # close the file
        f.close()
        
        # return the spectral object
        return s
        
    else:
        print("Unable to read file %s", infile)
        return -1
        
def readUSGS(infile):
    """
    reads the ascii format spectral data from the joint-fire-science-program
    and returns an object with the mean and +/- standard deviation
    reflectance data
    """
    if os.path.isfile(infile) and os.access(infile, os.R_OK):
        
        # open the file and read header info
        f = open(infile, 'r')
        x_start = 'gibberish'
        for line in f:
            if x_start in line:
                break
            if "Name:" in line:
                spectrum_name = line.strip().split("Name:")[-1].strip()
            if "X Units:" in line:
                band_unit = line.strip().split()
                band_unit = band_unit[-1].strip('()').capitalize()
            if "Y Units:" in line:
                refl_unit = line.strip().split()
                refl_unit = refl_unit[-1].strip('()').capitalize()
            if "First X Value:" in line:
                x_start = line.strip().split()[-1]
            if "Number of X Values:" in line:
                n_values = int(line.strip().split()[-1])
            
        # now that we got our header info, create the arrays
        #  necessary for output
        band_centers = np.empty(n_values)
        reflectance = np.empty(n_values)
        
        line = line.strip().split()
        band_centers[0] = float(line[0])
        reflectance[0] = float(line[1])
        
        # resume reading through file
        i = 1
        for line in f:
            line = line.strip().split()
            band_centers[i] = float(line[0])
            reflectance[i] = float(line[1])
            i += 1
            
        # some files read last -> first wavelength. resort as necessary
        if band_centers[0] > band_centers[-1]:
            band_centers = band_centers[::-1]
            reflectance = reflectance[::1]
            
        # convert units to nanometers and scale 0-1
        if band_unit == 'Micrometers':
            band_centers *= 1000.
            band_unit = 'Nanometers'
        if refl_unit == 'Percent':
            reflectance /= 100.
            
        # create the spectral object
        s = spectralObject(1, n_values, band_centers = band_centers,
                band_unit = band_unit, band_quantity = 'Wavelength')
                
        # assign relevant values
        s.spectra[0] = reflectance
        if spectrum_name:
            s.names[0] = spectrum_name
        
        # close the file
        f.close()
        
        # return the final object    
        return s
        
    else:
        print("Unable to read file %s", infile)
        return -1

###################
## build executable command wrappers
###################

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
        fparams.append(strJoin(['-rescale', strJoin(rescale)]))
    fparams.append(parseEtc(etc))
    fparams = strJoin(fparams)
    call = catCmd(cmd.lasmerge,['-i', inputs, '-o', output,
            fparams])
    ret=command.run(call)
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
    # return the docstring if user does not set inputs correctly
    if not (inputs or output):
        print(gdalwarp.__doc__)
        return -1
        
    # parse input params
    fparams = []
    if dstnodata:
        if not isinstance(dstnodata, (int, long, float, complex)):
            print('Unable to parse -dstnodata option %s' % (dstnodata))
            print('Must be scalar int, long or float')
            return -1
        else:
            fparams.append(strJoin(['-dstnodata', dstnodata]))
    if srcnodata:
        if not isinstance(srcnodata, (int, long, float, complex)):
            print('Unable to parse -srcnodata option %s' % (srcnodata))
            print('Must be scalar int, long or float')
            return -1
        else:
            fparams.append(strJoin(['-srcnodata', srcnodata]))
    if n_threads:
        if not isinstance(n_threads, (int, long, float, complex)):
            print('Unable to parse -n_threads option %s' % (n_threads))
            print('Must be scalar int, long or float')
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
            print('Unable to parse -of option %s' % (of))
            print('Must be scalar string')
            return -1
        else:
            fparams.append(strJoin(['-of', of]))
    if ot:
        ot=getOt(ot)
        fparams.append(strJoin(['-ot', ot]))
    if s_srs:
        if type(s_srs) is not str:
            print('Unable to parse -s_srs option %s' % (s_srs))
            print('Must be scalar string')
            return -1
        else:
            fparams.append(strJoin(['-s_srs', s_srs]))
    if t_srs:
        if type(t_srs) is not str:
            print('Unable to parse -t_srs option %s' % (t_srs))
            print('Must be a scalar string')
            return -1
        else:
            fparams.append(strJoin(['-t_srs', t_srs]))
    if te:
        if type(te) is not list and len(te) is not 4:
            print('Unable to parse -te option %s' % (te))
            print('Must be a 4-element list')
            return -1
        else:
            fparams.append(strJoin(['-te', strJoin(te)]))
    if tr:
        if type(te) is not list and len (tr) is not 2:
            print('Unable to parse -tr option %s' % (tr))
            print('Must be a 2-element list')
            return -1
        else:
            fparams.append(strjoin(['-tr', strJoin(tr)]))

    fparams.append(parseEtc(etc))
    fparams=strJoin(fparams)
    
    call = catCmd(cmd.gdalwarp,[fparams, inputs, output])
    #ret=command.run(cmd)
    print('PRINTING: %s' % (call))
    return(call)
#    try:
#        ret = call(cmd, stdout=STDOUT, stderr=STDOUT)
#    except:
#        print('Unable to call command %s' % (cmd))
#        print('Trying as list')
#        try:
#            ret=call([cmd], stdout=STDOUT, stderr=STDOUT)
#        except:
#            print('Still unable')
#        return -1
#    return ret

# set up the return functions from params. can be done by param
#  or with just returning all params
    
if __name__ == "__main__":
    main()
