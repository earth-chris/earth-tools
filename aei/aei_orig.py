#!/bin/env/python2.7
##
## The wonderous objects and functions for [AEI]
##
###################
##
from subprocess import call, STDOUT
import gdal as gdal
import osr as osr
import numpy as np
import spectral as spectral

###################
## set up system params and defaults
###################

# set up a file check function
def checkFile(infile, quiet=False):
    if os.path.isfile(infile) and os.access(infile, os.R_OK):
        return True
    else:
        if not quiet:
            print("[ ERROR ]: Unable to read file: %s" % infile)
        return False

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
            print('[ ERROR ]: Invalid -ot argument: %s, unable to output using String format' % (arg))
        elif (str(arg).lower() == '8'):
            print('[ ERROR ]: Invalid -ot argument: %s, unable to output using Structure format' % (arg))
        elif (str(arg).lower() == '9'):
            params.ot = 'CFloat64'
        elif (str(arg).lower() == '10'):
            print('[ ERROR ]: Invalid -ot argument: %s, unable to output using Pointer format' % (arg))
        elif (str(arg).lower() == '11'):
            print('[ ERROR ]: Invalid -ot argument: %s unable to output using Object format' % (arg))
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
        print('[ ERROR ]: Unable to parse -ot argument: %s' % (arg))
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
            print('[ ERROR ]: Unable to call command %s' % (self))
            return -1
        outfile.flush()
        if ret > 0:
            raise Exception("[ ERROR ]: Failed at command (return code %d): %s" % (ret, str(self)))
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
        
    def plot(self, inds = [], legend = False):
        """
        plots the spectra using a standard plot format
        can be set to only plot selected spectra
        
        usage: self.plot(inds = [], legend = False)
          where inds = optional 0-based indices for spectra to plot
                legend = set this to force a legend to be created
        """
        # import pyplot
        import matplotlib.pyplot as plt
        
        # set basic parameters
        plt.xlim((self.band_centers.min(), self.band_centers.max()))
        plt.xlabel('Wavelength (' + self.band_unit + ')')
        plt.ylabel('Reflectance (%)')
        plt.title('spectralObject plot')
        
        # check if indices were set and valid. if not, plot all items
        if inds:
            if max(inds) > len(self.names):
                inds = range(0, len(self.names))
                print("[ ERROR ]: invalid range set. using all spectra")
            if min(inds) < 0:
                inds = range(0, len(self.names))
                print("[ ERROR ]: invalid range set. using all spectra")
        else:
            inds = range(0, len(self.names))
            
        # turn on the legend if fewer than 10 entries
        if len(inds) < 10:
            legend = True
            
        # loop through each item to plot
        for i in inds:
            plt.plot(self.band_centers, self.spectra[i,:], 
                label = self.names[i])
            
        # add the legend with each spectrum's name
        if legend:
            plt.legend(fontsize = 'small', framealpha = 0.5, 
                fancybox = True)
        
        # display the plot
        plt.show()
        
    def bn(self, inds = []):
        """
        brightness normalizes the spectra
        
        usage: self.bn(inds = [])
          where inds = the indices to use for BN
        """
        # check if indices were set and valid. if not, plot all items
        if inds:
            if max(inds) > self.spectra.shape[-1]:
                inds = range(0, self.spectra.shape[-1])
                print("[ ERROR ]: invalid range set. using all spectra")
            if min(inds) < 0:
                inds = range(0, self.spectra.shape[-1])
                print("[ ERROR ]: invalid range set. using all spectra")
        else:
            inds = range(0, self.spectra.shape[-1])
            
        # perform the bn
        self.spectra = self.spectra[:,inds] / np.expand_dims(
            np.sqrt((self.spectra[:,inds]**2).sum(1)),1)
        
        # subset band centers to the indices selected, if they exist
        if self.band_centers.ndim != 0:
            self.band_centers = self.band_centers[inds]
        
# set series of landsat tools
class landsat:
    def __init__(self):
        pass
    
    @staticmethod
    def readMTL(mtl_file):
        """
        parses an MTL file and returns a dictionary with the available data
        
        usage: aei.landsat.readMTL(mtl_file)
        """
        # check that file exists
        if not checkFile(mtl_file):
            return -1
            
        # set up dictionary
        mtl_dict = {}
        with open(mtl_file, 'r') as f:
            for line in f:
                if "GROUP =" in line:
                    continue
                line = line.strip().split(" = ")
                if line[0] == "END":
                    break
                
                # convert to number if not string
                if not '"' in line[-1]:
                    try:
                        line[-1] = float(line[-1])
                    except ValueError:
                        pass
                
                # update the output dictionary    
                mtl_dict.update({line[0] : line[-1]})
        
        return mtl_dict
            
    # write function to calculate TOA reflectance
    @staticmethod
    def calcTOA(mtl_file, output_file, dn_file = '', of = 'GTiff'):
        """
        calculates top of atmosphere reflectance for a given
        landsat scene. Must specify MTL file.
        
        usage: aei.landsat.calcTOA(mtl_file, output_file, dn_file = '', of='GTiff')
          where dn_file is the optional path to a DN image stack.
          where of is the GDAL format output driver name
          typically gets the file paths from the MTL file and uses
          the unstacked tif data to create a reflectance stack.
        """
        # set creation options based on the output format
        if of == 'GTiff':
            options = ['INTERLEAVE=PIXEL', 'PHOTOMETRIC=RGB']
        else:
            options = []
        
        # check if the dn file is set
        dn_flag = False
        if dn_file:
            dn_flag = True
            dn = gdal.Open(dn_file)
            
            # check it is a 6-band raster
            if dn.RasterCount != 6:
                print("[ ERROR ]: Unable to calculate TOA reflectance")
                print("[ ERROR ]: Input reflectance file: %s", dn_file)
                print("[ ERROR ]: Is not a 6-band raster file")
                dn = None
                return -1
        
        # read the mtl file for parameters
        input_path = os.path.dirname(mtl_file)
        if input_path == '':
            input_path = '.'
        input_path = input_path + pathsep
        
        mtl = landsat.readMTL(mtl_file)
        if mtl == -1:
            print("[ ERROR ]: Unable to calculate TOA reflectance")
            return -1
        
        # set reflectance bands to use based on landsat sensor
        if mtl['SPACECRAFT_ID'] == '"LANDSAT_8"':
            bands = [2, 3, 4, 5, 6, 7]
        elif mtl['SPACECRAFT_ID'] == '"LANDSAT_5"':
            bands = [1, 2, 3, 4, 5, 7]
        elif mtl['SPACECRAFT_ID'] == '"LANDSAT_7"':
            bands = [1, 2, 3, 4, 5, 7]
        else:
            bands = [1, 2, 3, 4]
            
        # get the output file size based on the dn file or MTL file
        if dn_flag:
            dims = [dn.RasterXSize, dn.RasterYSize, dn.RasterCount]
        else:
            dims = [mtl['REFLECTIVE_SAMPLES'], mtl['REFLECTIVE_LINES'], 6.]
        
        # set up the output file
        print("[ STATUS ]: Creating output file %s" % output_file)
        outRaster = gdal.GetDriverByName(of).Create(
            output_file, int(dims[0]), int(dims[1]), int(dims[2]), 
            gdal.GDT_Float32, options = options)
        
        # read the input file(s) and apply gain/offset
        i = 0
        for j in bands:
            gain = mtl['REFLECTANCE_MULT_BAND_' + str(j)]
            offset = mtl['REFLECTANCE_ADD_BAND_' + str(j)]
            
            # read the data to an array
            if dn_flag:
                band = dn.GetRasterBand(j).ReadAsArray()
            else:
                dn = gdal.Open(input_path + mtl['FILE_NAME_BAND_' + str(j)].strip('"'))
                band = dn.GetRasterBand(1).ReadAsArray()
            
            # report status
            print("[ STATUS ]: Processing Band %d" % (i + 1))
            
            # apply the gain and offset    
            array = np.multiply(band,gain) + offset
            
            # write the output array to file and set metadata info
            outBand = outRaster.GetRasterBand(i+1)
            outBand.WriteArray(array)
            outBand.FlushCache()
            outBand.SetNoDataValue(offset)
            
            # set color interpretation so bands 5-4-3 are default rgb
            color_interp = [gdal.GCI_Undefined, gdal.GCI_Undefined,
                gdal.GCI_BlueBand, gdal.GCI_GreenBand, 
                gdal.GCI_RedBand, gdal.GCI_Undefined]
            outband.SetColorInterpetation(color_interp[i])
            
            # increment counter and move on to next band
            i += 1
            
        # get metadata info for the output file
        proj = dn.GetProjection()
        geot = dn.GetGeoTransform()
        
        # set output file info
        outRaster.SetProjection(proj)
        outRaster.SetGeoTransform(geot)
        outRaster.BuildOverviews(resampling = 'NEAREST',
            overviewlist = [2, 4, 8, 16, 32, 64, 128])
            
        # close up shop
        outBand = None
        outRaster = None
        
# set command merger
def catCmd(executable, args):
    parts = command(executable)
    parts.extend(args)
    return parts
    
# set function to join strings
def strJoin(list, join=' '):
    outstr = join.join(map(str,list))
    return outstr
    
# set function to return an array of random floats based on a min/max
def randomFloats(n_iterations, minval, maxval):
    return np.random.ranf(n_iterations) * (maxval - minval) + minval

# set function to join miscellaneous arguments
def parseEtc(etc): 
    if type(etc) is list:
        etc=strJoin(etc)
    if type(etc) is not str:
        print('[ ERROR ]: Unable to parse the etc argument.')
        type(etc)
        etc=''
    return etc
    
# set function to brightness normalize an array
def bn(array, axis = 1, inds = []):
    """
    performs a brightness normalization on a numpy array
    
    usage: output = aei.bn(array, axis = 1, inds = [])
      where: array = the input array to normalize
             axis  = the axis on which to normalize. default is the last dimension
             inds  = the indices 
    """
    # set the proper axis to use
    if axis > (array.ndim - 1):
        print("[ ERROR ]: invalid axis set: %s" % axis)
        print("[ ERROR ]: using dimenssion: %s" % (array.ndim-1))
        axis = array.ndim - 1
        
    # check that indices are set properly
    if inds:
        if max(inds) > array.shape[axis]:
            inds = range(0, array.shape[axis])
            print("[ ERROR ]: invalid range set. using all spectra")
        if min(inds) < 0:
            inds = range(0, array.shape[axis])
            print("[ ERROR ]: invalid range set. using all spectra")
    else:
        inds = range(0, array.shape[axis])
    
    # set up bn based on shape of array
    if array.ndim == 1:
        return array[inds] / np.sqrt((array[inds] ** 2).sum())
        
    elif array.ndim == 2:
        if axis == 0:
            return array[inds,:] / np.expand_dims(
                np.sqrt((array[inds,:] ** 2).sum(axis)), axis)
        elif axis == 1:
            return array[:,inds] / np.expand_dims(
                np.sqrt((array[:,inds] ** 2).sum(axis)), axis)
        
    elif array.ndim == 3:
        if axis == 0:
            return array[inds,:,:] / np.expand_dims(
                    np.sqrt((array[inds,:,:] ** 2).sum(axis)),axis)
        elif axis == 1:
            return array[:,inds,:] / np.expand_dims(
                np.sqrt((array[:,inds,:] ** 2).sum(axis)),axis)
        elif axis == 2:
            return array[:,:,inds] / np.expand_dims(
                np.sqrt((array[:,:,inds] ** 2).sum(axis)),axis)
    
    else:
        print("[ ERROR ]: unable to brightness normalize")
        print("[ ERROR ]: array must be 3 dimensions or smaller")
        print("[ ERROR ]: array provided has [%s] dimensions" % array.ndim)
        return -1
    
# set function to read the ascii spectra from the joint fire science
#  program (https://www.frames.gov/partner-sites/assessing-burn-severity/spectral/spectral-library-southern-california/)
def readJFSC(infile):
    """
    reads the ascii format spectral data from the joint-fire-science-program
    and returns an object with the mean and +/- standard deviation
    reflectance data
    """
    # check input file
    if not checkFile(infile):
        return -1
    
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
        
def readUSGS(infile):
    """
    reads the ascii format spectral data from USGS
    and returns an object with the mean and +/- standard deviation
    reflectance data
    """
    if not checkFile(infile):
        return -1
    
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
    
def readSpecLib(infile):
    """
    reads an envi format spectral library
    and returns an object with the mean 
    reflectance data
    """
    if not checkFile(infile):
        return -1
        
    # check for header file
    if checkFile(infile[:-4] + '.hdr', quiet=True):
        hdr = infile[:-4] + '.hdr'
    else:
        if checkFile(infile + '.hdr'):
            hdr = infile + '.hdr'
        else:
            return -1
        
    # read the spectral data
    slib = spectral.envi.open(hdr, infile)
    
    # create the spectral object
    s = spectralObject(slib.params.nrows, slib.params.ncols,
            band_centers = np.asarray(slib.bands.centers),
            band_unit = slib.bands.band_unit,
            band_quantity = slib.bands.band_quantity)
    
    # set the spectra and names        
    s.spectra = slib.spectra
    s.names = slib.names
    
    # return the final object
    return s

def tileRaster(infile, outfile, tiling=2, buff=0, of='GTiff'):
    """
    tiles a raster file into a set number of tiles based on
    a tiling parameter.
    
    syntax: tileRaster(infile, tiling=2, buff=0, of='GTiff')
    
    infile   [string] - the input raster file to tile
    outfile  [string] - the output base name to append each file ID to
    tiling [int/list] - the tiling parameter and style. setting as an
                           int should be the sqrt of the number of output
                           tiles (i.e. tiling = 2 means and output of 4 tiles
                           in a 2 x 2 orientation). setting as a list
                           sets the number of pixels to use for each dimension
                           (ie. tiling = [200, 200] will output 200x200 pixel tiles)
    buff        [int] - the number of pixels to buffer on each side
    of       [string] - the output format (in GDAL string format)
    """
    # set defaults for ensuring i/o is smooth
    if not checkFile(infile):
        return -1
    
    if not type(buff) is int:
        print("[ ERROR ]: Buffer specified (%s) is not an integer" % buff)
        return -1
        
    if of == "GTiff":
        etc = "-co COMPRESS=LZW"
        append = ".tif"
    else:
        append = ""
        etc = ""
        
    # get metadata from input file
    inf = gdal.Open(infile)
    srs = osr.SpatialReference()
    
    # get projection info
    wkt = inf.GetProjection()
    srs.ImportFromWkt(wkt)
    prj = srs.ExportToProj4()
    
    # get info on x/y sizes
    nx = inf.RasterXSize
    ny = inf.RasterYSize
    ulx, sx, ox, uly, oy, sy = inf.GetGeoTransform()
    
    # figure out tiling x/y min/max scheme based on int/list setup
    if type(tiling) is int:
        xmin = np.arange(tiling)
        xmax = np.arange(tiling)+1
        ymin = np.arange(tiling)+1
        ymax = np.arange(tiling)
        
        for i in range(tiling-1):
            xmin = np.row_stack((xmin, xmin))
            xmax = np.row_stack((xmax, xmax))
            ymin = np.row_stack((ymin, ymin))
            ymax = np.row_stack((ymax, ymax))
            
        ymin = ymin.transpose()
        ymax = ymax.transpose()
        
        nxt = tiling
        nyt = tiling
        
        nxp = int(math.floor(float(nx) / nxt))
        nyp = int(math.floor(float(ny) / nyt))
        
    elif type(tiling) is list:
        if len(tiling) != 2:
            print("[ ERROR ]: tiling parameter %s is not a 2-element list" % tiling)
            return -1
        
        nxp = tiling[0]
        nyp = tiling[1]
            
        nxt = int(math.ceil(float(nx) / nxp))
        nyt = int(math.ceil(float(ny) / nyp))
        
        xmin = np.arange(nxt)
        xmax = np.arange(nxt)+1
        ymin = np.arange(nyt)+1
        ymax = np.arange(nyt)
        
        for i in range(nxt-1):
            xmin = np.row_stack((xmin, xmin))
            xmax = np.row_stack((xmax, xmax))
            
        for i in range(nyt-1):
            ymin = np.row_stack((ymin, ymin))
            ymax = np.row_stack((ymax, ymax))
            
        ymin = ymin.transpose()
        ymax = ymax.transpose()
        
    else:
        print("[ ERROR ]: tiling parameter %s is not a 2-element list or int" % tiling)
        return -1
    
    # apply the geo transform info to get x and y coords
    xmin = xmin * nxp * sx - buff + ulx
    xmax = xmax * nxp * sx + buff + ulx
    ymin = ymin * nyp * sy - buff + uly
    ymax = ymax * nyp * sy + buff + uly
    
    print(xmax)
    print(ymin)
    
    # crop to the extent at the edges
    if xmax[0,nxt-1] > (ulx + (nx * sx) + buff):
        xmax[:,nxt-1] = ulx + (nx * sx) + buff
        
    if ymin[nyt-1,0] < (uly + (ny * sy) - buff):
        ymin[nyt-1,:] = uly + (ny * sy) - buff
        
    # report status
    print("[ STATUS ]: Tiling %s into %s output tiles") % (infile, nxt*nyt)
    
    # set up loop to run gdalwarp on tiles
    cnt = 0
    for i in range(nxt):
        for j in range(nyt):
            cnt += 1
            outtmp = (outfile + "_%03.f" + append) % cnt
            gdalwarp(infile, outtmp, etc, multi=True, of=of,
                te=[xmin[i,j], ymin[j,i], xmax[i,j], ymax[j,i]],
                overwrite=True)

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
        if not isinstance(n_threads, (int, long, float, complex)):
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
        ot=getOt(ot)
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

    fparams.append(parseEtc(etc))
    fparams=strJoin(fparams)
    
    # join and run the command
    ocmd = strJoin([cmd.gdalwarp[0], fparams, inputs, output])
    print('[ COMMAND ]: %s' % (ocmd))
    os.system(ocmd)
    
    return(ocmd)
    
def main():
    
    # import modules to namespace
    from aei import cmd
    from aei import landsat
    from aei import objects
    from aei import params
    from aei import read

# run main on load
if __name__ == "__main__":
    main()
