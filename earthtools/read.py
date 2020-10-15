#####
# the functions that read or parse frequently used data formats
#
# c. 2016-2017 Christopher Anderson
#####

# set function to return a standardized data type used by gdal, based on idl/python inputs
def ot(arg):
    """
    reads the argument passed by -ot and returns the output data type as a string.

    accepts both python style data type codes (e.g. Float32, Int16) and IDL
    style data type codes (e.g. 4 for Float32, 1 for Byte)

    syntax: read.ot(arg)
    """
    # Python style data type flags
    try:
        if str(arg).lower() == "byte":
            return "Byte"
        elif str(arg).lower() == "int16":
            return "Int16"
        elif str(arg).lower() == "int32":
            return "Int32"
        elif str(arg).lower() == "uint16":
            return "UInt16"
        elif str(arg).lower() == "uint32":
            return "UInt32"
        elif str(arg).lower() == "float32":
            return "Float32"
        elif str(arg).lower() == "float64":
            return "Float64"
        elif str(arg).lower() == "cint16":
            return "CInt16"
        elif str(arg).lower() == "cint32":
            return "CInt32"
        elif str(arg).lower() == "cfloat32":
            return "CFloat32"
        elif str(arg).lower() == "cfloat64":
            return "CFloat64"

        # IDL style data type flags
        elif str(arg).lower() == "1":
            return "Byte"
        elif str(arg).lower() == "2":
            return "Int16"
        elif str(arg).lower() == "3":
            return "Int32"
        elif str(arg).lower() == "4":
            return "Float32"
        elif str(arg).lower() == "5":
            return "Float64"
        elif str(arg).lower() == "6":
            return "CFloat32"
        elif str(arg).lower() == "7":
            print("[ ERROR ]: Invalid -ot argument: %s, unable to output using String format" % (arg))
            return -1
        elif str(arg).lower() == "8":
            print("[ ERROR ]: Invalid -ot argument: %s, unable to output using Structure format" % (arg))
            return -1
        elif str(arg).lower() == "9":
            return "CFloat64"
        elif str(arg).lower() == "10":
            print("[ ERROR ]: Invalid -ot argument: %s, unable to output using Pointer format" % (arg))
            return -1
        elif str(arg).lower() == "11":
            print("[ ERROR ]: Invalid -ot argument: %s unable to output using Object format" % (arg))
            return -1
        elif str(arg).lower() == "12":
            return "UInt16"
        elif str(arg).lower() == "13":
            return "UInt32"
        elif str(arg).lower() == "14":
            return "Int64"
        elif str(arg).lower() == "15":
            return "UInt64"

    # error message for unknown argument
    except ValueError:
        print("[ ERROR ]: Unable to parse -ot argument: %s" % (arg))
        return -1


# set function to join miscellaneous arguments
def etc(etcList):
    import earthtools as et

    etcStr = et.fn.strJoin(etcList)
    if type(etcStr) is not str:
        print("[ ERROR ]: Unable to parse the etc argument. Must be of type: list")
        print("[ ERROR ]: " + str(type(etcList)))
        etcStr = ""
    return etcStr


def ascii(infile):
    import earthtools as et

    # check input file
    if not et.fn.checkFile(infile):
        return -1

    # create the output list
    contents = []

    # open the file
    with open(infile, "r") as f:

        # read line by line
        for line in f.readlines():
            # and remove the \n or whitespace
            contents.append(line.strip())

    # return the list
    return contents


# set function to read the ascii spectra from the joint fire science
#  program (https://www.frames.gov/partner-sites/assessing-burn-severity/spectral/spectral-library-southern-california/)
def jfsc(infile):
    """
    reads the ascii format spectral data from the joint-fire-science-program
    and returns an object with the mean and +/- standard deviation
    reflectance data

    usage: et.read.jsfc(jsfcFile)
    """
    import numpy as np

    import earthtools as et

    # check input file
    if not et.fn.checkFile(infile):
        return -1

    # create the spectral object
    s = et.objects.spectral(1, type="asd")
    s.spectra_stdevm = np.zeros(s.spectra.shape)
    s.spectra_stdevp = np.zeros(s.spectra.shape)

    # open the file and read the data
    f = open(infile, "r")
    header = f.readline()
    i = 0
    for line in f:
        line = line.strip().split()
        s.spectra[0, i] = line[1]
        s.spectra_stdevp[0, i] = line[2]
        s.spectra_stdevm[0, i] = line[3]
        i += 1

    # close the file
    f.close()

    # return the spectral object
    return s


def usgs(infile):
    """
    reads the ascii format spectral data from USGS
    and returns an object with the mean and +/- standard deviation
    reflectance data
    """
    import numpy as np

    import earthtools as et

    if not et.fn.checkFile(infile):
        return -1

    # open the file and read header info
    f = open(infile, "r")
    x_start = "gibberish"
    for line in f:
        if x_start in line:
            break
        if "Name:" in line:
            spectrum_name = line.strip().split("Name:")[-1].strip()
        if "X Units:" in line:
            band_unit = line.strip().split()
            band_unit = band_unit[-1].strip("()").capitalize()
        if "Y Units:" in line:
            refl_unit = line.strip().split()
            refl_unit = refl_unit[-1].strip("()").capitalize()
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
    if band_unit == "Micrometers":
        band_centers *= 1000.0
        band_unit = "Nanometers"
    if refl_unit == "Percent":
        reflectance /= 100.0

    # create the spectral object
    s = et.objects.spectral(1, n_values, band_centers=band_centers, band_unit=band_unit, band_quantity="Wavelength")

    # assign relevant values
    s.spectra[0] = reflectance
    if spectrum_name:
        s.names[0] = spectrum_name

    # close the file
    f.close()

    # return the final object
    return s


def spectralLib(infile):
    """
    reads an envi format spectral library
    and returns an object with the mean
    reflectance data
    """
    import numpy as np
    import spectral as spectral

    import earthtools as et

    if not et.fn.checkFile(infile):
        return -1

    # check for header file
    if et.fn.checkFile(infile[:-4] + ".hdr", quiet=True):
        hdr = infile[:-4] + ".hdr"
    else:
        if et.fn.checkFile(infile + ".hdr"):
            hdr = infile + ".hdr"
        else:
            return -1

    # read the spectral data
    slib = spectral.envi.open(hdr, infile)

    # create the spectral object
    s = et.objects.spectral(
        slib.params.nrows,
        slib.params.ncols,
        band_centers=np.asarray(slib.bands.centers),
        band_unit=slib.bands.band_unit,
        band_quantity=slib.bands.band_quantity,
    )

    # set the spectra and names
    s.spectra = slib.spectra
    s.names = slib.names

    # return the final object
    return s


# a class to return image data through gdal
class raster:
    # the init class that stores the raster metadata as variables
    def __init__(self, input_file):
        """Reads metadata from a raster file and stores it in a shared object

        Args:
            input_file: a path to a raster file to read

        Returns:
            An object with the raster metadata as object variables (e.g.,
            ras = et.read.raster('file.tif') # file with dims x = 30, y = 50
            ras.nx will be 30, ras.ny will be 50, etc.)
        """
        import gdal

        # read the gdal reference as read-only
        ref = gdal.Open(input_file, 0)
        self.file_name = input_file

        # get file dimensions
        self.nx = ref.RasterXSize
        self.ny = ref.RasterYSize
        self.nb = ref.RasterCount

        # get georeferencing info
        self.prj = ref.GetProjection()
        geo = ref.GetGeoTransform()
        self.xmin = geo[0]
        self.xps = geo[1]
        self.xoff = geo[2]
        self.ymax = geo[3]
        self.yoff = geo[4]
        self.yps = geo[5]
        self.xmax = self.xmin + self.xoff + (self.nx * self.xps)
        self.ymin = self.ymax + self.yoff + (self.ny * self.yps)

        # get no-data info
        band = ref.GetRasterBand(1)
        self.no_data = band.GetNoDataValue()

        # get data type
        self.dt = band.DataType

        # create an empty 'data' variable to read into later
        self.data = None

        # get driver info
        self.driver_name = ref.GetDriver().ShortName

        # kill the gdal references
        ref = None
        band = None

    # a function to read raster data from a single band
    def read_band(self, band):
        """Reads the raster data from a user-specified band into the self.data variable

        Args:
            band: the 1-based index for the band to read

        Returns:
            the et.Raster object with the object.data variable updated with a
            numpy array of raster values
        """
        import gdal

        ref = gdal.Open(self.file_name, 0)
        band = ref.GetRasterBand(band)
        self.data = band.ReadAsArray()

    # a function to read raster data from all bands
    def read_all(self):
        """Reads all bands of raster data

        Args:
            None

        Returns:
            the et.raster object with the object.data variable updated with a
            numpy array of raster values
        """
        import gdal

        ref = gdal.Open(self.file_name, 0)
        self.data = ref.ReadAsArray()

    # a function to write raster data to a single band
    def write_band(self, band, data):
        """Writes new raster data to a user-specified band

        Args:
            band: a 1-based integer with the band to write to
            data: a numpy array with raster data to write

        Returns:
            None.
        """
        import gdal

        ref = gdal.Open(self.file_name, 1)
        band = ref.GetRasterBand(band)
        band.WriteArray(data)

    # a function to write raster data to all bands
    def write_all(self, data=None):
        """Writes new raster data to all bands

        Args:
            data: a numpy array with raster data to write. if not set, writes self.data

        Returns:
            None.
        """
        import gdal

        ref = gdal.Open(self.file_name, 1)

        # see which data to write
        if data:
            ref.WriteRaster(data)
        else:
            ref.WriteRaster(self.data)

    def write_metadata(self, ref=None):
        """Updates the metadata of a file if changed by the user

        Args:
            None

        Returns:
            None
        """
        import gdal

        if ref is None:
            ref = gdal.Open(self.file_name, 1)

        # update projection and geotransform at file-level
        ref.SetProjection(self.prj)
        ref.SetGeoTransform([self.xmin, self.xps, self.xoff, self.ymax, self.yoff, self.yps])

        # update no-data value band by band
        if self.no_data is not None:
            for band in range(1, self.nb + 1):
                b = ref.GetRasterBand(band)
                b.SetNoDataValue(self.no_data)

    def copy(self, file_name, nb=None, driver=None, dt=None, options=None):
        """Creates a new raster file and object with a copy of the raster metadata that can be
        used to write a new derived data product.

        Args:
            file_name: the name of the new file to create as a new reference
            nb       : the number of bands for the new file
            driver   : the name of the driver to use. default is the driver of the input reference
            dt       : the output data type
            options  : gdal raster creation options, passed as a list

        Returns:
            a new et.raster object with updated properties, and a new gdal file
            with those properties written to its header. no raster data is written to
            this file, only metadata.
        """
        import copy as cp

        import gdal

        # create a copy of the input object to manipulate
        new_obj = cp.copy(self)

        # update with the new parameters
        new_obj.file_name = file_name

        if nb:
            new_obj.nb = nb

        if driver:
            new_obj.driver_name = driver

        if dt:
            new_obj.dt = dt

        # create a new raster file with these parameters, but don't write any data
        ref = gdal.GetDriverByName(new_obj.driver_name).Create(
            new_obj.file_name, new_obj.nx, new_obj.ny, new_obj.nb, new_obj.dt, options=options
        )

        # set the projection and geotransform parameters
        new_obj.write_metadata(ref=ref)

        # kill the reference and return the new object
        ref = None

        return new_obj


# functions to convert decimal degrees and degrees-minutes-seconds
#  both functions written by Curtis Price, http://profile.usgs.gov/cprice
#  and pulled from ESRI forums (https://community.esri.com/thread/27279)
def dd_to_dms(dd1, dd2, ndec=6):
    """Convert a decimal degree coordinate pair to a six-tuple of degrees, minutes seconds.

    The returned values are not rounded.

    Arguments

    dd1, dd2 - coordinate pair, in decimal degrees

    Example

      >>> dd2dms(-74.25,32.1)
      (-74, 15, 6.9444444444444444e-05, 32, 6, 2.7777777777778172e-05)
    """

    def ToDMS(dd):
        dd1 = abs(float(dd))
        cdeg = int(dd1)
        minsec = dd1 - cdeg
        cmin = int(minsec * 60)
        csec = (minsec % 60) / float(3600)
        if dd < 0:
            cdeg = cdeg * -1
        return cdeg, cmin, csec

    try:
        # return a six-tuple
        return ToDMS(dd1) + ToDMS(dd2)
    except:
        raise Exception


def dms_to_dd(deg1, min1, sec1, deg2, min2, sec2):
    """Convert a degrees-minutes seconds coordinate pair to decimal degrees.

    The returned values are not rounded.

    Arguments

      deg1,min1,sec1,deg2,min2,sec2 - DMS coordinate pair (six values)

    Example

    >>> dms2deg(-74,45,0,34,10,20)
    (-74.75, 34.172222222222217)
    """

    def ToDD(deg, min=0, sec=0):
        dd = abs(deg) + min / 60.0 + sec / 3600.0
        if deg < 0:
            dd = dd * -1.0
        return dd

    try:
        return ToDD(deg1, min1, sec1), ToDD(deg2, min2, sec2)
    except Exception:
        raise Exception
