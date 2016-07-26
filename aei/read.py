#####
# the functions that read or parse various formats of data
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
        if (str(arg).lower() == 'byte'):
            return 'Byte'
        elif (str(arg).lower() == 'int16'):
            return 'Int16'
        elif (str(arg).lower() == 'int32'):
            return 'Int32'
        elif (str(arg).lower() == 'uint16'):
            return 'UInt16'
        elif (str(arg).lower() == 'uint32'):
            return 'UInt32'
        elif (str(arg).lower() == 'float32'):
            return 'Float32'
        elif (str(arg).lower() == 'float64'):
            return 'Float64'
        elif (str(arg).lower() == 'cint16'):
            return 'CInt16'
        elif (str(arg).lower() == 'cint32'):
            return 'CInt32'
        elif (str(arg).lower() == 'cfloat32'):
            return 'CFloat32'
        elif (str(arg).lower() == 'cfloat64'):
            return 'CFloat64'
        
        # IDL style data type flags
        elif (str(arg).lower() == '1'):
            return 'Byte'
        elif (str(arg).lower() == '2'):
            return 'Int16'
        elif (str(arg).lower() == '3'):
            return 'Int32'
        elif (str(arg).lower() == '4'):
            return 'Float32'
        elif (str(arg).lower() == '5'):
            return 'Float64'
        elif (str(arg).lower() == '6'):
            return 'CFloat32'
        elif (str(arg).lower() == '7'):
            print('[ ERROR ]: Invalid -ot argument: %s, unable to output using String format' % (arg))
            return -1
        elif (str(arg).lower() == '8'):
            print('[ ERROR ]: Invalid -ot argument: %s, unable to output using Structure format' % (arg))
            return -1
        elif (str(arg).lower() == '9'):
            return 'CFloat64'
        elif (str(arg).lower() == '10'):
            print('[ ERROR ]: Invalid -ot argument: %s, unable to output using Pointer format' % (arg))
            return -1
        elif (str(arg).lower() == '11'):
            print('[ ERROR ]: Invalid -ot argument: %s unable to output using Object format' % (arg))
            return -1
        elif (str(arg).lower() == '12'):
            return 'UInt16'
        elif (str(arg).lower() == '13'):
            return 'UInt32'
        elif (str(arg).lower() == '14'):
            return 'Int64'
        elif (str(arg).lower() == '15'):
            return 'UInt64'
        
    # error message for unknown argument
    except ValueError:
        print('[ ERROR ]: Unable to parse -ot argument: %s' % (arg))
        return -1
   
# set function to join miscellaneous arguments
def etc(etcList): 
    if type(etcList) is list:
        etcStr=strJoin(etcList)
    if type(etc) is not str:
        print('[ ERROR ]: Unable to parse the etc argument. Must be of type: list')
        print('[ ERROR ]: Type: ' + type(etcList))
        etcStr=''
    return etcStr

# set function to read the ascii spectra from the joint fire science
#  program (https://www.frames.gov/partner-sites/assessing-burn-severity/spectral/spectral-library-southern-california/)
def jfsc(infile):
    """
    reads the ascii format spectral data from the joint-fire-science-program
    and returns an object with the mean and +/- standard deviation
    reflectance data
    
    usage: aei.read.jsfc(jsfcFile)
    """
    import aei as aei
    import numpy as np
    
    # check input file
    if not aei.checkFile(infile):
        return -1
    
    # create the spectral object
    s = aei.objects.spectral(1,type='asd')
    s.spectra_stdevm = np.zeros(s.spectra.shape)
    s.spectra_stdevp = np.zeros(s.spectra.shape)
    
    # open the file and read the data
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
        
def usgs(infile):
    """
    reads the ascii format spectral data from USGS
    and returns an object with the mean and +/- standard deviation
    reflectance data
    """
    import aei as aei
    import numpy as np
    
    if not aei.checkFile(infile):
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
    s = aei.objects.spectral(1, n_values, band_centers = band_centers,
            band_unit = band_unit, band_quantity = 'Wavelength')
            
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
    import aei as aei
    import numpy as np
    import spectral as spectral
    
    if not aei.checkFile(infile):
        return -1
        
    # check for header file
    if aei.checkFile(infile[:-4] + '.hdr', quiet=True):
        hdr = infile[:-4] + '.hdr'
    else:
        if aei.checkFile(infile + '.hdr'):
            hdr = infile + '.hdr'
        else:
            return -1
        
    # read the spectral data
    slib = spectral.envi.open(hdr, infile)
    
    # create the spectral object
    s = aei.objects.spectral(slib.params.nrows, slib.params.ncols,
            band_centers = np.asarray(slib.bands.centers),
            band_unit = slib.bands.band_unit,
            band_quantity = slib.bands.band_quantity)
    
    # set the spectra and names        
    s.spectra = slib.spectra
    s.names = slib.names
    
    # return the final object
    return s

