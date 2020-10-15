####
# a series of basic tools for processing landsat data
#
# c. 2016-2017 Christopher Anderson
#####

# function to read a landsat-format MTL file and return a dictionary
def readMTL(mtlFile):
    """
    parses a landsat MTL file and returns a dictionary with the available data.

    usage: et.landsat.readMTL(mtl_file)
    """
    import earthtools as et

    # check that file exists
    if not et.fn.checkFile(mtlFile):
        return -1

    # set up dictionary
    mtlDict = {}
    with open(mtlFile, 'r') as f:
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
            mtlDict.update({line[0] : line[-1]})

    return mtlDict

# write function to calculate TOA reflectance
def calcTOA(mtlFile, outFile, dnFile = '', of = 'GTiff'):
    """
    calculates top of atmosphere reflectance for a given
    landsat scene. Must specify MTL file.

    usage: et.landsat.calcTOA(mtl_file, output_file, dn_file = '', of='GTiff')
      where dn_file is the optional path to a DN image stack.
      where of is the GDAL format output driver name
      typically gets the file paths from the MTL file and uses
      the unstacked tif data to create a reflectance stack.
    """
    import os
    import gdal
    import params

    # set creation options based on the output format
    if of == 'GTiff':
        options = ['INTERLEAVE=PIXEL', 'PHOTOMETRIC=RGB']
    else:
        options = []

    # check if the dn file is set
    dnFlag = False
    if dnFile:
        dnFlag = True
        dn = gdal.Open(dnFile)

        # check it is a 6-band raster
        if dn.RasterCount != 6:
            print("[ ERROR ]: Unable to calculate TOA reflectance")
            print("[ ERROR ]: Input reflectance file: %s", dnFile)
            print("[ ERROR ]: Is not a 6-band raster file")
            dn = None
            return -1

    # read the mtl file for parameters
    inputPath = os.path.dirname(mtlFile)
    if inputPath == '':
        inputPath = '.'
    inputPath = inputPath + params.pathsep

    mtl = readMTL(mtlFile)
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
    if dnFlag:
        dims = [dn.RasterXSize, dn.RasterYSize, dn.RasterCount]
    else:
        dims = [mtl['REFLECTIVE_SAMPLES'], mtl['REFLECTIVE_LINES'], 6.]

    # set up the output file
    print("[ STATUS ]: Creating output file %s" % output_file)
    outRaster = gdal.GetDriverByName(of).Create(
        outputFile, int(dims[0]), int(dims[1]), int(dims[2]),
        gdal.GDT_Float32, options = options)

    # set color interpretation so bands 5-4-3 are default rgb
    colorInterp = [gdal.GCI_Undefined, gdal.GCI_Undefined,
        gdal.GCI_BlueBand, gdal.GCI_GreenBand,
        gdal.GCI_RedBand, gdal.GCI_Undefined]

    # read the input file(s) and apply gain/offset
    for j in bands:
        gain = mtl['REFLECTANCE_MULT_BAND_' + str(j)]
        offset = mtl['REFLECTANCE_ADD_BAND_' + str(j)]

        # read the data to an array
        if dn_flag:
            band = dn.GetRasterBand(j).ReadAsArray()
        else:
            dn = gdal.Open(inputPath + mtl['FILE_NAME_BAND_' + str(j)].strip('"'))
            band = dn.GetRasterBand(1).ReadAsArray()

        # report status
        print("[ STATUS ]: Processing Band %d" % (i + 1))

        # apply the gain and offset
        array = np.multiply(band, gain) + offset

        # write the output array to file and set metadata info
        outBand = outRaster.GetRasterBand(i+1)
        outBand.WriteArray(array)
        outBand.FlushCache()
        outBand.SetNoDataValue(offset)

        # set color interpretation
        outband.SetColorInterpetation(colorInterp[j])

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
