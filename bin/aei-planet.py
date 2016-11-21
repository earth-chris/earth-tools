#!/usr/bin/python
##
## reads and applies standard processing parameters to planet rgb imagery
##
####################
##
import aei as aei
import numpy as np
import gdal as gdal
import osr as osr
import shutil
import json
import glob
import cv2
import sys

# define global parameters
class default_params:
    def __init__(self):
        
        # set base file parameters
        self.imageFile = None
        self.metadataFile = None
        self.geoinfoFile = None
        self.inputPath = None
        self.basename = None
        self.outputFile = None
        self.tempFile = None
        self.append = "_derived.tif"
        self.noData = -9999
        
        # set the number of output bands
        self.nOutputBands = 18
        
        # set some metadata parameters
        self.snr = None
        self.exposureTime = None
        self.altitude = None
        self.offNadir = None
        self.solarZenith = None
        self.solarAzimuth = None
        self.localTime = None
        self.cloudCover = None

# parse command line arguments
class parse_args:
    def __init__(self, arglist, params):
        
        # we'll be parsing the command line arguments passed
        #  via arglist, updating the default parameters,
        #  then returning the params list 
        
        # exit if no arguments passed
        if len(arglist) == 1:
            usage(exit=True)
    
        # read arguments from command line
        i = 1
        
        while i < len(arglist):
            arg = arglist[i]
        
            # check input flag    
            if arg.lower() == '-i':
                i += 1
                arg = arglist[i]
                
                # if input is a path, find the input files
                if aei.params.os.path.isdir(arg):
                    params.inputPath = arg
                    search = glob.glob(arg + aei.params.pathsep + "*")
                    
                    # reject if no files found
                    if len(search) == 0:
                        usage()
                        print("[ ERROR ]: no files found in input directory")
                        sys.exit(1)
                        
                    for file in search:
                        if ".tif" in file:
                            params.imageFile = file
                        if "_metadata.json" in file:
                            params.metadataFile = file
                        if "_geoinfo.json" in file:
                            params.geoinfoFile = file
                
                # if input is a .tif file            
                else:
                    dirname = aei.params.os.path.dirname(arg)
                    params.inputPath = dirname
                    params.imageFile = arg
                    
                    # find json files
                    loc = arg.find("_analytic.tif")
                    
                    if loc == -1:
                        usage()
                        print("[ ERROR ]: input tif file does not end with _analytic.tif")
                        print("[ ERROR ]: input tif file: %s" % arg)
                        sys.exit(1)
                        
                    else:
                        params.metadataFile = arg[:loc] + "_metadata.json"
                        params.geoinfoFile = arg[:loc] + "_geoinfo.json"
        
            # check if output is set                
            elif arg.lower() == '-o':
                i += 1
                arg = arglist[i]
                
                try:
                    params.outputFile = arg
                except:
                    print("[ ERROR ]: Unable to set output file: %s" % arg)
                    sys.exit(1)
                    
            else:
                print("[ ERROR ]: Unrecognized argument: %s" % arg)
                usage()
                sys.exit(1)
                
            i += 1

# read the metadata json file        
def readMetadata(params):
    
    # read the json metadata info into memory
    with open(params.metadataFile, 'r') as jsonf:
        jdict = json.load(jsonf)
    j = jdict[u'properties']
    
    # extract useful metadata info
    params.exposureTime = j[u'camera'][u'exposure_time']
    params.cloudCover = j[u'cloud_cover'][u'estimated']
    params.snr = j[u'image_statistics'][u'snr']
    params.altitude = j[u'sat'][u'alt']
    params.offNadir = j[u'sat'][u'off_nadir']
    params.solarZenith = 90-j[u'sun'][u'altitude']
    params.solarAzimuth = j[u'sun'][u'azimuth']
    params.localTime = j[u'sun'][u'local_time_of_day']

# make sure files exist    
def checkInputs(params):
    if not aei.fn.checkFile(params.imageFile, True):
        usage()
        aei.fn.checkFile(params.imageFile, False)
        sys.exit(1)
        
    if not aei.fn.checkFile(params.metadataFile, True):
        usage()
        aei.fn.checkFile(params.metadataFile, False)
        sys.exit(1)
           
    if not aei.fn.checkFile(params.geoinfoFile, True):
        usage()
        aei.fn.checkFile(params.geoinfoFile, False)
        sys.exit(1)
        
    # set up default output file based on the input basename
    if not params.outputFile:
        loc = params.imageFile.find("_analytic.tif")
        pathlen = len(params.inputPath)
        params.basename = params.imageFile[pathlen:loc]
        params.outputFile = params.inputPath + params.basename + params.append
        params.tempFile = params.inputPath + params.basename + "_temp" + params.append
        
# set up function to calculate simple ratios
def simpleRatio(band1, band2):
    return (band1.astype(np.float32)) / (band2.astype(np.float32))

# set up function to calculate normalized ratios    
def normalizedRatio(band1, band2):
    band1 = band1.astype(np.float32)
    band2 = band2.astype(np.float32)
    return (band1 - band2) / (band1 + band2)

# set up function to calculate distance from the edge of good data
def distanceFromEdge(array, gdalRef):
    """
    calculates the distance from the edge of an image.
     
    expects the input array to be the alpha band (a uint16 2-dimension array)
    """
    # find the good data values in the alpha band
    gd = np.where(array == array.max())
    
    # and find the converse bad data to mask later
    bd = np.where(array != array.max())
    
    # create a byte array and set good data to 255
    byte = np.zeros(array.shape, dtype = np.uint8)
    byte[gd[0], gd[1]] = 255
    
    # create a kernel to close up small chunks of no-data
    kernel = np.ones([50, 50], np.uint8)
    closing = cv2.morphologyEx(byte, cv2.MORPH_CLOSE, kernel)
    
    # detect edges on cleaned up data
    edges = cv2.Canny(closing, 0, 255)
    
    # set up a temp output file
    newFile = gdal.GetDriverByName("GTiff").Create("temp_edge.tif", array.shape[1], 
        array.shape[0], 1, gdal.GDT_Byte, options=["COMPRESS=LZW"])
    
    # set up projection info
    geo = gdalRef.GetGeoTransform()
    prj = gdalRef.GetProjection()
    newFile.SetGeoTransform(geo)
    newFile.SetProjection(prj)
    
    # write the data
    band = newFile.GetRasterBand(1)
    band.WriteArray(edges)
    band.FlushCache()
    
    # clear from memory
    band = None
    newFile = None
        
    # run gdal_proximity to get distance from edge
    aei.cmd.gdal_proximity("temp_edge.tif", "temp_proximity.tif", 
        of="GTiff", distunits="GEO")
        
    # read the output into memory
    dfeRef = gdal.Open("temp_proximity.tif")
    dfe = dfeRef.ReadAsArray()
    
    # clear the reference from memory
    dfeRef = None
    
    # delete the files 
    aei.params.os.remove("temp_edge.tif")
    aei.params.os.remove("temp_proximity.tif")
    
    # mask the bad data
    dfe[bd[0],bd[1]] = params.noData
    
    # return the final array
    return dfe

# set up a function to write an array to a raster band
def writeArrayToRaster(array, gdalRef, band, noData = None):
    """
    writes a 2-d array to an output raster. requires a gdal open file reference.
    """
    # set the gdal band reference
    bandRef = gdalRef.GetRasterBand(band)
    
    # set the No Data value if set
    if noData:
        bandRef.SetNoDataValue(noData)
    
    # write the array to memory
    bandRef.WriteArray(array)
    
    # clear the file cache
    #bandRef.FlushCache()
    
def usage(exit=False):
    """
    describes the aeiplanet.py procedure in case of incorrect parameter calls
    
    syntax: usage(exit=False)
    """

    print(
        """
$ aei-planet.py -i input_files
    
    output band list is:
    1. Red  2. Green  3. Blue  4. Brightness normalized (BN) Red
    5. BN Green 6. BN Blue 7. Brightness scalar 8. Distance from edge
    9. Simple Ratio (SR) red-green, 10. SR green-blue 11. SR red-blue
    12. Normalized difference (ND) red-green 13. ND green-blue 14. ND red-blue
    Exposure Time, Off-Nadir, SNR, Solar Zenith, Solar Azimuth, Sensor Altitude
        """)
    
    if exit:
        sys.exit(1)

def main ():
    """
    the main program for aei-planet.py
    
    this routine performs the following
    1) finds the indices with actual data (that passes planet's qa)
    2) performs a brightness normalization
    3) finds the distance from image edge
    4) calculates the simple ratios for all rgb bands
    5) calculates normalized ratios for all rgb bands
    6) assigns metadata (e.g. solar azimuth and zenith) to raster bands
    
    it writes data to file after each set of calculations.

    syntax: main()
    """
    import gdal
    import numpy as np
    import cv2

    # first, read the default parameters to get the processing object
    params = default_params()
    
    # then parse the arguments passed via command line
    args = aei.params.sys.argv
    parse_args(args, params)
    
    # check that inputs make sense
    checkInputs(params)
    
    # parse the json data
    readMetadata(params)
    
    # report starting
    print('[ STATUS ]: Running aei-planet.py')
    print('[ STATUS ]: Reading input file: %s' % params.imageFile)
    
    # get the input gdal reference
    gdalRef = gdal.Open(params.imageFile)
                                  
    ####################
    # create the output file
    print('[ STATUS ]: Creating output file: %s' % params.outputFile)
    outRef = gdal.GetDriverByName("GTiff").Create(params.tempFile, gdalRef.RasterXSize, 
        gdalRef.RasterYSize, params.nOutputBands, gdal.GDT_Float32, 
        options=["BIGTIFF=YES"])
      
    # add georeferencing info to output file
    outRef.SetGeoTransform(gdalRef.GetGeoTransform())
    outRef.SetProjection(gdalRef.GetProjection())
    
    ####################
    # read the rgb image into memory
    print("[ STATUS ]: Reading RGB data")
    rgba = gdalRef.ReadAsArray()
    
    # get the alpha band to find good data 
    alpha = rgba[3]
    
    # find the good data to limit analysis to
    gd = np.where(alpha != 0)
    
    # subset rgb to good data, and transform to floating point
    rgbf = rgba[0:3, gd[0], gd[1]] * 1.0
    
    # kill the rgba array to save memory
    rgba = None
    
    # set up a counter for the output bands
    bandCounter = 1
    
    # set up a temp array we'll use to write data band-by-band
    tmpArray = np.zeros((gdalRef.RasterYSize, gdalRef.RasterXSize), np.float32) + params.noData
    
    # write the raw rgb data to the output file
    for i in range(0,3):
        tmpArray[gd[0], gd[1]] = rgbf[i]
        writeArrayToRaster(tmpArray, outRef, bandCounter, noData = params.noData)
        bandCounter += 1
    
    ####################
    # perform the brightness normalization
    print('[ STATUS ]: Brightness normalizing')
    bn, bnScalar = aei.functions.bn(rgbf, returnScalar = True, axis = 0)
    
    # write to the output file
    for i in range(0,3):
        tmpArray[gd[0], gd[1]] = bn[i]
        writeArrayToRaster(tmpArray, outRef, bandCounter, noData = params.noData)
        bandCounter +=1
    
    # write just the scalar band
    tmpArray[gd[0], gd[1]] = bnScalar[0]
    writeArrayToRaster(tmpArray, outRef, bandCounter, noData = params.noData)
    bandCounter +=1
    
    # kill bn refs
    bn = None ; bnScalar = None
    
    ####################
    # find the distance from the edge
    print('[ STATUS ]: Finding distance from image edge')
    dfe = distanceFromEdge(alpha, gdalRef)
    
    # write it
    writeArrayToRaster(dfe, outRef, bandCounter, noData = params.noData)
    bandCounter += 1
    dfe = None
    
    ####################
    # run the simple ratios
    print('[ STATUS ]: Caculating simple ratios')
    
    # red/green
    tmpArray[gd[0], gd[1]] = simpleRatio(rgbf[0], rgbf[1])
    writeArrayToRaster(tmpArray, outRef, bandCounter, noData = params.noData)
    bandCounter +=1
    
    # green/blue
    tmpArray[gd[0], gd[1]] = simpleRatio(rgbf[1], rgbf[2])
    writeArrayToRaster(tmpArray, outRef, bandCounter, noData = params.noData)
    bandCounter += 1
    
    # red/blue
    tmpArray[gd[0], gd[1]] = simpleRatio(rgbf[0], rgbf[2])
    writeArrayToRaster(tmpArray, outRef, bandCounter, noData = params.noData)
    bandCounter += 1
    
    ####################
    # run the normalized ratios
    print('[ STATUS ]: Calculating normalized ratios')
    
    # red/green
    tmpArray[gd[0], gd[1]] = normalizedRatio(rgbf[0], rgbf[1])
    writeArrayToRaster(tmpArray, outRef, bandCounter, noData = params.noData)
    bandCounter += 1
    
    # green/blue
    tmpArray[gd[0], gd[1]] = normalizedRatio(rgbf[1], rgbf[2])
    writeArrayToRaster(tmpArray, outRef, bandCounter, noData = params.noData)
    bandCounter += 1
    
    # red/blue
    tmpArray[gd[0], gd[1]] = normalizedRatio(rgbf[0], rgbf[2])
    writeArrayToRaster(tmpArray, outRef, bandCounter, noData = params.noData)
    bandCounter += 1
    
    ####################
    # add bands from the metadata
    print('[ STATUS ]: Adding metadata bands')
    
    # exposure time
    tmpArray[gd[0], gd[1]] = params.exposureTime
    writeArrayToRaster(tmpArray, outRef, bandCounter, noData = params.noData)
    bandCounter += 1
    
    #md_offNadir = np.zeros([alpha.shape[0], alpha.shape[1]], dtype=np.float32) + params.offNadir
    
    # solar zenith angle
    tmpArray[gd[0], gd[1]] = params.solarZenith
    writeArrayToRaster(tmpArray, outRef, bandCounter, noData = params.noData)
    bandCounter += 1
    
    # solar azimuth
    tmpArray[gd[0], gd[1]] = params.solarAzimuth
    writeArrayToRaster(tmpArray, outRef, bandCounter, noData = params.noData)
    bandCounter += 1
    
    #md_altitude = np.zeros([alpha.shape[0], alpha.shape[1]], dtype=np.float32) + params.altitude
    
    # signal-to-noise ratio
    tmpArray[gd[0], gd[1]] = params.snr
    writeArrayToRaster(tmpArray, outRef, bandCounter, noData = params.noData)
    bandCounter += 1
    
    # add overviews
    #print('[ STATUS ]: Building overviews')
    #outRef.BuildOverviews("NEAREST", [2, 4, 8, 16])
        
    # kill the file references and placeholders
    outRef = None
    gdalRef = None
    tmpArray = None
    
    # compress the output file and delete the temp file
    print("[ STATUS ]: Finished performing calculations")
    print("[ STATUS ]: Compressing output file")
    aei.cmd.gdal_translate(params.tempFile, params.outputFile, etc=['-co','COMPRESS=LZW'])
    aei.params.os.remove(params.tempFile)
    
    # final report
    print('[ STATUS ]: aei-planet.py complete!')
    print('[ STATUS ]: Please see output file: %s' % params.outputFile)
    
    # high fives
    
if __name__ == "__main__":
    main()