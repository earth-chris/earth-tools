#!/usr/bin/python
##
## converts the bounding box of a raster file to a shape file
##
####################

import gdal
import ogr
import osr
import aei as aei

# define global parameters
class default_params:
    def __init__(self):
        
        import datetime
        self.now = datetime.datetime.now()
        
        ### OUTPUT / PROCESSING PARAMETERS
        
        # set the input file parameter
        self.inputFiles = None
        
        # set the output file parameters (if single file)
        # logic for output files is:
        #  if single input file: single output file
        #    if -o set: use -o parameter
        #    if not set: use base name + dir of input file
        #  if multi input file: 
        #    if -o set: single output file
        #    if -odir set: multi output file
        #    if not specified: multi output file in directory of input files
        self.outputFiles = None
        self.outputDir = None
        self.singleOutput = False
        
        # set output driver
        self.driver = ogr.GetDriverByName("ESRI Shapefile")
        
        # and the default file to append 
        self.append = ".shp"
        
        # set this to calculate raster statistics and include in attribute table.
        self.getStats = False

class parse_args:
    def __init__(self, arglist, params):
        
        # we'll be parsing the command line arguments passed
        #  via arglist, updating the default parameters,
        #  then returning the params list 
        
        # exit if no arguments passed
        if len(arglist) == 1:
            usage(exit=True)
            
        # get glob to search for files
        import glob
    
        # read arguments from command line
        i = 1
        
        while i < len(arglist):
            arg = arglist[i]
        
            # check input flag    
            if arg.lower() == '-i':
                i += 1
                arg = arglist[i]
                
                # if glob finds files, use what was set in command line
                if len(glob.glob(arg)) > 0:
                    params.inputFiles = glob.glob(arg)
                
                # if glob fails to find files, try splitting up if there
                #  are spaces between files
                else:
                    split = arg.split()
                    for files in split:
                        if not aei.checkFile(files):
                            aei.params.sys.exit(1)
                            
                    params.inputFiles = arg
                    
            # check output dir flag
            elif arg.lower() == '-odir':
                i += 1
                arg = arglist[i]
                
                params.outputDir = arg
                params.singleOutput = False
                
            # check single output file flag
            elif arg.lower() == '-o':
                i += 1
                arg = arglist[i]
                
                params.outputFiles = arg
                params.singleOutput = True
                
            elif arg.lower() == '-of':
                i += 1
                arg = arglist[i]
                
                # check that the of is legit
                if not arg in aei.params.ogrTypes:
                    print("[ ERROR ]: unrecognized -of argument: %s" % arg)
                    
                else:
                    params.driver = ogr.GetDriverByName(arg)
                
            elif arg.lower() == '-get_stats':
                params.getStats = True
                
            else:
                usage()
                print("[ ERROR ]: Unrecognized argument: %s" % arg)
                aei.params.sys.exit(1)
                
            i += 1

# updates params variable to reflect updates in processing options
def update_params(params):
    """
    this procedure takes the input parameters as set at 
      command line and updates within script
    """
    # go through logic of input and output files and update parameters
    #  check how many input files are specified. If more than one...
    if len(params.inputFiles) > 1:
        
        # if single output set, do nothing
        if params.singleOutput:
            pass
            
        # or multi output
        else:
            
            # check if output dir specified
            if params.outputDir:
                
                # if so, create an output file list of length (n input files)
                params.outputFiles = []
                for file in params.inputFiles:
                    file = aei.params.os.path.basename(file)
                    outFile = aei.fn.replaceExt(file, params.append)
                    
                    params.outputFiles.append(params.outputDir + \
                    aei.params.pathsep + outFile)
                
            # otherwise use directories from each input file
            else:
                params.outputFiles = []
                for file in params.inputFiles:
                    baseName = aei.params.os.path.basename(file)
                    outputDir = aei.params.os.path.dirname(file)
                    outFile = aei.fn.replaceExt(basename, params.append)
                    
                    params.outputFiles.append(outputDir + \
                    aei.params.pathsep + outFile)
                
    # if only one input file    
    else:
        
        # if -o flag set, nothing needs to be updated
        if params.singleOutput:
            pass
        
        # otherwise, set output file to match the input file
        else:
            baseName = aei.params.os.path.basename(params.inputFiles)
            outputDir = aei.params.os.path.dirname(params.inputFiles)
            outFile = aei.fn.replaceExt(baseName, params.append)
            
            params.outputFiles = outputDir + aei.params.pathsep + outFile

def singleOutput(params):
    
    # get spatial reference info
    srs = osr.SpatialReference()
    
    # open up the output file with shapefile driver
    dataSource = params.driver.CreateDataSource(params.outputFiles)
    
    # get projection info from first raster input file
    if type(params.inputFiles) is str:
        inf = gdal.Open(params.inputFiles)
        
        # convert to a list for future looping
        params.inputFiles = [params.inputFiles]
    
    else:
        print(params.inputFiles)
        print(type(params.inputFiles))
        inf = gdal.Open(params.inputFiles[0])
        
    wkt = inf.GetProjection()
    srs.ImportFromWkt(wkt)
    
    # create the output layer
    layer = dataSource.CreateLayer("BoundingBox", srs, ogr.wkbPolygon)
    
    # add field with file name
    fn = ogr.FieldDefn("FileName", ogr.OFTString)
    
    # find longest file name and use that for the length of the field
    length = []
    for file in params.inputFiles:
        length.append(len(file))
    
    # set the field    
    fn.SetWidth(max(length))
    layer.CreateField(fn)
    
    # add fields if -get_stats is set
    if params.getStats:
        rasterMin = ogr.FieldDefn("RasterMin", ogr.OFTReal)
        rasterMax = ogr.FieldDefn("RasterMax", ogr.OFTReal)
        rasterMean = ogr.FieldDefn("RasterMean", ogr.OFTReal)
        rasterSdev = ogr.FieldDefn("RasterSdev", ogr.OFTReal)
        layer.CreateField(rasterMin)
        layer.CreateField(rasterMax)
        layer.CreateField(rasterMean)
        layer.CreateField(rasterSdev)
    
    # loop through each input file, create features, and output
    for file in params.inputFiles:
        
        # get the bounding box for this file
        minx, miny, maxx, maxy = aei.fn.rasterBoundingBox(file)
        
        # create the input feature
        feature = ogr.Feature(layer.GetLayerDefn())
        feature.SetField("FileName", file)
    
        # create a polygon from points
        box = ogr.Geometry(ogr.wkbLinearRing)
        box.AddPoint(minx,maxy)
        box.AddPoint(maxx,maxy)
        box.AddPoint(maxx,miny)
        box.AddPoint(minx,miny)
        box.AddPoint(minx,maxy)
        
        # set the geometry
        poly = ogr.Geometry(ogr.wkbPolygon)
        poly.AddGeometry(box)
        feature.SetGeometry(poly)
        
        # add stats info if set
        if params.getStats:
            gdalRef = gdal.Open(file)
            gdalBand = gdalRef.GetRasterBand(1)
            rmin, rmax, rmean, rstdev = gdalBand.ComputeStatistics(0)
            feature.SetField("RasterMin", rmin)
            feature.SetField("RasterMax", rmax)
            feature.SetField("RasterMean", rmean)
            feature.SetField("RasterSdev", rstdev)
            gdalBand = None
            gdalRef = None
        
        # create the feature
        layer.CreateFeature(feature)
        
        # destroy the feature reference
        feature.Destroy()

    # destroy the shape file reference
    dataSource.Destroy()

def multiOutput(params):
    
    # get spatial reference info
    srs = osr.SpatialReference()
    
    # loop through each file and output a new file
    for i in range(len(params.inputFiles)):
    
        # open up the output file with shapefile driver
        dataSource = params.driver.CreateDataSource(params.outputFiles[i])
    
        # get projection info from first raster input file
        inf = gdal.Open(params.inputFiles[i])
        wkt = inf.GetProjection()
        srs.ImportFromWkt(wkt)
        
        # create the output layer
        layer = dataSource.CreateLayer("BoundingBox", srs, ogr.wkbPolygon)
        
        # add field with file name
        fn = ogr.FieldDefn("FileName", ogr.OFTString)
        
        # find longest file name and use that for the length of the field
        fn.SetWidth(len(params.inputFiles[i]))
        layer.CreateField(fn)
    
        # add fields if -get_stats is set
        if params.getStats:
            rasterMin = ogr.FieldDefn("RasterMin", ogr.OFTReal)
            rasterMax = ogr.FieldDefn("RasterMax", ogr.OFTReal)
            rasterMean = ogr.FieldDefn("RasterMean", ogr.OFTReal)
            rasterSdev = ogr.FieldDefn("RasterSdev", ogr.OFTReal)
            layer.CreateField(rasterMin)
            layer.CreateField(rasterMax)
            layer.CreateField(rasterMean)
            layer.CreateField(rasterSdev)
    
        # get the bounding box for this file
        minx, miny, maxx, maxy = aei.fn.rasterBoundingBox(params.inputFiles[i])
        
        # create the input feature
        feature = ogr.Feature(layer.GetLayerDefn())
        feature.SetField("FileName", params.inputFiles[i])
    
        # create a polygon from points
        box = ogr.Geometry(ogr.wkbLinearRing)
        box.AddPoint(minx,maxy)
        box.AddPoint(maxx,maxy)
        box.AddPoint(maxx,miny)
        box.AddPoint(minx,miny)
        box.AddPoint(minx,maxy)
        
        # set the geometry
        poly = ogr.Geometry(ogr.wkbPolygon)
        poly.AddGeometry(box)
        feature.SetGeometry(poly)
        
        # add stats info if set
        if params.getStats:
            gdalRef = gdal.Open(params.inputFiles[i])
            gdalBand = gdalRef.GetRasterBand(1)
            rmin, rmax, rmean, rstdev = gdalBand.ComputeStatistics(0)
            feature.SetField("RasterMin", rmin)
            feature.SetField("RasterMax", rmax)
            feature.SetField("RasterMean", rmean)
            feature.SetField("RasterSdev", rstdev)
            gdalBand = None
            gdalRef = None
        
        # create the feature
        layer.CreateFeature(feature)
        
        # destroy the feature reference
        feature.Destroy()

        # destroy the shape file reference
        dataSource.Destroy()
    
# set usage string    
def usage(exit=False):
    """
    describes the bb2shp.py procedure in case of incorrect parameter calls
    
    syntax: usage(exit=False)
    """

    print(
        """
$ bb2shp.py -i input_files [-o output file] [-odir output_directory]
    [-of output_format] [-get_stats]
        """
        )
    
    if exit:    
        aei.params.sys.exit(1)
        
def main ():
    """
    the main program for bb2shp.py
    
    syntax: main()
    """
    
    # first, read the default parameters to get the processing object
    params = default_params()
    
    # then parse the arguments passed via command line
    args = aei.params.sys.argv
    parse_args(args, params)
    
    # check that input files were specified
    if not params.inputFiles:
        print("[ ERROR ]: No input file(s) specified")
        aei.params.sys.exit(1)
    
    # update the params class for any dependencies set at runtime
    update_params(params)
    
    # report starting
    print("[ STATUS ]: Running bb2shp.py")
    print("[ STATUS ]: %04d-%02d-%02d  %02d:%02d" % (params.now.year, 
      params.now.month, params.now.day, params.now.hour, params.now.minute))
    
    # set up different paths if creating single vs. multiple output files
    if params.singleOutput:
        singleOutput(params)
        
    else:
        multiOutput(params)
    
    # report success
    print("[ STATUS ]: Completed running bb2shp.py")
    print("[ STATUS ]: See output directory for final files: %s" % params.outputDir)
    
if __name__ == "__main__":
    main()