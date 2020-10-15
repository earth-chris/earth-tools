#####
# sets up frequently used functions
#
# c. 2016-2017 Christopher Anderson
#####

# set up a file check function
def checkFile(infile, quiet=False):
    import os

    if os.path.isfile(infile) and os.access(infile, os.R_OK):
        return True
    else:
        if not quiet:
            print("[ ERROR ]: Unable to read file: %s" % infile)
        return False


# set function to join strings
def strJoin(list, join=" "):
    outStr = join.join(map(str, list))
    return outStr


# set function to replace a file extension
def replaceExt(inputFile, ext, fullPath=False):
    import params

    if inputFile.find(".") == -1:
        if fullPath:
            outputFile = inputFile + ext
        else:
            outputFile = params.os.path.basename(inputFile) + ext

    else:
        if fullPath:
            outputFile = inputFile[: inputFile.find(".")] + ext
        else:
            outputFile = params.os.path.basename(inputFile[: inputFile.find(".")]) + ext

    return outputFile


# set function to return an array of random floats based on a min/max
def randomFloats(nIterations, minVal, maxVal):
    import numpy as np

    return np.random.ranf(nIterations) * (maxVal - minVal) + minVal


# set function to brightness normalize an array
def bn(array, axis=None, inds=[], returnScalar=False):
    """
    performs a brightness normalization on a numpy array

    usage: output = aei.bn(array, axis = None, inds = [])
      where: array = the input array to normalize
             axis  = the axis on which to normalize. default is the last dimension
             inds  = the indices
             returnScalar = set to True to return a second array with the brighness scalar
    """
    import numpy as np

    # if axis isn't set, default to last axis in array
    if axis == None:
        axis = array.ndim - 1

    # correct axis if invalid
    if axis > (array.ndim - 1):
        print("[ ERROR ]: invalid axis set: %s" % axis)
        print("[ ERROR ]: using dimenssion: %s" % (array.ndim - 1))
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
        scalar = np.sqrt((array[inds] ** 2).sum())
        bn = array[inds] / scalar

    elif array.ndim == 2:
        if axis == 0:
            scalar = np.expand_dims(np.sqrt((array[inds, :] ** 2).sum(axis)), axis)
            bn = array[inds, :] / scalar

        elif axis == 1:
            scalar = np.expand_dims(np.sqrt((array[:, inds] ** 2).sum(axis)), axis)
            bn = array[:, inds] / scalar

    elif array.ndim == 3:
        if axis == 0:
            scalar = np.expand_dims(np.sqrt((array[inds, :, :] ** 2).sum(axis)), axis)
            bn = array[inds, :, :] / scalar

        elif axis == 1:
            scalar = np.expand_dims(np.sqrt((array[:, inds, :] ** 2).sum(axis)), axis)
            bn = array[:, inds, :] / scalar

        elif axis == 2:
            scalar = np.expand_dims(np.sqrt((array[:, :, inds] ** 2).sum(axis)), axis)
            bn = array[:, :, inds] / scalar

    else:
        print("[ ERROR ]: unable to brightness normalize")
        print("[ ERROR ]: array must be 3 dimensions or smaller")
        print("[ ERROR ]: array provided has [%s] dimensions" % array.ndim)
        return -1

    if returnScalar:
        return bn, scalar
    else:
        return bn


#####
# Raster functions
#####


def rasterTile(infile, outfile, tiling=2, buff=0, of="GTiff"):
    """
    tiles a raster file into a set number of tiles based on
    a tiling parameter.

    syntax: tileRaster(infile, outfile, tiling=2, buff=0, of='GTiff')

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
    import gdal as gdal
    import osr as osr
    import cmd as cmd
    import numpy as np
    import math

    # set defaults for ensuring i/o is smooth
    if not checkFile(infile):
        return -1

    if not type(buff) is int:
        print("[ ERROR ]: Buffer specified (%s) is not an integer" % buff)
        return -1

    if of == "GTiff":
        etc = ["-co COMPRESS=LZW"]
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
        xmax = np.arange(tiling) + 1
        ymin = np.arange(tiling) + 1
        ymax = np.arange(tiling)

        for i in range(tiling - 1):
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
        xmax = np.arange(nxt) + 1
        ymin = np.arange(nyt) + 1
        ymax = np.arange(nyt)

        for i in range(nxt - 1):
            xmin = np.row_stack((xmin, xmin))
            xmax = np.row_stack((xmax, xmax))

        for i in range(nyt - 1):
            ymin = np.row_stack((ymin, ymin))
            ymax = np.row_stack((ymax, ymax))

        # ymin = ymin.transpose()
        # ymax = ymax.transpose()

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
    if xmax[0, nxt - 1] > (ulx + (nx * sx) + buff):
        xmax[:, nxt - 1] = ulx + (nx * sx) + buff

    if ymin[nyt - 1, 0] < (uly + (ny * sy) - buff):
        ymin[nyt - 1, :] = uly + (ny * sy) - buff

    # report status
    print("[ STATUS ]: Tiling %s into %s output tiles") % (infile, nxt * nyt)

    # set up loop to run gdalwarp on tiles
    cnt = 0
    for i in range(nxt):
        for j in range(nyt):
            cnt += 1
            outtmp = (outfile + "_%03.f" + append) % cnt
            cmd.gdalwarp(
                infile,
                outtmp,
                etc,
                multi=True,
                of=of,
                te=[xmin[i, j], ymin[j, i], xmax[i, j], ymax[j, i]],
                overwrite=True,
            )


def rasterBoundingBox(infile):
    """
    creates a shape file from the bounding box of a raster file
    """
    import ogr as ogr
    import osr as osr
    import gdal as gdal

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

    # calculate minx, miny
    lrx = ulx + (sx * nx)
    lry = uly + (sy * ny)

    # close the raster file
    inf = None

    # return in order of minx, miny, maxx, maxy
    return [ulx, lry, lrx, uly]


# function to return all possible pair-wise combinations of a list
def pairwise(inlist, n_combos=2, all_combos=False):
    """Returns a list of all unique pairwise combinations of an input list

    Args:
        inlist: the list of items to generate pairwise combinations from
        n_combos: the number of items to match (2 for pairwise, 3 for unique groups of 3, etc.)
        all_combos: set if all combinations should be returned instead of only unique combos

    Returns:
        outlist: a list of pairwise combinations for each item in inlist
    """
    import itertools

    # if all combos are specified, use a separate iterator
    if all_combos:
        iterator = itertools.combinations_with_replacement
    else:
        iterator = itertools.combinations

    # populate the iterator
    combos = iterator(inlist, n_combos)

    # create the output list to return, and add all the combinations
    outlist = []
    for unique in combos:
        outlist.append(unique)

    return outlist
