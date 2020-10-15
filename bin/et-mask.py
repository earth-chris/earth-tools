#!/usr/bin/python
#####
# et-mask.py
#
#  masks a raster file using other raster files of the same dimensions.
#  primary use is to ensure no-data pixels from one file is aligned with another.
#
# c. 2016 Christopher Anderson
#####

import sys

import gdal as gdal
import numpy as np

import earthtools as et


class parse_args:
    def __init__(self, arglist):

        # set up main variables and defaults to parse
        self.refFile = None
        self.maskFiles = []
        self.nMaskFiles = 0
        self.bands = []
        self.subsetBands = False
        self.refNoData = None
        self.maskNoData = None
        self.useMaskNoData = False
        self.useRefNoData = False
        self.fileChanged = False

        # exit if no arguments passed
        if len(arglist) == 1:
            usage(exit=True)

        # read arguments from command line
        i = 1
        while i < len(arglist):
            arg = arglist[i]

            # parse the base reference file flag
            if arg.lower() == "-ref":
                i += 1
                arg = arglist[i]

                if not et.fn.checkFile(arg, quiet=True):
                    usage()
                    et.fn.checkFile(arg)
                    et.params.sys.exit(1)

                self.refFile = arg

            # check mask flag
            elif arg.lower() == "-mask":
                i += 1
                arg = arglist[i]

                # loop through masks until we find a new parameter
                newArg = False
                while not newArg:
                    if arg[0] == "-":
                        newArg = True
                        i -= 1
                        continue

                    self.maskFiles.append(arg)

                    # increment counter and move on
                    i += 1
                    if i >= len(arglist):
                        newArg = True
                        i -= 1
                        continue
                    arg = arglist[i]

            # check -b option
            elif arg.lower() == "-b" or arg.lower() == "-band":
                i += 1
                arg = arglist[i]

                # set the band flag
                self.subsetBands = True

                # loop through band args until we find a new parameter
                newArg = False
                while not newArg:
                    if arg[0] == "-":
                        newArg = True
                        i -= 1
                        continue
                    # ensure integers set
                    try:
                        arg = int(arg)
                    except:
                        usage()
                        print("[ ERROR ]: Invalid argument passed to -b: %s" % arg)
                        sys.exit(1)

                    self.bands.append(arg)

                    # increment counter and move on
                    i += 1
                    if i >= len(arglist):
                        newArg = True
                        i -= 1
                        continue
                    arg = arglist[i]

            # use a user-defined no-data value if set
            elif arg.lower() == "-masknodata":
                i += 1
                arg = arglist[i]

                # set the parameter
                self.maskNoData = float(arg)
                self.useMaskNoData = True

            # use a user-defined ref no-data value if set
            elif arg.lower() == "-refnodata":
                i += 1
                arg = arglist[i]

                # set the parameter
                self.refNoData = float(arg)
                self.useRefNoData = True

            # set up catch-all for incorrect parameter call
            else:
                usage()
                print("[ ERROR ]: Unrecognized argument: %s" % arg)
                print("%s" % i)
                et.params.sys.exit(1)

            # increment the counter for next argument
            i += 1


# set up a function to write an array to a raster band
def writeArrayToRaster(array, gdalRef, band, noData=None):
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


# set up function to check that arguments have been properly specified
def check_args(args):
    # check there is a reference file
    if len(args.refFile) == 0:
        usage()
        print("[ ERROR ]: No base reference file specified.")
        sys.exit(1)

    # check number of mask files
    args.nMaskFiles = len(args.maskFiles)
    if args.nMaskFiles == 0:
        usage()
        print("[ ERROR ]: No mask files specified")
        sys.exit(1)

    # return fixed arguments
    return args


def usage(exit=False):
    """
    describes the et-mask.py procedure in case of incorrect parameter calls

    syntax: usage()
    """
    print(
        """
$ et-mask.py -ref fileToMask -mask inputMask
        [-b bandList] [-masknodata noDataVal] [-refnodata noDataVal]
        """
    )
    if exit:
        sys.exit(1)


def main():
    """
    the main program for et-mask.py

    the order of operations for this procedure is to
    1) ensure the masks are the same dimensions
    2) loop through each band, or the specified bands, of the mask files and
       mask the input reference file.

    syntax: main()
    """
    # parse the argument list
    args = parse_args(sys.argv)

    # check the argument list to ensure consistent arguments are set
    args = check_args(args)

    # report starting
    print("[ STATUS ]: Starting et-mask.py")
    print("[ STATUS ]: Reference file : %s" % args.refFile)
    print("[ STATUS ]: Mask file(s)   : %s" % et.fn.strJoin(args.maskFiles))
    print("[ STATUS ]: ----------")
    print("[ STATUS ]: Reading reference data")

    # open the reference file
    refFile = gdal.Open(args.refFile, 1)

    # get nx, ny info
    refnx = refFile.RasterXSize
    refny = refFile.RasterYSize
    refnb = refFile.RasterCount

    # read the raster data into an array
    refArr = refFile.ReadAsArray()

    # determine the no-data value to use if not set
    if not args.useRefNoData:
        refBand = refFile.GetRasterBand(1)
        args.refNoData = refBand.GetNoDataValue()

        # check that it exists
        if args.refNoData is None:
            print("[ ERROR ]: No data value not specified in reference file")
            print("[ ERROR ]: Please specify using -refnodata, or update ref file")
            sys.exit(1)

        # kill the reference
        refBand = None

    # loop through each mask file
    for i in range(args.nMaskFiles):

        # report
        print("[ STATUS ]: Reading mask file %02d: %s" % (i + 1, args.maskFiles[i]))

        # open the mask file reference
        maskRef = gdal.Open(args.maskFiles[i])

        # get x, y, nb info
        masknx = maskRef.RasterXSize
        maskny = maskRef.RasterYSize
        masknb = maskRef.RasterCount

        # ignore this file if x and y dims don't match
        if (refnx != masknx) or (refny != maskny):
            print("[ ERROR ]: Dimensions of mask and reference files do not match")
            print("[ ERROR ]: Reference X, Y : %s, %s" % (refnx, refny))
            print("[ ERROR ]: Mask file X, Y : %s, %s" % (masknx, maskny))
            continue

        # check the flag to see if -b was set
        if args.subsetBands == True:

            # ensure that the bands set make sense
            if max(args.bands) > masknb:
                print("[ ERROR ]: Invalid bands set to subset: %s" % max(args.bands))
                print("[ ERROR ]: Total bands in mask file   : %s" % masknb)
                bands = range(1, masknb + 1)
            else:
                bands = args.bands

        else:
            bands = range(1, masknb + 1)

        # loop through each band and get the indices
        for b in bands:

            # report
            print("[ STATUS ]: Reading band: %03d" % b)

            # read the band reference
            bandRef = maskRef.GetRasterBand(b)

            # get the noData val if not set
            if not args.useMaskNoData:
                args.maskNoData = bandRef.GetNoDataValue()

                # make sure it exists
                if args.maskNoData is None:
                    print("[ ERROR ]: No data value not specified in mask band: %s" % b)
                    print("[ ERROR ]: Please specify using -masknodata, or update mask file")
                    continue

            # read the data as an array
            bandArr = bandRef.ReadAsArray()

            # find the no-data indices
            nd = np.where(bandArr == args.maskNoData)

            # check that it is not empty
            if not nd[0].any():
                print("[ ERROR ]: No no-data values found in band: %03d" % b)
                bandArr = None
                continue

            # mask the reference data
            if refnb > 1:
                refArr[:, nd[0], nd[1]] = args.refNoData
            else:
                refArr[nd[0], nd[1]] = args.refNoData

            # kill some references
            bandRef = None
            bandArr = None

            # use this flag if the ref file needs to be updated
            args.fileChanged = True

        # kill some references
        maskRef = None

    # now that each mask band has been read, write the output file
    if args.fileChanged:
        print("[ STATUS ]: Writing output file")
        if refnb > 1:
            for i in range(refnb):
                writeArrayToRaster(refArr[i], refFile, i + 1, noData=args.refNoData)
        else:
            writeArrayToRaster(refArr, refFile, 1, noData=args.refNoData)

    else:
        print("[ ERROR ]: No updates were made to input refernce file")

    # kill some references
    refArr = None
    refFile = None

    # report finished
    print("[ STATUS ]: ----------")
    print("[ STATUS ]: Finished masking raster data!")
    print("[ STATUS ]: Please see output file : %s" % args.refFile)


# call the aain routine when run from command lne
if __name__ == "__main__":
    main()
