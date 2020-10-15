#!/usr/bin/python
#####
# calculates the statistical moments for each band in a raster file
#
# c. 2017 Christopher Anderson
#####

import argparse
import gdal
import sys
import numpy as np
from scipy import stats

import earthtools as et


# function to parse the input arguments
def parse_args():
    # create the parser object
    parser = argparse.ArgumentParser(description="Processes ALOS-PALSAR data")

    # create options for input files
    parser.add_argument(
        "-i", metavar="input_file", type=str, nargs=1, help="path to the input raster data", required=True
    )

    # create option for output file
    parser.add_argument("-o", metavar="output_file", type=str, nargs=1, help="path to the output file", required=True)

    # create option for setting output resolution
    parser.add_argument(
        "--tr", metavar="resolution", type=float, nargs=2, help="the target spatial resolution", required=False
    )

    # create option for setting output window size
    parser.add_argument(
        "--ws",
        metavar="window_size",
        type=int,
        nargs=2,
        help="the window size to calculate moments for",
        required=False,
    )

    # create option for the number of moments to calculate
    parser.add_argument(
        "-m", metavar="moments", type=int, nargs=1, help="the number of moments to calculate", default=4, required=False
    )

    # create option for output no-data value
    parser.add_argument(
        "--dstnodata",
        metavar="nodata",
        nargs=1,
        type=float,
        help="the output nodata value",
        default=-9999.0,
        required=False,
    )

    # parse the arguments and return the object
    args = parser.parse_args()
    return args


# function to determine output raster size based on target resolution
def parse_tr(args, inref):
    # first, determine the window size to calculate moments for
    xsize = args.tr[0] / inref.xps
    ysize = abs(args.tr[1] / inref.yps)

    # next, calculate dimensions for the output raster
    xres = int(np.ceil(inref.nx / xsize))
    yres = int(np.ceil(inref.ny / ysize))

    # create the output array of size (n_moments, output_yres, output_xres)
    arr = np.zeros((args.m, yres, xres))

    # initialize to the nodata value
    arr += args.dstnodata

    # create lists of the x and y indices
    xstart = np.round(xsize * np.arange(xres)).astype(long)
    xend = np.round(xsize * (np.arange(xres) + 1)).astype(long)
    xend[-1] = inref.nx

    ystart = np.round(ysize * np.arange(yres)).astype(long)
    yend = np.round(ysize * (np.arange(yres) + 1)).astype(long)
    yend[-1] = inref.ny

    # update the raster reference with new metadata
    inref.xps = args.tr[0]
    inref.yps = -args.tr[1]
    inref.nx = xres
    inref.ny = yres
    inref.dt = gdal.GDT_Float32

    # return the output array and the x/y indices
    return [arr, xstart, xend, ystart, yend]


# function to determine output raster size based on target resolution
def parse_ws(args, inref):
    # first, set the window size parameters
    xsize = args.ws[0]
    ysize = args.ws[1]

    # next, calculate dimensions for the output raster
    xres = int(np.ceil(inref.nx / xsize))
    yres = int(np.ceil(inref.ny / ysize))

    # create the output array of size (n_moments, output_yres, output_xres)
    arr = np.zeros((args.m, yres, xres))

    # initialize to the nodata value
    arr += args.dstnodata

    # create lists of the x and y indices
    xstart = np.round(xsize * np.arange(xres)).astype(long)
    xend = np.round(xsize * (np.arange(xres) + 1)).astype(long)
    xend[-1] = inref.nx

    ystart = np.round(ysize * np.arange(yres)).astype(long)
    yend = np.round(ysize * (np.arange(yres) + 1)).astype(long)
    yend[-1] = inref.ny

    # update the raster reference with new metadata
    inref.xps = xsize * inref.xps
    inref.yps = ysize * inref.yps
    inref.nx = xres
    inref.ny = yres
    inref.dt = gdal.GDT_Float32

    # return the output array and the x/y indices
    return [arr, xstart, xend, ystart, yend]


# the main program function
def main():
    # parse the input bash arguments
    args = parse_args()

    # read the raster metadata into memory
    inref = et.read.raster(args.i[0])

    # first, we'll determine the output file parameters and how we loop through the raster
    # check that either window size or output resolution are set.
    #  if not: return an error and quit
    #  if so: return the parameters for the output raster
    if args.tr[0] is not None:
        output_array, xstart, xend, ystart, yend = parse_tr(args, inref)
    elif args.ws[0] is not None:
        output_array, xstart, xend, ystart, yend = parse_ws(args, inref)
    else:
        print("[ ERROR ]: Neither --tr nor --ws was set")
        print("[ ERROR ]: Set one of these and re-run")
        sys.exit(1)

    # now that we have parsed the data to establish final x/y resolutions, let's create
    #  an output file with our new metadata and a new variable to work with
    outref = inref.copy(args.o[0], nb=inref.nb * args.m, dt=gdal.GDT_Float32)

    # we'll loop through each band in the input raster and calculate the moments for each pixel
    oband = 1
    for i in range(inref.nb):
        # read the data for this band into memory
        print("[ STATUS ]: Reading band: {:03d}".format(i + 1))
        inref.read_band(i + 1)

        # if there is a no-data value, set the no-data pixels to NaN to ignore in calculations
        if inref.no_data is not None:
            inref.data = inref.data.astype(float)
            inref.data[inref.data == inref.no_data] = np.nan

        # and we'll loop through each window (by row, column)
        for j in range(len(ystart)):
            for k in range(len(xstart)):
                window = inref.data[ystart[j] : yend[j], xstart[k] : xend[k]]

                # if there is no good data, skip this window
                if np.isnan(window).all():
                    continue

                # otherwise, compute the moments and store in the output array
                for l in range(args.m):
                    if l == 0:
                        moment = np.nanmean(window)
                    else:
                        moment = stats.moment(window, l + 1, nan_policy="omit", axis=None)
                    output_array[l, j, k] = moment

        # write the output to disk
        for l in range(args.m):
            outref.write_band(oband, output_array[l, :, :])
            oband += 1

    # update the no-data values
    outref.no_data = args.dstnodata
    outref.write_metadata()

    # report finished
    print("[ STATUS ]: Finished writing output file: {}".format(args.o[0]))
    print("[ STATUS ]: et-moment.py complete!")


# call the main routine when run from command lne
if __name__ == "__main__":
    main()
