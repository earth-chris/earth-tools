#!/usr/bin/python
#####
# et-palsar.py
#  
#  computes standard data products for ALOS-PALSAR imagery
#
# c. 2017 Christopher Anderson
#####

import os
import earthtools as et
import argparse


# function to parse the input arguments
def parse_args():
    # create the parser object
    parser = argparse.ArgumentParser(description="Processes ALOS-PALSAR data")

    # create options for input files
    parser.add_argument('--hh', metavar='hh', type=str,
                        nargs=1, help='path to the HH polarization data',
                        required=True)

    parser.add_argument('--hv', metavar='hv', type=str,
                        nargs=1, help='path to the HV polarization data',
                        required=True)

    # create option for output file
    parser.add_argument('-o', metavar='output_file', type=str,
                        nargs=1, help='path to the output file',
                        required=True)

    # create options for the despeckling algorithm to use
    parser.add_argument('--filter', metavar='filter', type=str,
                        nargs=1, help='the filter type',
                        choices=['lee', 'frost', 'gammamap', 'kuan'],
                        default='lee', required=False)

    # options for the despeckling radius    
    parser.add_argument('--rad', metavar='radius', type=int,
                        nargs=1, help='the despeckling radius',
                        default=3, required=False)

    # option to set the output nodata value    
    parser.add_argument('--nodata', metavar='no data val', nargs=1,
                        help='the output nodata value', default=-9999,
                        required=False)

    # option to not mask the output file
    parser.add_argument('--nomask', action='store_true',
                        help='set to prevent masking', required=False)

    # parse the arguments and return the object
    args = parser.parse_args()
    return args


# function to calculate the raw data to decibel
def raw_to_decibel(infile, outfile):
    # add compression to output file
    outfile = '"{}?&gdal:co:COMPRESS=LZW"'.format(outfile)

    # run the command
    et.cmd.otb.BandMath(infile, outfile, '"10*log10(im1b1^2)-83"')


# function to despeckle the data
def despeckle(infile, outfile, filtername, radius):
    # add compression to output file
    outfile = '"{}?&gdal:co:COMPRESS=LZW"'.format(outfile)

    # run the command
    cmd = "otbcli_Despeckle -in {} -out {} -filter {} -filter.{}.rad {}".format(
        infile, outfile, filtername, filtername, radius)
    print(cmd)
    et.cmd.run(cmd)


# function to calculate the polarimetric entropy
def calc_polarimetry(hh, hv, outfile):
    # add compression to output file
    outfile = '"{}?&gdal:co:COMPRESS=LZW"'.format(outfile)

    # build the expression string
    exp = '"(-1)/(im1b1+im2b1)*(im1b1*log2(im1b1/(im1b1+im2b1))+im2b1*log2(im2b1/(im1b1+im2b1)))"'

    # run the command
    et.cmd.otb.BandMath([hh, hv], outfile, exp)


# function to stack all three images
def stack_outputs(hh, hv, pol, outfile):
    # add compression to output file
    outfile = '"{}?&gdal:co:COMPRESS=LZW"'.format(outfile)

    # run the command
    et.cmd.otb.ConcatenateImages([hh, hv, pol], outfile)


# function to build a mask from the input data
def get_mask(infile, outfile):
    # add compression to output file
    outfile = '"{}?&gdal:co:COMPRESS=LZW"'.format(outfile)

    # run the command
    cmd = "otbcli_ManageNoData -in {} -out {} uint8".format(
        infile, outfile)
    print(cmd)
    et.cmd.run(cmd)


# function to apply the mask
def apply_mask(infile, outfile, maskfile, ndval):
    # add compression to output file
    outfile = '"{}?&gdal:co:COMPRESS=LZW"'.format(outfile)

    # run the command
    cmd = "otbcli_ManageNoData -in {} -out {} -mode apply -mode.apply.mask {} -mode.apply.ndval {}".format(
        infile, outfile, maskfile, ndval)
    print(cmd)
    et.cmd.run(cmd)


# the main program function
def main():
    # parse the input bash arguments
    args = parse_args()

    # create three temp files for hh, hv, polarimetry data
    hh_dec = 'temp_HH_decibel.tif'
    hh_msk = 'temp_HH_decibel_masked.tif'
    hh_des = 'temp_HH_despeckle.tif'
    hv_dec = 'temp_HV_decibel.tif'
    hv_msk = 'temp_HV_decibel_masked.tif'
    hv_des = 'temp_HV_despeckle.tif'
    polari = 'temp_polarimetry.tif'
    tmpmsk = 'temp_mask.tif'
    premsk = 'temp_stack.tif'

    # build a mask file, if set
    if not args.nomask:
        get_mask(args.hh[0], tmpmsk)

    # calculate the decibel data
    raw_to_decibel(args.hh[0], hh_dec)
    raw_to_decibel(args.hv[0], hv_dec)

    # mask the decibel data prior to despeckling
    if not args.nomask:
        apply_mask(hh_dec, hh_msk, tmpmsk, args.nodata)
        apply_mask(hv_dec, hv_msk, tmpmsk, args.nodata)
        os.remove(hh_dec)
        os.remove(hv_dec)
        hh_dec = hh_msk
        hv_dec = hv_msk

    # despeckle the decible data
    despeckle(hh_dec, hh_des, args.filter, args.rad)
    despeckle(hv_dec, hv_des, args.filter, args.rad)

    # calculate the polarimetry
    calc_polarimetry(hh_des, hv_des, polari)

    # stack the outputs in different order if masking set
    if args.nomask:

        # concatenate all files into a single stack
        stack_outputs(hh_des, hv_des, polari, args.o[0])

    else:

        # concatenate all files into a single stack
        stack_outputs(hh_des, hv_des, polari, premsk)

        # then apply the mask
        apply_mask(premsk, args.o[0], tmpmsk, args.nodata)

    # delete the temp files
    # os.remove(hh_dec)
    # os.remove(hv_dec)
    # os.remove(hh_des)
    # os.remove(hv_des)
    # os.remove(polari)
    # os.remove(tmpmsk)
    # os.remove(premsk)

    # report finished
    print("[ STATUS ]: Finished ALOS-PALSAR processing!")


# call the main routine when run from command lne
if __name__ == "__main__":
    main()
