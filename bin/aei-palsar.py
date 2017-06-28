#!/usr/bin/python
#####
# aei-palsar.py
#  
#  computes standard data products for ALOS-PALSAR imagery
#
# c. 2017 Christopher Anderson
#####

import os
import sys
import aei
import gdal as gdal
import osr as osr
import numpy as np
import argparse

# function to parse the input arguments
class parse_args:
    def __init__(self):
        
        # create the parser object
        parser = argparse.ArgumentParser(description = "Processes ALOS-PALSAR data")
        
        # create options for input files
        parser.add_argument('--hh', metavar = 'hh', type = str,
            nargs = 1, help = 'path to the HH polarization data',
            required = True)
        
        parser.add_argument('--hv', metavar = 'hv', type = str,
            nargs = 1, help = 'path to the HV polarization data',
            required = True)
            
        # create option for output file
        parser.add_argument('-o', metavar = 'output_file', type = str,
            nargs = 1, help = 'path to the output file',
            required = True)
        
        # parse the arguments and return the object
        args = parser.parse_args()
        return args
        
# function to calculate the raw data to decibel
def raw_to_decibel(infile, outfile):
    
    # add compression to output file
    outfile = '"{}?&gdal:co:COMPRESS=LZW"'.format(outfile)
    
    # run the command
    aei.cmd.otb.BandMath(infile, outfile, "10*log10(im1b1^2)-83")
    
# function to despeckle the data
def despeckle(infile, outfile, filtername):
    
    # add compression to output file
    outfile = '"{}?&gdal:co:COMPRESS=LZW"'.format(outfile)
    
    # run the command
    cmd = "otbcli_Despeckle -in {} -out {} -filter {} -filter.{}.rad".format(
        infile, outfile, filtername)
    aei.cmd.run(cmd)
    
# function to calculate the polarimetric entropy
def calc_polarimetry(hh, hv, outfile):
    
    # add compression to output file
    outfile = '"{}?&gdal:co:COMPRESS=LZW"'.format(outfile)
    
    # build the expression string
    exp = "(-1)/(im1b1+im2b1)*(im1b1*log2(im1b1/(im1b1+im2b1))+im2b1*log2(im2b1/(im1b1+im2b1)))"
    
    # run the command
    aei.cmd.otb.BandMath([hh, hv], outfile, exp)
    
# function to stack all three images
def stack_outputs(hh, hv, pol, outfile):
    
    # add compression to output file
    outfile = '"{}?&gdal:co:COMPRESS=LZW"'.format(outfile)
    
    # run the command
    aei.cmd.otb.ConcatenateImages([hh, hv, pol], outfile)
    
# the main program function
def main():
    
    # parse the input arguments
    args = parse_args()
    
    # create three temp files for hh, hv, polarimetry data
    hh_dec = 'temp_HH_decibel.tif'
    hh_des = 'temp_HH_despeckle.tif'
    hv_dec = 'temp_HV_decibel.tif'
    hv_des = 'temp_HV_despeclkle.tif'
    polari = 'temp_polarimetry.tif'
    
    # calculate the decibel and polarimetery data
    raw_to_decibel(args.hh, hh_dec)
    raw_to_decibel(args.hv, hv_dec)
    calc_polarimetry(hh_dec, hv_dec, polari)
    
    # concatenate all files into a single stack
    stack_outputs(hh_dec, hv_dec, polari, args.o)
    
    # delete the temp files
    os.remove(hh_dec)
    os.remove(hv_dec)
    os.remove(polari)
    
    # report finished
    print("[ STATUS ]: Finished ALOS-PALSAR processing!")
    
# call the aain routine when run from command lne
if __name__ == "__main__":
    main()