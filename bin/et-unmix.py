#!/usr/bin/python
#####
# aei_unmix.py performs least-squares 
#  spectral ummixing on an image file using the 
#  HyperspectralUnmixing module from the Orfeo Toolbox
#
# c. 2016 Christopher Anderson
#####

import os
import sys
import aei
import random
import gdal as gdal
import numpy as np

# create a class to parse out the arguments passed to the main function
class parse_args:
    def __init__(self, arglist):
        
        # set up main variables and defaults to parse
        self.infile = ''
        self.outfile = ''
        self.outpath = ''
        self.spectral_libs = []
        self.n = 20
        self.ua = 'ucls'
        self.bands = []
        self.subset_bands = False
        self.numpy_dtype = np.float32
        self.gdal_dtype = gdal.GDT_Float32
        self.of = 'GTiff'
        self.normalize = False
        self.rescale = False
        self.rescale_val = 1.0
        self.temp_libs = []
        self.temp_mixtures = []
        self.temp_averages = []
        self.temp_stdevs = []
        self.temp_vrt = 'temp_vrt.vrt'
        self.include_shade = False
        
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
                
                if type(arg) is str:
                    self.infile = arg
                    if not aei.fn.checkFile(self.infile, quiet = True):
                        usage()
                        aei.fn.checkFile(self.infile)
                        aei.params.sys.exit(1)
            
            # check output flag
            elif arg.lower() == '-o':
                i += 1
                arg = arglist[i]
                
                if type(arg) is str:
                    self.outfile = arg
                    outpath = os.path.dirname(self.outfile)
                    if outpath == '':
                        outpath = '.'
                    if not os.access(outpath, os.W_OK):
                        usage()
                        print("[ ERROR ]: unable to write to output path: %s" % outpath)
                        aei.params.sys.exit(1)
                        
                    self.outpath = outpath
            
            # check spectral library paths            
            elif arg.lower() == "-lib":
                i += 1
                arg = arglist[i]
                libs = arg.split(" ")
                
                # throw an error if only one lib specified
                if len(libs) == 1:
                    usage()
                    print("[ ERROR ]: unable to unmix with one spectral library: %s" % libs[0])
                    aei.params.sys.exit(1)
                
                # loop through each lib and update spec_lib list
                for j in range(len(libs)):
                    if not aei.fn.checkFile(libs[j], quiet=True):
                        usage()
                        aei.fn.checkFile(libs[j])
                        aei.params.sys.exit()
                        
                    self.spectral_libs.append(libs[j])
                    
            # check number of iterations
            elif arg.lower() == "-n":
                i += 1
                arg = arglist[i]
                
                try:
                    self.n = int(arg)
                except ValueError:
                    usage()
                    print("[ ERROR ]: -n argument is not an integer: %s" % arg)
                    aei.params.sys.exit(1)
                    
            # check indices
            elif arg.lower() == "-bands":
                i += 1
                arg = arglist[i]
                band = arg.split(" ")
                
                # loop through and make sure each is a number
                for j in range(len(band)):
                    
                    try:
                        int(band[j])
                    except ValueError:
                        usage()
                        print("[ ERROR ]: invalid index set: %s" % ind[j])
                        aei.params.sys.exit(1)
                    
                    self.bands.append(int(band[j]))
                    
                # set flag to indicate to subset the bands
                self.subset_bands = True
                    
            # check normalize flag
            elif arg.lower() == "-normalize":
                self.normalize = True
                
            # check output format
            elif arg.lower() == "-of":
                i += 1
                arg = arglist[i]
                
                if not arg in aei.params.gdalTypes:
                    usage()
                    print("[ ERROR }: invalid output format: %s" % arg)
                    aei.params.sys.exit(1)
                    
                self.of = arg
                
            # check if we should include photometric shade (i.e. a vector of 0s)
            elif arg.lower() == '-include_shade':
                self.include_shade = True
                
            # check if rescaling parameter is set to scale spec libs to image
            elif arg.lower() == '-rescale':
                i += 1
                arg = arglist[i]
                
                self.rescale = True
                self.rescale_val = arglist[i]
                
            # check the unmixing algorithm to use
            elif arg.lower() == '-ua':
                i += 1
                arg = arglist[i]
                
                # check that algorithm set is legit
                if not arg.lower() in ['ucls', 'ncls', 'isra', 'mdmdnmf']:
                    usage()
                    print("[ ERROR ]: Unsupported unmixing method: %s" % arg)
                    print("[ ERROR ]: Options are: ucls / ncls / isra / mdmdnmf")
                    aei.params.sys.exit(1)
                
                self.ua = arg
                
            # set up catch-all for incorrect parameter call
            else:
                usage()
                print("[ ERROR ]: Unrecognized argument: %s" % arg)
                print("%s" % i)
                aei.params.sys.exit(1)
            
            i += 1

# set up function to convert between gdal and numpy data types for file reading
def get_numpy_data_type(args):
    if args.gdal_dtype == 1:
        args.numpy_dtype = np.byte
    elif args.gdal_dtype == 2:
        args.numpy_dtype = np.uint16
    elif args.gdal_dtype == 3:
        args.numpy_dtype = np.int16
    elif args.gdal_dtype == 4:
        args.numpy_dtype = np.uint32
    elif args.gdal_dtype == 5:
        args.numpy_dtype = np.int32
    elif args.gdal_dtype == 6:
        args.numpy_dtype = np.float32
    elif args.gdal_dtype == 7:
        args.numpy_dtype = np.float64
    else:
        print("[ ERROR ]: Unrecognized gdal data type: %s" % args.gdal_dtype)

# set up function to write the temporary endmember libraries for otbcli
def write_temp_libraries(args,lib,n_libs,lnb):
    
    # loop through each library, write temp files, and store info for later deletion
    for i in range(0, args.n):
        
        # set up the output file name
        temp_file = args.outpath + aei.params.pathsep + str('temp_lib_%02d.tif' % (i+1))
        args.temp_libs.append(temp_file)
        
        # create array to hold the endmembers. if shade is included, add an
        #  additional endmember
        if args.include_shade:
            bundles = np.zeros((n_libs + 1, lnb[0]))
            nl = n_libs + 1
        else:
            bundles = np.zeros((n_libs, lnb[0]))
            nl = n_libs
        
        # load random spectra from each lib into the output array
        for j in range(n_libs):
            bundles[j,:] = lib['lib_%s' % j].spectra[lib['rnd_%s' % j][i]]
            
        # rescale the spectral libraries if set
        if args.rescale:
            bundles *= float(args.rescale_val)
        
        # open the temp file
        file_ref = gdal.GetDriverByName("GTiff").Create(temp_file, 1, nl, \
          int(lnb[0]), args.gdal_dtype)
        
        # write the arrays band by band  
        for j in range(int(lnb[0])):
            band = file_ref.GetRasterBand(j + 1)
            band.WriteArray(np.expand_dims(bundles[:,j], 1).astype(args.numpy_dtype))
            band.FlushCache()
            
        # and clear references from memory
        band = None
        file_ref = None

def usage(exit=False):
    """
    describes the aei-unmix.py procedure in case of incorrect parameter calls
    
    syntax: usage()
    """
    print(
        """
$ aei-unmix.py -lib "lib1 lib2 ... libx" [-n n_random_selections]
      [-bands "band1 band2 ... bandx"] [-normalize] [-of output_format]
      [-rescale rescale_parameter] [-include_shade] [-ua unmixing_algorithm]
      [-i] input_file [-o] output_file
        """
        )
    if exit:
        sys.exit(1)

def main():
    """
    the main program for aei-unmix.py
    
    syntax: main()
    """
    # parse the argument list
    args = parse_args(sys.argv)
    
    # load the spectral libraries
    lib = {}
    n_libs = len(args.spectral_libs)
    lnb = np.zeros(n_libs, dtype=np.int16)
    
    # load the input image and get parameters
    inf = gdal.Open(args.infile)
    ns = inf.RasterXSize
    nl = inf.RasterYSize
    nb = inf.RasterCount
    geo = inf.GetGeoTransform()
    prj = inf.GetProjection()
    n_bundles = lnb.shape[0]
    b1 = inf.GetRasterBand(1)
    
    # get raster data type
    args.gdal_dtype = b1.DataType
    
    # determine data types for gdal and numpy
    get_numpy_data_type(args)
    
    # destroy the input file reference
    inf = None
    b1 = None
    
    # get info from each library
    for i in range(n_libs):
        
        # assign libraries to a dictionary
        lib['lib_%s' % i] = aei.read.spectralLib(args.spectral_libs[i])
        lnb[i] = (lib['lib_%s' % i].spectra.shape[-1])
        
        # get indices to randomly sample each library
        lib['rnd_%s' % i] = random.sample(range(
            lib['lib_%s' % i].spectra.shape[0]), args.n)
            
        # and brightness normalize if set
        if args.normalize:
            lib['lib_%s' % i].bn(inds=args.bands)
            
    # write randomly sampled libraries to new files for otbcli input
    write_temp_libraries(args,lib,n_libs,lnb)
    
    # if bands were not set in command line, use all bands
    if not args.subset_bands:
        args.bands = range(nb)
    
    # check that the spectral libraries match the bands
    if nb != int(lnb[0]):
        print(nb, int(lnb[0]))
        print("[ ERROR ]: Number of image bands does not match number of spectral library bands")
        sys.exit(1)
    
    # report a little
    print("[ STATUS ]: Beginning aei-unmix.py")
    print("[ STATUS ]: Input file : %s" % args.infile)
    print("[ STATUS ]: Output file: %s" % args.outfile)
        
    # normalize the data if set
    #if args.normalize:
    #    img = aei.bn(img, inds = args.bands)
    #    args.bands = range(len(args.bands))
    
    # set up a dictionary with the arguments to be passed to the otb band math command
    bm_args = {}
    for i in range(n_libs):
        bm_args['b%s' % i] = []
        
    # loop through each random index and unmix
    for i in range(args.n):
        
        # report status
        print("[ STATUS ]: Iteration [%02d] of [%02d]" % ((i + 1), args.n))
        
        # set up the output temporary file name for each mixture iteration
        temp_outfile = args.outpath + aei.params.pathsep + str('temp_mixture_%02d.tif' % (i+1))
        args.temp_mixtures.append(temp_outfile)
        
        # perform the unmixing
        aei.cmd.otb.HyperspectralUnmixing(args.infile, args.temp_mixtures[i],
          args.temp_libs[i], ua = args.ua)
          
        # update variables for averaging
        for j in range(n_libs):
            bm_args['b%s' % j].append('im%sb%s' % (i+1, j+1))
        
    # report completion of unmixing
    print("[ STATUS ]: Completed unmixing iterations!")
    print("[ STATUS ]: Calculating average mixture")
    
    # calculate average and stdev for the bands
    for i in range(n_libs):
        
        # report status
        print("[ STATUS ]: Iteration [%02d] of [%02d]" % ((i + 1), n_libs))
        
        # define the output files
        temp_avgfile = args.outpath + aei.params.pathsep + str('temp_average_%02d.tif' % (i+1))
        temp_sdvfile = args.outpath + aei.params.pathsep + str('temp_stdev_%02d.tif' % (i+1))
        
        args.temp_averages.append(temp_avgfile)
        args.temp_stdevs.append(temp_sdvfile)
        
        # define the expression to compute
        avgexp = '"avg(' + aei.fn.strJoin(bm_args['b%s' % i], ', ') + ')"'
        avg_otb = "im%sb%s" % (args.n + 1, 1)
        
        # its a little more complicated for stdev since stdev isn't supported
        sdvsub = []
        for j in range(args.n):
            sdvsub.append('(' + bm_args['b%s' % i][j] + '-' + avg_otb + ')^2')
        sdvexp = '"sqrt(avg(' + aei.fn.strJoin(sdvsub, ', ') + '2))"'
        
        # run the commands
        aei.cmd.otb.BandMath(args.temp_mixtures, args.temp_averages[i], avgexp)
        aei.cmd.otb.BandMath([aei.fn.strJoin(args.temp_mixtures), args.temp_averages[i]], 
          args.temp_stdevs[i], sdvexp)
    
    # report completion of averaging
    print("[ STATUS ]: Completed averaging mixture bundles!")
    print("[ STATUS ]: Stacking into single output file")
    
    temp_vrtfile = args.outpath + aei.params.pathsep + args.temp_vrt
    
    vrt_input = [aei.fn.strJoin(args.temp_averages), aei.fn.strJoin(args.temp_stdevs)]
    
    # use gdalbuildvrt and gdal_translate to build stack
    aei.cmd.gdalbuildvrt(aei.fn.strJoin(vrt_input), temp_vrtfile, separate=True)
    aei.cmd.gdal_translate(temp_vrtfile, args.outfile, etc=['-co', 'COMPRESS=LZW'])
    
    # clean up temporary files
    print("[ STATUS ]: Removing temporary files")
    
    # remove the spectra libraries and intermediate mixtures
    #for i in range(args.n):
    #    os.remove(args.temp_libs[i])
    #    os.remove(args.temp_mixtures[i])
    
    # remove the averaged files
    #for i in range(n_libs):
    #    os.remove(args.temp_averages[i])
    #    os.remove(args.temp_stdevs[i])
    
    # remove the vrt file
    #os.remove(temp_vrtfile)
    
    # report finished
    print("[ STATUS ]: Completed aei-unmix.py!")
    print("[ STAUTS ]: Please see output file: %s" % args.outfile)

if __name__ == "__main__":
    main()