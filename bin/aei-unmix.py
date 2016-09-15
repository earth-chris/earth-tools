#####
# aei_unmix.py performs unconstrained least-squares 
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
        self.bands = []
        self.subset_bands = False
        self.of = 'GTiff'
        self.normalize = False
        self.temp_libs = []
        self.temp_mixtures = []
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
                    if not aei.checkFile(self.infile, quiet = True):
                        usage()
                        aei.checkFile(self.infile)
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
                    if not aei.checkFile(libs[j], quiet=True):
                        usage()
                        aei.checkFile(libs[j])
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
                
            # set up catch-all for incorrect parameter call
            else:
                usage()
                print("[ ERROR ]: Unrecognized argument: %s" % arg)
                print("%s" % i)
                aei.params.sys.exit(1)
            
            i += 1

# set up function to write the temporary endmember libraries for otbcli
def write_temp_libraries(args,lib,n_libs):
    
    # loop through each library, write temp files, and store info for later deletion
    for i in range(0, n_libs):
        
        # set up the output file name
        temp_file = args.outpath + aei.params.pathsep + str('temp_lib_%02d.tif' % i)
        args.temp_libs.append(temp_file)
        
        # create array to hold the endmembers. if shade is included, add an
        #  additional endmember
        if args.include_shade:
            bundles = np.zeros((n_libs + 1, len(args.bands)))
            nl = n_libs + 1
        else:
            bundles = np.zeros((n_libs, len(args.bands)))
            nl = n_libs
        
        # load random spectra from each lib into the output array
        for j in range(n_libs):
            bundles[j,:] = lib['lib_%s' % j].spectra[lib['rnd_%s' % j][i],args.bands]
            
        # open the temp file
        file_ref = gdal.GetDriverByName("GTiff").Create(temp_file, 1, nl, 
          len(args.bands), gdal.GDT_Float32)
        
        # write the arrays band by band  
        for j in range(len(args.bands)):
            band = file_ref.GetRasterBand(j + 1)
            band.WriteArray(np.expand_dims(bundles[:,j], 1))
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
      [-rescale rescale_parameter] [-include_shade]
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
    lnb = np.zeros(n_libs)
    
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
    write_temp_libraries(args,lib,n_libs)
    
    # load the input image and get parameters
    inf = gdal.Open(args.infile)
    ns = inf.RasterXSize
    nl = inf.RasterYSize
    nb = inf.RasterCount
    geo = inf.GetGeoTransform()
    prj = inf.GetProjection()
    n_bundles = lnb.shape[0]
    b1 = inf.GetRasterBand(1)
    
    # if bands were not set in command line, use all bands
    if not args.subset_bands:
        args.bands = range(nb)
    
    # check that the spectral libraries match the bands
    if not 0 in np.where((lnb-nb) != 0)[0].shape:
        print("[ ERROR ]: number of image bands does not match number of spectral library bands")
        inf = None
        sys.exit(1)
    
    # report a little
    print("[ STATUS ]: Beginning aei-unmix.py")
    print("[ STATUS ]: Input file : %s" % args.infile)
    print("[ STATUS ]: Output file: %s" % args.outfile)
        
    # normalize the data if set
    #if args.normalize:
    #    img = aei.bn(img, inds = args.bands)
    #    args.bands = range(len(args.bands))
    
    b1_arg = []
    b2_arg = []
    b3_arg = []
        
    # loop through each random index and unmix
    for i in range(args.n):
        
        # report status
        print("[ STATUS ]: Iteration [%02d] of [%02d]" % ((i + 1), args.n))
        
        # set up the output temporary file name for each mixture iteration
        temp_outfile = args.outpath + aei.params.pathsep + str('temp_mixture_%02d.tif' % i)
        args.temp_mixtures.append(temp_outfile)
        
        # perform the unmixing
        aei.cmd.otb.HyperspectralUnmixing(args.infile, args.temp_mixtures[i],
          args.temp_libs[i], ua = 'ucls')
          
        # update variables for averaging
        b1_arg.append('im%sb1' % i)
        b2_arg.append('im%sb2' % i)
        b3_arg.append('im%sb3' % i)
        
    # report completion of unmixing
    print("[ STATUS ]: Completed unmixing iterations")
    print("[ STATUS ]: Calculating average mixture")
    
    # average the bands

if __name__ == "__main__":
    main()