#!/usr/bin/python
#####
# aei_sample-raster.py 
#  takes an input training raster data set (e.g. a classification map)
#  and extracts a number of random samples from each class.
#  the default behavior is to treat 0 as the class to randomly sample,
#  and sample all values of 1 or higher. the training stack should explicitly
#  contain a no data value if anything is to be ignored; otherwise it is
#  included in the sampling.
#
# the default output is a csv file, so titles can be added to the predictor
#  columns, if necessary, and data can be manipulated more easily
#
# c. 2016 Christopher Anderson
#####

import os
import sys
import aei
import random
import gdal as gdal
import numpy as np

class parse_args:
    def __init__(self, arglist):
        
        # set up main variables and defaults to parse
        self.trainingFile = ''
        self.predictorFiles = []
        self.outfile = ''
        self.nSamples = 10000
        self.balance = False
        self.includeXY = False
        self.trainingNoData = False

        # exit if no arguments passed
        if len(arglist) == 1:
            usage(exit=True)
    
        # read arguments from command line
        i = 1
        while i < len(arglist):
            arg = arglist[i]
        
            # check predictor data paths            
            if arg.lower() == '-predictors':
                i += 1
                arg = arglist[i]
                #libs = arg.split(" ")
                # loop through predictor args until we find a new parameter
                newArg = False
                while not newArg:
                    if arg[0] == '-':
                        newArg = True
                        i -= 1
                        continue
                    if not aei.fn.checkFile(arg, quiet=True):
                        usage()
                        aei.fn.checkFile(arg)
                        aei.params.sys.exit()
                        
                    self.predictorFiles.append(arg)
                    
                    # increment counter and move on
                    i += 1
                    if i >= len(arglist): 
                        newArg = True
                        i -= 1
                        continue
                    arg = arglist[i]
                    
            # parse the training flag    
            elif arg.lower() == '-training':
                i += 1
                arg = arglist[i]
                
                self.trainingFile = arg
                if not aei.fn.checkFile(self.trainingFile, quiet = True):
                    usage()
                    aei.fn.checkFile(self.trainingFile)
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
            
            # check number of iterations
            elif arg.lower() == "-n":
                i += 1
                arg = arglist[i]
                
                try:
                    self.nSamples = int(arg)
                except ValueError:
                    usage()
                    print("[ ERROR ]: -n argument is not an integer: %s" % arg)
                    aei.params.sys.exit(1)
                    
            elif arg.lower() == '-balance':
                self.balance = True
            
            elif arg.lower() == '-include_xy':
                self.includeXY = True
        
            # set up catch-all for incorrect parameter call
            else:
                usage()
                print("[ ERROR ]: Unrecognized argument: %s" % arg)
                print("%s" % i)
                aei.params.sys.exit(1)
            
            i += 1
            
def usage(exit=False):
    """
    describes the aei-sample-raster.py procedure in case of incorrect parameter calls
    
    syntax: usage()
    """
    print(
        """
$ aei-sample-raster.py -training training_raster -predictors 
      predictor1, predictor2 ... predictorx -o output_file
      [-n n_random_samples] [-balance] [-include_xy]
        """
        )
    if exit:
        sys.exit(1)

def main():
    """
    the main program for aei-sample-raster.py
    
    the order of operations for this procedure is to
    1) read in a training data raster file
    2) find the unique class numbers after accounting for no-data vals
    3) generate random samples from each of the classes 
       (based on balanced/unbalanced)
    4) close the training data file
    5) loop through the predictor raster data sets
    6) loop through each band in the predictor data sets and 
       extract the same random samples from each class
    7) output these data to a csv file
    8) profit
    
    syntax: main()
    """
    # parse the argument list
    args = parse_args(sys.argv)
    
    # report starting
    print("[ STATUS ]: Starting aei-sample-raster.py")
    print("[ STATUS ]: Training stack : %s" % args.trainingFile)
    print("[ STATUS ]: Random samples : %s" % args.nSamples)
    print("[ STATUS ]: Balanced samples : %s" % args.balance)
    print("[ STATUS ]: ----------")
    print("[ STATUS ]: Reading training data")
    
    # open the training data set and open the first band
    tdRef = gdal.Open(args.trainingFile)
    tdb1 = tdRef.GetRasterBand(1)
    
    # get metadata from the training file
    tdns = tdRef.RasterXSize
    tdnl = tdRef.RasterYSize
    tdulx, tdxps, tdxoff, tduly, tdyoff, tdyps = tdRef.GetGeoTransform()
    tdnd = tdb1.GetNoDataValue()
    
    # if there is a nodata value for the training band, be sure to ignore it
    if tdnd is not None:
        args.trainingNoData = True
        
    # read the training data into memory
    td = tdb1.ReadAsArray()
    
    # remove no data if it exists
    if args.trainingNoData:
        nd = np.where(td != tdnd)
        td = td[nd[0], nd[1]]
    else:
        nd = [range(tdns), range(tdnl)]
        td = td[nd[0], nd[1]]
    
    # get the number of classes to loop through
    startClass = td.min()
    endClass = td.max()
    
    # spit out an error if there are only 0s
    if startClass == endClass:
        print("[ ERROR ]: Training data raster has only one value: %s" % startClass)
        print("[ ERROR ]: Please check the input file: %s" % args.trainingFile)
    
    # compute a histogram to get the number of samples for each class
    tdHisto = np.histogram(td, range(startClass, endClass + 2))
    
    # calculate the number of classes with good data
    classNumbers = np.where(tdHisto[0] > 0)[0]
    nClasses = classNumbers.shape[0]
    
    # set a variable to store which class has the fewest samples, 
    #  a 2-element array for [class index, n_samples]
    nSamples = np.zeros((nClasses, 2), np.int64)
    nSamples[:,0] = classNumbers
    nSamples[:,1] = tdHisto[0][classNumbers]
    
    # now we have the number of samples for each class, we can randomly subset them
    # set up a dictionary to add the indices for each class
    indDict = {}
    
    # if we are balancing the training data sets, set the number of samples to 
    #  the smaller number of a) the smallest training class, or b) the minimum set
    if args.balance:
        args.nSamples = min([args.nSamples, nSamples[:,1].min()])
    
    # set up a counter to measure total number of values to collect
    nTotalSamples = []
    
    # set up a counter for cumulative n samples
    nIterSamples = []
    nIterSamples.append(0)
    
    # set up a loop to run through each class and get the 
    for i in range(nClasses):
        gd = np.where(td == classNumbers[i])
        
        # get number of good data values
        ngd = gd[0].shape[0]
        
        # create random indices for this class if there are more 
        #  points than the number of random samples set. othewise, get all.
        if ngd <= args.nSamples:
            finalInds = gd
        else:
            rnd = np.random.randint(0, ngd+1, args.nSamples)
            finalInds = []
            for j in range(gd[0].ndim):
                finalInds.append(gd[j][rnd])
        
        # add the indices from this class to the dictionary
        indDict['class_%03d' % classNumbers[i]] = finalInds
        
        # and add the number of samples to the running total
        nTotalSamples.append(finalInds[0].shape[0])
        
        # add cumulative n samples
        nIterSamples.append(nIterSamples[i] + nTotalSamples[i])
        
    # now that we have our final indices for each class, close the training data set
    tdRef = None
    tdb1 = None
    td = None
    
    # report finished
    print("[ STATUS ]: Finished collecting indices for training data")
    print("[ STATUS ]: Number of classes found: %s" % nClasses)
    for i in range(nClasses):
        print("[ STATUS ]: Random samples for class %s: %s" % (classNumbers[i], 
            indDict['class_%03d' % classNumbers[i]][0].shape[0]))
    print("[ STATUS ]: ----------")
    print("[ STATUS ]: Beginning extraction of predictor data")
    
    # create an array that will grow as more predictor data is read in
    nTotalValues = sum(nTotalSamples)
    outputArray = np.zeros((nTotalValues, 1))
    outputLabels = ['Class']
    predictorCounter = 0
    
    # fill the array with the class numbers
    outputArray[0:nTotalSamples[0]] = classNumbers[0]
    for i in range(nClasses):
        outputArray[nIterSamples[i]:nIterSamples[i+1]] = classNumbers[i]
    
    # we'll loop through each predictor data set, read the data, and extract
    #  the random indices
    for i in range(len(args.predictorFiles)):
        print("[ STATUS ]: Extracting from predictor set %s: %s" % (i+1, args.predictorFiles[i]))
        
        # read the file metadata
        pdRef = gdal.Open(args.predictorFiles[i])
        
        # get the dimensions from the file
        pdns = pdRef.RasterXSize
        pdnl = pdRef.RasterYSize
        pdulx, pdxps, pdxoff, pduly, pdyoff, pdyps = pdRef.GetGeoTransform()
        
        # check that the dimensions are the same
        if (pdns != tdns) or (pdnl != tdnl):
            print("[ ERROR ]: File dimensions do not match for training and predictor data")
            print("[ ERROR ]: Training data file : %s" % args.trainingFile)
            print("[ ERROR ]: Predictor data file: %s" % args.predictorFiles[i])
            print("[ ERROR ]: Training ns, nl : %s, %s" % (tdns, tdnl))
            print("[ ERROR ]: Predictor ns, nl: %s, %s" % (pdns, pdnl))
            print("[ ERROR ]: skipping...")
            pdRef = None
            continue
        
        # now that we know the data are good, we'll loop through each band and extract the data
        for j in range(pdRef.RasterCount):
            
            # add a new label for this predictor band
            predictorCounter += 1
            outputLabels.append("Predictor %s" % predictorCounter)
            
            # add a new column to the output array
            outputArray = np.concatenate([outputArray, np.zeros((nTotalValues, 1))],1)
            
            # set up a reference to the band we use
            pdBand = pdRef.GetRasterBand(j+1)
            
            # read the data into memory
            pdData = pdBand.ReadAsArray()
            
            # subset to the extent of the training NoData mask
            pdData = pdData[nd[0], nd[1]]
            
            # extract the data from each class and assign it to the output array
            for k in range(nClasses):
                outputArray[nIterSamples[k]:nIterSamples[k+1], predictorCounter] = \
                    pdData[indDict['class_%03d' % classNumbers[k]][0]]
            
    # write the output file
    np.savetxt(args.outfile, outputArray, delimiter = ",", header = aei.fn.strJoin(outputLabels, ','))
        
    # report finished
    print("[ STATUS ]: Completed random sampling!")
    print("[ STATUS ]: See output file: %s" % args.outfile)
        
# call the aain routine when run from command lne
if __name__ == "__main__":
    main()