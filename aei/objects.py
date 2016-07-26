#
# contains the classes and functions of frequently used objects
#  to support aei
###################

# set up class for spectral objects
class spectral:
    
    def __init__(self, n_spectra=1, n_wl=2151, type='', band_unit = '',
            band_quantity = '', band_centers = []):
        
        import numpy as np
        
        # set to asd type if no params set to change n_wl
        if n_wl == 2151:
            type = 'asd'
        
        # set up pre-defined types
        if (str(type).lower()) == 'asd':
            n_wl = 2151
            band_unit = 'Nanometers'
            band_quantity = 'Wavelength'
            band_centers = np.arange(350.,2501.)
        
        # return a list same size as number of spectra
        names = []
        for i in range(n_spectra):
            names.append('Spectrum ' + str(i))
        self.names = names
        
        # set up the band definitions
        try:
            self.band_unit = band_unit
        except NameError:
            pass
        try:
            self.band_quantity = band_quantity
        except NameError:
            pass
        try:
            self.band_centers = band_centers
        except NameError:
            pass
        
        # return an np array size of n spectra x n wavelengths
        self.spectra = np.zeros([n_spectra, n_wl])
        
    def remove_water_bands(self):
        """ 
        a function to set data from water absorption bands,
        i.e. 1350 - 1460 nm and 1790 - 1960 nm, to zero
        
        usage: self.remove_water_bands()
        """
        import numpy as np
        
        water_bands = [[1350.0, 1460.0], [1790.0, 1960.0]]
        
        # start with nir-swir1 transition
        gt = np.where(self.band_centers > water_bands[0][0])
        lt = np.where(self.band_centers < water_bands[0][1])
        nd = np.intersect1d(gt[0], lt[0])
        self.spectra[:,nd] = 0.0
        
        # then swir1-swir2 transition
        gt = np.where(self.band_centers > water_bands[1][0])
        lt = np.where(self.band_centers < water_bands[1][1])
        nd = np.intersect1d(gt[0], lt[0])
        self.spectra[:,nd] = 0.0
    
    def get_shortwave_bands(self, bands=[]):
        """
        a function that returns an index of the bands that encompass
        the shortwave range (350 - 2500 nm)
        
        usage: self.get_shortwave_bands(bands=[])
        
        returns: an index of bands to subset to the shortwave range
        """
        import numpy as np
        
        # set range to return in nanometers
        shortwave_range = [350., 2500.]
        
        # normalize if wavelength units are different
        if self.band_unit == 'Micrometers':
            shortwave_range /= 1000.
            
        # find overlapping range
        gt = np.where(self.band_centers > shortwave_range[0])
        lt = np.where(self.band_centers < shortwave_range[1])
        overlap = np.intersect1d(gt[0], lt[0])
        
        # return output
        return overlap
        
    def plot(self, inds = [], legend = False):
        """
        plots the spectra using a standard plot format
        can be set to only plot selected spectra
        
        usage: self.plot(inds = [], legend = False)
          where inds = optional 0-based indices for spectra to plot
                legend = set this to force a legend to be created
        """
        # import pyplot
        import matplotlib.pyplot as plt
        
        # set basic parameters
        plt.xlim((self.band_centers.min(), self.band_centers.max()))
        plt.xlabel('Wavelength (' + self.band_unit + ')')
        plt.ylabel('Reflectance (%)')
        plt.title('spectralObject plot')
        
        # check if indices were set and valid. if not, plot all items
        if inds:
            if max(inds) > len(self.names):
                inds = range(0, len(self.names))
                print("[ ERROR ]: invalid range set. using all spectra")
            if min(inds) < 0:
                inds = range(0, len(self.names))
                print("[ ERROR ]: invalid range set. using all spectra")
        else:
            inds = range(0, len(self.names))
            
        # turn on the legend if fewer than 10 entries
        if len(inds) < 10:
            legend = True
            
        # loop through each item to plot
        for i in inds:
            plt.plot(self.band_centers, self.spectra[i,:], 
                label = self.names[i])
            
        # add the legend with each spectrum's name
        if legend:
            plt.legend(fontsize = 'small', framealpha = 0.5, 
                fancybox = True)
        
        # display the plot
        plt.show()
        
    def bn(self, inds = []):
        """
        brightness normalizes the spectra
        
        usage: self.bn(inds = [])
          where inds = the indices to use for BN
        """
        # check if indices were set and valid. if not, plot all items
        if inds:
            if max(inds) > self.spectra.shape[-1]:
                inds = range(0, self.spectra.shape[-1])
                print("[ ERROR ]: invalid range set. using all spectra")
            if min(inds) < 0:
                inds = range(0, self.spectra.shape[-1])
                print("[ ERROR ]: invalid range set. using all spectra")
        else:
            inds = range(0, self.spectra.shape[-1])
            
        # perform the bn
        self.spectra = self.spectra[:,inds] / np.expand_dims(
            np.sqrt((self.spectra[:,inds]**2).sum(1)),1)
        
        # subset band centers to the indices selected, if they exist
        if self.band_centers.ndim != 0:
            self.band_centers = self.band_centers[inds]
  
# class for las/laz data  
class las:
    
    import numpy as np
    
    def __init__(self, n_pts=1):
    
        # the stuff
        self.x = np.zeros(n_pts)
        self.y = np.zeros(n_pts)
        self.z = np.zeros(n_pts)