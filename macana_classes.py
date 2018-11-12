import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import netCDF4
from matplotlib.patches import Ellipse

class macana_plotter:
    """Macana Plotter:
    
    This python class includes a number of methods for various visual 
    representations of AzTEC beammap data produced by the macana c++ data 
    reduction pipeline.  The purpose of this code is to provide quick, but robust
    ways to ascertain the properties and qualtiy of the beammmaps.
    
    The following python modules are required:
        numpy
        matplotlib
        netCDF4
        
    The following methods are currently implemented:
        
        load_nc: 
            Load an nc file into memory, storing it in the variable 'nc.'
            The nc variable names are stored in the variable 'nc_variables.'
            The number of detectors found in the nc file is stored in
            'workingDetectors.'
        
        save_setup: 
            Specify the path to where figures will be saved if 'saveon' is True. The
            current options include:
                beam_loc
                gauss_loc 
                resid_loc
                corner_loc 
                mosaic_loc
                array_plot_loc
        
        beammap: 
            Get's the beammap Signal or Weight array for the specified detector
            and stores it in memory.  Can plot if 'plotting' is True.  If saveon is True
            the image is saved to the specified location.  A filename can be specified
            by the save_name parameter.
        
        get_gauss_params: 
            Get's the fit parameters for a given beammap and stores 
            them in memory as separate variables (i.e. .Amplitude, etc).  Defaults to 
            stored detector if none provided.
        
        make_gauss: 
            Creates a model of the Gaussian specified by the fit paramters
            for a given detector.  Defaults to the stored detector if none provided.
            The dimensions of the model array match those of the beammap. Can plot if 
            'plotting' is True.  If saveon is True the image is saved to the specified 
            location.  A filename can be specified by the save_name parameter.
        
        gauss_resid: 
            Subtracts the the Gaussian model from the Signal beammap.  
            Defaults to the stored detector if none provided.  Will call get_gauss_params
            and make_gauss if not called previously.  Can plot if  'plotting' is True.  
            If saveon is True the image is saved to the specified location.  A filename 
            can be specified by the save_name parameter.
        
        plot_corner:  
            Creates a 'corner plot' for the given detector, which includes
            the beammap and 1 dimensional distributions along each axis (azimuth and 
            elevation).  Defaults to the stored detector if none provided.  If saveon is
            True the image is saved to the specified location.  A filename can be s
            pecified by the save_name parameter.
        
        make_mosaic_image:  
            Generates a mosaic of the beammaps (signal or weights)
            for all detectors.  Calls beammap iteratively.  If saveon is True the image 
            is saved to the specified location.  A filename can be specified by the 
            save_name parameter.
        
        plot_detectors:  
            Creates a plot where each detector is represented by an
            ellipse whose axes are the FWHM of the Gaussian fit parameters.  Each
            ellipse is centered at the xoffset and yoffset fit values.  A fit variable
            can be specified when plot_detectors is called and will be used as the
            colormap of the ellipses.  Calls get_gauss_params iteratively.  If saveon is
            True the image is saved to the specified location.  A filename can be 
            specified by the save_name parameter.
        
        fit_hists: 
            Makes histograms of the fit parameters and their errors using all
            detectors.  Calls get_gauss_params iteratively.  If saveon is 
            True the images are saved to the specified location.  A filename can be 
            specified by the save_name parameter.  The parameter name is appended to the
            end of the input file name to prevent overwriting.
        
    Some methods are only used internally.  These are:
        
        get_rows_cols:
            Gets the rows and columns scales and lengths for the nc file.
            
        zoom:
            Handles the optional region parameter (see below).
        
        get_1D_dist:
            Used by plot_corner to calculate the 1-dimensional count distribution
            along each axis.  This are found by integrating (trapezoidal method)
            along one axis, for each point in the other axis.
        
    
    Some additional info:
        All axes are in arcseconds.
        
        The methods 'beammap,' 'make_gauss,' 'gauss_resid,' 'plot_corner,' and 
        'make_mosaic_image' have an optional 'region' parameter that can be
        specified when the methods are called.  If 'an int or float is specified for
        region,' the beammap will be a map of the square of length 2xregion centered
        at the xoffset (azimuth) and yoffset (elevation) of the Gaussian fit.
        
    Below is an example of each method run in sequence:
        n = 99
        obs53701 = macana_plotter('53701')
        obs53701.load_nc('/Users/quirkyneku/Downloads/Toltec/53701pixel.nc')
        obs53701.save_setup(beam_loc = '/Users/quirkyneku/Downloads/Toltec/')
        obs53701.beammap(n, 'Signal', plotting = True, saveon = True, save_name = 'test-beammap')
        obs53701.make_gauss(plotting = True, region = None)
        obs53701.gauss_resid(plotting = True, region = None)
        obs53701.plot_corner()
        obs53701.make_mosaic_image('Signal')
        obs53701.plot_array('Amplitude')
        
    This code is a work-in-progress and will be updated frequently with new 
    methods for different visual representations and changes to improve 
    and usability add the support for TolTEC beammap plotting."""
    
    def __init__(self, name):
        self.name = name
        print('Welcome to macana_plotter!')

        #Some useful variables/controls
        self.rad_to_arcsec = 3600*180/np.pi
        self.close_fig = False
        
        print('''------------------------------------------------''')
    
    def save_setup(self, beam_loc = None, gauss_loc = None, resid_loc = None,
                         corner_loc = None, mosaic_loc = None, 
                         array_plot_loc = None):
                         
        if type(beam_loc) == str:
            if 'beam_loc' not in dir(self):
                    self.beam_loc = beam_loc
        if type(gauss_loc) == str:
            if 'gauss_loc' not in dir(self):
                self.gauss_loc = gauss_loc
        if type(resid_loc) == str:
            if 'resid_loc' not in dir(self):
                self.resid_loc = resid_loc
        if type(corner_loc) == str:
            if 'corner_loc' not in dir(self):
                self.corner_loc = corner_loc
        if type(mosaic_loc) == str:
            if 'mosaic_loc' not in dir(self):
                self.mosaic_loc = mosaic_loc
        if type(array_plot_loc) == str:
            if 'array_plot_loc' not in dir(self):
                self.array_plot_loc = array_plot_loc
    
    def load_nc(self, ncfile):
        print('Loading nc file %s' % (ncfile))
        self.nc = netCDF4.Dataset(ncfile)
        self.nc_methods = dir(self.nc)

        if 'variables' in self.nc_methods:
            self.nc_variables = self.nc.variables
            self.get_rows_cols()
            
            #if 'boloData' in self.nc_variables:
            self.workingDetectors = 113#self.nc_variables['boloData'].shape[1]
        else:
            print('    Cannot find variables')
                
        print('''------------------------------------------------''')
        
    def get_rows_cols(self):
        self.row = self.nc.variables['rowCoordsPhys'][:]*self.rad_to_arcsec
        self.col = self.nc.variables['colCoordsPhys'][:]*self.rad_to_arcsec
        
        self.nrows = len(self.row)
        self.ncols = len(self.col)
        
        self.extent = [self.row.min(), self.row.max(), self.col.min(), self.col.max()]

    def zoom(self, region):
        if type(region) == int or type(region) == float:
            print('Zooming into the %i x %i region around the source' % (region, region))
            if 'gauss' not in dir(self):
                print('    >No gaussian found')
                print('    >Making gaussian for detector %i' % (self.detector))
                self.make_gauss(self.detector)
            
            xlower = self.xoffset - region
            xupper = self.xoffset + region
            ylower = self.yoffset - region
            yupper = self.yoffset + region

            self.xlimits = [xlower, xupper]
            self.ylimits = [ylower, yupper]
            
            region = None
    
    def beammap(self, detector, maptype = None, region=None, plotting = False, saveon = False, save_name = None):
        self.maptype = maptype
        self.detector = detector
        print('Getting beammap ' + self.maptype + ' for detector %i' % (detector))
        self.var = 'beammap' + self.maptype + str(detector)
        self.matrix = np.array(self.nc.variables[self.var])
        
        for i in range(self.nrows):
            for j in range(self.ncols):
                if self.matrix[i,j] != self.matrix[i,j]:
                    self.matrix[i,j]  = 0.0
        self.matrix_rot = np.rot90(self.matrix[:])
                
        if plotting == True:
            print('    >Plotting beammap ' + self.maptype + ' for detector %i' % (detector))
            
            self.zoom(region)
            
            plt.figure()
            plt.imshow(self.matrix_rot, extent=self.extent, cmap = 'hot')
            if type(region) == int or type(region) == float:
                plt.xlim(self.xlimits)
                plt.ylim(self.ylimits)
            plt.xlabel('Azimuth (arcseconds)')
            plt.ylabel('Elevation (arcseconds)')
            plt.title('Beammap' + self.maptype + ' for detector %i' % (detector))
            plt.colorbar()
            if saveon == True:
                if 'beam_loc' in dir(self):
                    if type(save_name) == str:
                        plt.savefig(self.beam_loc + save_name)
                        print('    >Figure saved to ' + self.beam_loc + save_name + '.png')
                    else:
                        print('    >invalid filename')
                else:
                    print('    >No file location given.  Cannot save figure.')
                if self.close_fig == True:
                    plt.close()
        #print('''------------------------------------------------''')

    def get_gauss_params(self, detector = None):
        if detector == None:
            if 'detector' in dir(self):
                print('yeah')
                print('Getting fit parameters for detector %i' % (self.detector))
                d = self.detector
            else:
                print('    >No beammap specified and no beammap stored')
        else:
            if type(detector) == int:
                print('Getting fit parameters for detector %i' % (detector))
                d = detector
        if 'nc_variables' in dir(self):
            self.amp = self.nc_variables['beammapSignal' + str(d)].getncattr('amplitude')
            self.azfwhm = self.nc_variables['beammapSignal' + str(d)].getncattr('FWHM_x')
            self.elfwhm = self.nc_variables['beammapSignal' + str(d)].getncattr('FWHM_y')
            self.xoffset = self.nc_variables['beammapSignal' + str(d)].getncattr('offset_x')
            self.yoffset = self.nc_variables['beammapSignal' + str(d)].getncattr('offset_y')
            
            self.amp_err = self.nc_variables['beammapSignal' + str(d)].getncattr('amplitude_err')
            self.azfwhm_err = self.nc_variables['beammapSignal' + str(d)].getncattr('FWHM_x_err')
            self.elfwhm_err = self.nc_variables['beammapSignal' + str(d)].getncattr('FWHM_y_err')
            self.xoffset_err = self.nc_variables['beammapSignal' + str(d)].getncattr('offset_x_err')
            self.yoffset_err = self.nc_variables['beammapSignal' + str(d)].getncattr('offset_y_err')
        
        return d
            
    def make_gauss(self, region=None, plotting = False, saveon = False, save_name = None):
        detector = self.get_gauss_params()
        print('Making map of gauss fit parameters for detector %i' % (detector))

        self.gauss = np.zeros([self.nrows, self.ncols])
                    
        for i in range(self.nrows):
            for j in range(self.ncols):
                self.gauss[i,j] = self.amp*np.exp(-4*np.log(2) * ((self.col[j]-self.yoffset)**2/self.elfwhm**2 + (self.row[i] - self.xoffset)**2/self.azfwhm**2))
        self.gauss = np.rot90(self.gauss)
        
        if plotting == True:
            self.zoom(region)
            
            plt.figure()
            plt.imshow(self.gauss, extent = self.extent, cmap = 'hot')
            if type(region) == int or type(region) == float:
                plt.xlim(self.xlimits)
                plt.ylim(self.ylimits)
            plt.xlabel('Azimuth (arcseconds)')
            plt.ylabel('Elevation (arcseconds)')
            plt.title('Gaussian Fit Reconstruction for detector %i' % (detector))
            plt.colorbar()
            if saveon == True:
                if 'gauss_loc' in dir(self):
                    if type(save_name) == str:
                        plt.savefig(self.gauss_loc + save_name)
                        print(' >Figure saved to ' + self.gauss_loc + save_name)
                    else:
                        print('    >invalid filename')
                else:
                    print('    >No file location given.  Cannot save figure.')
                if self.close_fig == True:
                    plt.close()
        print('''------------------------------------------------''')
    

    def gauss_resid(self, region=None, plotting = False, saveon = False, save_name = None):
        print('Calculating gauss fit residual for detector %i' % (self.detector))
        if self.maptype == 'Signal':
            if 'matrix' in dir(self):
                if 'gauss' not in dir(self):
                    print('    >No gaussian found')
                    print('    >Making gaussian for detector %i' % (self.detector))
                    d = self.get_gauss_params(self.detector)
                    self.make_gauss(self.detector)
                self.resid = self.matrix_rot - self.gauss
                self.chi2 = np.sum(np.sqrt(self.resid**2.))
                print('    >Residual chi squared: %f: ' % (self.chi2))
                    
                if plotting == True:
                    self.zoom(region)
                    
                    plt.figure()
                    plt.imshow(self.resid, extent = self.extent, cmap = 'hot')
                    if type(region) == int or type(region) == float:
                        plt.xlim(self.xlimits)
                        plt.ylim(self.ylimits)
                    plt.xlabel('Azimuth (arcseconds)')
                    plt.ylabel('Elevation (arcseconds)')
                    plt.xlabel('Azimuth (arcseconds)')
                    plt.ylabel('Elevation (arcseconds)')
                    plt.title('Residual Map for detector %i' % (self.detector))
                    plt.colorbar()
                    if saveon == True:
                        if 'resid_loc' in dir(self):
                            if type(save_name) == str:
                                plt.savefig(self.resid_loc + save_name)
                                print(' >Figure saved to ' + self.resid_loc + save_name)
                            else:
                                print('    >invalid filename')
                        else:
                            print('    >No file location given.  Cannot save figure.')
                        if self.close_fig == True:
                            plt.close()
        else:
            print('    >This is not a signal beammap')
        print('''------------------------------------------------''')
        
    def get_1D_dist(self):
        #Integration for corner plot.
        el_integ = np.zeros(self.nrows)
        az_integ = np.zeros(self.ncols)
        
        for i in range(self.nrows):
            el_integ[i] = np.trapz(np.array(self.matrix[i,:]))
        
        for j in range(self.ncols):
            az_integ[j] = np.trapz(np.array(self.matrix[:,j]))
    
        return az_integ, el_integ
        
    def plot_corner(self, bins = 25, region=None, saveon = False, save_name = None):
        print('Making a corner plot for detector %i' % (self.detector))
        
        az_integ, el_integ = self.get_1D_dist()

        f, a = plt.subplots(2,2, gridspec_kw = {'width_ratios':[1, 1], 'height_ratios':[1, 1]}, figsize = (7,7))
        a[0,0].step(self.row, el_integ, 'k')
        a[0,0].xaxis.set_ticks_position('none')
        a[0,0].axvline(self.xoffset, color="g",  linewidth=1)
        a[0,0].set_xticklabels([])
        
        a[1,0].contourf(self.row, -self.col, self.matrix_rot,cmap = 'hot', origin = 'upper')
        a[1,0].axvline(self.xoffset, color="g",  linewidth=1)
        a[1,0].axhline(self.yoffset, color="g",  linewidth=1)

        a[1,0].set_xlabel('Azimuth')
        a[1,0].set_ylabel('Elevation')
        
        a[1,1].step(az_integ, self.col,'k')
        a[1,1].yaxis.set_ticks_position('none')
        a[1,1].axhline(self.yoffset, color="g",  linewidth=1)
        a[1,1].set_yticklabels([])
        
        a[0,0].axis('off')
        a[0,0].axis('tight')
        a[0,1].axis('off')
        a[1,1].axis('off')

        f.subplots_adjust(wspace=0, hspace=0)
        
        if saveon == True:
            if 'corner_loc' in dir(self):
                if type(save_name) == str:
                    plt.savefig(self.corner_loc + save_name)
                    print(' >Figure saved to ' + self.corner_loc + save_name)
                else:
                    print('    >invalid filename')
            else:
                print('    >No file location given.  Cannot save figure.')
            if self.close_fig == True:
                plt.close()        
        
        print('''------------------------------------------------''')

    def make_mosaic_image(self, maptype = None, region = None, saveon = False, save_name = None):
        print('Making a mosaic plot of all detectors')
        tmp_region = region
        fig = plt.figure(figsize=(8,8))
        nplots = int(np.sqrt(self.workingDetectors) + 1)
        ax = [fig.add_subplot(nplots,nplots,i+1) for i in range(nplots**2)]
        
        fig.subplots_adjust(wspace=0, hspace=0)

        for a in ax:
            a.set_xticklabels([])
            a.set_yticklabels([])
            a.set_aspect('equal')
        
        for i in range(nplots**2):
            if i < self.workingDetectors:
                self.beammap(i,maptype)
                self.zoom(region)
                region = tmp_region
                ax[i].imshow(self.matrix_rot, extent=self.extent,cmap='hot')
                if type(region) == int or type(region) == float:
                    ax[i].set_xlim(self.xlimits)
                    ax[i].set_ylim(self.ylimits) 
                    del self.gauss
       
            else:
                ax[i].axis('off')
    
            fig.text(0.5, 0.1, 'Azimuth', ha='center', va='center')
            fig.text(0.1, 0.5, 'Elevation', ha='center', va='center', rotation='vertical')
            ax[int(nplots/2)].set_title('Signal Beammaps for AzTEC Detectors')
            
            if saveon == True:
                if 'mosaic_loc' in dir(self):
                    if type(save_name) == str:
                        plt.savefig(self.mosaic_loc + save_name)
                        print(' >Figure saved to ' + self.mosaic_loc + save_name)
                    else:
                        print('    >invalid filename')
                else:
                    print('    >No file location given.  Cannot save figure.')
                if self.close_fig == True:
                    plt.close()
        print('''------------------------------------------------''')
    def plot_array(self, color_param, saveon = False, save_name = None):
        print('Making a plot of the array')
        param_array = {}
        param_array['Amplitude'] = np.zeros(self.workingDetectors)
        param_array['xoffset'] = np.zeros(self.workingDetectors)
        param_array['yoffset'] = np.zeros(self.workingDetectors)
        param_array['azfwhm'] = np.zeros(self.workingDetectors)
        param_array['elfwhm'] = np.zeros(self.workingDetectors)
        
        param_array['Amplitude Error'] = np.zeros(self.workingDetectors)
        param_array['xoffset Error'] = np.zeros(self.workingDetectors)
        param_array['yoffset Error'] = np.zeros(self.workingDetectors)
        param_array['azfwhm Error'] = np.zeros(self.workingDetectors)
        param_array['elfwhm Error'] = np.zeros(self.workingDetectors)
        
        f, ax = plt.subplots()        
        for i in range(self.workingDetectors):
            #self.beammap(i,'Signal')
            d = self.get_gauss_params(i)
            param_array['Amplitude'][i] = self.amp
            param_array['xoffset'][i] = self.xoffset
            param_array['yoffset'][i] = self.yoffset
            param_array['azfwhm'][i] = self.azfwhm
            param_array['elfwhm'][i] = self.elfwhm

            param_array['Amplitude Error'][i] = self.amp_err
            param_array['xoffset Error'][i] = self.xoffset_err
            param_array['yoffset Error'][i] = self.yoffset_err
            param_array['azfwhm Error'][i] = self.azfwhm_err
            param_array['elfwhm Error'][i] = self.elfwhm_err
                
        cmap = plt.cm.hot
        norm = matplotlib.colors.Normalize(vmin=param_array[color_param].min(), vmax=param_array[color_param].max())  
        for i in range(self.workingDetectors):
            ax.add_artist(Ellipse((param_array['xoffset'][i], param_array['yoffset'][i]), 
            height = param_array['azfwhm'][i], width = param_array['elfwhm'][i], color=cmap(norm(param_array[color_param][i]))))
        ax.set_aspect('equal')
        ax.set_xlim([self.extent[0], self.extent[1]])
        ax.set_ylim([self.extent[2], self.extent[3]])
        
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        plt.colorbar(sm, label=color_param)
        plt.xlabel('Azimuth')
        plt.ylabel('Elevation')
        
        if saveon == True:
            if 'array_plot_loc' in dir(self):
                if type(save_name) == str:
                    plt.savefig(self.array_plot_loc + save_name)
                    print(' >Figure saved to ' + self.array_plot_loc + save_name)
                else:
                    print('    >invalid filename')
            else:
                print('    >No file location given.  Cannot save figure.')
            if self.close_fig == True:
                plt.close()
        print('''------------------------------------------------''')
        
    def fit_hists(self, bins = 25, saveon = None, save_name = False):
        print('Plotting histograms of all fit parameters')
        param_array = np.zeros([self.workingDetectors, 10])
        param_names = ['Amplitude', 'xoffset', 'yoffset', 'Azimuth FWHM', 
                    'Elevation FWHM', 'Amplitude error', 'xoffset error', 'yoffset error', 'Azimuth FWHM error', 
                    'Elevation FWHM error']

        for i in range(self.workingDetectors):
            #self.beammap(i,'Signal')
            self.get_gauss_params(i)
            param_array[i,0] = self.amp
            param_array[i,1] = self.xoffset
            param_array[i,2] = self.yoffset
            param_array[i,3] = self.azfwhm
            param_array[i,4] = self.elfwhm
            
            param_array[i,5] = self.amp_err
            param_array[i,6] = self.xoffset_err
            param_array[i,7] = self.yoffset_err
            param_array[i,8] = self.azfwhm_err
            param_array[i,9] = self.elfwhm_err
                    
        for i in range(len(param_array[0,:])):
            plt.figure()
            plt.hist(param_array[:,i], bins = bins, histtype = 'stepfilled', facecolor = 'w',
                     edgecolor = 'k')
            plt.xlabel(param_names[i])
            plt.ylabel('N')
        
            if saveon == True:
                if 'array_plot_loc' in dir(self):
                    if type(save_name) == str:
                        plt.savefig(self.hist_plot_loc + save_name)
                        print(' >Figure saved to ' + self.array_plot_loc + save_name + '_' + param_names[i])
                    else:
                        print('    >invalid filename')
                else:
                    print('    >No file location given.  Cannot save figure.')
                if self.close_fig == True:
                    plt.close()
        print('''------------------------------------------------''')
