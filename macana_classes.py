import numpy as np

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse, Rectangle
from mpl_toolkits.mplot3d import Axes3D
import copy

import netCDF4
import glob
import re
import sys

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
        obs53701 = macana_plotter()
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
    
    def __init__(self, name = None):
        self.name = name
        print('Welcome to macana_plotter!')

        #Some useful variables/controls
        self.rad_to_arcsec = 3600*180/np.pi
        self.close_fig = False
        
        self.param_names = ['Amplitude', 'xoffset', 'yoffset', 'Azimuth FWHM', 
            'Elevation FWHM', 'Amplitude error', 'xoffset error', 'yoffset error', 'Azimuth FWHM error', 
            'Elevation FWHM error']
        
        print('''------------------------------------------------''')
    
    def save_setup(self, beam_loc = None, gauss_loc = None, resid_loc = None,
                         corner_loc = None, mosaic_loc = None, 
                         array_plot_loc = None):
                         
        if type(beam_loc) == str:
            #if 'beam_loc' not in dir(self):
                    self.beam_loc = beam_loc
        if type(gauss_loc) == str:
            #if 'gauss_loc' not in dir(self):
                self.gauss_loc = gauss_loc
        if type(resid_loc) == str:
            #if 'resid_loc' not in dir(self):
                self.resid_loc = resid_loc
        if type(corner_loc) == str:
            #if 'corner_loc' not in dir(self):
                self.corner_loc = corner_loc
        if type(mosaic_loc) == str:
            #if 'mosaic_loc' not in dir(self):
                self.mosaic_loc = mosaic_loc
        if type(array_plot_loc) == str:
            #if 'array_plot_loc' not in dir(self):
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
            self.boloname = self.nc_variables['beammapSignal' + str(d)].getncattr('bolo_name')
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
        
        
    def array_click(self, event, param_array):
        ix, iy = event.xdata, event.ydata
        dist2 = (ix - param_array['xoffset'])**2 + (iy - param_array['yoffset'])**2
        index = np.where(dist2 == min(dist2))
        if min(dist2) <10:
            print('This is detector %s ' % (param_array['Bolo Name'][index[0][0]]))
            #print('This is detector %d ' % (index[0]))
    
    def plot_array(self, color_param, saveon = False, save_name = None):
        print('Making a plot of the array')
        param_array = {}
        param_array['Bolo Name'] = []
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
        f.canvas.mpl_connect('button_press_event', lambda event: self.array_click(event, param_array))
        #ax.patch.set_facecolor('gray')
        plt.grid()
        
        for i in range(self.workingDetectors):
            #self.beammap(i,'Signal')
            d = self.get_gauss_params(i)
            param_array['Bolo Name'].append(self.boloname)
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
                
        cmap = plt.cm.plasma
        norm = mpl.colors.Normalize(vmin=param_array[color_param].min(), vmax=param_array[color_param].max())  
        for i in range(self.workingDetectors):
            ax.add_artist(Ellipse((param_array['xoffset'][i], param_array['yoffset'][i]), 
            height = param_array['azfwhm'][i], width = param_array['elfwhm'][i], color=cmap(norm(param_array[color_param][i]))))
        ax.set_aspect('equal')
        scale = 0
        ax.set_xlim([self.extent[0]+scale, self.extent[1]-scale])
        ax.set_ylim([self.extent[2]+scale, self.extent[3]-scale])
        
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
            plt.xlabel(self.param_names[i])
            plt.ylabel('N')
        
            if saveon == True:
                if 'array_plot_loc' in dir(self):
                    if type(save_name) == str:
                        plt.savefig(self.hist_plot_loc + save_name)
                        print(' >Figure saved to ' + self.array_plot_loc + save_name + '_' + self.param_names[i])
                    else:
                        print('    >invalid filename')
                else:
                    print('    >No file location given.  Cannot save figure.')
                if self.close_fig == True:
                    plt.close()
        print('''------------------------------------------------''')


class beammap_analyzer:
    """
    This code is an early work in progress.

    This code takes a list of beammap nc files and pulls the fit values and 
    errors from each detector for each observation and stores it in a python 
    dictionary.  For now, this code plots the values and errors against the 
    observation number.  Each observaton point can be clicked on for any of the 
    figures, which will open new figures of the fit values and errors plotted
    against the detector number for that observation.  This let's you find what 
    detector in a particular observation is bad.
    
    The following methods are included:
        
        get_files:
            This takes as input a path to a directory containing nc files and
            get's the name of all the files within that directory.
        
        load:
            This method uses macana_plotter to get the fit values and errors for
            each of the nc files stored using 'load' and stores them into a
            python dictionary.
        
        plot_obs:
            This plots the fit values and errors (10 figures) against the
            observation number (arbitrary based on your input list of files).
            There are two options for plotting currently specifying the 
            'fig_type' input to 'plot_obs:'
                
                fig_type = scatter: draw points as a scatter plot, where each 
                point is a detector.
                
                fig_type = bar: make bars for each observation where the height 
                is the min or max of the fit value or error.
            
            If you click on an observation point, it will display the name of
            the nc file in the list of files.
            
            If you double click, this will open two new figures of the fit 
            values and errors plotted against the detector number for that
            observation only:
                
                Clicking on a detector gives that detector number (and name 
                soon).
                
                Double clicking will call macana_plotter's beammap method and
                make a beammap for that detector.
                
                Right clicking (two finger on Mac) will call macana_plotter's
                plot_array function.  The colormap will be for picked based on
                the subplot clicked on.
    """
    
    def __init__(self):
        print('Welcome to beammap_analyzer!')
        
        #Number of detectors is hardcoded for now.
        self.ndetectors = 113
        
        #These are the "reasonable" limits for the fit errors.  If the fit errors are larger, the observation number and detector will
        #output in the loop below.
        self.amp_lim = 0.01
        self.xoffset_lim = 100.0
        self.yoffset_lim = 100.0
        self.azfwhm_lim = 10.0
        self.elfwhm_lim = 10.0
        
        self.amp_err_lim = 0.0001
        self.xoffset_err_lim = 1.0
        self.yoffset_err_lim = 1.0
        self.azfwhm_err_lim = 1.0
        self.elfwhm_err_lim = 1.0

        #Names for plotting.
        self.param_names = ['Amplitude', 'xoffset', 'yoffset', 'azfwhm', 
        'elfwhm', 'Amplitude Error', 'xoffset Error', 'yoffset Error',
              'azfwhm Error', 'elfwhm Error']
        
        self.nparams = len(self.param_names)
        
        self.selected = False
        self.found = False
    
    def setup_detector_figs(self, fig_type):
        #Values
        f1 = plt.figure()
        #plt.title('Fit Values vs Detector number for observation %i' % (int(round(ix))))
        
        ax1 = plt.subplot2grid(shape=(2,6), loc=(0,0), colspan=2)
        ax1.aname = self.param_names[0]
        
        ax2 = plt.subplot2grid((2,6), (0,2), colspan=2)
        ax2.aname = self.param_names[1]
        
        ax3 = plt.subplot2grid((2,6), (0,4), colspan=2)
        ax3.aname = self.param_names[2]
        
        ax4 = plt.subplot2grid((2,6), (1,1), colspan=2)
        ax4.aname = self.param_names[3]
        
        ax5 = plt.subplot2grid((2,6), (1,3), colspan=2)
        ax5.aname = self.param_names[4]
        
        #Errors
        f2 = plt.figure()
        #plt.title('Fit Errors vs Detector number for observation %i' % (int(ix)))
        
        ax6 = plt.subplot2grid(shape=(2,6), loc=(0,0), colspan=2)
        ax6.aname = self.param_names[5]
        
        ax7 = plt.subplot2grid((2,6), (0,2), colspan=2)
        ax7.aname = self.param_names[6]
        
        ax8 = plt.subplot2grid((2,6), (0,4), colspan=2)
        ax8.aname = self.param_names[7]
        
        ax9 = plt.subplot2grid((2,6), (1,1), colspan=2)
        ax9.aname = self.param_names[8]
        
        ax10 = plt.subplot2grid((2,6), (1,3), colspan=2)
        ax10.aname = self.param_names[9]
        
        if fig_type == 'obs_line':
            ax1.set_ylabel('Amplitude')
            ax1.set_xlabel('Detector number')
            
            ax2.set_ylabel('xoffset')
            ax2.set_xlabel('Detector number')
            
            ax3.set_ylabel('yoffset')
            ax3.set_xlabel('Detector number')
            
            ax4.set_ylabel('azfwhm')
            ax4.set_xlabel('Detector number')
            
            ax5.set_ylabel('elfwhm')
            ax5.set_xlabel('Detector number')
            
            ax6.set_ylabel('Amplitude error')
            ax6.set_xlabel('Detector number')
            
            ax7.set_ylabel('xoffset error')
            ax7.set_xlabel('Detector number')
            
            ax8.set_ylabel('yoffset error')
            ax8.set_xlabel('Detector number')
            
            ax9.set_ylabel('azfwhm error')
            ax9.set_xlabel('Detector number')
            
            ax10.set_ylabel('elfwhm error')
            ax10.set_xlabel('Detector number')
            
        if fig_type == 'detec_line':
            ax1.set_ylabel('Amplitude')
            ax1.set_xlabel('Observation number')
            
            ax2.set_ylabel('xoffset')
            ax2.set_xlabel('Observation number')
            
            ax3.set_ylabel('yoffset')
            ax3.set_xlabel('Observation number')
            
            ax4.set_ylabel('azfwhm')
            ax4.set_xlabel('Observation number')
            
            ax5.set_ylabel('elfwhm')
            ax5.set_xlabel('Observation number')
            
            ax6.set_ylabel('Amplitude error')
            ax6.set_xlabel('Observation number')
            
            ax7.set_ylabel('xoffset error')
            ax7.set_xlabel('Observation number')
            
            ax8.set_ylabel('yoffset error')
            ax8.set_xlabel('Observation number')
            
            ax9.set_ylabel('azfwhm error')
            ax9.set_xlabel('Observation number')
            
            ax10.set_ylabel('elfwhm error')
            ax10.set_xlabel('Observation number')
            
        if fig_type == 'hist':
            ax1.set_ylabel('N')
            ax1.set_xlabel('Amplitude')
            
            ax2.set_ylabel('N')
            ax2.set_xlabel('xoffset')
            
            ax3.set_ylabel('N')
            ax3.set_xlabel('yoffset')
            
            ax4.set_ylabel('N')
            ax4.set_xlabel('azfwhm')
            
            ax5.set_ylabel('N')
            ax5.set_xlabel('elfwhm')
            
            ax6.set_ylabel('N')
            ax6.set_xlabel('Amplitude Error')
            
            ax7.set_ylabel('N')
            ax7.set_xlabel('xoffset Error')
            
            ax8.set_ylabel('N')
            ax8.set_xlabel('yoffset Error')
            
            ax9.set_ylabel('N')
            ax9.set_xlabel('azfwhm error')
            
            ax10.set_ylabel('N')
            ax10.set_xlabel('elfwhm error')
    
        return f1, f2, ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9, ax10
        
    def obs_click(self, event):
        global ix, iy
        ix, iy = event.xdata, event.ydata
        
        print('You selected observation %i which corresponds to %s' %(ix, self.beammaps[int(round(ix))]))
        
        if event.button == 1:
            if event.dblclick:
                global plot_type
                plot_type = 'obs_line'
                f1, f2, ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9, ax10 = self.setup_detector_figs(fig_type=plot_type)
                f1.canvas.mpl_connect('button_press_event', self.detector_click)
                f2.canvas.mpl_connect('button_press_event', self.detector_click)
                
                #Values 
                ax1.plot(range(self.ndetectors), self.plot_array['Amplitude'][:,int(ix)], c='k')
                ax2.plot(range(self.ndetectors), self.plot_array['xoffset'][:,int(ix)], c='k')
                ax3.plot(range(self.ndetectors), self.plot_array['yoffset'][:,int(ix)], c='k')
                ax4.plot(range(self.ndetectors), self.plot_array['azfwhm'][:,int(ix)], c='k')
                ax5.plot(range(self.ndetectors), self.plot_array['elfwhm'][:,int(ix)], c='k')
                
                f1.tight_layout()
            
                #Errors            
                ax6.plot(range(self.ndetectors), self.plot_array['Amplitude Error'][:,int(ix)], c='k')
                ax7.plot(range(self.ndetectors), self.plot_array['xoffset Error'][:,int(ix)], c='k')
                ax8.plot(range(self.ndetectors), self.plot_array['yoffset Error'][:,int(ix)], c='k')
                ax9.plot(range(self.ndetectors), self.plot_array['azfwhm Error'][:,int(ix)], c='k')
                ax10.plot(range(self.ndetectors), self.plot_array['elfwhm Error'][:,int(ix)], c='k')
            
                f2.tight_layout()
        
        if event.button == 3:
            if event.dblclick:
                global plot_type
                plot_type = 'hist'
                f1, f2, ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9, ax10 = self.setup_detector_figs(fig_type=plot_type)
                f1.canvas.mpl_connect('button_press_event', self.hist_click)
                f2.canvas.mpl_connect('button_press_event', self.hist_click)
                
                #Values
                indices = np.where(self.plot_array['Amplitude'] != np.inf)
                n1, b1, _ = ax1.hist(self.plot_array['Amplitude'][indices[0],int(ix)], bins = 25, color='k', histtype = 'step')
                
                indices = np.where(self.plot_array['xoffset'] != np.inf)
                n2, b2, _ = ax2.hist(self.plot_array['xoffset'][indices[0],int(ix)], bins = 25, color='k', histtype = 'step')
                
                indices = np.where(self.plot_array['yoffset'] != np.inf)
                n3, b3, _ = ax3.hist(self.plot_array['yoffset'][indices[0],int(ix)], bins = 25, color='k', histtype = 'step')
                
                indices = np.where(self.plot_array['azfwhm'] != np.inf)
                n4, b4, _ = ax4.hist(self.plot_array['azfwhm'][indices[0],int(ix)], bins = 25, color='k', histtype = 'step')
                
                indices = np.where(self.plot_array['elfwhm'] != np.inf)
                n5, b5, _ = ax5.hist(self.plot_array['elfwhm'][indices[0],int(ix)], bins = 25, color='k', histtype = 'step')
    
                f1.tight_layout()
                
                #Errors            
                indices = np.where(self.plot_array['Amplitude Error'] != np.inf)
                n6, b6, _ = ax6.hist(self.plot_array['Amplitude Error'][indices[0],int(ix)], bins = 25, color='k', histtype = 'step')
                
                indices = np.where(self.plot_array['xoffset Error'] != np.inf)
                n7, b7, _ = ax7.hist(self.plot_array['xoffset Error'][indices[0],int(ix)], bins = 25, color='k', histtype = 'step')
                
                indices = np.where(self.plot_array['yoffset Error'] != np.inf)
                n8, b8, _ = ax8.hist(self.plot_array['yoffset Error'][indices[0],int(ix)], bins = 25, color='k', histtype = 'step')
                
                indices = np.where(self.plot_array['azfwhm Error'] != np.inf)
                n9, b9, _ = ax9.hist(self.plot_array['azfwhm Error'][indices[0],int(ix)], bins = 25, color='k', histtype = 'step')
                
                indices = np.where(self.plot_array['elfwhm Error'] != np.inf)
                n10, b10, _ = ax10.hist(self.plot_array['elfwhm Error'][indices[0],int(ix)], bins = 25, color='k', histtype = 'step')
            
                f2.tight_layout()
                
                global bins_dict 
                
                bins_dict = {}
                bins_dict[self.param_names[0]] = b1
                bins_dict[self.param_names[1]] = b2
                bins_dict[self.param_names[2]] = b3
                bins_dict[self.param_names[3]] = b4
                bins_dict[self.param_names[4]] = b5
                bins_dict[self.param_names[5]] = b6
                bins_dict[self.param_names[6]] = b7
                bins_dict[self.param_names[7]] = b8
                bins_dict[self.param_names[8]] = b9
                bins_dict[self.param_names[9]] = b10
        
        plt.show()
    
    def detec_click(self, event):
        global ix, iy
        ix, iy = event.xdata, event.ydata
        
        print('You selected detector %i' %(ix))
        
        if event.button == 1:
            if event.dblclick:
                global plot_type
                plot_type = 'detec_line'
                f1, f2, ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9, ax10 = self.setup_detector_figs(fig_type=plot_type)
                f1.canvas.mpl_connect('button_press_event', self.observation_click)
                f2.canvas.mpl_connect('button_press_event', self.observation_click)
                
                #Values 
                ax1.plot(range(self.beam_num), self.plot_array['Amplitude'][int(ix),:], c='k')
                ax2.plot(range(self.beam_num), self.plot_array['xoffset'][int(ix),:], c='k')
                ax3.plot(range(self.beam_num), self.plot_array['yoffset'][int(ix),:], c='k')
                ax4.plot(range(self.beam_num), self.plot_array['azfwhm'][int(ix),:], c='k')
                ax5.plot(range(self.beam_num), self.plot_array['elfwhm'][int(ix),:], c='k')
                
                f1.tight_layout()
            
                #Errors            
                ax6.plot(range(self.beam_num), self.plot_array['Amplitude Error'][int(ix),:], c='k')
                ax7.plot(range(self.beam_num), self.plot_array['xoffset Error'][int(ix),:], c='k')
                ax8.plot(range(self.beam_num), self.plot_array['yoffset Error'][int(ix),:], c='k')
                ax9.plot(range(self.beam_num), self.plot_array['azfwhm Error'][int(ix),:], c='k')
                ax10.plot(range(self.beam_num), self.plot_array['elfwhm Error'][int(ix),:], c='k')
            
                f2.tight_layout()
        
        if event.button == 3:
            if event.dblclick:
                global plot_type
                plot_type = 'hist'
                f1, f2, ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9, ax10 = self.setup_detector_figs(fig_type=plot_type)
                f1.canvas.mpl_connect('button_press_event', self.hist_click)
                f2.canvas.mpl_connect('button_press_event', self.hist_click)
                
                #Values 
                indices = np.where(self.plot_array['Amplitude'][int(ix),:] != np.inf)
                n1, b1, _ = ax1.hist(self.plot_array['Amplitude'][int(ix),indices[0]], bins = 25, color='k', histtype = 'step')
                
                indices = np.where(self.plot_array['xoffset'][int(ix),:] != np.inf)
                n2, b2, _ = ax2.hist(self.plot_array['xoffset'][int(ix),indices[0]], bins = 25, color='k', histtype = 'step')
                
                indices = np.where(self.plot_array['yoffset'][int(ix),:] != np.inf)
                n3, b3, _ = ax3.hist(self.plot_array['yoffset'][int(ix),indices[0]], bins = 25, color='k', histtype = 'step')
                
                indices = np.where(self.plot_array['azfwhm'][int(ix),:] != np.inf)
                n4, b4, _ = ax4.hist(self.plot_array['azfwhm'][int(ix),indices[0]], bins = 25, color='k', histtype = 'step')
                
                indices = np.where(self.plot_array['elfwhm'][int(ix),:] != np.inf)
                n5, b5, _ = ax5.hist(self.plot_array['elfwhm'][int(ix),indices[0]], bins = 25, color='k', histtype = 'step')
    
                f1.tight_layout()
                
                #Errors            
                indices = np.where(self.plot_array['Amplitude Error'][int(ix),:] != np.inf)
                n6, b6, _ = ax6.hist(self.plot_array['Amplitude Error'][int(ix),indices[0]], bins = 25, color='k', histtype = 'step')
                
                indices = np.where(self.plot_array['xoffset Error'][int(ix),:] != np.inf)
                n7, b7, _ = ax7.hist(self.plot_array['xoffset Error'][int(ix),indices[0]], bins = 25, color='k', histtype = 'step')
                
                indices = np.where(self.plot_array['yoffset Error'][int(ix),:] != np.inf)
                n8, b8, _ = ax8.hist(self.plot_array['yoffset Error'][int(ix),indices[0]], bins = 25, color='k', histtype = 'step')
                
                indices = np.where(self.plot_array['azfwhm Error'][int(ix),:] != np.inf)
                n9, b9, _ = ax9.hist(self.plot_array['azfwhm Error'][int(ix),indices[0]], bins = 25, color='k', histtype = 'step')
                
                indices = np.where(self.plot_array['elfwhm Error'][int(ix),:] != np.inf)
                n10, b10, _ = ax10.hist(self.plot_array['elfwhm Error'][int(ix),indices[0]], bins = 25, color='k', histtype = 'step')
            
                f2.tight_layout()

                global bins_dict 
                
                bins_dict = {}
                bins_dict[self.param_names[0]] = b1
                bins_dict[self.param_names[1]] = b2
                bins_dict[self.param_names[2]] = b3
                bins_dict[self.param_names[3]] = b4
                bins_dict[self.param_names[4]] = b5
                bins_dict[self.param_names[5]] = b6
                bins_dict[self.param_names[6]] = b7
                bins_dict[self.param_names[7]] = b8
                bins_dict[self.param_names[8]] = b9
                bins_dict[self.param_names[9]] = b10
                    
        plt.show()
        
    def hist_click(self,event):
        global ix3, iy3, iinaxes3
        ix3, iy3, iinaxes3 = event.xdata, event.ydata, event.inaxes
        
        print('''------------------------------------------------''')
        print('You selected %s = %f' % (iinaxes3.aname, ix3))

        if 'ix' in globals():
            if event.dblclick:
                rows, bin_width, lower_bin, upper_bin = self.find_on_hist()
                rows = rows[0]
                if self.selected:
                    self.R.remove()
                self.R = iinaxes3.add_artist(Rectangle((lower_bin,0),width=bin_width,
                    height=len(rows), color='r'))
                self.selected = True
                plt.draw()
                print('Beammap: %s ' % (self.beammaps[int(round(ix))]))
                for i in range(len(rows)):
                    print('   >Detector: %s ' % (rows[i]))
        else:
            if event.dblclick:
                rows, cols, bin_width, lower_bin, upper_bin = self.find_on_hist()
                if self.selected:
                    self.R.remove()
                self.R = iinaxes3.add_artist(Rectangle((lower_bin,0),width=bin_width,
                    height=len(rows), color='r'))
                self.selected = True
                plt.draw()
                cols, rows = zip(*sorted(zip(cols, rows)))
                for i in range(len(rows)):
                    print('''------------------------------------------------''')
                    print('Beammap: %s ' % (self.beammaps[cols[i]]))
                    print('   >Detector: %s ' % (rows[i]))
    
    def find_on_hist(self):
        bin_edges = bins_dict[iinaxes3.aname]
        bin_width = bin_edges[1] - bin_edges[0]
        index = np.where(abs(bin_edges - ix3) == min(abs(bin_edges - ix3)))
        closest_bin = bin_edges[index[0]]
        if closest_bin > ix3:
            upper_bin = closest_bin
            lower_bin = bin_edges[index[0]-1]
            
        elif closest_bin < ix3:
            upper_bin = bin_edges[index[0]+1]
            lower_bin = closest_bin
                
        bin_mid = np.mean([lower_bin, upper_bin])
                
        if 'ix' in globals():
            rows = np.where(abs(self.plot_array[iinaxes3.aname][:,int(round(ix))] - bin_mid) < bin_width/2.0)
            return rows, bin_width, lower_bin, upper_bin
        else:
            rows, cols = np.where(abs(self.plot_array[iinaxes3.aname] - bin_mid) < bin_width/2.0)
            return rows, cols, bin_width, lower_bin, upper_bin
    
    def imshow_click(self,event):
        global ix4, iy4
        ix4, iy4 = event.xdata, event.ydata
        
        print('You selected observation %i which corresponds to %s' %(iy4, self.beammaps[int(round(iy4))]))
        print('This is detector %i ' % (int(round(ix4))))

        
    
    #For the figures produced by "obs_click", this outputs the nearest detector
    #number to the clicked point.
    def detector_click(self, event):
        global ix2, iy2, iinaxes2
        ix2, iy2, iinaxes2 = event.xdata, event.ydata, event.inaxes
        print('This is detector %i ' % (int(round(ix2))))
        
        if event.button == 1:
            if event.dblclick:
                obs = macana_plotter(str(ix))
                obs.load_nc(self.beammaps[int(ix)])
                obs.beammap(int(ix2), 'Signal', plotting=True)
            
        if event.button == 3:
            obs = macana_plotter(str(ix))
            obs.load_nc(self.beammaps[int(ix)])
            obs.plot_array(iinaxes2.aname)
    
    def observation_click(self, event):
        global ix2, iy2, iinaxes2
        ix2, iy2, iinaxes2 = event.xdata, event.ydata, event.inaxes
        print('You selected observation %i which corresponds to %s' %(ix2, self.beammaps[int(round(ix2))]))
    
    def get_files(self,path):
        self.beammaps = glob.glob(path + '/*')
        self.beammaps.sort()
        self.beam_num = len(self.beammaps)
    
    def load(self, npy_file = None, find_bad = False):
            if type(npy_file) == str:
                self.param_array = np.load(npy_file, encoding = 'latin1').item()
                self.beam_num = len(self.param_array['Amplitude'][0,:])

            else:
                #Setting up the dictionary.  Dimensions are number of detectors by number of obseravations.
                self.param_array = {}
                self.param_array['Amplitude'] = np.zeros([self.ndetectors, self.beam_num])
                self.param_array['xoffset'] = np.zeros([self.ndetectors, self.beam_num])
                self.param_array['yoffset'] = np.zeros([self.ndetectors, self.beam_num])
                self.param_array['azfwhm'] = np.zeros([self.ndetectors, self.beam_num])
                self.param_array['elfwhm'] = np.zeros([self.ndetectors, self.beam_num])
                
                self.param_array['Amplitude Error'] = np.zeros([self.ndetectors, self.beam_num])
                self.param_array['xoffset Error'] = np.zeros([self.ndetectors, self.beam_num])
                self.param_array['yoffset Error'] = np.zeros([self.ndetectors, self.beam_num])
                self.param_array['azfwhm Error'] = np.zeros([self.ndetectors, self.beam_num])
                self.param_array['elfwhm Error'] = np.zeros([self.ndetectors, self.beam_num])
                
                #Loop through the observations putting the fit values and errors for each detector into the dictionary.
                for j in range(self.beam_num):
                    obs = macana_plotter(j)
                    obs.load_nc(self.beammaps[j])
                    for i in range(self.ndetectors):
                        obs.get_gauss_params(detector = i)
                        
                        self.param_array['Amplitude'][i,j] = obs.amp
                        self.param_array['xoffset'][i,j] = obs.xoffset
                        self.param_array['yoffset'][i,j] = obs.yoffset
                        self.param_array['azfwhm'][i,j] = obs.azfwhm
                        self.param_array['elfwhm'][i,j] = obs.elfwhm
                        
                        self.param_array['Amplitude Error'][i,j] = obs.amp_err
                        self.param_array['xoffset Error'][i,j] = obs.xoffset_err
                        self.param_array['yoffset Error'][i,j] = obs.yoffset_err
                        self.param_array['azfwhm Error'][i,j] = obs.azfwhm_err
                        self.param_array['elfwhm Error'][i,j] = obs.elfwhm_err
        
    def finder(self, values = True, errors = True, exclude = 'bad', recursive = False, print_bad = False):
        print('''------------------------------------------------''')
        print('Don''t forget to set the limits!  Current limits are set to: \n')
        print('%s %.3f' % ('Amplitude limit:', self.amp_lim))
        print('%s %7.3f' % ('Xoffset limit:', self.xoffset_lim))
        print('%s %7.3f' % ('Yoffset limit:', self.yoffset_lim))
        print('%s %8.3f' % ('Azfwhm limit:', self.azfwhm_lim))
        print('%s %8.3f \n' % ('Elfwhm limit:', self.elfwhm_lim))
        
        print('%s %.3f' % ('Amplitude error limit:', self.amp_err_lim))
        print('%s %7.3f' % ('Xoffset error limit:', self.xoffset_err_lim))
        print('%s %7.3f' % ('Yoffset error limit:', self.yoffset_err_lim))
        print('%s %8.3f' % ('Azfwhm error limit:', self.azfwhm_err_lim))
        print('%s %8.3f' % ('Elfwhm error limit:', self.elfwhm_err_lim))
        
        lim_array = np.array([self.amp_lim, self.xoffset_lim, self.yoffset_lim, 
                              self.azfwhm_lim, self.elfwhm_lim, self.amp_err_lim, 
                              self.xoffset_err_lim, self.yoffset_err_lim,
                              self.azfwhm_err_lim, self.elfwhm_err_lim])
            
        if recursive == False:
            self.bad_array = copy.deepcopy(self.param_array)
            if 'bad_array' not in dir(self):
                print('cannot do recursive. Have you run finder once yet?')
        
        for i in range(self.nparams):
            if exclude == 'bad':
                rows, cols = np.where(abs(self.param_array[self.param_names[i]]) > lim_array[i])
            elif exclude == 'good':
                rows, cols = np.where(abs(self.param_array[self.param_names[i]]) <= lim_array[i])
            for j in range(len(rows)):
                self.bad_array[self.param_names[i]][rows[j], cols[j]] = np.inf
                if print_bad == True:
                    print('Beammap: %s ' % (self.beammaps[cols[j]]))
                    print('   >Detector: %s ' % (rows[j]))
        
        return self.bad_array           
                    
    def excluder(self,exclude):
        if exclude == 'None':
            return self.param_array
            
        elif exclude == 'bad':
            return self.finder(values = True, errors = True, exclude = 'bad', print_bad = False)      
        
        elif exclude == 'good':
            return self.finder(values = True, errors = True, exclude = 'good', print_bad = False)      
     
    def plot_obs(self, fig_type = 'obs', plot_types = 'scatter', exclude = None, bins = 25):
        #Plots the fit values and errors against the observation number.  
        #Each figure is slightly interactive in that one can click on an 
        #observation to open up additional figures for that observation.  See intro.
        
        cmap = plt.cm.plasma
        
        self.plot_array = self.excluder(exclude)
                
        if fig_type == 'obs':
            norm = mpl.colors.Normalize(0.0, self.ndetectors) 
            for j in range(len(self.param_names)):
                f = plt.figure(j)
                f.canvas.mpl_connect('button_press_event', self.obs_click)
                if plot_types == '3d':
                    ax = f.add_subplot(111, projection='3d')
                else:
                    ax = f.add_subplot(111)
                
                if plot_types == 'scatter':
                    plt.grid(True)
                    for i in range(self.ndetectors):
                        ax.scatter(range(self.beam_num), self.plot_array[self.param_names[j]][i,:], s=20,
                        marker = 's', color=cmap(norm(i)))
                        
                elif plot_types == 'bar':
                    max_bars = np.zeros(self.beam_num)
                    min_bars = np.zeros(self.beam_num)
                    for i in range(self.beam_num):
                        max_bars[i] = max(self.plot_array[self.param_names[j]][:,i])
                        min_bars[i] = min(self.plot_array[self.param_names[j]][:,i])
                    ax.bar(range(self.beam_num), height=max_bars)
                    ax.bar(range(self.beam_num), height=min_bars)
                    
                elif plot_types == '3d':
                    plt.grid(True)
                    for i in range(self.beam_num):
                        ax.scatter(range(self.beam_num), self.plot_array[self.param_names[j]][i,:], i, s=20,
                        marker = 's', color=cmap(norm(i)))
                
                ax.set_xlabel('observation')
                ax.set_ylabel(self.param_names[j])
                plt.axis('tight')
        
        elif fig_type == 'detec':
            norm = mpl.colors.Normalize(0.0, self.beam_num) 
            
            for j in range(len(self.param_names)):
                f = plt.figure(j)
                f.canvas.mpl_connect('button_press_event', self.detec_click)
                if plot_types == '3d':
                    ax = f.add_subplot(111, projection='3d')
                else:
                    ax = f.add_subplot(111)
                    
                if plot_types == 'scatter':
                    plt.grid(True)
                    for i in range(self.beam_num):
                        ax.scatter(range(self.ndetectors), self.plot_array[self.param_names[j]][:,i], s=20,
                        marker = 's', color=cmap(norm(i)))
                    
                elif plot_types == 'bar':
                    max_bars = np.zeros(self.ndetectors)
                    min_bars = np.zeros(self.ndetectors)
                    for i in range(self.ndetectors):
                        max_bars[i] = max(self.plot_array[self.param_names[j]][i,:])
                        min_bars[i] = min(self.plot_array[self.param_names[j]][i,:])
                    ax.bar(range(self.ndetectors), height=max_bars)
                    ax.bar(range(self.ndetectors), height=min_bars)
                
                elif plot_types == '3d':
                    plt.grid(True)
                    for i in range(self.beam_num):
                        ax.scatter(range(self.ndetectors), self.plot_array[self.param_names[j]][:,i], zs = i, s=20,
                        marker = 's', color=cmap(norm(i)))
                        
                ax.set_xlabel('Detector Number')
                ax.set_ylabel(self.param_names[j])
                plt.axis('tight')
        
        elif fig_type == 'hist':
            params_ravel = {}
            for i in range(self.nparams):
                params_ravel[self.param_names[i]] = np.ravel(self.plot_array[self.param_names[i]])
                params_ravel[self.param_names[i]] = params_ravel[self.param_names[i]][np.where(params_ravel[self.param_names[i]]!=np.inf)]
            f1, f2, ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9, ax10 = self.setup_detector_figs(fig_type=fig_type)
            f1.canvas.mpl_connect('button_press_event', self.hist_click)
            f2.canvas.mpl_connect('button_press_event', self.hist_click)
            
            #Values 
            n1, b1, _ = ax1.hist(params_ravel['Amplitude'], bins = bins, color='k', histtype = 'step')
            n2, b2, _ = ax2.hist(params_ravel['xoffset'], bins = bins, color='k', histtype = 'step')
            n3, b3, _ = ax3.hist(params_ravel['yoffset'], bins = bins, color='k', histtype = 'step')
            n4, b4, _ = ax4.hist(params_ravel['azfwhm'], bins = bins, color='k', histtype = 'step')
            n5, b5, _ = ax5.hist(params_ravel['elfwhm'], bins = bins, color='k', histtype = 'step')

            f1.tight_layout()
            
            #Errors            
            n6, b6, _ = ax6.hist(params_ravel['Amplitude Error'], bins = bins, color='k', histtype = 'step')
            n7, b7, _ = ax7.hist(params_ravel['xoffset Error'], bins = bins, color='k', histtype = 'step')
            n8, b8, _ = ax8.hist(params_ravel['yoffset Error'], bins = bins, color='k', histtype = 'step')
            n9, b9, _ = ax9.hist(params_ravel['azfwhm Error'], bins = bins, color='k', histtype = 'step')
            n10, b10, _ = ax10.hist(params_ravel['elfwhm Error'], bins = bins, color='k', histtype = 'step')
        
            f2.tight_layout()

            global bins_dict 
            bins_dict = {}
            bins_dict[self.param_names[0]] = b1
            bins_dict[self.param_names[1]] = b2
            bins_dict[self.param_names[2]] = b3
            bins_dict[self.param_names[3]] = b4
            bins_dict[self.param_names[4]] = b5
            bins_dict[self.param_names[5]] = b6
            bins_dict[self.param_names[6]] = b7
            bins_dict[self.param_names[7]] = b8
            bins_dict[self.param_names[8]] = b9
            bins_dict[self.param_names[9]] = b10
        
        elif fig_type == 'imshow':            
            for i in range(self.nparams):
                f = plt.figure(i)
                f.canvas.mpl_connect('button_press_event', self.imshow_click)
                ax = f.add_subplot(111)
                img = ax.imshow(self.plot_array[self.param_names[i]], origin = 'lower')
                ax.set_xlabel('Observation')
                ax.set_ylabel('Detector')
                f.colorbar(img, ax = ax, label = self.param_names[i])
            
        plt.show()

plt.close('all')
beam = beammap_analyzer()
beam.get_files('/Users/quirkyneku/Documents/TolTEC-Project/beammaps_out/')
beam.load()#'/Users/quirkyneku/Documents/macana_python_plotting/fit_array.npy')
beam.finder(print_bad = False)
beam.plot_obs(fig_type = 'imshow', plot_types = 'scatter', bins=50, exclude = 'bad')

plt.show()