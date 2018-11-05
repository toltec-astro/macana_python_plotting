Macana Plotter:
    
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
            True the image is saved to the specified location.  A filename can be specified 
            by the save_name parameter.
        
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
    methods for different visual representations and changes to improve usability and to
    add the support for TolTEC beammap plotting.
