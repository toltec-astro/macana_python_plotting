'''
This code is an early work in progress.

This code takes a list of beammap nc files and pulls the fit values and errors from each detector for each observation and
stores it in a python dictionary.  For now, this code plots the values and errors against the observation number.  Each 
observaton point can be clicked on for any of the figures, which will open new figures of the fit values and errors plotted
against the detector number for that observation.  This let's you find what detector in a particular observation is bad.
'''

#Import needed modules
import numpy as np
import glob
import re
import sys
import matplotlib.pyplot as plt
import matplotlib as mpl

#This code uses macana_plotter, so this line specifies the path to it.
sys.path.insert(0, '/home/mmccrackan/macana_python_plotting/')

#We need this for opening the nc files easily.
from macana_classes import macana_plotter

plt.close('all')

'''
Function to handle clicking.  This makes two figures when a point is *double* clicked on in the value/error vs observation number
Figures (see below).  Each of the generated figures have 5 subplots of the fit values or errors plotted against the detector
number.
'''
def obs_click(event):
    if event.dblclick:
        global ix, iy
        ix, iy = event.xdata, event.ydata
        print('You selected observation %i which corresponds to %s' %(ix, beammaps[int(ix)]))
        
        #Values   
        f1 = plt.figure()
        plt.title('Fit Values vs Detector number for observation %i' % (int(ix)))
    
        cid = f1.canvas.mpl_connect('button_press_event', detector_click)
    
        ax1 = plt.subplot2grid(shape=(2,6), loc=(0,0), colspan=2)
        ax2 = plt.subplot2grid((2,6), (0,2), colspan=2)
        ax3 = plt.subplot2grid((2,6), (0,4), colspan=2)
        ax4 = plt.subplot2grid((2,6), (1,1), colspan=2)
        ax5 = plt.subplot2grid((2,6), (1,3), colspan=2)
        
        ax1.plot(range(ndetectors), param_array['Amplitude'][:,int(ix)], c='k')
        ax1.set_ylabel('Amplitude')
        ax1.set_xlabel('Detector number')
    
        ax2.plot(range(ndetectors), param_array['xoffset'][:,int(ix)], c='k')
        ax2.set_ylabel('xoffset')
        ax2.set_xlabel('Detector number')
    
        ax3.plot(range(ndetectors), param_array['yoffset'][:,int(ix)], c='k')
        ax3.set_ylabel('yoffset')
        ax3.set_xlabel('Detector number')
    
        ax4.plot(range(ndetectors), param_array['azfwhm'][:,int(ix)], c='k')
        ax4.set_ylabel('azfwhm')
        ax4.set_xlabel('Detector number')
        
        ax5.plot(range(ndetectors), param_array['elfwhm'][:,int(ix)], c='k')
        ax5.set_ylabel('elfwhm')
        ax5.set_xlabel('Detector number')
        
        plt.tight_layout()
        
        #Errors
        f2 = plt.figure()
        plt.title('Fit Errors vs Detector number for observation %i' % (int(ix)))
    
        cid = f2.canvas.mpl_connect('button_press_event', detector_click)
    
        ax6 = plt.subplot2grid(shape=(2,6), loc=(0,0), colspan=2)
        ax7 = plt.subplot2grid((2,6), (0,2), colspan=2)
        ax8 = plt.subplot2grid((2,6), (0,4), colspan=2)
        ax9 = plt.subplot2grid((2,6), (1,1), colspan=2)
        ax10 = plt.subplot2grid((2,6), (1,3), colspan=2)
        
        ax6.plot(range(ndetectors), param_array['Amplitude Error'][:,int(ix)], c='k')
        ax6.set_ylabel('Amplitude error')
        ax6.set_xlabel('Detector number')
    
        ax7.plot(range(ndetectors), param_array['xoffset Error'][:,int(ix)], c='k')
        ax7.set_ylabel('xoffset error')
        ax7.set_xlabel('Detector number')
    
        ax8.plot(range(ndetectors), param_array['yoffset Error'][:,int(ix)], c='k')
        ax8.set_ylabel('yoffset error')
        ax8.set_xlabel('Detector number')
    
        ax9.plot(range(ndetectors), param_array['azfwhm Error'][:,int(ix)], c='k')
        ax9.set_ylabel('azfwhm error')
        ax9.set_xlabel('Detector number')
        
        ax10.plot(range(ndetectors), param_array['elfwhm Error'][:,int(ix)], c='k')
        ax10.set_ylabel('elfwhm error')
        ax10.set_xlabel('Detector number')
    
        plt.tight_layout()
        plt.show()

#For the figures produced by "obs_click", this outputs the nearest detector
#number to the clicked point.
def detector_click(event):
    ix2, iy2 = event.xdata, event.ydata
    print('This is detector %i ' % (int(ix2)))
    
    if event.dblclick:
        obs = macana_plotter(str(ix))
        obs.load_nc(beammaps[int(ix)])
        obs.beammap(int(ix2), 'Signal', plotting=True)
    
#Parameters for clicking function.
nplots = 3
nparams = 5

#Path of directory including all the beammap nc files.
beammap_path = '/home/mmccrackan/Documents/TolTEC-Project/beammaps_out/'

#Gets the name of the files in the directory.
beammaps = glob.glob(beammap_path + '*')
beam_num = len(beammaps)

#Number of detectors is hardcoded for now.
ndetectors = 113

#These are the "reasonable" limits for the fit errors.  If the fit errors are larger, the observation number and detector will
#output in the loop below.
amp_lim = 0.01
xoffset_lim = 1.0
yoffset_lim = 1.0
azfwhm_lim = 1.0
elfwhm_lim = 1.0

#Setting up the dictionary.  Dimensions are number of detectors by number of obseravations.
param_array = {}
param_array['Amplitude'] = np.zeros([ndetectors, beam_num])
param_array['xoffset'] = np.zeros([ndetectors, beam_num])
param_array['yoffset'] = np.zeros([ndetectors, beam_num])
param_array['azfwhm'] = np.zeros([ndetectors, beam_num])
param_array['elfwhm'] = np.zeros([ndetectors, beam_num])

param_array['Amplitude Error'] = np.zeros([ndetectors, beam_num])
param_array['xoffset Error'] = np.zeros([ndetectors, beam_num])
param_array['yoffset Error'] = np.zeros([ndetectors, beam_num])
param_array['azfwhm Error'] = np.zeros([ndetectors, beam_num])
param_array['elfwhm Error'] = np.zeros([ndetectors, beam_num])

#Loop through the observations putting the fit values and errors for each detector into the dictionary.
for j in range(beam_num):
    obs = macana_plotter(j)
    obs.load_nc(beammaps[j])
    for i in range(ndetectors):
        obs.get_gauss_params(detector = i)
        
        param_array['Amplitude'][i,j] = obs.amp
        param_array['xoffset'][i,j] = obs.xoffset
        param_array['yoffset'][i,j] = obs.yoffset
        param_array['azfwhm'][i,j] = obs.azfwhm
        param_array['elfwhm'][i,j] = obs.elfwhm
        
        param_array['Amplitude Error'][i,j] = obs.amp_err
        param_array['xoffset Error'][i,j] = obs.xoffset_err
        param_array['yoffset Error'][i,j] = obs.yoffset_err
        param_array['azfwhm Error'][i,j] = obs.azfwhm_err
        param_array['elfwhm Error'][i,j] = obs.elfwhm_err
        
        if abs(obs.azfwhm_err) >=azfwhm_lim:
            print('        Az FWHM too large for beam %i detector %i with value of %f' % (j,i, obs.azfwhm_err))
        if abs(obs.elfwhm_err) >=elfwhm_lim:
            print('        El FWHM too large for beam %i detector %i with value of %f' % (j,i, obs.elfwhm_err))

 
#Names for plotting.
param_names = ['Amplitude', 'xoffset', 'yoffset', 'azfwhm', 'elfwhm', 
               'Amplitude Error', 'xoffset Error', 'yoffset Error',
              'azfwhm Error', 'elfwhm Error']

#Plots the fit values and errors against the observation number.  Each figure is slightly interactive in that one can click on
#an observation to open up additional figures for that observation.  See intro.
for j in range(len(param_names)):
    f = plt.figure(j)
    cid = f.canvas.mpl_connect('button_press_event', obs_click)
    ax = f.add_subplot(111)
    plt.grid(True)
    for i in range(ndetectors):
        ax.scatter(range(beam_num), param_array[param_names[j]][i,:], s=5,
        marker = 's', c='k')
        ax.set_xlabel('observation')
        ax.set_ylabel(param_names[j])
    plt.axis('tight')
    
plt.show()
