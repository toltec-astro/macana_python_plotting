import numpy as np
import glob
import re
import sys
import matplotlib.pyplot as plt
import matplotlib as mpl

plt.close('all')

def onclick(event):
    global ix, iy
    ix, iy = event.xdata, event.ydata
    print 'x = %d, y = %f'%(
        ix, iy)
 
    #Values   
    plt.figure()   
    ax1 = plt.subplot2grid(shape=(2,6), loc=(0,0), colspan=2)
    ax2 = plt.subplot2grid((2,6), (0,2), colspan=2)
    ax3 = plt.subplot2grid((2,6), (0,4), colspan=2)
    ax4 = plt.subplot2grid((2,6), (1,1), colspan=2)
    ax5 = plt.subplot2grid((2,6), (1,3), colspan=2)
    
    ax1.plot(range(ndetectors), param_array['Amplitude'][:,int(ix)])
    ax1.set_ylabel('Amplitude')
    ax1.set_xlabel('Detector number')

    ax2.plot(range(ndetectors), param_array['xoffset'][:,int(ix)])
    ax2.set_ylabel('xoffset')
    ax2.set_xlabel('Detector number')

    ax3.plot(range(ndetectors), param_array['yoffset'][:,int(ix)])
    ax3.set_ylabel('yoffset')
    ax3.set_xlabel('Detector number')

    ax4.plot(range(ndetectors), param_array['azfwhm'][:,int(ix)])
    ax4.set_ylabel('azfwhm')
    ax4.set_xlabel('Detector number')
    
    ax5.plot(range(ndetectors), param_array['elfwhm'][:,int(ix)])
    ax5.set_ylabel('elfwhm')
    ax5.set_xlabel('Detector number')

    plt.tight_layout()
    
    #Errors
    plt.figure()   
    ax6 = plt.subplot2grid(shape=(2,6), loc=(0,0), colspan=2)
    ax7 = plt.subplot2grid((2,6), (0,2), colspan=2)
    ax8 = plt.subplot2grid((2,6), (0,4), colspan=2)
    ax9 = plt.subplot2grid((2,6), (1,1), colspan=2)
    ax10 = plt.subplot2grid((2,6), (1,3), colspan=2)
    
    ax6.plot(range(ndetectors), param_array['Amplitude Error'][:,int(ix)])
    ax6.set_ylabel('Amplitude error')
    ax6.set_xlabel('Detector number')

    ax7.plot(range(ndetectors), param_array['xoffset Error'][:,int(ix)])
    ax7.set_ylabel('xoffset error')
    ax7.set_xlabel('Detector number')

    ax8.plot(range(ndetectors), param_array['yoffset Error'][:,int(ix)])
    ax8.set_ylabel('yoffset error')
    ax8.set_xlabel('Detector number')

    ax9.plot(range(ndetectors), param_array['azfwhm Error'][:,int(ix)])
    ax9.set_ylabel('azfwhm error')
    ax9.set_xlabel('Detector number')
    
    ax10.plot(range(ndetectors), param_array['elfwhm Error'][:,int(ix)])
    ax10.set_ylabel('elfwhm error')
    ax10.set_xlabel('Detector number')

    plt.tight_layout()
    plt.show()

coords = []
nplots = 3
nparams = 5

#plt.ioff()

sys.path.insert(0, '/home/mmccrackan/macana_python_plotting/')

from macana_classes import macana_plotter

beammap_path = '/home/mmccrackan/Documents/TolTEC-Project/beammaps_out/'

beammaps = glob.glob(beammap_path + '*')
beam_num = len(beammaps)

ndetectors = 113

amp_lim = 0.01
xoffset_lim = 1.0
yoffset_lim = 1.0
azfwhm_lim = 1.0
elfwhm_lim = 1.0

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
        
        if obs.azfwhm_err >=azfwhm_lim:
            print('        Az FWHM too large for beam %i detector %i with value of %f' % (j,i, obs.azfwhm_err))
        if obs.elfwhm_err >=azfwhm_lim:
            print('        El FWHM too large for beam %i detector %i with value of %f' % (j,i, obs.elfwhm_err))

#np.save('/home/mmccrackan/Documents/TolTEC-Project/fit_array', param_array)

param_names = ['Amplitude', 'xoffset', 'yoffset', 'azfwhm', 'elfwhm', 
               'Amplitude Error', 'xoffset Error', 'yoffset Error',
               'azfwhm Error', 'elfwhm Error']
'''
for i in range(ndetectors):
    for j in range(len(param_names)):
        print('On detector %i, plotting %s' % (i, param_names[j]))
        plt.figure(0)
        plt.hist(param_array[param_names[j]][i,:], 25, histtype='step', edgecolor='k')
        plt.xlabel(param_names[j])
        plt.title(param_names[j] + ' for detector ' + str(i))
        plt.savefig('/home/mmccrackan/Documents/TolTEC-Project/detector_figs/detector' + str(i) + '_' + param_names[j])
        plt.close('all')
        
for i in range(beam_num):
    for j in range(len(param_names)):
        print('On beam %i, plotting %s' % (i, param_names[j]))
        plt.figure()
        plt.hist(param_array[param_names[j]][:,i], 25, histtype='step', edgecolor='k')
        plt.xlabel(param_names[j])
        plt.title(param_names[j] + ' for beam ' + str(i))
        plt.savefig('/home/mmccrackan/Documents/TolTEC-Project/beam_figs/beam' + str(i) + '_' + param_names[j])
        plt.close('all')

for j in range(len(param_names)):
    f = plt.figure()
    cid = f.canvas.mpl_connect('button_press_event', onclick)
    ax = f.add_subplot(111)
    plt.grid(True)
    for i in range(beam_num):
        ax.scatter(range(ndetectors), param_array[param_names[j]][:,i])
        ax.set_xlabel('Detector')
        ax.set_ylabel(param_names[j])
'''
for j in range(len(param_names)):
    f = plt.figure()
    cid = f.canvas.mpl_connect('button_press_event', onclick)
    ax = f.add_subplot(111)
    plt.grid(True)
    for i in range(ndetectors):
        ax.scatter(range(beam_num), param_array[param_names[j]][i,:])
        ax.set_xlabel('observation')
        ax.set_ylabel(param_names[j])
    plt.axis('tight')

#cid = plt.gcf().canvas.mpl_connect('button_press_event', onclick)

plt.show()