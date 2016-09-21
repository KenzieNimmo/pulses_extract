from waterfaller import waterfall, plot_waterfall
from matplotlib import pyplot as plt
import matplotlib.cm
import numpy as np
import sys
import optparse
import copy
import psr_utils
import rfifind
import psrfits
import spectra

def main(fits, time, DM, directory=''):# (fits, time, DM, time_resolution=81.9e-6, )
    if directory != '':
    	directory += '/'
    else:
    	pass

	rawdata = psrfits.PsrfitsFile(fits)
	#time_resolution = 81.9e-6 #only need this if need downsamp
	#Downsamp = time[0]/sample[0]/time_resolution
	count = 0

	for i in time: 
		start_time = i - 0.01
		data, nbinsextra, nbins, start = waterfall(rawdata, start_time, 0.03, DM[0],\
				nbins=None, nsub=None, subdm = DM, zerodm=False, downsamp=1,\
				scaleindep=False, width_bins=1, mask=False, maskfn=None,\
				bandpass_corr=False, ref_freq=None)

		plot_waterfall(data, start, 0.03, 
	                   integrate_ts=True, integrate_spec=False, show_cb=False, 
	                   cmap_str="gist_yarg", sweep_dms=[], sweep_posns=[],
	                   ax_im=None, ax_ts=None, ax_spec=None, interactive=False)

		count += 1
		#plt.title()
		plt.savefig('%stest_burst_%d_zoomed.png'%(directory, count))

		#Zoomed out version
		start_time = i - 0.05

		data, nbinsextra, nbins, start = waterfall(rawdata, start_time, 0.1, DM[0],\
				nbins=None, nsub=None, subdm = DM, zerodm=False, downsamp=1,\
				scaleindep=False, width_bins=1, mask=False, maskfn=None,\
				bandpass_corr=False, ref_freq=None)

		plot_waterfall(data, start, 0.1, 
	                   integrate_ts=True, integrate_spec=False, show_cb=False, 
	                   cmap_str="gist_yarg", sweep_dms=[], sweep_posns=[],
	                   ax_im=None, ax_ts=None, ax_spec=None, interactive=False)

		plt.savefig('%stest_burst_%d.png'%(directory, count))
		#plt.savefig('test_burst_%d.png'%count)






if __name__ == '__main__':
	DM, time, sample = np.loadtxt(sys.argv[2], usecols=(0,2,3), unpack=True) #can simplify if don't need downsamp info

	main(sys.argv[1],time,DM) #
