"""
auto_waterfaller.py

Kelly Gourdji - Sept. 2016

"""

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
import os

def main(fits, time, DM, directory='.', FRB_name='FRB121102'):

	rawdata = psrfits.PsrfitsFile(fits)
	observation = str(fits)[:-5]
	observation = os.path.basename(observation)

	#Open header of the fits file
	with psrfits.pyfits.open(fits, memmap=True) as fn:
  		header = fn['SUBINT'].header + fn['PRIMARY'].header

	#Start MJD (days) of the beginning of the observation
	IMJD = header['STT_IMJD'] 

	#Seconds past UTC 00h  of the pulses
	SMJD = header['STT_SMJD'] + time

	for i, t in enumerate(time): 
		start_time = t - 0.01
		data, nbinsextra, nbins, start = waterfall(rawdata, start_time, 0.03, DM[i],\
				nbins=None, nsub=None, subdm = DM, zerodm=False, downsamp=1,\
				scaleindep=False, width_bins=1, mask=False, maskfn=None,\
				bandpass_corr=False, ref_freq=None)

		plt.figure()
		plot_waterfall(data, start, 0.03, 
	                   integrate_ts=True, integrate_spec=False, show_cb=False, 
	                   cmap_str="gist_yarg", sweep_dms=[], sweep_posns=[],
	                   ax_im=None, ax_ts=None, ax_spec=None, interactive=False)

		plt.suptitle('%s %i.%i [close up]\n %s'%(FRB_name, IMJD, SMJD[i], observation), y=1.05)
		plt.savefig('%s/%s_%i.%i_zoomed.png'%(directory, FRB_name, IMJD, SMJD[i]), bbox_inches='tight', pad_inches=0.1)

		#Zoomed out version
		start_time = t - 0.05

		data, nbinsextra, nbins, start = waterfall(rawdata, start_time, 0.1, DM[i],\
				nbins=None, nsub=None, subdm = DM, zerodm=False, downsamp=1,\
				scaleindep=False, width_bins=1, mask=False, maskfn=None,\
				bandpass_corr=False, ref_freq=None)

		plt.figure()
		plot_waterfall(data, start, 0.1, 
	                   integrate_ts=True, integrate_spec=False, show_cb=False, 
	                   cmap_str="gist_yarg", sweep_dms=[], sweep_posns=[],
	                   ax_im=None, ax_ts=None, ax_spec=None, interactive=False)
		
		plt.suptitle('%s %i.%i \n %s'%(FRB_name, IMJD, SMJD[i],observation), y=1.05)
                plt.savefig('%s/%s_%i.%i.png'%(directory, FRB_name, IMJD, SMJD[i]), bbox_inches='tight', pad_inches=0.1)
		plt.close('all')

if __name__ == '__main__':
	DM, time, sample = np.loadtxt(sys.argv[2], usecols=(0,2,3), unpack=True)

	main(sys.argv[1],time,DM)
