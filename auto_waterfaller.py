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

def main(fits, time, DM, top_freq=0.0, sigma=0.0, duration=0.0, directory='.', FRB_name='FRB121102'):

	rawdata = psrfits.PsrfitsFile(fits)
	observation = str(fits)[:-5]
	observation = os.path.basename(observation)

	#Open header of the fits file
	with psrfits.pyfits.open(fits, memmap=True) as fn:
  		header = fn['SUBINT'].header + fn['PRIMARY'].header

	#MJD of the beginning of the observation
	start_MJD = header['STT_IMJD'] + header['STT_SMJD'] / 86400. 

	#MJD of the pulses
	pulse_MJD = start_MJD + time / 86400.

	for i, t in enumerate(time[:1]): 
		start_time = t - 0.01
		data, nbinsextra, nbins, start = waterfall(rawdata, start_time, 0.03, DM[i],\
				nbins=None, nsub=None, subdm = DM, zerodm=False, downsamp=1,\
				scaleindep=False, width_bins=1, mask=False, maskfn=None,\
				bandpass_corr=False, ref_freq=None)

		#plot it
		fig = plt.figure(figsize=(16, 7))
		ax1 = plt.subplot2grid((3,4), (1,1), rowspan=3, colspan=3)
		ax2 = plt.subplot2grid((3,4), (0,0), rowspan=3)
		ax3 = plt.subplot2grid((3,4), (0,1), colspan=3)
		ax2.axis([0,10,0,10])
		ax2.annotate('DM = %0.2f'%DM[i], xy=(0,9))
		ax2.annotate('Time (s) = %0.4f'%t, xy=(0,6))
		ax2.annotate('MJD = %d'%start_MJD, xy=(0,8))
		ax2.annotate('MJD time = %0.8f'%(t/86400), xy=(0,7))
		ax2.annotate('Duration = %0.4f'%duration, xy=(0,5))
		ax2.annotate('Top frequency = %0.2f'%top_freq, xy=(0,4))
		ax2.annotate('Sigma = %0.2f'%sigma, xy=(0,3))

		plot_waterfall(data, start, 0.03, 
	                   integrate_ts=True, integrate_spec=False, show_cb=False, 
	                   cmap_str="gist_yarg", sweep_dms=[], sweep_posns=[],
	                   ax_im=ax1, ax_ts=ax3, ax_spec=None, interactive=False)

		ax3.axvline(t, c='r')
		ax2.axis('off')
		fig.tight_layout(w_pad = 10, h_pad = 0.5)
		plt.suptitle('%s %0.8f [close up]\n %s'%(FRB_name, pulse_MJD[i], observation), y=1.05)
		plt.savefig('%s/%s_%0.8f_zoomed.png'%(directory, FRB_name, pulse_MJD[i]), bbox_inches='tight', pad_inches=0.1)

		#Zoomed out version
		start_time = t - 0.05

		data, nbinsextra, nbins, start = waterfall(rawdata, start_time, 0.1, DM[i],\
				nbins=None, nsub=None, subdm = DM, zerodm=False, downsamp=1,\
				scaleindep=False, width_bins=1, mask=False, maskfn=None,\
				bandpass_corr=False, ref_freq=None)

		fig = plt.figure()
		ax1 = plt.subplot2grid((3,4), (1,1), rowspan=3, colspan=3)
		ax2 = plt.subplot2grid((3,4), (0,0), rowspan=3)
		ax3 = plt.subplot2grid((3,4), (0,1), colspan=3)
		ax2.axis([0,10,0,10])
		ax2.annotate('DM = %0.2f'%DM[i], xy=(0,9))
		ax2.annotate('Time (s) = %0.4f'%t, xy=(0,6))
		ax2.annotate('MJD = %d'%start_MJD, xy=(0,8))
		ax2.annotate('MJD time = %0.8f'%(t/86400), xy=(0,7))
		ax2.annotate('Duration = %0.4f'%duration, xy=(0,5))
		ax2.annotate('Top frequency = %0.2f'%top_freq, xy=(0,4))
		ax2.annotate('Sigma = %0.2f'%sigma, xy=(0,3))

		plot_waterfall(data, start, 0.1, 
	                   integrate_ts=True, integrate_spec=False, show_cb=False, 
	                   cmap_str="gist_yarg", sweep_dms=[], sweep_posns=[],
	                   ax_im=None, ax_ts=None, ax_spec=None, interactive=False)

		ax3.axvline(t, c='r')
		ax2.axis('off')
		fig.tight_layout(w_pad = 10, h_pad = 0.5)		
		plt.suptitle('FRB121102 %0.8f \n %s'%(pulse_MJD[i],observation), y=1.05)
		plt.savefig('%s/FRB121102_%0.8f.png'%(directory, pulse_MJD[i]), bbox_inches='tight', pad_inches=0.1)
		plt.close('all')

if __name__ == '__main__':
	DM, time, sample = np.loadtxt(sys.argv[2], usecols=(0,2,3), unpack=True)

	main(sys.argv[1],time,DM)
