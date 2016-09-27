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

def plotter(data, start, plot_duration, t, DM, IMJD, SMJD, duration, top_freq, sigma, directory, FRB_name, observation, pulse_id, zoom=True, idx='', downsamp=True):
	fig = plt.figure(figsize=(8,5))
	ax1 = plt.subplot2grid((3,4), (1,1), rowspan=3, colspan=3)
	ax2 = plt.subplot2grid((3,4), (0,0), rowspan=3)
	ax3 = plt.subplot2grid((3,4), (0,1), colspan=3)

	ax2.axis([0,7,0,7])
	ax2.annotate('Pulse ID \n%i'%pulse_id, xy=(0,6))
	ax2.annotate('DM \n%d'%DM, xy=(0,5))
	ax2.annotate('MJD \n%.8f'%(IMJD+SMJD), xy=(0,4))
	ax2.annotate('Time (s) \n%0.3f'%t, xy=(0,3))
	ax2.annotate('Duration (ms) \n%0.2f'%(duration*1000.), xy=(0,2))
	ax2.annotate('Top frequency \n%0.2f'%top_freq, xy=(0,1))
	ax2.annotate('Sigma \n%0.2f'%sigma, xy=(0,0))

	plot_waterfall(data, start, plot_duration, 
	               integrate_ts=True, integrate_spec=False, show_cb=False, 
	               cmap_str="gist_yarg", sweep_dms=[], sweep_posns=[],
	               ax_im=ax1, ax_ts=ax3, ax_spec=None, interactive=False)


	ax3.axvline(t, c='r')
	ax2.axis('off')
	fig.tight_layout(w_pad = 8, h_pad = 0.5)

	if zoom and downsamp:
		title = (' [close up with downsamp = %d]'%downsamp)
		name = '_zoomed_downsamped'
	elif zoom:
		title = (' [close up]')
		name = '_zoomed'
	else:
		title = name = ''
          
	plt.suptitle('%s %id %.8fs %s\n %s'%(FRB_name, IMJD, SMJD, title, observation), y=1.05)
	plt.savefig('%s/%s_%i_%.8f_%s%s.png'%(directory, FRB_name, IMJD, SMJD, idx, name), bbox_inches='tight', pad_inches=0.2)
	fig.clf()
	plt.close('all')

def main(fits, time, DM=560., sigma=0., duration=0., pulse_id=0, top_freq=0., directory='.', FRB_name='FRB121102', downsamp=1.):
        num_elements = time.size
        if isinstance(DM, float) or isinstance(DM, int): DM = np.zeros(num_elements) + DM
        if isinstance(sigma, float) or isinstance(sigma, int): sigma = np.zeros(num_elements) + sigma
        if isinstance(duration, float) or isinstance(duration, int): duration = np.zeros(num_elements) + duration
        if isinstance(pulse_id, float) or isinstance(pulse_id, int): pulse_id = np.zeros(num_elements) + pulse_id
        
	rawdata = psrfits.PsrfitsFile(fits)
	observation = str(fits)[:-5]
	observation = os.path.basename(observation)

	#Open header of the fits file
	with psrfits.pyfits.open(fits, memmap=True) as fn:
  		header = fn['SUBINT'].header + fn['PRIMARY'].header

	#Start MJD (days) of the beginning of the observation
	IMJD = header['STT_IMJD'] 

	#Fractional day
	SMJD = (header['STT_SMJD'] + time) / 86400.

	for i, t in enumerate(time[:1]): 
		start_time = t - 0.05
		plot_duration = 0.1

		data, nbinsextra, nbins, start = waterfall(rawdata, start_time, plot_duration, DM[i],\
				nbins=None, nsub=None, subdm = DM, zerodm=False, downsamp=1,\
				scaleindep=False, width_bins=1, mask=False, maskfn=None,\
				bandpass_corr=False, ref_freq=None)

		plotter(data, start, plot_duration, t, DM[i], IMJD, SMJD[i], duration[i], top_freq,\
			sigma[i], directory, FRB_name, observation, zoom=False, idx=i, pulse_id=pulse_id[i], downsamp=False)
                
        #Zoomed version
 		start_time = t - 0.01
		plot_duration = 0.03

		data, nbinsextra, nbins, start = waterfall(rawdata, start_time, plot_duration, DM[i],\
				nbins=None, nsub=None, subdm = DM, zerodm=False, downsamp=1,\
				scaleindep=False, width_bins=1, mask=False, maskfn=None,\
				bandpass_corr=False, ref_freq=None)

		plotter(data, start, plot_duration, t, DM[i], IMJD, SMJD[i], duration[i], top_freq,\
			sigma[i], directory, FRB_name, observation, zoom=True, idx=i, pulse_id=pulse_id[i], downsamp=False)               
        
        #downsamped version (zoomed)
		#print downsamp
		data, nbinsextra, nbins, start = waterfall(rawdata, start_time, plot_duration, DM[i],\
				nbins=None, nsub=None, subdm = DM, zerodm=False, downsamp=downsamp,\
				scaleindep=False, width_bins=1, mask=False, maskfn=None,\
				bandpass_corr=False, ref_freq=None)

		plotter(data, start, plot_duration, t, DM[i], IMJD, SMJD[i], duration[i], top_freq,\
			sigma[i], directory, FRB_name, observation, zoom=True, idx=i, pulse_id=pulse_id[i], downsamp=downsamp)
                

def histogram(data, title='', xlabel='', color='k'):

  #...

  pass





if __name__ == '__main__':
	DM, time, sample = np.loadtxt(sys.argv[2], usecols=(0,2,3), unpack=True)

	main(sys.argv[1],time,DM)
