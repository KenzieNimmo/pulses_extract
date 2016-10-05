"""
auto_waterfaller.py

Kelly Gourdji - Sept. 2016

"""
from waterfaller import waterfall, plot_waterfall
from matplotlib import pyplot as plt
from matplotlib import cm
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
import pandas as pd
import matplotlib.lines as mlines

def plotter(data, start, plot_duration, t, DM, IMJD, SMJD, duration, top_freq, sigma, 
			directory, FRB_name, observation, pulse_id, zoom=True, idx='', downsamp=True):
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
          
	plt.suptitle('%s %.8f %s\n %s'%(FRB_name, IMJD + SMJD, title, observation), y=1.05)
	if not os.path.isdir('%s/%s'%(directory, pulse_id)): os.makedirs('%s/%s'%(directory, pulse_id))
	plt.savefig('%s/%s/%s_%.8f_%s%s.png'%(directory, pulse_id, FRB_name, IMJD + SMJD, idx, name),\
										   bbox_inches='tight', pad_inches=0.2)
	fig.clf()
	plt.close('all')

def histogram(data, ax, title='', xlabel='', color='', name='', bins=None, stacked=False): #optionally choose number of bins. Currently: 2*sqrt(counts)
	"""
	Creates a histogram of a given burst property. 

	Inputs:

		name: burst property being plotted

	"""
	if not bins:
		bins = 2 * np.sqrt(len(data))
	ax.hist(data, bins=bins, color=color, histtype='step', lw=2, stacked=stacked)
	ax.set_xlabel(xlabel, fontsize=8)
	ax.set_ylabel('Counts', fontsize=8)
	ax.set_title(title, fontsize=8)
	ax.tick_params(axis='x', labelsize=8)
	ax.tick_params(axis='y', labelsize=8)

def toa_plotter(time, SN, duration, Rank, observation='', ax=None):
	"""
	Plots a bar at each candidate time. Bar width corresponds to pulse duration, while its
	height corresponds to the signal to noise ratio of the burst.

	Inputs:

		time: numpy array of arrival times
		SN: numpy array of signal-to-noise ratios
		duration: numpy array of burst durations

	Optional Input:
		Rank: numpy array of pulse rankings for color-mapping. Default is no coloring.
		observation: observation name (str). Default is none.


	"""
	rank_colors = cm.colors.LinearSegmentedColormap.from_list('rank_colors', [(0,'green'), (0.5,'#D4AC0D'), (1,'red')]) 
	norm=cm.colors.Normalize(vmin=0, vmax=2)
	ranks=norm(Rank)
	ax.bar(time, SN, duration, color=rank_colors(ranks), edgecolor=rank_colors(ranks))
	ax.set_xlabel('Time (s)',fontsize=10)
	ax.set_ylabel('S/N', fontsize=10)
	ax.set_title('Times of Arrival v. Signal to Noise Ratio\n%s'%observation, fontsize=12)
	ax.tick_params(axis='x', labelsize=9)
	ax.tick_params(axis='y', labelsize=9)

def plot_statistics(dm, SNR, duration, Rank=False):
	fig = plt.figure(figsize=(8,5))
	ax1 = plt.subplot2grid((2,3), (0,0))
	ax2 = plt.subplot2grid((2,3), (0,1))
	ax3 = plt.subplot2grid((2,3), (0,2))
	ax4 = plt.subplot2grid((2,3), (1,0), colspan=3)

        colors = ['green', '#D4AC0D', 'red']
        toa_plotter(time, SNR, duration, Rank, ax=ax4)
    
        DMs = [dm[np.where(Rank==0)], dm[np.where(Rank==1)], dm[np.where(Rank==2)]]
        histogram(DMs, ax = ax1, title='Distribution of Dispersion Measures \n%s'%observation,\
                                                    xlabel=(r'DM (pc cm$^{-3}$)'), color=colors, name='DM')
    
        SNRs = [SNR[np.where(Rank==0)], SNR[np.where(Rank==1)], SNR[np.where(Rank==2)]]
        histogram(SNRs, ax= ax2, title='Distribution of Signal to Noise Ratios\n%s'%observation,\
                                                    xlabel='S/N', color=colors, name='SN')

        durations = [duration[np.where(Rank==0)]*1000., duration[np.where(Rank==1)]*1000., duration[np.where(Rank==2)]*1000.]
        histogram(durations, ax=ax3, title='Distribution of Burst Durations\n%s'%observation,\
                                                    xlabel='Duration (ms)', color=colors, name='width')
        fig.tight_layout(w_pad = 0.3, h_pad = 0.5)
        if ranked: rank = '_ranked'
        else: rank = ''
        plt.savefig('%s/statistical_plots_%s%s.png'%(folder,observation,rank), bbox_inches='tight')


def main(fits, time, DM=560., sigma=0., duration=0.01, pulse_id=0, top_freq=0., directory='.',\
		  FRB_name='FRB121102', downsamp=1.):
        num_elements = time.size
        if isinstance(DM, float) or isinstance(DM, int): DM = np.zeros(num_elements) + DM
        if isinstance(sigma, float) or isinstance(sigma, int): sigma = np.zeros(num_elements) + sigma
        if isinstance(duration, float) or isinstance(duration, int): duration = np.zeros(num_elements) + duration
        if isinstance(pulse_id, float) or isinstance(pulse_id, int): pulse_id = np.zeros(num_elements) + pulse_id
        if isinstance(downsamp, float) or isinstance(downsamp, int): downsamp = np.zeros(num_elements) + downsamp
        
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

	for i, t in enumerate(time): 
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
		data, nbinsextra, nbins, start = waterfall(rawdata, start_time, plot_duration, DM[i],\
				nbins=None, nsub=None, subdm = DM, zerodm=False, downsamp=downsamp[i],\
				scaleindep=False, width_bins=1, mask=False, maskfn=None,\
				bandpass_corr=False, ref_freq=None)

		plotter(data, start, plot_duration, t, DM[i], IMJD, SMJD[i], duration[i], top_freq,\

				sigma[i], directory, FRB_name, observation, zoom=True, idx=i, pulse_id=pulse_id[i],\
			 	downsamp=downsamp[i])

if __name__ == '__main__':
	#DM, sigma, time, downfact = np.loadtxt(sys.argv[2], usecols=(0,1,2,4), unpack=True)
	#downsamp = np.zeros(len(downfact)) + 1. #just a place holder so my code runs upon testing.
	#main(sys.argv[1],time, DM, sigma, downsamp = downsamp)
	database = '/psr_temp/hessels/AO-FRB/P3054/FRB_pipeline/output/puppi_57614_C0531+33_0803/OLD/SinglePulses.hdf5'

	pulses = pd.read_hdf(database,'pulses')
	dm = np.array(pulses.DM)
	SNR = np.array(pulses.Sigma)
	time = np.array(pulses.Time)
	duration = np.array(pulses.Duration)
	observation = "test_observation"
	Rank = np.random.randint(0,3,len(time))
	plot_statistics(dm, SNR, duration)
