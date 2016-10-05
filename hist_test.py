import pandas as pd
import auto_waterfaller
import numpy as np

database = '/psr_temp/hessels/AO-FRB/P3054/FRB_pipeline/output/puppi_57614_C0531+33_0803/OLD/SinglePulses.hdf5'

pulses = pd.read_hdf(database,'pulses')

dm = np.array(pulses.DM)
SNR = np.array(pulses.Sigma)
time = np.array(pulses.Time)
duration = np.array(pulses.Duration)
observation = "test_observation"
Rank = np.random.randint(0,3,len(time))

auto_waterfaller.plot_statistics()

		histogram(dm[np.where(Rank==0)], ax = ax1, title='Distribution of Dispersion Measures \n%s'%observation,\
								xlabel=(r'DM (pc cm$^{-3}$)'), color='r', name='DM', stacked=True)
		histogram(dm[np.where(Rank==1)], ax = ax1, title='Distribution of Dispersion Measures \n%s'%observation,\
								xlabel=(r'DM (pc cm$^{-3}$)'), color='g', name='DM', stacked=True)
		histogram(dm[np.where(Rank==2)], ax = ax1, title='Distribution of Dispersion Measures \n%s'%observation,\
								xlabel=(r'DM (pc cm$^{-3}$)'), color='k', name='DM', stacked=True)
		#ranked SNR hist
		histogram(SNR[np.where(Rank==0)], ax= ax2, title='Distribution of Signal to Noise Ratios\n%s'%observation,\
								xlabel='S/N', color='r', name='SN', stacked=True)
		histogram(SNR[np.where(Rank==1)], ax= ax2, title='Distribution of Signal to Noise Ratios\n%s'%observation,\
								xlabel='S/N', color='g', name='SN', stacked=True)
		histogram(SNR[np.where(Rank==2)], ax= ax2, title='Distribution of Signal to Noise Ratios\n%s'%observation,\
								xlabel='S/N', color='k', name='SN', stacked=True)

		#ranked duration hist
		histogram((duration[np.where(Rank==0)])*1000., ax=ax3, title='Distribution of Burst Durations\n%s'%observation,\
								xlabel='Duration (ms)', color='r', name='width', stacked=True)
		histogram((duration[np.where(Rank==1)])*1000., ax=ax3, title='Distribution of Burst Durations\n%s'%observation,\
								xlabel='Duration (ms)', color='g', name='width', stacked=True)
		histogram((duration[np.where(Rank==2)])*1000., ax=ax3, title='Distribution of Burst Durations\n%s'%observation,\
								xlabel='Duration (ms)', color='k', name='width', stacked=True)