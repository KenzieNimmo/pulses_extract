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
auto_waterfaller.histogram(dm, title='Distribution of Dispersion Measures \n%s'%observation,\
							xlabel=(r'DM (pc cm$^{-3}$)'), color='r', name='DM')
auto_waterfaller.histogram(SNR, title='Distribution of Signal to Noise Ratios\n%s'%observation,\
							xlabel='S/N', color='b', name='SN')
auto_waterfaller.histogram((duration*1000.), title='Distribution of Burst Durations\n%s'%observation,\
							xlabel='Duration (ms)', color='g', name='width')
auto_waterfaller.toa_plotter(time, SNR, duration) #no rankings passed

