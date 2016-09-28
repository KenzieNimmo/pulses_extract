import pandas as pd
import auto_waterfaller
import numpy as np

database = '/psr_temp/hessels/AO-FRB/P3054/FRB_pipeline/output/puppi_57614_C0531+33_0803/OLD/SinglePulses.hdf5'

pulses = pd.read_hdf(database,'pulses')

dm = np.array(pulses.DM)
SNR = np.array(pulses.Sigma)
time = np.array(pulses.Time)
duration = np.array(pulses.Duration)




auto_waterfaller.histogram(dm)