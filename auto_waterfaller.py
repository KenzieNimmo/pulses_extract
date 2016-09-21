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

def main(fits, time, DM, directory='.', FRB_name='FRB121102'):
      rawdata = psrfits.PsrfitsFile(fits)
      count = 0
      observation = str(fits)[:-5]

      #Open header of the fits file
      with psrfits.pyfits.open(fits, memmap=True) as fn:
              header = fn['SUBINT'].header + fn['PRIMARY'].header

      #MJD of the beginning of the observation
      start_MJD = header['STT_IMJD'] + header['STT_SMJD'] / 86400. 

      #MJD of the pulses
      pulse_MJD = start_MJD + time / 86400.

      for i in time: 
              start_time = i - 0.01
              data, nbinsextra, nbins, start = waterfall(rawdata, start_time, 0.03, DM[0],\
                              nbins=None, nsub=None, subdm = DM, zerodm=False, downsamp=1,\
                              scaleindep=False, width_bins=1, mask=False, maskfn=None,\
                              bandpass_corr=False, ref_freq=None)

              plt.figure()
              plot_waterfall(data, start, 0.03, 
                          integrate_ts=True, integrate_spec=False, show_cb=False, 
                          cmap_str="gist_yarg", sweep_dms=[], sweep_posns=[],
                          ax_im=None, ax_ts=None, ax_spec=None, interactive=False)

              plt.suptitle('%s %0.8f [close up]\n %s'%(FRB_name, pulse_MJD[count], observation), y=1.05)
              plt.savefig('%s/%s_%0.8f_zoomed.png'%(directory, FRB_name, pulse_MJD[count]), bbox_inches='tight', pad_inches=0.1)

              #Zoomed out version
              start_time = i - 0.05

              data, nbinsextra, nbins, start = waterfall(rawdata, start_time, 0.1, DM[0],\
                              nbins=None, nsub=None, subdm = DM, zerodm=False, downsamp=1,\
                              scaleindep=False, width_bins=1, mask=False, maskfn=None,\
                              bandpass_corr=False, ref_freq=None)
              plt.figure()
              plot_waterfall(data, start, 0.1, 
                          integrate_ts=True, integrate_spec=False, show_cb=False, 
                          cmap_str="gist_yarg", sweep_dms=[], sweep_posns=[],
                          ax_im=None, ax_ts=None, ax_spec=None, interactive=False)
              
              plt.suptitle('FRB121102 %0.8f \n %s'%(pulse_MJD[count],observation), y=1.05)
              plt.savefig('%s/FRB121102_%0.8f.png'%(directory, pulse_MJD[count]), bbox_inches='tight', pad_inches=0.1)
              plt.close('all')
              count += 1





if __name__ == '__main__':
	DM, time, sample = np.loadtxt(sys.argv[2], usecols=(0,2,3), unpack=True)

	main(sys.argv[1], time, DM)
