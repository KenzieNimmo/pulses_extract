import os
import argparse

import matplotlib.pyplot as plt
import numpy as np
import psrchive
import scipy.misc

def parser():
    # Command-line options
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                    description="Plot dynamic spectrum from single archives with 1 polarisation and 1 subintegration.")
    parser.add_argument('archive_name', help="Name of the psrchive file to plot.")
    parser.add_argument('-show', help="Show the plot.", action='store_false')
    parser.add_argument('-save_fig', help="Save the plot.", action='store_true')
    parser.add_argument('-zap', help="Plot to manually zap bins out.", action='store_true')
    return parser.parse_args()

def main():
  args = parser()
  DS, extent = load_DS(args.archive_name)
  if args.zap: extent = None
  zap(args.archive_name, DS)
  plot_DS(DS, args.archive_name, extent=extent, show=args.show, save=args.save_fig)
  
def plot_DS(DS, archive_name, extent=None, show=True, save=False):
  fig = plt.figure()
  
  #Dynamic spectrum
  ax1 = plt.subplot2grid((5,5), (1,0), rowspan=4, colspan=4)
  smooth_DS = scipy.misc.imresize(DS, 0.5, interp='cubic').astype(np.float)
  smooth_DS -= np.median(smooth_DS)
  smooth_DS /= smooth_DS.max()
  ax1.imshow(smooth_DS, cmap='RdGy_r', origin='upper', aspect='auto', interpolation='nearest', extent=extent)
  ax1.set_xlabel("Time (ms)")
  ax1.set_ylabel("Frequency (MHz)")
  
  if not extent: extent = [0, smooth_DS.shape[1], 0, smooth_DS.shape[0]]
  
  #Pulse profile
  ax2 = plt.subplot2grid((5,5), (0,0), colspan=4, sharex=ax1)
  prof = np.mean(smooth_DS, axis=0)
  x = np.linspace(extent[0], extent[1], prof.size)
  ax2.plot(x, prof, 'k-')
  ax2.tick_params(axis='y', which='both', left='off', right='off', labelleft='off')
  ax2.tick_params(axis='x', labelbottom='off')
  ax2.set_xlim(extent[0:2])
  
  #Baseline
  ax3 = plt.subplot2grid((5,5), (1,4), rowspan=4, sharey=ax1)
  bl = np.mean(smooth_DS, axis=1)
  y = np.linspace(extent[3], extent[2], bl.size)
  ax3.plot(bl, y, 'k-')
  ax3.tick_params(axis='x', which='both', top='off', bottom='off', labelbottom='off')
  ax3.tick_params(axis='y', labelleft='off')
  ax3.set_ylim(extent[2:4])
  
  #General plot settings
  fig.tight_layout()
  fig.subplots_adjust(hspace=0)
  title = os.path.split(os.path.basename(archive_name))[0]
  plt.title(title)
  
  if show: plt.show()
  if save: fig.savefig(title)
  return 

def load_DS(archive_name):
  load_archive = psrchive.Archive_load(archive_name)
  load_archive.remove_baseline()
  DS = load_archive.get_data().squeeze()
  
  freq = load_archive.get_centre_frequency()
  bw = abs(load_archive.get_bandwidth())
  duration = load_archive.integration_length() * 1000
  
  return DS, [0., duration, freq-bw/2, freq+bw/2]
    
def zap(archive_name, DS):
  zap_list = load_zap_list(archive_name)
  
  med = np.median(DS)
  for chan in zap_list:
    DS[chan[0], chan[1]:chan[2]] = med
  DS -= med
  DS /= DS.max()
  return
  
def load_zap_list(archive_name):
  '''
  List of bins to zap in the archive.
  Single list per archive where first column is the frequency channel, second column is the starting time and third column is the ending time.
  None value in time mean all values.
  '''
  
  if os.path.basename(archive_name) == 'puppi_57638_C0531+33_1218_2797.Tp':
    zap_list = [\
[0,   None, None],
[27,  2800, 3000],
[99,  2200, 3200],
[100, 2000, 3000],
[159, 200 , 500 ],
[168, 2200, 3200],
[169, 2200, 3200],
[190, 0   , 100 ],
[191, 0   , 100 ],
[201, None, None],
[202, None, None],
[275, None, None],
[276, None, None],
[286, None, None],
[308, None, None],
[309, None, None],
[310, None, None],
[311, None, None],
[320, None, None],
[321, None, None],
[323, None, None],
[324, None, None],
[334, None, None],
[335, None, None],
[343, None, None],
[344, None, None],
[353, None, None]]
    
  elif os.path.basename(archive_name) == 'puppi_.Tp':
    zap_list = [\
[],
[]]
    
  else:
    print "Archive not known. It will not be zapped. Select bins to zap out if you wish."
    zap_list = []
    
  return zap_list 


if __name__ == '__main__':
  main()
