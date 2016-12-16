import os
import argparse
from glob import glob

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import psrchive
import scipy.ndimage

mpl.rc('font',size=20)
mpl.rc('lines', linewidth=2)

def parser():
    # Command-line options
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                    description="Plot dynamic spectrum from multiple archives with 1 polarisation and 1 subintegration.")
    parser.add_argument('archives_list', help="Name of the psrchive files to plot.", nargs='+')
    parser.add_argument('-o', help="Plot name.", default='bursts.png')
    parser.add_argument('-show', help="Show the plot.", action='store_false')
    parser.add_argument('-save_fig', help="Save the plot.", action='store_true')
    parser.add_argument('-zap', help="Plot to manually zap bins out.", action='store_true')
    parser.add_argument('-ncols', help="Number of columns in the general plot.", default=1, type=int)
    parser.add_argument('-nrows', help="Number of rows in the general plot.", default=1, type=int)
    parser.add_argument('-t_scrunch', help="Time scrunch archives by this factor.", default=1., type=float)
    parser.add_argument('-f_scrunch', help="Frequency scrunch archives by this factor.", default=1., type=float)
    parser.add_argument('-time_window', help="Time window around the burst in ms.", default=False, type=float)
    parser.add_argument('-f_min', help="Minimum frequency to plot in GHz.", default=None, type=float)
    parser.add_argument('-f_max', help="Maximum frequency to plot in GHz.", default=None, type=float)
    parser.add_argument('-cmap', help="Colormap to use in the plots. Other useful: RdGy_r", default='Greys')
    parser.add_argument('-log_scale', help="Logaritmic color scale.", action='store_true')
    return parser.parse_args()



def main():
  #Define general variables
  args = parser()
  plot_grid = gridspec.GridSpec(args.nrows, args.ncols)  #Grid of burst plots
  fig = plt.figure(figsize=(4*8.27, 4*11.69))  #A4
  
  #Load archive list
  if len(args.archives_list) == 1: ar_list = glob(args.archives_list[0])
  else: ar_list = args.archives_list
  
  #Loop on each archive
  skip = 0
  for idx, archive_name in enumerate(ar_list):
    #Skip plots in the first row
    if idx / args.ncols == 0:
      plots_to_skip = args.nrows * args.ncols - len(ar_list) - 1
      if args.ncols - idx == plots_to_skip: skip += plots_to_skip
      
    #Load archive
    DS, extent = load_DS(archive_name)

    #Zap the archive
    if args.zap: extent = None
    zap(archive_name, DS)
    
    #Invert the band if flipped
    if extent[3] < extent[2]:
      temp = extent[3]
      extent[3] = extent[2]
      extent[2] = temp
      DS = np.flipud(DS)
    
    #Plot the archive
    if idx == args.nrows * args.ncols / 2:
      skip += 1
      plot_DM_curves(extent, plot_grid[idx], fig, fmin=args.f_min, fmax=args.f_max)

    idx += skip
    plot(DS, plot_grid[idx], fig, extent=extent, ncols=args.ncols, nrows=args.nrows, t_scrunch=args.t_scrunch, f_scrunch=args.f_scrunch,\
         index=idx, width=args.time_window, fmin=args.f_min, fmax=args.f_max, cmap=args.cmap, log_scale=args.log_scale)

  
  #General plot settings
  fig.tight_layout()
  fig.subplots_adjust(hspace=0.1, wspace=0.05)
  if args.show: plt.show()
  if args.save_fig: fig.savefig(args.o, papertype = 'a4', orientation = 'portrait', format = 'png')
  return
  
  
  
def plot_DM_curves(extent, subplot_spec, fig, fmin=None, fmax=None, width=False):
  plot_grid = gridspec.GridSpecFromSubplotSpec(2, 2, subplot_spec, wspace=0., hspace=0., height_ratios=[1,5], width_ratios=[5,1])
  ax = plt.Subplot(fig, plot_grid[2])

  if not fmin: fmin = extent[2]
  if not fmax: fmax = extent[3]
  f = np.linspace(fmin, fmax, 1000)
  
  def dt(f, dDM, fmin):
    return 4.14881e6 * ((f*1000)**-2 - (fmax*1000)**-2) * dDM
  
  t = dt(f, 0., fmin)
  ax.plot(t, f, 'k-')
  ax.annotate("0" , xy=(t[0] + 0.1, f[0] + 5e-3), horizontalalignment='left', verticalalignment='baseline')
  
  t = dt(f, +1., fmin)
  ax.plot(t, f, 'k--')
  ax.annotate("+1", xy=(t[0] + 0.1, f[0] + 5e-3), horizontalalignment='left', verticalalignment='baseline')
  
  t = dt(f, -1., fmin)
  ax.plot(t, f, 'k--')
  ax.annotate("-1", xy=(t[0] + 0.1, f[0] + 5e-3), horizontalalignment='left', verticalalignment='baseline')
  
  t = dt(f, +2., fmin)
  ax.plot(t, f, 'k-.')
  ax.annotate("+2", xy=(t[0] + 0.1, f[0] + 5e-3), horizontalalignment='left', verticalalignment='baseline')
  
  t = dt(f, -2., fmin)
  ax.plot(t, f, 'k-.')
  ax.annotate("-2", xy=(t[0] + 0.1, f[0] + 5e-3), horizontalalignment='left', verticalalignment='baseline')
  
  if width: ax.set_xlim(-width/2., width/2.)
  ax.set_ylim(fmin, fmax)
  fig.add_subplot(ax)
  return
  
  
def plot(DS, subplot_spec, fig, extent=None, ncols=1, nrows=1, t_scrunch=1., f_scrunch=1., index=None, width=False, fmin=None, fmax=None, cmap='Greys', log_scale=False):
  #Define subplots
  plot_grid = gridspec.GridSpecFromSubplotSpec(2, 2, subplot_spec, wspace=0., hspace=0., height_ratios=[1,5], width_ratios=[5,1])
  ax1 = plt.Subplot(fig, plot_grid[2])
  ax2 = plt.Subplot(fig, plot_grid[0], sharex=ax1)
  ax3 = plt.Subplot(fig, plot_grid[3], sharey=ax1)

  #Dynamic spectrum
  if extent:
    smooth_DS = scipy.ndimage.zoom(DS, (1./f_scrunch,1./t_scrunch), order=1)
    smooth_DS -= np.median(smooth_DS)
    smooth_DS /= smooth_DS.max()
    units = ("GHz", "ms")
  else: 
    smooth_DS = DS
    units = ("chan", "bin")

  if extent and not fmin: fmin = extent[2]
  if extent and not fmax: fmax = extent[3]
    
  #Applies temporal and frequency windows
  if extent:
    prof = np.mean(smooth_DS, axis=0)
    res_t = extent[1] / prof.size 
    peak_ms = float(prof.argmax()) * res_t
    extent[0] = - width / 2.
    extent[1] = width / 2.
    if width: smooth_DS = smooth_DS[:, int(prof.argmax() - np.ceil(width / 2. / res_t)) : int(prof.argmax() + np.ceil(width / 2. / res_t))]
    
    fmin_bin = int(np.floor((fmin - extent[2]) / (extent[3] - extent[2]) * smooth_DS.shape[0]))
    fmax_bin = int(np.ceil((fmax - extent[2]) / (extent[3] - extent[2]) * smooth_DS.shape[0]))
    smooth_DS = smooth_DS[fmin_bin:fmax_bin]
    extent[2] = fmin
    extent[3] = fmax
        
  else: extent = [0, smooth_DS.shape[1]-1, 0, smooth_DS.shape[0]-1]
  
  if log_scale:
    smooth_DS -= smooth_DS.min() + 1e-6
    smooth_DS /= smooth_DS.max()
    smooth_DS = np.log(smooth_DS)
    
  ax1.imshow(smooth_DS, cmap=cmap, origin='lower', aspect='auto', interpolation='nearest', extent=extent)
  
  if width: ax1.set_xlim(-width/2., width/2.)
  ax1.set_ylim(fmin, fmax)
  
  #Give labels only to edge plots
  if index % ncols == 0: ax1.set_ylabel("Frequency ({})".format(units[0]))
  else: ax1.tick_params(axis='y', labelleft='off')
  if (index < ncols * (nrows - 1)) and not width: ax1.tick_params(axis='x', labelbottom='off')
  else: ax1.set_xlabel("Time ({})".format(units[1]))
  
  #Pulse profile
  prof = np.mean(smooth_DS, axis=0)
  x = np.linspace(extent[0], extent[1], prof.size)
  ax2.plot(x, prof, 'k-')
  ax2.tick_params(axis='y', which='both', left='off', right='off', labelleft='off')
  ax2.tick_params(axis='x', labelbottom='off')
  if width: ax2.set_xlim(-width/2., width/2.)
  else: ax2.set_xlim(extent[0:2])
  
  #Baseline
  bl = np.mean(smooth_DS, axis=1)
  y = np.linspace(extent[3], extent[2], bl.size)
  ax3.plot(bl, y, 'k-')
  ax3.tick_params(axis='x', which='both', top='off', bottom='off', labelbottom='off')
  ax3.tick_params(axis='y', labelleft='off')
  ax3.set_ylim(fmin, fmax)
  
  fig.add_subplot(ax1)
  fig.add_subplot(ax2)
  fig.add_subplot(ax3)
  return 



def load_DS(archive_name):
  load_archive = psrchive.Archive_load(archive_name)
  load_archive.remove_baseline()
  DS = load_archive.get_data().squeeze()
  
  freq = load_archive.get_centre_frequency()
  bw = load_archive.get_bandwidth()
  f0 = (freq - bw / 2) / 1000
  f1 = (freq + bw / 2) / 1000
  duration = load_archive.integration_length() * 1000
  
  return DS, [0., duration, f0, f1]
   
   
   
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
  
  if os.path.basename(archive_name) == 'puppi_57649_C0531+33_0106_413.Tp':
    zap_list = [\
[99,  2700, 3700],
[103, 1700, 2700],
[275, None, None],
[276, None, None],
[323, None, None],
[324, None, None],
[334, None, None],
[0, None, None]
]
    
  elif os.path.basename(archive_name) == 'puppi_57644_C0531+33_0021_2461.Tp':
    zap_list = [\
[0, None, None],
[286, None, None],
[310, None, None],
[311, None, None],
[320,765, 1500],
[321, 0, 1500],
[332, None, None],
[341, None, None],
[342, None, None],
[353, None, None],
[323, None, None],
[324, None, None],
[343, None, None],
[344, None, None],
[287, None, None],
[510, None, None]
]

  elif os.path.basename(archive_name) == 'puppi_57638_C0531+33_1218_280.Tp':
    zap_list = [\
[102, 3750, None],
[170, None, None],
[180, None, None],
[275, None, None],
[276, None, None],
[286, None, None],
[310, 3750, None],
[311, 3750, None],
[320, 700, 1000],
[332, None, None],
[334, None, None],
[335, None, None],
[341, None, None],
[342, None, None],
[353, None, None],
[410, 2500, None],
[431, 2500, None],
[0, None, None]
]

  elif os.path.basename(archive_name) == 'puppi_57648_C0531+33_0048_821.Tp':
    zap_list = [\
[0  , None, None],
[1  , None, None],
[2  , None, None],
[3  , None, None],
[4  , None, None],
[5  , None, None],
[6  , None, None],
[7  , None, None],
[8  , None, None],
[9  , None, None],
[10 , None, None],
[11 , None, None],
[12 , None, None],
[13 , None, None],
[14 , None, None],
[15 , None, None],
[16 , None, None],
[17 , None, None],
[18 , None, None],
[19 , None, None],
[20 , None, None],
[27 , None, None],
[28 , None, None],
[29 , None, None],
[102, None, 2000],
[168, None, None],
[169, None, None],
[170, None, None],
[171, None, None],
[172, None, None],
[180, None, None],
[168, None, None],
[168, None, None],
[275, None, None],
[276, None, None],
[286, None, None],
[287, None, None],
[323, None, None],
[324, None, None],
[332, None, None],
[334, None, None],
[341, None, None],
[342, None, None],
[100, 1300, 2300],
[410, None, None]
]
  
  elif os.path.basename(archive_name) == 'puppi_57640_C0531+33_1274_1421.Tp':
    zap_list = [\
[100, None, 1000],
[275, None, None],
[276, None, None],
[310, None, None],
[311, None, None],
[320, None, None],
[321, None, None],
[334, None, None],
[335, None, None],
[343, None, 3000],
[344, None, 3000],
[353, None, None],
[0, None, None]
]
  
  elif os.path.basename(archive_name) == 'puppi_57641_C0531+33_1312_185.Tp':
    zap_list = [\
[178, None, None],
[275, None, None],
[276, None, None],
[286, None, None],
[287, None, None],
[311, 1000, 1750],
[320, 800, 1350],
[321, 800, 1200],
[334, None, None],
[335, None, None],
[343, 0, 1000],
[344, 0, 1000],
[0, None, None]
]

  elif os.path.basename(archive_name) == 'puppi_57641_C0531+33_1312_521.Tp':
    zap_list = [\
[102, None, 1500],
[103, None, None],
[236, None, None],
[237, None, None],
[275, None, None],
[276, None, None],
[310, None, 1000],
[311, None, 1000],
[334, None, None],
[335, None, None],
[343, None, None],
[344, None, 2000],
[353, None, None],
[0, None, None]
]

  elif os.path.basename(archive_name) == 'puppi_57642_C0531+33_1322_1965.Tp':
    zap_list = [\
[0, None, None],
[99, 3700, None],
[103, None, 500],
[113, None, None],
[114, None, None],
[115, None, None],
[116, None, None],
[159, None, None],
[170, None, None],
[171, None, None],
[180, None, None],
[275, None, None],
[276, None, None],
[286, None, None],
[287, None, None],
[311, 1250, 1700],
[323, None, None],
[324, None, None],
[327, None, None],
[332, None, None],
[334, None, None],
[335, None, None],
[341, None, None],
[342, None, None],
[353, None, None]
]

  elif os.path.basename(archive_name) == 'puppi_57648_C0531+33_0048_378.Tp':
    zap_list = [\
[0, None, None],
[28 , None, None],
[83 , None, None],
[123, None, None],
[124, None, None],
[125, None, None],
[152, None, None],
[170, None, None],
[171, None, None],
[180, None, None],
[207, None, None],
[208, None, None],
[235, None, None],
[236, None, None],
[237, None, None],
[275, None, None],
[276, None, None],
[277, None, None],
[286, None, None],
[332, None, None],
[333, None, None],
[334, None, None],
[335, None, None],
[341, None, None],
[342, None, None],
[323, None, None],
[343, None, None],
[113, None, None],
[324, None, None]
]

  elif os.path.basename(archive_name) == 'puppi_57642_C0531+33_1322_7699.Tp':
    zap_list = [\
[0, None, None],
[170, None, None],
[171, None, None],
[180, None, None],
[275, None, None],
[276, None, None],
[286, None, None],
[287, None, None],
[320, None, None],
[323, None, None],
[324, None, None],
[334, None, None],
[335, None, None],
[342, 3500, None],
[343, 3200, None],
[406, 3600, None],
[406, 3000, None],
[407, 2600, None],
[408, 2100, None],
[409, 1700, None],
[410, 1200, None],
[411, 800, None],
[412, None, None],
[417, 3200, None],
[418, 2800, None],
[419, 2500, None],
[420, 1800, None],
[421, 1200, None],
[494, None, None],
[495, None, None],
[496, None, None],
[509, None, None],
[510, None, None],
[511, None, None],
[115, None, None],
[422, 800, None]
]

  elif os.path.basename(archive_name) == 'puppi_57646_C0531+33_0085_2476.Tp':
    zap_list = [\
[0, None, None],
[320, None, None],
[343, None, None],
[344, None, None],
[320, None, None],
[343, None, None],
[344, None, None],
[428, None, None]
]

  elif os.path.basename(archive_name) == 'puppi_57646_C0531+33_0085_4275.Tp':
    zap_list = [\
[0, None, None],
[275, None, None],
[310, None, None],
[311, None, None],
[311, None, 600],
[343, 3200, 3600],
[344, 3000, 3400],
[353, None, None],
[320, None, None],
[321, None, 500]
]

  elif os.path.basename(archive_name) == 'puppi_57638_C0531+33_1218_2797.Tp':
    zap_list = [\
[0, None, None],
[27, 2400, 3000],
[99, 2000, 3100],
[100, 1850, 2900] ,
[159, 100, 300],
[168, 2000, 3200],
[169, 2000, 3200],
[201, 1200, 3000],
[202, 1200, 3000],
[275, None, None],
[276, None, None],
[286, None, None],
[310, None, 3000],
[311, None, 3500],
[320, None, None],
[321, None, None],
[323, None, None],
[324, None, None],
[334, None, None],
[343, None, None],
[344, 1500, 3000],
[353, None, None]
]

  else:
    print "Archive not known. It will not be zapped. Select bins to zap out if you wish."
    zap_list = []
    
  return zap_list 


if __name__ == '__main__':
  main()
