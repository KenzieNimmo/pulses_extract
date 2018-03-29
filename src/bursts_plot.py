import os
import argparse
from glob import glob

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import psrchive
import scipy.ndimage
from matplotlib.ticker import MultipleLocator
import astropy.constants as cc

mpl.rcParams['font.size'] = 7
mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['axes.linewidth'] = 1
mpl.rcParams['ps.useafm'] = True
mpl.rcParams['pdf.use14corefonts'] = True
mpl.rcParams['text.usetex'] = True

mm_to_in = 0.0393701

def parser():
  # Command-line options
  parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                  description="Plot dynamic spectrum from multiple archives with 1 polarisation and 1 subintegration.")
  parser.add_argument('archives_list', help="Name of the psrchive files to plot, scrunched in time and polarisation.", nargs='+')
  parser.add_argument('-o', help="Plot name.", default='bursts.png')
  parser.add_argument('-show', help="Show the plot.", action='store_false')
  parser.add_argument('-save_fig', help="Save the plot.", action='store_true')
  parser.add_argument('-zap', help="Plot to manually zap bins out.", action='store_true')
  parser.add_argument('-ncols', help="Number of columns in the general plot.", default=1, type=int)
  parser.add_argument('-nrows', help="Number of rows in the general plot.", default=1, type=int)
  parser.add_argument('-t_scrunch', help="Time scrunch archives by this factor.", default=1, type=int)
  parser.add_argument('-f_scrunch', help="Frequency scrunch archives by this factor.", default=1, type=int)
  parser.add_argument('-time_window', help="Time window around the burst in ms.", default=20., type=float)
  parser.add_argument('-f_min', help="Minimum frequency to plot in GHz (Arecibo L-band: 1.15).", default=None, type=float)
  parser.add_argument('-f_max', help="Maximum frequency to plot in GHz (Arecibo L-band: 1.73.", default=None, type=float)
  parser.add_argument('-cmap', help="Colormap to use in the plots. Other useful: RdGy_r", default='Greys')
  parser.add_argument('-log_scale', help="Logaritmic color scale.", action='store_true')
  parser.add_argument('-pol', help="Plot polarisation information.", action='store_true')
  parser.add_argument('-plot_DM_curves', help="Plot curves of DM sweeps.", action='store_true')
  parser.add_argument('-DM', help="DM value to de-disperse the bursts.", default=False, type=float)
  parser.add_argument('-no_spectra', help="Plot burst spectra.", action='store_false')
  parser.add_argument('-plot_PA', help="Plot Position Angle.", action='store_true')

  return parser.parse_args()



def main():
  #Define general variables
  args = parser()
  plot_grid = gridspec.GridSpec(args.nrows, args.ncols, wspace=0.1)  #Grid of burst plots
  #fig = plt.figure(figsize=[120*mm_to_in,120*mm_to_in])  #Nature  
  fig = plt.figure(figsize=[7,6])  #Nature extended 

  #Zapping mode
  if args.zap:
    args.t_scrunch = args.f_scrunch = False
    args.pol = args.time_window = False
  
  #Load archive list
  if len(args.archives_list) == 1: ar_list = glob(args.archives_list[0])
  else: ar_list = args.archives_list

  #labels = [0,6,13,'GBT-1','GBT-2']
  labels = np.arange(len(ar_list)+1)
  #Loop on each archive
  skip = 0
  for idx, archive_name in enumerate(ar_list):
    i = idx + 1
    #Skip plots in the first row
    if idx / args.ncols == 0:
      plots_to_skip = args.nrows * args.ncols - len(ar_list) - 1
      if args.ncols - idx == plots_to_skip: skip += plots_to_skip

    
    if os.path.basename(archive_name).startswith('puppi_57364_C0531+33_4998_129.02'): 
      t_scrunch = 16 * args.t_scrunch
      f_scrunch = 4 * args.f_scrunch
      pol = False
      DM_curve = (-2.5, -0.2, 2.2)
      width = 2 * args.time_window
    else: 
      t_scrunch = args.t_scrunch
      f_scrunch = args.f_scrunch
      pol = args.pol
      DM_curve = False
      width = args.time_window

    #Load archive
    DS, spectrum, ts, extent = load_DS(archive_name, t_scrunch=t_scrunch, f_scrunch=f_scrunch, DM=args.DM)
    components = burst_components(archive_name)
    if args.zap: extent = None
    
    #Invert the band if flipped
    if extent and (extent[3] < extent[2]):
      temp = extent[3]
      extent[3] = extent[2]
      extent[2] = temp
      DS = np.flipud(DS)
      spectrum = np.flipud(spectrum)
    
    #Plot the archive
    idx += skip
    
    if args.plot_DM_curves:
      if (idx > 0) and (idx == args.nrows * args.ncols / 2):
        plot_DM_curves(extent, plot_grid[idx], fig, fmin=args.f_min, fmax=args.f_max, width=args.time_window)
        skip += 1
        idx += 1
    
    plot(DS, spectrum, ts, extent, plot_grid[idx], fig, ncols=args.ncols, nrows=args.nrows,\
         index=idx, width=width, fmin=args.f_min, fmax=args.f_max, cmap=args.cmap, log_scale=args.log_scale, components=components,\
         zap=args.zap, pol=pol, t_scrunch=t_scrunch, f_scrunch=f_scrunch, burst_n=labels[i], DM_curve=DM_curve, plot_spectra=args.no_spectra,\
         plot_PA=args.plot_PA)

  
  #General plot settings
  #fig.tight_layout()
  #fig.subplots_adjust(hspace=0.1, wspace=0.05, left=0.1,right=.98,bottom=.1,top=.98)
  fig.subplots_adjust(hspace=0.1, wspace=0.05, left=0.07,right=.99,bottom=.05,top=.99)
  if args.show: plt.show()
  if args.save_fig: fig.savefig(args.o, format = 'eps', dpi=150)
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
  ax.tick_params(axis='x', labelbottom='off')
  ax.tick_params(axis='y', labelleft='off')
  fig.add_subplot(ax)
  return
  
  
def plot(DS, spectrum, ts, extent, subplot_spec, fig, ncols=1, nrows=1, t_scrunch=1., f_scrunch=1., index=None,\
    width=False, fmin=None, fmax=None, cmap='Greys', log_scale=False, components=None, zap=False, pol=False,\
    burst_n=False, DM_curve=False, plot_spectra=True, plot_PA=False):

  rows = 2
  if plot_PA: rows += 1
  cols = 1
  if plot_spectra: cols += 1

  plot_grid = gridspec.GridSpecFromSubplotSpec(rows, cols, subplot_spec, wspace=0., hspace=0., height_ratios=[1,]*(rows-1)+[4,], width_ratios=[5,]+[1,]*(cols-1))
  ax1 = plt.Subplot(fig, plot_grid[rows-1,0])
  ax2 = plt.Subplot(fig, plot_grid[rows-2,0], sharex=ax1)
  if plot_spectra: ax3 = plt.Subplot(fig, plot_grid[rows-1,1], sharey=ax1)
  if plot_PA: ax4 = plt.Subplot(fig, plot_grid[0,0], sharex=ax1)


  #Applies temporal and frequency windows
  if zap:
    fmin = fmax = None
    units = ("chan", "bin")
  else: units = ("GHz", "ms")

  if not zap and not fmin: fmin = extent[2]
  if not zap and not fmax: fmax = extent[3]
  
  if zap: extent = [0, DS.shape[1]-1, 0, DS.shape[0]-1]
  else:
    res_t = extent[1] / ts.shape[1]
    peak = ts[0].argmax()
    peak_ms = float(peak) * res_t
    extent[0] = - width / 2.
    extent[1] = width / 2.
    components_ms = components / t_scrunch * res_t
    components = (components / t_scrunch).astype(int)
    if width:
      center = np.ceil(width / 2. / res_t)
      t0 = int(np.clip(peak - center, 0, np.inf))
      t1 = int(peak + center)
      components_ms -= peak_ms
      components -= int(peak - center)
    else: t0 = t1 = None
    fmin_bin = int(np.floor((fmin - extent[2]) / (extent[3] - extent[2]) * spectrum.size))
    fmax_bin = int(np.ceil((fmax - extent[2]) / (extent[3] - extent[2]) * spectrum.size))
    extent[2] = fmin
    extent[3] = fmax
  
    DS = DS[fmin_bin : fmax_bin, t0 : t1]
    spectrum = spectrum[fmin_bin : fmax_bin]
    ts = ts[:, t0 : t1]
  
  if log_scale:
    DS -= DS.min() + 1e-6
    DS /= DS.max()
    DS = np.log(DS)
    
  ax1.imshow(DS, cmap=cmap, origin='lower', aspect='auto', interpolation='nearest', extent=extent)
  if DM_curve and not zap:
    f = np.linspace(extent[2], extent[3], 1000)
    t = 4.14881e6 * ((f*1000)**-2 - (fmax*1000)**-2) * DM_curve[0] + DM_curve[1]
    ax1.plot(t, f, 'w-')
    t = 4.14881e6 * ((f*1000)**-2 - (fmax*1000)**-2) * DM_curve[0] + DM_curve[2]
    ax1.plot(t, f, 'w-')
    
  if width: ax1.set_xlim(-width/2.-0.001, width/2.+0.001)
  ax1.set_ylim(fmin, fmax)
  
  #Give labels only to edge plots
  if index % ncols == 0: ax1.set_ylabel("Frequency ({})".format(units[0]))
  else: ax1.tick_params(axis='y', labelleft='off')
  if (index < ncols * (nrows - 1)) and width: ax1.tick_params(axis='x', labelbottom='off')
  else: ax1.set_xlabel("Time ({})".format(units[1]))
  ax1.yaxis.set_major_locator(MultipleLocator(.2))
  ax1.xaxis.set_major_locator(MultipleLocator(1))
 
  #Pulse profile
  x = np.linspace(extent[0], extent[1], ts.shape[1])
  ax2.plot(x, ts[0], 'k-')
  if pol:
    ax2.plot(x, ts[1], 'r-')
    ax2.plot(x, ts[2], 'b-')
  ax2.tick_params(axis='y', which='both', left='off', right='off', labelleft='off')
  ax2.tick_params(axis='x', labelbottom='off', top='off')
  #if width: ax2.set_xlim(-width/2., width/2.)
  #else: ax2.set_xlim(extent[0:2])
  y_range = ts[0].max() - ts[0].min()
  ax2.set_ylim(-y_range/4., y_range*6./5.)
  #ax2.set_yticks([ts[0].max(),])
  #ax2.annotate("{:.0f} mJy".format(ts[0].max()), xy=(0.05,0.5), xycoords='axes fraction')
  if burst_n: ax2.annotate("\#{}".format(burst_n), xy=(0.98,0.5), xycoords='axes fraction', ha='right')

  #Spectrum
  if plot_spectra:
    y = np.linspace(extent[2], extent[3], spectrum.size)
    ax3.plot(spectrum, y, 'k-')
    ax3.tick_params(axis='x', which='both', top='off', bottom='off', labelbottom='off')
    ax3.tick_params(axis='y', labelleft='off')
    ax3.set_ylim(fmin, fmax)
  
    #Plot components
    try:
      colors = ['b', 'r', 'g', 'c', 'm', 'y', 'lime', 'darkblue']
      for i in range(len(components)-1):
        ax2.axvspan(components_ms[i], components_ms[i+1], color=colors[i], ls='-', ymax=.1)
        bl_c = np.sum(DS[:, components[i] : components[i+1]], axis=1)
        ax3.plot(bl_c, y, c=colors[i], lw=1.)
    except UnboundLocalError: pass

  if plot_PA:
    PA = ts[3]
    #PA -= np.nanmean(PA)
    err_PA = ts[4]
    ax4.errorbar(x, np.rad2deg(PA), yerr=np.rad2deg(err_PA), fmt='None', ecolor='k', capsize=0, zorder=1)
    ax4.axhline(np.rad2deg(np.nanmedian(PA)), color='grey', zorder=0)
    ax4.tick_params(axis='x', which='both', top='on', bottom='on', labelbottom='off')
    #ax4.set_ylim([np.rad2deg(np.nanmean(PA)) - 25, np.rad2deg(np.nanmean(PA)) + 25])
    if index % ncols == 0: ax4.set_ylabel('PA (deg)')
    else: ax4.tick_params(axis='y', labelleft='off')
    ax4.yaxis.set_major_locator(MultipleLocator(15))

    ax4.set_ylim([40,90])

  fig.add_subplot(ax1)
  fig.add_subplot(ax2)
  if plot_spectra: fig.add_subplot(ax3)
  if plot_PA: fig.add_subplot(ax4)
  return 


def load_DS(archive_name, zap=False, t_scrunch=False, f_scrunch=False, DM=False):
  #Load archive
  load_archive = psrchive.Archive_load(archive_name)
  load_archive.tscrunch()
  if DM:
    load_archive.dededisperse()
    load_archive.update()
    load_archive.set_dispersion_measure(DM)
    load_archive.dedisperse()
  else: load_archive.dedisperse()

  #load_archive.fscrunch(2)
  if archive_name.startswith('BL'): load_archive.fscrunch(32)
  #load_archive.bscrunch(2)
  
  if load_archive.get_npol() > 1:
    pol_info = True
  else: pol_info = False

  for ch in load_zap_list(archive_name):
    load_archive.get_first_Integration().set_weight(ch,0)
  
  if pol_info: load_archive.convert_state('Stokes')

  if t_scrunch > 1: load_archive.bscrunch(t_scrunch)
  if f_scrunch > 1: load_archive.fscrunch(f_scrunch)
  load_archive.remove_baseline()

  w = load_archive.get_weights().squeeze()
  archive = load_archive.get_data().squeeze()


  if pol_info:
    def RM_model(x, *p):
      exp = 2.
      ph0_deg, RM = p
      ph0_rad = np.deg2rad(ph0_deg)
      y = np.exp(2j*(RM*cc.c.value**2/(x*1e6)**exp + ph0_rad))
      return y

    #Reference PA to infinity
    Q = archive[1, :, :]
    U = archive[2, :, :]
    QU = Q+1.j*U
    QU /= RM_model(np.array([load_archive.get_centre_frequency()]), 0, load_archive.get_rotation_measure())[0]
    archive[1, :, :] = QU.real
    archive[2, :, :] = QU.imag
  
  if pol_info:
    archive = np.multiply(w, np.rollaxis(archive,2))
    archive = np.rollaxis(np.rollaxis(archive, 1, 0), 2, 1)
  else: 
    archive = (w*archive.T).T

  #Load dynamic spectrum
  if pol_info: DS = archive[0]
  else: DS = archive
  
  #Load frequency spectrum
  spectrum = np.sum(DS, axis=1)
  
  #Load timeseries
  I = np.sum(DS, axis=0)
  if pol_info:
    err_ds = archive.copy()
    idx = np.where(w==0)[0]
    err_ds[:,idx,:] = np.nan
    mean = np.nanmedian(np.reshape(err_ds, [err_ds.shape[0], err_ds.shape[1]*err_ds.shape[2]]), axis=1)
    err_ds -= mean[:, np.newaxis, np.newaxis]
    for i,n in enumerate(err_ds):
      if np.abs(n.min()) > n.max():
        err_ds[i] *= -1
    err_I = err_ds[0].copy()
    std = np.nanstd(err_I)
    for i,n in enumerate(err_ds):
      if np.abs(n.min()) > n.max():
        err_ds[i, err_I > 3*std] = np.nan
    err_ds = np.nanstd(np.reshape(err_ds, [err_ds.shape[0], err_ds.shape[1]*err_ds.shape[2]]), axis=1)

    I = np.sum(archive[0], axis=0)
    err_I = err_ds[0] * np.sqrt(archive.shape[1])
    Q = np.sum(archive[1], axis=0)
    err_Q = err_ds[1] * np.sqrt(archive.shape[1])
    U = np.sum(archive[2], axis=0)
    err_U = err_ds[2] * np.sqrt(archive.shape[1])
    L = np.sqrt(U**2 + Q**2)
    #err_L = np.sqrt((err_Q*Q.T)**2 + (err_U*U.T)**2).T / L
    L -= np.median(L)
    V = np.sum(archive[3], axis=0)
    PA = np.arctan2(U, Q) / 2.
    PA[I < 5*err_I] = np.nan
    err_PA = np.sqrt((err_Q*U)**2+(err_U*Q)**2)/2./(Q**2+U**2)
  else: L = V = PA = np.zeros_like(I)
  ts = np.vstack((I, L, V, PA, err_PA))
  
  #Load extensions
  if zap: extent = None
  else:
    freq = load_archive.get_centre_frequency()
    bw = load_archive.get_bandwidth()
    f0 = (freq - bw / 2) / 1000
    f1 = (freq + bw / 2) / 1000
    duration = load_archive.integration_length() * 1000
    extent = [0., duration, f0, f1]
    
  return DS, spectrum, ts, extent
   

def von_mises(x, mu, k, A):
  #To open .m files from paas
  #mu, k and A are the three parameters in the file
  vm = np.exp(k*np.cos((x/x.size-mu)*2.*np.pi))/2./np.pi/np.i0(k)
  return vm / vm.max() * A
 
   
def zap_ar(archive_name, DS):
  zap_list = load_zap_list(archive_name)
  
  med = np.median(DS)
  for chan in zap_list:
    DS[chan[0], chan[1]:chan[2]] = med
  DS -= med
  DS /= DS.max()
  
  return
  
  

def burst_components(archive_name):
  if os.path.basename(archive_name).startswith('puppi_57649_C0531+33_0106_413'):
    components = [1910, 2012, 2120, 2307]
  
  elif os.path.basename(archive_name).startswith('puppi_57644_C0531+33_0021_2461'):
    components = [2100, 2700]
  
  elif os.path.basename(archive_name).startswith('puppi_57638_C0531+33_1218_280'):
    components = [1100, 1600]
  
  elif os.path.basename(archive_name).startswith('puppi_57648_C0531+33_0048_821'):
    components = [2295, 2434, 2522, 2800]
  
  elif os.path.basename(archive_name).startswith('puppi_57640_C0531+33_1274_1421'):
    components = [2100, 2600]
  
  elif os.path.basename(archive_name).startswith('puppi_57641_C0531+33_1312_185'):
    components = [1570, 1834, 1969, 2211, 2264, 2361, 2510, 2685]
  
  elif os.path.basename(archive_name).startswith('puppi_57641_C0531+33_1312_521'):
    components = [2368, 2471, 2530, 2750]
  
  elif os.path.basename(archive_name).startswith('puppi_57642_C0531+33_1322_1965'):
    components = [1847, 1944, 2032, 2124, 2240, 2360]
  
  elif os.path.basename(archive_name).startswith('puppi_57648_C0531+33_0048_378'):
    components = [1187, 1347, 1482]
  
  elif os.path.basename(archive_name).startswith('puppi_57642_C0531+33_1322_7699'):
    components = [2300, 2385, 2474, 2560, 2825]
  
  elif os.path.basename(archive_name).startswith('puppi_57646_C0531+33_0085_2476'):
    components = [1714, 1797, 1910, 2070, 2222]
  
  elif os.path.basename(archive_name).startswith('puppi_57646_C0531+33_0085_4275'):
    components = [626, 753, 890, 1194]
  
  elif os.path.basename(archive_name).startswith('puppi_57638_C0531+33_1218_2797'):
    components = [1601, 1659, 1735, 1860, 2075]
  
  else:
    print "Archive not known, components will not be identified."
    components = []
  
  return np.array(components)





  
  
def load_zap_list(archive_name):
  '''
  List of bins to zap in the archive.
  Single list per archive where first column is the frequency channel, second column is the starting time and third column is the ending time.
  None value in time mean all values.
  '''
  
  if os.path.basename(archive_name).startswith('puppi_57649_C0531+33_0106_413'):
    zap_list = [\
99,
103,
275,
276,
323,
324,
334,
0,
100,

]
    
  elif os.path.basename(archive_name).startswith('puppi_57644_C0531+33_0021_2461'):
    zap_list = [\
0,
286, 
310,
311,
320,
321,
332,
341,
342,
353,
323,
324,
343,
344,
287,
510,
146,
]

  elif os.path.basename(archive_name).startswith('puppi_57638_C0531+33_1218_280'):
    zap_list = [\
102,
170,
180,
275,
276,
286,
310,
311,
320,
332,
334,
335,
341,
342,
353,
410,
431,
0,
323,
]

  elif os.path.basename(archive_name).startswith('puppi_57648_C0531+33_0048_821'):
    zap_list = [\
0  ,
1  ,
2  ,
3  ,
4  ,
5  ,
6  ,
7  ,
8  ,
9  ,
10 ,
11 ,
12 ,
13 ,
14 ,
15 ,
16 ,
17 ,
18 ,
19 ,
20 ,
27 ,
28 ,
29 ,
102,
168,
169,
170,
171,
172,
180,
168,
168,
275,
276,
286,
287,
323,
324,
332,
334,
341,
342,
100,
410,
107,
320,
]
  
  elif os.path.basename(archive_name).startswith('puppi_57640_C0531+33_1274_1421'):
    zap_list = [\
100,
275,
276,
310,
311,
320,
321,
334,
335,
343,
344,
353,
0,
161, 
]
  
  elif os.path.basename(archive_name).startswith('puppi_57641_C0531+33_1312_185'):
    zap_list = [\
178,
275,
276,
286,
287,
311,
320,
321,
334,
335,
343,
344,
0,
]

  elif os.path.basename(archive_name).startswith('puppi_57641_C0531+33_1312_521'):
    zap_list = [\
102,
103,
236,
237,
275,
276,
310,
311,
334,
335,
343,
344,
353,
0,
161, 
]

  elif os.path.basename(archive_name).startswith('puppi_57642_C0531+33_1322_1965'):
    zap_list = [\
0,
99,
103,
113,
114,
115,
116,
159,
170,
171,
180,
275,
276,
286,
287,
311,
323,
324,
327,
332,
334,
335,
341,
342,
353,
338,
]

  elif os.path.basename(archive_name).startswith('puppi_57648_C0531+33_0048_378'):
    zap_list = [\
0, 
28 ,
83 ,
123,
124,
125,
152,
170,
171,
180,
207,
208,
235,
236,
237,
275,
276,
277,
286,
332,
333,
334,
335,
341,
342,
323,
343,
113,
324,
327,
56,
285,
291,
292,
]

  elif os.path.basename(archive_name).startswith('puppi_57642_C0531+33_1322_7699'):
    zap_list = [\
0, 
170,
171,
180,
275,
276,
286,
287,
320,
323,
324,
334,
335,
342,
343,
406,
406,
407,
408,
409,
410,
411,
412,
413,
417,
418,
419,
420,
421,
493,
494,
495,
496,
497,
509,
510,
511,
115,
422,
423,
424,
100,
279,
]

  elif os.path.basename(archive_name).startswith('puppi_57646_C0531+33_0085_2476'):
    zap_list = [\
0,
320,
343,
344,
320,
343,
344,
428,
]

  elif os.path.basename(archive_name).startswith('puppi_57646_C0531+33_0085_4275'):
    zap_list = [\
0, 
275,
310,
311,
311,
343,
344,
353,
320,
321,
286,
]

  elif os.path.basename(archive_name).startswith('puppi_57638_C0531+33_1218_2797'):
    zap_list = [\
0,
27,
99, 
100,
159,
168,
169,
201,
202,
275,
276,
286,
310,
311,
320,
321,
323,
324,
334,
343,
344,
353,
]
    
  elif os.path.basename(archive_name).startswith('puppi_57364_C0531+33_4998_129.02'):
    zap_list = [\
335,
334,
333,
323,
286,
275,
188,
324,
276,
278,
]

  else:
    print "Archive {} not known. It will not be zapped. Select bins to zap out if you wish.".format(os.path.basename(archive_name))
    zap_list = []
    
  return zap_list 


if __name__ == '__main__':
  main()
