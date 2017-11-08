import psrchive
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from matplotlib.ticker import MultipleLocator
from glob import glob
from os import path
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.patches as patches
import astropy.constants as cc
from matplotlib.colors import hsv_to_rgb
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec
from matplotlib.backends.backend_pdf import PdfPages

mpl.rcParams['font.size'] = 7
mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['axes.linewidth'] = 1
mpl.rcParams['ps.useafm'] = True
mpl.rcParams['pdf.use14corefonts'] = True
mpl.rcParams['text.usetex'] = True
mm_to_in = 0.0393701


ar_pars = {
  '01_puppi_57747_C0531+33_0558_5':    [1672,1866],
  '02_puppi_57747_C0531+33_0558_1183': [2030,2190],
  '03_puppi_57747_C0531+33_0558_1202': [2195,2258],
  #'04_puppi_57747_C0531+33_0558_25437':[1490,1600],
  #'05_puppi_57747_C0531+33_0558_3683': [1550,1760],
  '06_puppi_57747_C0531+33_0558_3687': [1383,1424],
  '07_puppi_57747_C0531+33_0558_3688': [1760,1886],
  '08_puppi_57747_C0531+33_0558_3689': [2112,2400],
  #'09_puppi_57747_C0531+33_0558_3690': [1750,1975],
  #'10_puppi_57747_C0531+33_0558_12568':[970,1090],
  #'11_puppi_57748_C0531+33_0594_2':    [1970,2100],
  '12_puppi_57748_C0531+33_0594_48':   [1200,1292],
  '13_puppi_57748_C0531+33_0594_49':   [2060,2238],
  '14_puppi_57748_C0531+33_0594_50':   [1790,1840],
  '15_puppi_57748_C0531+33_0594_1269': [1750,1917],
  '16_puppi_57772_C0531+33_0007_2695': [870,1050],
  'puppi_57747_B0525+21_0555_0001': [600,650],
  'puppi_57748_B0525+21_0591_0001': [680,730],
  'puppi_57772_B0525+21_0005_0001': [220,270],
  '11A_16sec':   [200, 360],
  '11C_284sec':  [1345, 1410],
  '11D_323sec':  [1260, 1350],
  '11E_344sec':  [700, 800],
  '11F_356sec':  [1900, 2040],
  #'11G_580sec':  [0, 0],
  '11H_597sec':  [620, 680],
  '11I_691sec':  [1680, 1740],
  '11J_704sec':  [1960, 2035],
  '11K_769sec':  [660, 745],
  '11M_993sec':  [1020, 1040],
  '11N_1036sec': [1540, 1570],
  '11O_1142sec': [1572, 1625],
  '11Q_1454sec': [140, 170],
  '12A_104sec':  [1190, 1160],
  '12C_1526sec': [10, 120],
}



ar_prop = {
  '01_puppi_57747_C0531+33_0558_5':    [57747.1295648440, 102703, 4, 59.54, 0.9, 567.5],
  '02_puppi_57747_C0531+33_0558_1183': [57747.1308531433, 102703, 4, 59.54, 0.9, 572.0],
  '03_puppi_57747_C0531+33_0558_1202': [57747.1399375575, 102703, 4, 59.54, 0.9, 563.5],
  '04_puppi_57747_C0531+33_0558_25437':[57747.1515739801, 102703, 4, 59.54, 0.9, 567.0],
  '05_puppi_57747_C0531+33_0558_3683': [57747.1544675163, 102703, 4, 59.54, 0.9, 564.0],
  '06_puppi_57747_C0531+33_0558_3687': [57747.1602892694, 102703, 4, 59.54, 0.9, 559.7],
  '07_puppi_57747_C0531+33_0558_3688': [57747.1603436884, 102703, 4, 59.54, 0.9, 560.2],
  '08_puppi_57747_C0531+33_0558_3689': [57747.1658277078, 102703, 4, 59.54, 0.9, 580.0],
  '09_puppi_57747_C0531+33_0558_3690': [57747.1663749556, 102703, 4, 59.54, 0.9, 565.0],
  '10_puppi_57747_C0531+33_0558_12568':[57747.1759673770, 102703, 4, 59.54, 0.9, 561.0],
  '11_puppi_57748_C0531+33_0594_2':    [57748.1256436151, 102516, 4, 59.54, 0.9, 560.2],
  '12_puppi_57748_C0531+33_0594_48':   [57748.1535244755, 102516, 4, 59.54, 0.9, 563.0],
  '13_puppi_57748_C0531+33_0594_49':   [57748.1552149312, 102516, 4, 59.54, 0.9, 562.0],
  '14_puppi_57748_C0531+33_0594_50':   [57748.1576076747, 102516, 4, 59.54, 0.9, 562.0],
  '15_puppi_57748_C0531+33_0594_1269': [57748.1756968496, 102516, 4, 59.54, 0.9, 560.2],
  '16_puppi_57772_C0531+33_0007_2695': [57772.1290302205, 103035, 3, 59.54, 0.9, 565.0],
  '11D_323sec':  [57991.5801413263],
  '11H_597sec':  [57991.5833162103],
}


def load_ar(ar_name, RM=None, use_psrchive_baseline=True, fscrunch=None):
    load_archive = psrchive.Archive_load(ar_name)
    if RM is not None:
      I = load_archive.get_Integration(0)
      load_archive.set_rotation_measure(RM / I.get_doppler_factor()**2)
      load_archive.defaraday()
    load_archive.tscrunch()
    load_archive.convert_state('Stokes')
    load_archive.dedisperse()
    if fscrunch is not None: load_archive.fscrunch(fscrunch)
    if use_psrchive_baseline: load_archive.remove_baseline()
    band_inverted = load_archive.get_bandwidth() < 0
    ds = load_archive.get_data().squeeze()

    if RM is not None:
      Q = ds[1, :, :]
      U = ds[2, :, :]
      QU = Q+1.j*U
      QU /= RM_model(np.array([load_archive.get_centre_frequency()]), 0, RM)[0]
      ds[1, :, :] = QU.real
      ds[2, :, :] = QU.imag

    if band_inverted: ds = ds[:,::-1,:]
    lim = ar_pars[path.basename(ar_name).split('.')[0]]
    w = load_archive.get_weights().flatten()
    if band_inverted: w = w[::-1]
    idx = np.where(w==0)[0]
    ds = np.multiply(ds, w[np.newaxis,:,np.newaxis])
    ds[:,idx,:] = np.nan
    if lim[0] > ds.shape[-1]/2: err_lim = [ds.shape[-1]/20, lim[0]-ds.shape[-1]/10]
    else: err_lim = [lim[1]+ds.shape[-1]/10, ds.shape[-1]-ds.shape[-1]/20]
    #plt.plot(np.nansum(ds[0],axis=0), 'k')
    #plt.axvspan(err_lim[0], err_lim[1], color='r', alpha=.5)
    #plt.axvspan(lim[0], lim[1], color='g', alpha=.5)
    #plt.show()
    err_ds = ds[:,:,err_lim[0]:err_lim[1]].copy()
    mean = np.nanmean(err_ds, axis=2)
    if not use_psrchive_baseline: ds -= mean[:,:,np.newaxis]
    err_ds -= mean[:,:,np.newaxis]
    err_ds = np.nanstd(err_ds, axis=2)
    return ds, err_ds

def load_UQ_freq(ar_name, RM=False, use_psrchive_baseline=True, rms_level=3., fscrunch=None):
    ds, err_ds = load_ar(ar_name, RM=RM, fscrunch=fscrunch)
    lim = ar_pars[path.basename(ar_name).split('.')[0]]
    ds = np.sum(ds[:, :, lim[0]:lim[1]], axis=2)
    err_ds *= np.sqrt(lim[1] - lim[0])
    L = np.sqrt(ds[1]**2+ds[2]**2)
    err_L = np.sqrt((err_ds[1] * ds[1] / L)**2 + (err_ds[2] * ds[2] / L)**2)
    ds[:, L < rms_level * err_L] = np.nan
    return ds[0], ds[1], ds[2], err_ds[0], err_ds[1], err_ds[2]

def load_UQ_time(ar_name, RM=False, rms_level=3.):
    ds, err_ds = load_ar(ar_name, RM=RM)
    lim = ar_pars[path.basename(ar_name).split('.')[0]]
    ds = np.nansum(ds[:, :, lim[0]:lim[1]], axis=1)
    err_ds = np.nanmean(err_ds, axis=-1)
    err_ds *= np.sqrt(ds.shape[1])
    ds[:, ds[0] < rms_level * err_ds[0]] = np.nan
    return ds[0], ds[1], ds[2], err_ds[0], err_ds[1], err_ds[2]

def load_freq(ar_name, fscrunch=None):
    load_archive = psrchive.Archive_load(ar_name)
    if fscrunch is not None: load_archive.fscrunch(fscrunch)
    I = load_archive.get_Integration(0)
    band_inverted = load_archive.get_bandwidth() < 0
    x = np.array([I.get_centre_frequency(i) for i in range(I.get_nchan())])
    x *= I.get_doppler_factor()

    #x *= (1 + 0.19273)

    if band_inverted: x = x[::-1]
    return x


def load_QU(x, Q, U, err_Q, err_U, check_NaN=True, return_idx=False):
    L = np.sqrt(U**2+Q**2)
    err_L = np.sqrt((err_Q * Q / L)**2 + (err_U * U / L)**2)
    QL = Q / L
    UL = U / L
    err_QL = err_Q / L
    err_UL = err_U / L
    QU = np.hstack([QL, UL])
    err_QU = np.hstack([err_QL, err_UL])
    if check_NaN:
      idx = np.where(~np.isnan(QU))[0]
      idx_x = idx[idx<x.size]
      x = x[idx_x]
      QU = QU[idx]
      err_QU = err_QU[idx]
      if return_idx: return idx_x, x, QU, err_QU
    return x, QU, err_QU


def cos(x, *p):
  return RM_model(x, *p).real

def sin(x, *p):
  return RM_model(x, *p).imag

def RM_model(x, *p):
    exp = 2.
    ph0_deg, RM = p
    ph0_rad = np.deg2rad(ph0_deg)
    y = np.exp(2j*(RM*cc.c.value**2/(x*1e6)**exp + ph0_rad))
    return y

def model_func(x, *p):
  return complex_to_real(RM_model(x, *p))

def complex_to_real(A):
  return np.hstack([A.real, A.imag])

def real_to_complex(A):
  if A.shape[0]%2:
    raise ValueError("A has shape %s, first dimension is not even" % (A.shape,))
  return A[:A.shape[0]//2]+1.j*A[A.shape[0]//2:]

def wrap_ang(ang):
  return np.mod(ang + 90., 180.) - 90.

def multiple_plots(plot=False):
  #ar_list = glob('*_puppi_*.DM2.clean')
  #ar_list = [val for val in ar_list if path.basename(val).split('.')[0] in ar_pars.keys()]
  #RMinit = [[52.2,102727.1], [63.5, 102702.7], [60.1, 102701.9], [30.7, 102796.3], [41.1, 102135.2], [51.3, 102727.4], [76.6, 102641.3], [66.6, 102665.4], [69.2, 102513.7], [77.0, 102639.3], [88.7, 102406.3], [74.4, 102465.2], [62.3, 102502.4], [67.9, 102481.5], [52.2, 102531.2], [63.0, 103011.2]]


  #RMfit = np.load('LI_RM_Arecibo.npy')
  #RMrange = RMfit[-1]
  #RMfit = RMfit[:-1]
  #RMinit = RMrange[RMfit.argmax(axis=1)]

  #BL
  ar_list = glob('BL/*.4p.clean.f8')
  ar_list = [val for val in ar_list if path.basename(val).split('.')[0] in ar_pars.keys()]
  #RMinit = [[50.,93673.], [50.,93533.]]

  for i,ar in enumerate(ar_list):
    if plot: fig, (ax1, ax2) = plt.subplots(2, sharey=True, sharex=True, figsize=(20,4))
 
    _, Q, U, _, err_Q, err_U = load_UQ_freq(ar, rms_level=5)

    if Q[~np.isnan(Q)].size < 20: continue

    xQU = load_freq(ar)

    """Partial channels
    xQU = xQU[1::2]
    Q = Q[1::2]
    U = U[1::2]
    err_Q = err_Q[1::2]
    err_U = err_U[1::2]
    """

    idx, xQU, QU, err_QU = load_QU(xQU, Q, U, err_Q, err_U, return_idx=True)
    Q = Q[idx]
    U = U[idx]
    err_Q = err_Q[idx]
    err_U = err_U[idx]

    if plot:
      ax1.errorbar(xQU, QU[:QU.size/2], yerr=err_QU[:QU.size/2], fmt='r.')
      ax2.errorbar(xQU, QU[QU.size/2:], yerr=err_QU[QU.size/2:], fmt='r.')

    #p0 = [RMinit[i][0], RMinit[i][1]*(1 + 0.19273)**2]
    #p0 = RMinit[i]
    p0 = [60., 93500.]
    try:
      par, pcov = curve_fit(model_func, xQU, QU, p0=p0, sigma=err_QU)#, bounds=[[-np.inf, p0[1]-1],[np.inf,p0[1]+1]])
      err_par = np.sqrt(np.diag(pcov))
      model = model_func(xQU, *par)
      chi2 = np.sum(((model-QU)/err_QU)**2)
      red_chi2 = np.mean(((model-QU)/err_QU)**2)
    except (TypeError, ValueError) as e: continue

    x = np.linspace(xQU.min(), xQU.max(), 1e4)
    if plot:
      ax1.plot(x, cos(x, *par), 'k-')
      ax2.plot(x, sin(x, *par), 'k-')

    if plot:
      print "\nBurst {} - {}".format(i+1, ar)
      print "RM: {:.0f} +- {:.0f}, PA: {:.3f} +- {:.3f}".format(par[1], err_par[1], wrap_ang(par[0]), err_par[0])
      if len(p0) == 3: 
        print "exp: {} +- {} ".format(par[2], err_par[2])
        print par[2] - 3*err_par[2], par[2], par[2] + 3*err_par[2]
      print "red chi2:", red_chi2
      print "dof:", xQU.size
      print "offset initial conditions:", p0[1]-par[1]

    else:
      #print par[2] - 3*err_par[2], par[2], par[2] + 3*err_par[2]
      print "Burst {}: RM={:.0f}+-{:.0f}, PA={:.3f}+-{:.3f} - chi2={:.2f}, dof={}, offset={:.0f}".format(path.basename(ar)[:3],par[1], err_par[1], wrap_ang(par[0]), err_par[0], red_chi2, xQU.size, p0[1]-par[1])


    if plot:
      #ax1.set_xlim([4140, 4860])
      ax1.set_ylim([-2,2])
      ax1.set_ylabel("Q/L")
      ax2.set_ylabel("U/L")
      ax2.set_xlabel("Frequency (MHz)")
      ax1.yaxis.set_major_locator(MultipleLocator(1))
      ax2.yaxis.set_major_locator(MultipleLocator(1))
      fig.tight_layout()
      plt.show()
    
    
    if plot:
      #PA flatness 
      fig, ax1 = plt.subplots(1, sharey=True, sharex=True, figsize=(8,8))
      _, Q, U, _, err_Q, err_U = load_UQ_freq(ar, RM=par[1])
      #xQU = load_freq(ar)
      #idx = np.where(~np.logical_or(np.isnan(Q), np.isnan(U)))[0]
      #xQU = xQU[idx]
      U = U[idx]
      Q = Q[idx]
      err_U = err_U[idx]
      err_Q = err_Q[idx]
      ang = np.rad2deg(np.arctan2(U, Q) / 2.)
      ang = np.mod(ang, 180.)
      #ang = np.mod(ang + np.median(ang), 180)
      err_ang = np.rad2deg(np.sqrt((err_Q*U)**2+(err_U*Q)**2)/2./(Q**2+U**2))
      ax1.errorbar(xQU, ang, yerr=err_ang, ecolor='k', marker='None', capsize=0, linestyle='None')
      par, pcov = np.polyfit(xQU, ang, 1, w=1/err_ang, cov=True)
      err_par = np.sqrt(np.diag(pcov))
      print "Slope PA: {:.3f} +- {:.3f}".format(par[0], err_par[0])
      lin = np.poly1d(par)
      ax1.plot(xQU, lin(xQU), 'k-')
      #ax1.set_xlim([4140, 4860])
      ax1.set_ylabel("PA (deg)")
      ax1.set_xlabel("Frequency (MHz)")
      fig.tight_layout()
      plt.show()
    
  return



def multidim_model(x, *p):
    if len(p) == 2:
      exp = 2.
      ph0_deg, RM = p
    if len(p) == 4:
      exp = 2.
      ph0_deg, RM1, RM2, RM3 = p
    elif (len(p) == 5) & (p[-1] < 10):
      #exp
      ph0_deg, RM1, RM2, RM3, exp = p
      #exp = np.array([exp1,]*6+[exp2,]*4)
    elif (len(p) == 5) & (p[-1] >= 10):
      exp = 2.
      ph0_deg, RM1, RM2, RM3, RM4 = p

    x[x==0] = np.nan
    #RM = np.array([RM1,]*10+[RM2,]*5+[RM3,])
    #RM = np.array([RM1,]*10+[RM2,]*5+[RM3,]+[RM4,]*2)
    ph0_rad = np.deg2rad(ph0_deg)
    #y = np.exp(2j*(RM[:,np.newaxis]*cc.c.value**2/(x*1e6)**exp + ph0_rad))
    y = np.exp(2j*(RM*cc.c.value**2/(x*1e6)**exp + ph0_rad))
    y = np.dstack([y.real, y.imag]).flatten()
    return y[~np.isnan(y)]

def global_fit():
  #ar_list = glob('BL/*.4p')
  #ar_list = ['BL/11D_323sec.calib.4p', 'BL/11H_597sec.calib.4p']
  #ar_list = glob('*_puppi_*.DM2')
  #ar_list = glob('*_puppi_*.DM2.clean') + ['BL/11D_323sec.calib.4p.clean.f8', 'BL/11H_597sec.calib.4p.clean.f8']
  ar_list = glob('*.calib.4p.clean.f8')
  ar_list = [val for val in ar_list if path.basename(val).split('.')[0] in ar_pars.keys()]

  x_size = load_freq(ar_list[-1])
  xsize = x_size.size
  n_ar = len(ar_list)
  QU = np.zeros([n_ar, xsize, 2])
  x = np.zeros([n_ar, xsize])
  err_QU = np.zeros([n_ar, xsize, 2])
  for i,ar in enumerate(ar_list):
    x_t = load_freq(ar)
    xsize = x_t.size
    x[i, :xsize] = x_t
    _, Q, U, _, err_Q, err_U = load_UQ_freq(ar, rms_level=5)
    _ , QU_t, err_QU_t = load_QU(x[i], Q, U, err_Q, err_U, check_NaN=False)
    QU[i, :xsize, :] = QU_t.reshape([2, -1]).T
    err_QU[i, :xsize, :] = err_QU_t.reshape([2, -1]).T

  y = QU.flatten()
  yerr = err_QU.flatten()
  y[y==0] = np.nan
  idx = np.where(~np.isnan(y))[0]
  y = y[idx]
  yerr = yerr[idx]
  x = np.where(~np.isnan(QU[:,:,0]), x, 0)

  #p0 = [60., 102702., 102515., 103027., 93500.]
  #p0 = [60., 102702., 102515., 103027.]
  #p0 = [60., 102702., 102515., 103027., 2.]  #exp
  p0 = [40., 93600.] #BL
  par, pcov = curve_fit(multidim_model, x, y, p0=p0, sigma=yerr)
  err_par = np.sqrt(np.diag(pcov))
  model = multidim_model(x, *par)
  chi2 = np.sum(((model-y)/yerr)**2)
  red_chi2 = np.mean(((model-y)/yerr)**2)

  print "PA: {:.2f} +- {:.2f}".format(par[0], err_par[0])
  print "RM1: {:.0f} +- {:.0f}".format(par[1], err_par[1])
  #print "RM2: {:.0f} +- {:.0f}".format(par[2], err_par[2])
  #print "RM3: {:.0f} +- {:.0f}".format(par[3], err_par[3])
  #print "RM4: {:.0f} +- {:.0f}".format(par[4], err_par[4])

  print "chi2red:", red_chi2



  #""" Plot BL
  fig, (ax1, ax2, ax3) = plt.subplots(3, figsize=[183*mm_to_in,183/3.*mm_to_in], sharex=True)
  c = cm.rainbow(np.resize(np.linspace(0, 1, len(ar_list)), [2,(len(ar_list)+1)/2]).T.flatten())
  p = [par[0], par[1]]
  for i, [xn, yn, err_yn] in enumerate(zip(x, QU, err_QU)):
    ax1.plot(xn, yn[:,0], '.', zorder=5, markeredgewidth=0, ms=7, color=c[i])
    ax2.plot(xn, yn[:,1], '.', zorder=5, markeredgewidth=0, ms=7, color=c[i])

    model = RM_model(xn, *p)
    #model = np.vstack([model.real, model.imag]).T
    yn_c = yn[:,0] + 1j*yn[:,1]
    err_c = np.angle(yn_c / model, deg=True) / 2.
    err_a = np.rad2deg(np.sqrt(err_yn[:,0]**2 + err_yn[:,1]**2)) / 2.
    ax3.errorbar(xn, err_c, yerr=err_a, ecolor=c[i], marker='None', ms=3., color=c[i], capsize=0, linestyle='None')
    ax3.scatter(xn, err_c, marker='.', s=3., color='k', zorder=10)

  x_f = np.linspace(np.nanmin(x_size), np.nanmax(x_size), 1e4)
  y_f = RM_model(x_f, *p)
  ax1.plot(x_f, y_f.real, 'k-')
  ax2.plot(x_f, y_f.imag, 'k-')

  ax3.axhline(0, color='k')
  ax1.set_ylabel('Q/L')
  ax2.set_ylabel('U/L')
  ax3.set_ylabel('$\Delta$ (deg)')
  ax3.set_xlabel('Frequency (MHz)')
  ax1.set_xlim([4800, 8250])
  ax1.set_ylim([-1.5, 1.5])
  ax2.set_ylim([-1.5, 1.5])
  ax3.set_ylim([-30, 30])
  ax1.tick_params(axis='x', labelbottom='off')
  ax2.tick_params(axis='x', labelbottom='off')
  ax1.yaxis.set_major_locator(MultipleLocator(1))
  ax2.yaxis.set_major_locator(MultipleLocator(1))
  ax3.yaxis.set_major_locator(MultipleLocator(25))
  #"""


  """Plot MJD57747
  ar_list = ar_list[:10]
  x_size = load_freq(ar_list[0])
  xsize = x_size.size
  n_ar = len(ar_list)
  QU = np.zeros([n_ar, xsize, 2])
  x = np.zeros([n_ar, xsize])
  err_QU = np.zeros([n_ar, xsize, 2])

  fig, (ax1, ax2, ax3) = plt.subplots(3, figsize=[183*mm_to_in,183/3.*mm_to_in], sharex=True)
  c = cm.rainbow(np.resize(np.linspace(0, 1, 10), [2,10/2]).T.flatten())
  p = [par[0], par[1]]
  for i,ar in enumerate(ar_list):
    x = load_freq(ar)
    _, Q, U, _, err_Q, err_U = load_UQ_freq(ar, rms_level=5)
    x , QU, err_QU = load_QU(x, Q, U, err_Q, err_U)

    ax1.plot(x, QU[:QU.size/2], '.', zorder=5, markeredgewidth=0, ms=7, color=c[i])
    ax2.plot(x, QU[QU.size/2:], '.', zorder=5, markeredgewidth=0, ms=7, color=c[i])

    model = RM_model(x, *p)
    #model = np.vstack([model.real, model.imag]).T
    yn_c = QU[:QU.size/2] + 1j*QU[QU.size/2:]
    err_c = np.angle(yn_c / model, deg=True) / 2.
    err_a = np.rad2deg(np.sqrt(err_QU[:err_QU.size/2]**2 + err_QU[err_QU.size/2:]**2)) / 2.
    ax3.errorbar(x, err_c, yerr=err_a, ecolor=c[i], marker='None', ms=3., color=c[i], capsize=0, linestyle='None')
    ax3.scatter(x, err_c, marker='.', s=3., color='k', zorder=10)

  x_f = np.linspace(np.nanmin(x_size), np.nanmax(x_size), 1e4)
  y_f = RM_model(x_f, *p)
  ax1.plot(x_f, y_f.real, 'k-')
  ax2.plot(x_f, y_f.imag, 'k-')

  ax1.set_ylabel('Q/L')
  ax2.set_ylabel('U/L')
  ax3.set_ylabel('$\Delta$ (deg)')
  ax3.set_xlabel('Frequency (MHz)')
  ax1.set_xlim([4140, 4860])
  ax1.set_ylim([-1.5, 1.5])
  ax2.set_ylim([-1.5, 1.5])
  ax3.set_ylim([-30, 30])
  ax1.tick_params(axis='x', labelbottom='off')
  ax2.tick_params(axis='x', labelbottom='off')
  ax1.yaxis.set_major_locator(MultipleLocator(1))
  ax2.yaxis.set_major_locator(MultipleLocator(1))
  ax3.yaxis.set_major_locator(MultipleLocator(25))
  fig.subplots_adjust(left=0.07,right=.99,bottom=.15,top=.99, hspace=0.)
  """


  """temp
  err_y = np.angle(QU / model) / 2.
  ax3.erro

  x_f = np.linspace(np.nanmin(x_size), np.nanmax(x_size), 1e5)
  y_f = RM_model(x_f, *p)
  ax1.plot(x_f, y_f.real, 'grey')
  ax2.plot(x_f, y_f.imag, 'grey')
  """


  #Plots
  """pdf_plot
  store = 'UQ_global.pdf'
  with PdfPages(store) as pdf:
    for i, [xn, yn, err_yn] in enumerate(zip(x, QU, err_QU)):
      fig, (ax1, ax2, ax3) = plt.subplots(3, figsize=[20,6], sharex=True, sharey=True)

      x_f = np.linspace(np.nanmin(xn), np.nanmax(xn), 1e4)

      ax1.plot(xn, yn[:,0], 'ro', zorder=5, markeredgewidth=0)
      ax2.plot(xn, yn[:,1], 'bo', zorder=5, markeredgewidth=0)

      #ax1.errorbar(xn, yn[:,0], yerr=err_yn[:,0], ecolor='r', marker='None', capsize=0, linestyle='None')
      #ax2.errorbar(xn, yn[:,1], yerr=err_yn[:,1], ecolor='b', marker='None', capsize=0, linestyle='None')

      if i < 10: p = [par[0], par[1]]
      elif i < 15: p = [par[0], par[2]]
      else: p = [par[0], par[3]]

      y_f = RM_model(x_f, *p)
      ax1.plot(x_f, y_f.real, 'k-')
      ax2.plot(x_f, y_f.imag, 'k-')

      model = RM_model(xn, *p)
      ax3.errorbar(xn, yn[:,0]-model.real, yerr=err_yn[:,0], ecolor='r', marker='None', capsize=0, linestyle='None')
      ax3.errorbar(xn, yn[:,1]-model.imag, yerr=err_yn[:,1], ecolor='b', marker='None', capsize=0, linestyle='None')

      ax1.set_ylabel('Q/L')
      ax2.set_ylabel('U/L')
      ax3.set_ylabel('res.')
      ax3.set_xlabel('Frequency (MHz)')
      ax1.set_xlim([4140, 4860])
      ax1.set_ylim([-1.5, 1.5])
      ax1.tick_params(axis='x', labelbottom='off')
      ax2.tick_params(axis='x', labelbottom='off')
      ax1.yaxis.set_major_locator(MultipleLocator(.5))
      ax2.yaxis.set_major_locator(MultipleLocator(.5))
      ax3.yaxis.set_major_locator(MultipleLocator(.5))
      fig.subplots_adjust(hspace=0)
      pdf.savefig(bbox_inches='tight',dpi=200)
      plt.close()
  """

  """ gridspec
  fig = plt.figure(figsize=[183*mm_to_in,183*mm_to_in])
  plot_grid = gridspec.GridSpec(16, 1, wspace=0.1)

  for i, [xn, yn, err_yn] in enumerate(zip(x, QU, err_QU)):
    sub_grid = gridspec.GridSpecFromSubplotSpec(3, 1, plot_grid[i], hspace=0.)
    ax1 = plt.Subplot(fig, sub_grid[0])
    ax2 = plt.Subplot(fig, sub_grid[1], sharex=ax1, sharey=ax1)
    ax3 = plt.Subplot(fig, sub_grid[2], sharex=ax1)

    ax1.errorbar(xn, yn[:,0], yerr=err_yn[:,0], ecolor='r', marker='None', capsize=0, linestyle='None')
    ax2.errorbar(xn, yn[:,1], yerr=err_yn[:,1], ecolor='b', marker='None', capsize=0, linestyle='None')

    if i < 10: p = [par[0], par[1]]
    elif i < 15: p = [par[0], par[2]]
    else: p = [par[0], par[3]]

    x_f = np.linspace(np.nanmin(xn), np.nanmax(xn), 1e4)
    y_f = RM_model(x_f, *p)
    ax1.plot(x_f, y_f.real, 'k-')
    ax2.plot(x_f, y_f.imag, 'k-')

    model = RM_model(xn, *p)
    ax3.errorbar(xn, yn[:,0]-model.real, yerr=err_yn[:,0], ecolor='r', marker='None', capsize=0, linestyle='None')
    ax3.errorbar(xn, yn[:,1]-model.imag, yerr=err_yn[:,1], ecolor='b', marker='None', capsize=0, linestyle='None')

    fig.add_subplot(ax1)
    fig.add_subplot(ax2)
    fig.add_subplot(ax3)

  """


  """test
  x_f = np.linspace(np.nanmin(xn), np.nanmax(xn), 1e4)
  ph0 = par[0]
  RM = par[3]
  y_f = np.exp(1j*(RM*cc.c.value**2/(x_f*1e6)**2. *2 + ph0))
  Q = y_f.real
  U = y_f.imag
  ax1.plot(x_f, Q, 'k-')
  ax2.plot(x_f, U, 'k-')
  
  _, Q, U, _, err_Q, err_U = load_UQ_freq(ar)
  xQU = load_freq(ar)
  xQU, QU, err_QU = load_QU(xQU, Q, U, err_Q, err_U)
  ax1.errorbar(xQU, QU[:QU.size/2], yerr=err_QU[:QU.size/2], fmt='r.')
  ax2.errorbar(xQU, QU[QU.size/2:], yerr=err_QU[QU.size/2:], fmt='r.')
  p = [ph0, RM]
  model = model_func(xQU, *p)
  chi2 = np.sum(((model-QU)/err_QU)**2)
  red_chi2 = np.mean(((model-QU)/err_QU)**2)
  print red_chi2
  """


  

  #plt.savefig('UQ.pdf', format='pdf', dpi=300)

  plt.show()
  return



def fit_BL(x, *p):
    exp = 2.
    ph0_deg, RM1, RM2, RM3, RM4 = p
    x[x==0] = np.nan
    RM = np.array([RM1,]*10+[RM2,]*5+[RM3,]+[RM4,]*15)
    ph0_rad = np.deg2rad(ph0_deg)
    y = np.exp(2j*(RM[:,np.newaxis]*cc.c.value**2/(x*1e6)**exp + ph0_rad))
    y = np.dstack([y.real, y.imag]).flatten()
    return y[~np.isnan(y)]

def RM_BL():
  fig, (ax1, ax2) = plt.subplots(2, sharey=True, sharex=True, figsize=[183*mm_to_in,183/3.*mm_to_in])

  #Arecibo
  ar_list = glob('*_puppi_577*.DM2')
  ar_list = [val for val in ar_list if path.splitext(path.basename(val))[0] in ar_pars.keys()]

  x_size = load_freq(ar_list[0])
  xsize = x_size.size
  n_ar = len(ar_list)
  QU = np.zeros([n_ar, xsize, 2])
  err_QU = np.zeros([n_ar, xsize, 2])
  x1 = np.zeros([n_ar, xsize])
  for i,ar in enumerate(ar_list):
    x1[i] = load_freq(ar)
    _, Q, U, _, err_Q, err_U = load_UQ_freq(ar, rms_level=4)
    _ , QU_t, err_QU_t = load_QU(x1[i], Q, U, err_Q, err_U, check_NaN=False)
    #ax1.plot(xQU, QU_t[:QU_t.size/2], '.', zorder=2, ms=5)
    #ax2.plot(xQU, QU_t[QU_t.size/2:], '.', zorder=2, ms=5)
    QU[i, :, :] = QU_t.reshape([2, -1]).T
    err_QU[i, :, :] = err_QU_t.reshape([2, -1]).T

  y1 = QU.flatten()
  yerr1 = err_QU.flatten()
  x1 = np.where(~np.isnan(QU[:,:,0]), x1, 0)


  #BL
  ar_list = glob('BL/*.4p')
  ar_list = [val for val in ar_list if path.splitext(path.basename(val))[0] in ar_pars.keys()]

  x_size = load_freq(ar_list[0])
  xsize = x_size.size
  n_ar = len(ar_list)
  QU = np.zeros([n_ar, xsize, 2])
  err_QU = np.zeros([n_ar, xsize, 2])
  x2 = np.zeros([n_ar, xsize])
  for i,ar in enumerate(ar_list):
    x2[i] = load_freq(ar)
    _, Q, U, _, err_Q, err_U = load_UQ_freq(ar, rms_level=4)
    _ , QU_t, err_QU_t = load_QU(x2[i], Q, U, err_Q, err_U, check_NaN=False)
    #ax1.plot(x[i], QU_t[:QU_t.size/2], '.', zorder=2, ms=5)
    #ax2.plot(x[i], QU_t[QU_t.size/2:], '.', zorder=2, ms=5)
    QU[i, :, :] = QU_t.reshape([2, -1]).T
    err_QU[i, :, :] = err_QU_t.reshape([2, -1]).T

  y2 = QU.flatten()
  yerr2 = err_QU.flatten()
  x2 = np.where(~np.isnan(QU[:,:,0]), x2, 0)

  #Fit
  x1 = np.pad(x1, ((0, 0), (0, x2.shape[1] - x1.shape[1])), 'constant')
  x = np.vstack([x1, x2])
  y = np.hstack([y1, y2])
  yerr = np.hstack([yerr1, yerr2])
  y = np.hstack([y1, y2])
  y = np.hstack([y1, y2])

  idx = np.where(~np.isnan(y))[0]
  y = y[idx]
  yerr = yerr[idx]

  p0 = [0., 102702., 102515., 103027., 93500.]
  par, pcov = curve_fit(fit_BL, x, y, p0=p0, sigma=yerr)
  err_par = np.sqrt(np.diag(pcov))
  model = fit_BL(x, *par)
  chi2 = np.sum(((model-y)/yerr)**2)
  red_chi2 = np.mean(((model-y)/yerr)**2)

  print "PA: {:.2f} +- {:.2f}".format(wrap_ang(par[0]), err_par[0])
  print "RM1: {:.0f} +- {:.0f}".format(par[1], err_par[1])
  print "RM2: {:.0f} +- {:.0f}".format(par[2], err_par[2])
  print "RM3: {:.0f} +- {:.0f}".format(par[3], err_par[3])
  print "RM4: {:.0f} +- {:.0f}".format(par[4], err_par[4])

  print "chi2red:", red_chi2

  return



def average(val, sig):
  val = np.array(val, dtype=float)
  sig = np.array(sig, dtype=float)
  avg = np.sum(val*sig**-2) / np.sum(sig**-2)
  std = np.sqrt(1. / np.sum(sig**-2))
  return avg, std


def RMevolution():
  fig, axarr = plt.subplots(3, 4, sharey='row', sharex='col', figsize=[183*mm_to_in,183.*mm_to_in], gridspec_kw={'height_ratios': [1.,.5,1.], 'wspace': .1, 'hspace': .05})
  c = cm.rainbow(np.resize(np.linspace(0, 1, 14), [2,14/2]).T.flatten())

  """Load mjd from archive
  def load_x(ar_list):
    x = np.zeros(len(ar_list))
    for i,ar in enumerate(ar_list):
      load_archive = psrchive.Archive_load(ar)
      I = load_archive.get_first_Integration()

      x[i] = I.get_epoch().get_secs() + I.find_max_phase() * I.get_duration()
    #x -= x.min()
    x /= 60.
    return x
  """

  def load_x(ar_list):
    x = np.zeros(len(ar_list))
    for i,ar in enumerate(ar_list):
      mjd = ar_prop[path.basename(ar).split('.')[0]][0]
      mins = (mjd % 1) * 24*60
      x[i] = mins
    return x


  def plot(x, RM, err_RM, c, RM_glob, err_RM_glob, PA, err_PA, PA_glob, err_PA_glob, ax0, ax1):
    ax0.errorbar(x, RM, yerr=err_RM, ecolor=c, fmt='k.', ms=5., capsize=0, linestyle='None', lw=3., zorder=1)
    ax0.axhspan(RM_glob-err_RM_glob, RM_glob+err_RM_glob, color='r', zorder=0)
    ax1.errorbar(x, PA, yerr=err_PA, ecolor=c, fmt='k.', ms=5., capsize=0, linestyle='None', lw=3., zorder=1)
    ax1.axhspan(PA_glob-err_PA_glob, PA_glob+err_PA_glob, color='r', zorder=0)
    #ax0.set_xlim([-5, x.max()+5])
    return

  PA_glob = 59.5
  err_PA_glob = 0.9

  #MJD 57747
  ar_list = glob('*_puppi_57747*.DM2.clean')
  ar_list = [val for val in ar_list if path.basename(val).split('.')[0] in ar_pars.keys()]
  x = load_x(ar_list)
  RM = [102739, 102750, 102694, 102713, 102666, 102656,]
  err_RM = [7, 24, 12, 22, 23, 15,]
  RM_glob = 102703
  err_RM_glob = 4
  PA = [49.866, 50.719, 63.130, 56.300, 70.552, 69.624,]
  err_PA = [6.996, 24.457, 12.371, 22.097, 22.985, 14.609,]

  plot(x, RM, err_RM, c[:6], RM_glob, err_RM_glob, PA, err_PA, PA_glob, err_PA_glob, axarr[0,0], axarr[2,0])

  #MJD 57748
  ar_list = glob('*_puppi_57748*.DM2.clean')
  ar_list = [val for val in ar_list if path.basename(val).split('.')[0] in ar_pars.keys()]
  x = load_x(ar_list)
  RM = [102515, 102532, 102486, 102493,]
  err_RM = [25, 12, 13, 31,]
  RM_glob = 102516
  err_RM_glob = 4
  PA = [61.468, 55.974, 67.556, 63.751,]
  err_PA = [25.376, 12.108, 13.078, 31.402,]

  plot(x, RM, err_RM, c[6:10], RM_glob, err_RM_glob, PA, err_PA, PA_glob, err_PA_glob, axarr[0,1], axarr[2,1])

  #MJD 57772
  ar_list = glob('*_puppi_57772*.DM2.clean')
  ar_list = [val for val in ar_list if path.basename(val).split('.')[0] in ar_pars.keys()]
  x = load_x(ar_list)
  RM = [103016,]
  err_RM = [10,]
  RM_glob = 103035
  err_RM_glob = 3
  PA = [64.459,]
  err_PA = [9.715,]

  plot(x, RM, err_RM, c[10], RM_glob, err_RM_glob, PA, err_PA, PA_glob, err_PA_glob, axarr[0,2], axarr[2,2])

  #axarr[0,0].set_ylabel('RM (rad$\,$m$^{-2}$)')
  #axarr[2,0].set_ylabel('PA (deg)')
  fig.text(0.01, 0.65, 'RM (rad$\,$m$^{-2}$)', va='center', rotation='vertical')
  fig.text(0.01, 0.231, 'PA (deg)', va='center', rotation='vertical')
  for ax in axarr[2]:
    ax.set_xlabel('Minutes')
  axarr[0,0].set_title('MJD 57747')
  axarr[0,1].set_title('MJD 57748')
  axarr[0,2].set_title('MJD 57772')
  axarr[0,3].set_title('MJD 57991 (GBT)')
  axarr[0,0].ticklabel_format(useOffset=False)

  for ax in axarr[0]:
    # hide the spines between ax and ax2
    ax.spines['bottom'].set_visible(False)
    ax.xaxis.tick_top()
    ax.tick_params(labeltop='off')  # don't put tick labels at the top
  for ax2 in axarr[1]:
    ax2.spines['top'].set_visible(False)
    ax2.xaxis.tick_bottom()

  axarr[0,0].set_ylim([102750-300, 102750+350])
  axarr[1,0].set_ylim([93550-750/4, 93550+750/4])
  axarr[2,0].set_ylim([20, 100])
  axarr[0,0].yaxis.set_major_locator(MultipleLocator(100))
  axarr[1,0].yaxis.set_major_locator(MultipleLocator(100))
  axarr[2,0].set_xlim([181, 249])
  axarr[2,0].xaxis.set_major_locator(MultipleLocator(10))
  axarr[2,1].set_xlim([215, 259])
  axarr[2,1].xaxis.set_major_locator(MultipleLocator(10))
  axarr[2,2].set_xlim([175, 194])
  axarr[2,2].xaxis.set_major_locator(MultipleLocator(5))
  axarr[2,3].set_xlim([830, 845])
  axarr[2,3].xaxis.set_major_locator(MultipleLocator(5))


  d = .04  # how big to make the diagonal lines in axes coordinates
  for ax in axarr[0]:
    # arguments to pass to plot, just so we don't keep repeating them
    kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
    ax.plot((-d, +d), (-d/2., +d/2.), **kwargs)        # top-left diagonal
    ax.plot((1 - d, 1 + d), (-d/2., +d/2.), **kwargs)  # top-right diagonal
  for ax2 in axarr[1]:
    kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
    ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
    ax2.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal

  #GBT
  ar_list = ['BL/11D_323sec.calib.4p.clean.f8', 'BL/11H_597sec.calib.4p.clean.f8']
  ar_list = [val for val in ar_list if path.basename(val).split('.')[0] in ar_pars.keys()]
  x = load_x(ar_list)
  RM = [93557, 93504,]
  err_RM = [45, 31,]
  RM_glob = 93656
  err_RM_glob = 9
  PA = [70.028, 74.180,]
  err_PA = [4.698, 3.052]


  plot(x, RM, err_RM, c[11:], RM_glob, err_RM_glob, PA, err_PA, PA_glob, err_PA_glob, axarr[1,3], axarr[2,3])


  """
  y = np.zeros(len(ar_prop.keys()))
  y_err = np.zeros_like(y)
  for i,n in enumerate(ar_prop.values()):
    y[i], y_err[i] = n[1], n[2]

  x = np.zeros_like(y)
  for i,ar in enumerate(ar_prop.keys()):
    load_archive = psrchive.Archive_load(ar)
    x[i] = load_archive.get_first_Integration().get_epoch().in_days()

  x -= x.min()
  idx = np.argsort(x)
  x = x[idx]
  y = y[idx]
  y_err = y_err[idx]

  fig, ax = plt.subplots(figsize=(4,4))
  ax.errorbar(x.astype(int), y, yerr=y_err, fmt='None', ecolor='k', capsize=0, lw=3)
  ax.set_xlim([-1.5,27])
  ax.set_xlabel('Days')
  ax.set_ylabel('RM (rad$\,$m$^{-2}$)')
  ax.ticklabel_format(useOffset=False)
  """

  """
  mjd = 57747
  idx = [i for i,n in enumerate(ar_prop.keys()) if n.find(str(mjd))>-1]
  y_mjd, y_err_mjd = average(y[idx], y_err[idx], returned=True)
  ax.errorbar(0, y_mjd, yerr=y_err_mjd, fmt='r+', zorder=10)
  mjd = 57748
  idx = [i for i,n in enumerate(ar_prop.keys()) if n.find(str(mjd))>-1]
  y_mjd, y_err_mjd = average(y[idx], y_err[idx], returned=True)
  ax.errorbar(1, y_mjd, yerr=y_err_mjd, fmt='r+', zorder=10)
  """
  """
  idx = np.where(x.astype(int) == 0)[0]
  y_mjd, _ = average(y[idx], y_err[idx])
  y_err_mjd = np.std(y[idx])
  ax.errorbar(0, y_mjd, yerr=y_err_mjd, fmt='None', zorder=0, lw=6, capsize=0, ecolor='r')

  idx = np.where(x.astype(int) == 1)[0]
  y_mjd, _ = average(y[idx], y_err[idx])
  y_err_mjd = np.std(y[idx])
  ax.errorbar(1, y_mjd, yerr=y_err_mjd, fmt='None', zorder=0, lw=6, capsize=0, ecolor='r')

  
  #Inset
  axins = inset_axes(ax, width="70%", height="60%", loc=4)
  ax.add_patch(patches.Rectangle((-0.5, 102250), 1, 200, fill=False )) 
  idx = np.where(x.astype(int) == 0)[0]
  x = x[idx] * 24*60
  y = y[idx]
  y_err = y_err[idx]
  y_mjd, _ = average(y[idx], y_err[idx])
  y_err_mjd = np.std(y[idx])
  axins.errorbar(x, y, yerr=y_err, fmt='None', ecolor='k', capsize=0, lw=3)
  axins.errorbar(np.mean(x), y_mjd, yerr=y_err_mjd, fmt='None', zorder=0, lw=6, capsize=0, ecolor='r')
  axins.set_xlim(x.min()-5, x.max()+5)
  axins.set_ylim(y.min()-y_err.max()-10, y.max()+y_err.max()+10)
  axins.tick_params(axis='both', labelleft='off', left='off', right='off', labelbottom='off', top='on', bottom='off', labeltop='on')
  axins.set_xlabel('Minutes')
  axins.xaxis.set_label_position('top') 
  axins.ticklabel_format(useOffset=False)

  #mark_inset(ax, axins, loc1a=2, loc1b=1, loc2a=4, loc2b=4, fc="w", ec="k")
  

  #ticks = axins.get_xticks()*24*60
  #axins.set_xticklabels(ticks)
  """

  fig.subplots_adjust(left=0.09,right=.98, bottom=.05,top=.95, hspace=0.1, wspace=0.1)
  plt.savefig('RMevolution.pdf', format='pdf', dpi=300)
  plt.show()



  """
  x2 = x[:6]*24*60
  y2 = y[:6]
  y_err2 = y_err[:6]

  fig, ax = plt.subplots(figsize=(4,4))
  ax.errorbar(x2, y2, yerr=y_err2, fmt='k.')
  ax.set_xlim([-4,59])
  #ax.set_ylim([102200, 102500])
  ax.set_xlabel('Minutes')
  ax.set_ylabel('RM (rad$\,$m$^{-2}$)')
  ax.ticklabel_format(useOffset=False)
  p = np.polyfit(x2, y2, 1, w=1./y_err2)
  x_f = np.array([x2.min()-10, x2.max()+10])
  y_f = p[0] * x_f + p[1]
  ax.plot(x_f, y_f, 'r-')

  fig.subplots_adjust(left=0.15,right=.98,bottom=.1,top=.98)
  plt.savefig('RMevolution_zoom.pdf', format='pdf', dpi=300)
  plt.show()
  """

  return


def LI_plot():
  ar_list = glob('/data/FRB121102/analyses/C-band/*_puppi_*.calibP')
  ar_list = [val for val in ar_list if path.splitext(path.basename(val))[0] in ar_pars.keys()]

  fig, ax = plt.subplots(1, figsize=[89*mm_to_in,89*mm_to_in])

  n_ar = len(ar_list)
  I = np.zeros([n_ar,512])
  U = np.zeros([n_ar,512])
  Q = np.zeros([n_ar,512])
  err_I = np.zeros(n_ar)
  err_U = np.zeros(n_ar)
  err_Q = np.zeros(n_ar)
  i = 0
  for ar in ar_list:
    I[i], Q[i], U[i], err_I[i], err_Q[i], err_U[i] = load_UQ(ar)
    i += 1

  L = np.sqrt(U**2+Q**2)
  err_L = np.sqrt((err_Q*Q.T)**2 + (err_U*U.T)**2).T / L

  LI = L/I
  err_LI = np.sqrt((err_L/I)**2+(err_I*L.T/I.T**2).T**2) 

  x_plot = np.linspace(4500.78125-400, 4500.78125+400, LI.shape[1])
  x_fit = np.linspace(-400., +400., LI.shape[1])

  for i,n in enumerate(LI):
    ax.errorbar(x_plot, n, yerr=err_LI[i], fmt='None', capsize=0, zorder=2)

  #x = np.tile(np.linspace(4500.78125-400, 4500.78125+400, LI.shape[1]), LI.shape[1])
  x = np.tile(x_fit, LI.shape[1])
  LI = LI.flatten()
  err_LI = err_LI.flatten()
  idx = np.where(~np.isnan(LI))[0]
  x = x[idx]
  LI = LI[idx]
  err_LI = err_LI[idx]

  #lin = np.poly1d(np.polyfit(x, LI, 1, w=1/err_LI))
  p,pcov = np.polyfit(x, LI, 1, w=1/err_LI, cov=True)
  perr = np.sqrt(np.diag(pcov))
  lin = np.poly1d(p)

  print 'm:', p[0], '+-', perr[0]
  print 'q:', p[1], '+-', perr[1]

  #x_f = np.array([x_sing[0], x_sing[-1]])
  ax.plot(x_plot, lin(x_fit), color='k', zorder=3, lw=2)
  #ax.plot(x_plot, x_fit*(p[0]+perr[0])+(p[1]+perr[1]), 'k--', zorder=3)
  #ax.plot(x_plot, x_fit*(p[0]-perr[0])+(p[1]-perr[1]), 'k--', zorder=3)

  #squ = np.poly1d(np.polyfit(x, LI, 2, w=1/err_LI))
  #ax.plot(x_plot, squ(x_fit), color='k', zorder=3)

  ax.set_xlim([4140, 4860])
  ax.set_ylim([0,2])
  ax.set_ylabel("L/I")
  ax.set_xlabel("Frequency (MHz)")
  #ax.yaxis.set_major_locator(MultipleLocator(1))
  #ax.yaxis.set_major_locator(MultipleLocator(1))
  

  fig.subplots_adjust(left=.11,right=.98,bottom=.1,top=.98)
  plt.savefig('LI_f.pdf', format='pdf', dpi=300)
  plt.show()
  return


def Gauss_model(x, *p):
  a, b, c = p
  return a * np.exp(-(x - b)**2 / 2. / c**2)

def LI_RM():
  fig, ax = plt.subplots(1, figsize=[89*mm_to_in,120*mm_to_in])

  LI_AO = np.load('/data/FRB121102/analyses/C-band/LI_RM_Arecibo.npy')
  x = LI_AO[-1]
  LI_AO = LI_AO[:-1]

  LI_BL = np.load('/data/FRB121102/analyses/C-band/LI_RM_GBT.npy')
  LI_BL = LI_BL[:-1]

  y = np.vstack([LI_AO, LI_BL])

  y -= np.median(y, axis=1)[:, np.newaxis]
  y /= np.max(y, axis=1)[:, np.newaxis]

  c = iter(cm.rainbow(np.resize(np.linspace(0, 1, 6), [2, 6/2]).T.flatten()))
  l = iter(['57747 (Arecibo)', '57748 (Arecibo)', '57772 (Arecibo)', '57991 (GBT)'])
  for i,n in enumerate(y):
    if i == 0 or i == 10 or i == 15 or i == 16: 
      col = next(c)
      lab = next(l)
      ax.fill_between(x, -10, n - i/5., facecolor='w', edgecolor=col, label=lab)
    else:
      ax.fill_between(x, -10, n - i/5., facecolor='w', edgecolor=col)


  RM_avg = x[np.argmax(y[:10], axis=1)].mean()
  ax.axvline(RM_avg, c='grey',lw=2, zorder=10)

  ax.set_xlabel('RM (rad$\,$m$^{-2}$)')
  ax.set_ylabel('L (norm.)')
  ax.set_ylim([-i/5+.3, 1.05])
  ax.set_xlim([85000, 115000])
  ax.tick_params(axis='y', labelleft='off', left='off', right='off')
  ax.ticklabel_format(useOffset=False)

  ax.legend(loc=2)

  fig.subplots_adjust(left=0.1,right=.98,bottom=.1,top=.98)
  plt.savefig('LI_RM.pdf', format='pdf', dpi=300)
  plt.show()

  return


def Faraday_spectrum():
  ar_list = glob('/data/FRB121102/analyses/C-band/*.DM2')

  #ar_list = ['/data/FRB121102/analyses/C-band/BL/11D_323sec.calib.4p', '/data/FRB121102/analyses/C-band/BL/11H_597sec.calib.4p']


  ar_list = [val for val in ar_list if path.splitext(path.basename(val))[0] in ar_pars.keys()]

  #fig, ax = plt.subplots(1, figsize=[89*mm_to_in,89*mm_to_in])

  d = len(ar_list) + 1
  #RM_range = np.linspace(52500, 152500, 1e3)
  RM_range = np.linspace(102000, 104000, 1e3)
  LI = np.zeros([d, RM_range.size])
  LI[-1] = RM_range
  err_LI = np.zeros_like(LI)
  for j,ar_name in enumerate(ar_list):
    print "Archive processing n.",j

    def load_LI(archive, rm):
      I, Q, U, err_I, err_Q, err_U = load_UQ_time(ar_name, RM=rm)
      L = np.sqrt(U**2 + Q**2)
      err_L = np.sqrt((err_Q*Q)**2 + (err_U*U)**2) / L
      L /= err_L
      #lim = ar_pars[path.splitext(path.basename(ar_name))[0]]
      #err_L = np.std(L[err_lim])
      #L = np.nansum(L[lim[0]:lim[1]])
      #I = np.nansum(I[lim[0]:lim[1]])
      #err_L = np.sqrt(np.sum(err_L[lim[0]:lim[1]]**2))
      #err_I *= np.sqrt(I.size)
      #I = np.nansum(I)
      L = np.nansum(L)
      err_L = np.sqrt(np.nansum(err_L**2))
      #LI = L / I
      #err_LI = np.sqrt((err_L/I)**2+(err_I*L/I**2)**2)
      return L, err_L

    for i,rm in enumerate(RM_range):
      if (i*len(RM_range))%len(RM_range) == 0: print "Finished: {:.2f}%".format(float(i)/len(RM_range)*100)
      LI[j,i], err_LI[j,i] = load_LI(ar_name, rm)

    #ax.fill_between(RM_range, 0, LI[j]-np.median(LI[j])+d, facecolor='w', edgecolor='k')
    #d -= 0.5

  np.save('LI_RM_AO_zoom', LI)

  """
  #ax.plot(RM_range, LI, 'k.')
  #ax.errorbar(RM_range, LI, yerr=err_LI, fmt='k.', ecolor='k', capsize=0)
  ax.set_xlabel('RM (rad$\,$m$^{-2}$)')
  ax.set_ylabel('L')
  ax.set_ylim([d, len(ar_list)+1.1])
  ax.tick_params(axis='y', labelleft='off', left='off', right='off')
  ax.set_xlim([RM_range[0], RM_range[-1]])
  ax.ticklabel_format(useOffset=False)

  fig.subplots_adjust(left=0.1,right=.98,bottom=.1,top=.98)
  plt.savefig('LI_RM.pdf', format='pdf', dpi=300)
  #plt.show()
  """

  return



def DM_RM_all():
  fig, ax = plt.subplots(1, figsize=[89*mm_to_in,89*mm_to_in])
  DM_RM = np.loadtxt("ATNF.dat", usecols=[2,3], delimiter=';', skiprows=1).T
  PSR_list = np.loadtxt("ATNF.dat", usecols=[1], delimiter=';', skiprows=1, dtype='str')
  ax.plot(DM_RM[0], np.abs(DM_RM[1]), 'k.')

  magnetars = ['J1550-5418', 'J1809-1943', 'J1745-2900', '1622-4950']
  sorter = np.argsort(PSR_list)
  idx = sorter[np.searchsorted(PSR_list, magnetars, sorter=sorter)]
  for i in idx:
    ax.plot(DM_RM[0,i], np.abs(DM_RM[1,i]), 'ro', zorder=10)

  
  ax.plot([57,227], 1.1927**2*np.array([102677.45,102677.45]), 'g-', lw=3)

  ax.annotate('J1746-2849', xy=(1456-100,10104), ha='right', va='top')
  ax.annotate('J1746-2856', xy=(1168-100,13253), ha='right', va='bottom')
  ax.annotate('J1745-2900', xy=(1778-300,66080), ha='right', va='center')
  ax.annotate('FRB 121102', xy=(115,1.1927**2*102677.45-30000), ha='center', va='top')

  ax.set_xscale("log", nonposx='clip')
  ax.set_yscale("log", nonposx='clip')
  ax.set_xlim([2,2400])
  ax.set_ylim([1,4e5])
  ax.set_xlabel('DM (pc$\,$cm$^{-3}$)')

  ax.set_ylabel('RM (rad$\,$m$^{-2}$)')
  fig.subplots_adjust(left=0.15,right=.98,bottom=.15,top=.98)
  plt.savefig('DM_RM_all.pdf', format='pdf', dpi=300)
  plt.show()
  return


def PA_f():
  ar_list = glob('*_puppi_*.DM2.clean')
  ar_list = [val for val in ar_list if path.basename(val).split('.')[0] in ar_pars.keys()]

  fig, (ax1, ax2, ax3) = plt.subplots(3, figsize=[120*mm_to_in,120*mm_to_in])
  xt_flag = 0
  obs_end = []
  c = cm.rainbow(np.resize(np.linspace(0, 1, len(ar_list)), [2,(len(ar_list)+1)/2]).T.flatten())
  LI_tot = np.zeros(512/16)
  err_LI_tot = np.zeros(512/16)
  ang_tot = np.zeros(512/16)
  err_ang_tot = np.zeros(512/16)
  for i,ar in enumerate(ar_list):
    """
    I, Q, U, err_I, err_Q, err_U = load_UQ_freq(ar, rms_level=3)
    xQU = load_freq(ar)
    xQU, QU, err_QU = load_QU(xQU, Q, U, err_Q, err_U)
    p0 = [50.,102700.]
    par, pcov = curve_fit(model_func, xQU, QU, p0=p0, sigma=err_QU)
    """

    RM = ar_prop[path.basename(ar).split('.')[0]][1]

    x = load_freq(ar, fscrunch=16)
    I, Q, U, err_I, err_Q, err_U = load_UQ_freq(ar, RM=RM, rms_level=5*4, fscrunch=16)
    L = np.sqrt(U**2+Q**2)
    err_L = np.sqrt((err_Q*Q.T)**2 + (err_U*U.T)**2).T / L
    LI = L/I
    err_LI = np.sqrt((err_L/I)**2+(err_I*L.T/I.T**2).T**2)
    LI_tot += np.nan_to_num(LI / err_LI**2)
    err_LI_tot += np.nan_to_num(err_LI**-2)
    ax1.errorbar(x, LI, yerr=err_LI, ecolor=c[i], marker='None', capsize=0, linestyle='None', lw=1.5)

    ang = np.rad2deg(np.arctan2(U, Q) / 2.)
    ang = np.mod(ang, 180.)
    err_ang = np.rad2deg(np.sqrt((err_Q*U)**2+(err_U*Q)**2)/2./(Q**2+U**2))
    ang_tot += np.nan_to_num(ang / err_ang**2)
    err_ang_tot += np.nan_to_num(err_ang**-2)
    ax2.errorbar(x, ang, yerr=err_ang, ecolor=c[i], marker='None', capsize=0, linestyle='None', lw=1.5)

    I, Q, U, err_I, err_Q, err_U = load_UQ_time(ar, RM=RM, rms_level=5)
    ang = np.rad2deg(np.arctan2(U, Q) / 2.)
    err_ang = np.rad2deg(np.sqrt((err_Q*U)**2+(err_U*Q)**2)/2./(Q**2+U**2))
    xt = np.arange(ang.size) + xt_flag
    ax3.errorbar(xt*10.24, ang, yerr=err_ang, marker='None', ecolor=c[i], capsize=0, linestyle='None')

    xt_flag += ang.size
    #if (i==9) or (i==14): obs_end.append(xt_flag)
    if (i==5) or (i==9): obs_end.append(xt_flag)

  ax1.scatter(x, LI_tot/err_LI_tot, facecolors='none', edgecolors='k', zorder=10)
  ax2.scatter(x, ang_tot/err_ang_tot, facecolors='none', edgecolors='k', zorder=10)

  x_f = np.linspace(3500., 5000., 1e4)
  theta = 102650.*cc.c.value**2*1.5625*1e6/(x_f*1e6)**3
  ax1.plot(x_f, np.sinc(theta), 'k-', zorder=5)

  ax3.axvline(obs_end[0]*10.24, color='k', lw=2., ls='--')
  ax3.axvline(obs_end[1]*10.24, color='k', lw=2., ls='--')
  ax3.axhspan(59.9-0.7, 59.9+0.7, color='y', zorder=0, alpha=0.3)
  ax2.axhspan(59.9-0.7, 59.9+0.7, color='y', zorder=0, alpha=0.3)

  ax1.set_xlim([4500.78125-400+1.5625/2, 4500.78125+400-1.5625/2])
  ax1.set_ylim([0.7,1.2])
  ax1.set_ylabel("L/I")
  ax1.set_xlabel("Frequency (MHz)")
  ax2.set_ylim([55,65])
  ax2.set_ylabel("PA (deg)")
  ax2.set_xlabel("Frequency (MHz)")
  ax3.set_ylim([25,100])
  ax3.set_ylabel("PA (deg)")
  ax3.set_xlabel("Time ($\mu$s)")
  #ax3.tick_params(axis='x', labelbottom='off', bottom='off')

  #ax.yaxis.set_major_locator(MultipleLocator(1))
  #ax.yaxis.set_major_locator(MultipleLocator(1))

  fig.tight_layout()
  #fig.subplots_adjust(left=.11,right=.98,bottom=.1,top=.98)
  plt.savefig('PA.pdf', format='pdf', dpi=300)
  plt.show()
  
  """
  fig, ax1 = plt.subplots(1, figsize=[89*mm_to_in,89*mm_to_in])

  RM_PA = np.zeros([11,4])
  i = 0
  for n in ar_prop.values():
    if n[1] > 0:
      RM_PA[i] = n[1:5]
      i += 1

  RM_PA = RM_PA.T
  RM_PA[2] = np.mod(RM_PA[2], np.pi)
  
  #RM_PA[2, RM_PA[2] < np.pi/4] += np.pi / 2.
  #RM_PA[2, RM_PA[0] < 102555] += np.pi / 2.
  ##RM_PA[2, RM_PA[2] < 60] += np.pi / 2.
  #RM_PA[2, RM_PA[0] > 102800] -= np.pi / 2.
  
  ax1.errorbar(RM_PA[0], np.rad2deg(RM_PA[2]), xerr=RM_PA[1], yerr=np.rad2deg(RM_PA[3]), fmt='ko')
  ax1.set_ylabel("PA (deg)")
  ax1.set_xlabel("RM (rad/m2)")
  plt.savefig('PA_RM.pdf', format='pdf', dpi=300)
  plt.show()
  """

  return


def fill_table():
  ar_pars_all = {
    '01_puppi_57747_C0531+33_0558_5':    [2200,4096],
    '02_puppi_57747_C0531+33_0558_1183': [0,1700],
    '03_puppi_57747_C0531+33_0558_1202': [0,1700],
    '04_puppi_57747_C0531+33_0558_25437':[2200,4096],
    '05_puppi_57747_C0531+33_0558_3683': [2200,4096],
    '06_puppi_57747_C0531+33_0558_3687': [1800,4096],
    '07_puppi_57747_C0531+33_0558_3688': [2300,4096],
    '08_puppi_57747_C0531+33_0558_3689': [0,1700],
    '09_puppi_57747_C0531+33_0558_3690': [2500,4096],
    '10_puppi_57747_C0531+33_0558_12568':[1800,4096],
    '11_puppi_57748_C0531+33_0594_2':    [0,1600],
    '12_puppi_57748_C0531+33_0594_48':   [1800,4096],
    '13_puppi_57748_C0531+33_0594_49':   [0,1600],
    '14_puppi_57748_C0531+33_0594_50':   [2200,4096],
    '15_puppi_57748_C0531+33_0594_1269': [2400,4096],
    '16_puppi_57772_C0531+33_0007_2695': [1700,4096],
    '11D_323sec':  [0, 1100],
    '11H_597sec':  [900, 2048],
  }

  def load_ar_table(ar_name, ar_list):
    load_archive = psrchive.Archive_load(ar_name)
    load_archive.dedisperse()
    load_archive.pscrunch()
    load_archive.remove_baseline()
    #ds = load_archive.get_data().squeeze()
    #w = load_archive.get_weights().squeeze()
    #ds = ds * w[:, np.newaxis]
    #prof = ds.sum(axis=0)
    load_archive.fscrunch()
    prof = load_archive.get_data().squeeze()
    err_lim = ar_list[path.basename(ar_name).split('.')[0]]
    err_prof = prof[err_lim[0]:err_lim[1]].copy()
    mean = err_prof.mean()
    err_prof -= mean
    prof -= mean
    std = err_prof.std()
    err_prof /= std
    prof /= std
    return prof


  ar_list = glob('/data/FRB121102/analyses/C-band/*_puppi_*.DM2')
  prof = np.zeros([len(ar_list), 4096])

  #ar_list = [val for val in ar_list if path.splitext(path.basename(val))[0] in ar_pars_all.keys()]

  ar_list = glob('/data/FRB121102/analyses/C-band/BL/*.clean.f8')
  prof = np.zeros([len(ar_list), 2048])

  print ''
  print 'Archives'
  for i,ar_name in enumerate(ar_list):
    print ar_name
    prof[i] = load_ar_table(ar_name, ar_pars_all)

  print ''
  print 'S/N'
  for n in prof:
    print n.sum()

  

  print ''
  print 'S'
  for n in prof:
    print n.max() * 30./7./(2*800*10.24)**.5

  print ''
  print 'F'
  for n in prof:
    print n.sum() * 30./7./(2*800*10.24)**.5 *1e-3


def IQUV_plot():
  ar_name = '/data/FRB121102/analyses/C-band/13_puppi_57748_C0531+33_0594_49.DM2'
  #ar_name = '/data/FRB121102/analyses/C-band/01_puppi_57747_C0531+33_0558_5.DM2'

  ds, err_ds = load_ar(ar_name)
  #lim = ar_pars[path.splitext(path.basename(ar_name))[0]]
  argmax = np.nanargmax(np.nansum(ds[0], axis=0))
  ds_bord = 200
  ds = ds[:, :, argmax - ds_bord : argmax + ds_bord]

  Q,U = ds[1:3]
  #ang = (np.arctan2(Q, U) / np.pi + 1) / 2.
  ang = np.mod(np.arctan2(U, Q)/2., 2*np.pi) / np.pi 
  L = np.sqrt(Q**2 + U**2)
  cost = np.zeros_like(L) + 0.5
   
  idx = np.where(~np.isnan(np.max(L,axis=1)))[0]
  L = L[idx]
  cost = cost[idx]
  ang = ang[idx]
  

  L -= L.min()
  L /= L.max()
  #ang -= ang.min()
  #ang /= ang.max()

  HSV = np.dstack((ang, L, cost))
  RGB = hsv_to_rgb(HSV)

  plt.imshow(RGB, origin='lower', aspect='auto', interpolation='nearest', extent=[-ds_bord*10.24e-3, ds_bord*10.24e-3, 4500.78125-400, 4500.78125+400])
  plt.xlabel('Time (ms)')
  plt.ylabel('Frequency (MHz)')

  plt.savefig('QU_img.pdf', format='pdf', dpi=300)
  plt.show()

  return


if __name__ == '__main__':
  multiple_plots() 
  #RMevolution()
  #LI_plot()
  #Faraday_spectrum()
  #LI_RM()
  #DM_RM_all()
  #PA_f()
  #fill_table()
  #IQUV_plot()
  #global_fit()
  #RM_BL() 







