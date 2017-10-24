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

mpl.rcParams['font.size'] = 7
mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['axes.linewidth'] = 1
mpl.rcParams['ps.useafm'] = True
mpl.rcParams['pdf.use14corefonts'] = True
mpl.rcParams['text.usetex'] = True
mm_to_in = 0.0393701


ar_pars = {
  '02_puppi_57747_C0531+33_0558_1183': [2030,2190],
  '03_puppi_57747_C0531+33_0558_1202': [2195,2258],
  '06_puppi_57747_C0531+33_0558_3687': [1383,1424],
  '07_puppi_57747_C0531+33_0558_3688': [1760,1886],
  '08_puppi_57747_C0531+33_0558_3689': [2112,2400],
  '01_puppi_57747_C0531+33_0558_5':    [1672,1866],
  '15_puppi_57748_C0531+33_0594_1269': [1750,1917],
  '12_puppi_57748_C0531+33_0594_48':   [1200,1292],
  '13_puppi_57748_C0531+33_0594_49':   [2060,2238],
  '14_puppi_57748_C0531+33_0594_50':   [1790,1840],
  '16_puppi_57772_C0531+33_0007_2695': [870,1050],

  'pulse_81845749.calibP': [680,730],
  'pulse_82398320.calibP': [225,270],
  'pulse_81822537.calibP': [680,730],
  'pulse_82398320.calibP.RM': [225,270],
}



ar_prop = {
  '01_puppi_57747_C0531+33_0558_5':    [57747.1232313932, 102729, 6, 1.833, 0.059, 567.5],
  '02_puppi_57747_C0531+33_0558_1183': [57747.1308531025, 102759, 27, 1.702, 0.244, 572.0],
  '03_puppi_57747_C0531+33_0558_1202': [57747.1399375004, 102694, 15, 2.208, 0.137, 563.5],
  #'04_puppi_57747_C0531+33_0558_25437':[57747.1452406758,      0,  0,     0,     0, 567.0],
  #'05_puppi_57747_C0531+33_0558_3683': [57747.1481344927,      0,  0,     0,     0, 564.0],
  '06_puppi_57747_C0531+33_0558_3687': [57747.1539561231, 102721, 25, 1.885, 0.226, 559.7],
  '07_puppi_57747_C0531+33_0558_3688': [57747.1540104922, 102656, 22, 2.549, 0.218, 560.2],
  '08_puppi_57747_C0531+33_0558_3689': [57747.1594948018, 102671, 16, 2.308, 0.138, 580.0],
  #'09_puppi_57747_C0531+33_0558_369':  [57747.1600418535,      0,  0,     0,     0, 565.0],
  #'10_puppi_57747_C0531+33_0558_125':  [57747.1696345353,      0,  0,     0,     0, 561.0],
  #'11_puppi_57748_C0531+33_0594_2':    [57748.1193274896,      0,  0,     0,     0, 560.2],
  '12_puppi_57748_C0531+33_0594_48':   [57748.1472089520, 102469, 29, 2.599, 0.277, 563.0],
  '13_puppi_57748_C0531+33_0594_49':   [57748.1488993263, 102529, 11, 1.982, 0.088, 562.0],
  '14_puppi_57748_C0531+33_0594_50':   [57748.1512921576, 102500, 14, 2.223, 0.126, 562.0],
  '15_puppi_57748_C0531+33_0594_1269': [57748.1693814308, 102659, 18, 0.606, 0.177, 560.2],
  '16_puppi_57772_C0531+33_0007_2695': [57772.1236149924, 103028, 7, 2.160, 0.067, 565.0],
}


def load_ar(ar_name, RM=False):
    load_archive = psrchive.Archive_load(ar_name)
    if RM:
      load_archive.set_rotation_measure(RM)
      load_archive.defaraday()
    load_archive.tscrunch()
    load_archive.dedisperse()
    #load_archive.remove_baseline()
    band_inverted = load_archive.get_bandwidth() < 0
    ds = load_archive.get_data().squeeze()
    if band_inverted: ds = ds[:,::-1,:]
    lim = ar_pars[path.splitext(path.basename(ar_name))[0]]
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
    ds -= mean[:,:,np.newaxis]
    err_ds -= mean[:,:,np.newaxis]
    err_ds = np.nanstd(err_ds, axis=2)
    return ds, err_ds

def load_UQ_freq(ar_name, RM=False):
    ds, err_ds = load_ar(ar_name, RM=RM)
    lim = ar_pars[path.splitext(path.basename(ar_name))[0]]
    ds = np.sum(ds[:, :, lim[0]:lim[1]], axis=2)
    err_ds *= np.sqrt(lim[1] - lim[0])
    ds[:, ds[0] < 3 * err_ds[0]] = np.nan
    return ds[0], ds[1], ds[2], err_ds[0], err_ds[1], err_ds[2]

def load_UQ_time(ar_name, RM=False):
    ds, err_ds = load_ar(ar_name, RM=RM)
    lim = ar_pars[path.splitext(path.basename(ar_name))[0]]
    ds = np.nansum(ds[:, :, lim[0]:lim[1]], axis=1)
    err_ds = np.nanmean(err_ds, axis=-1)
    err_ds *= np.sqrt(ds.shape[1])
    ds[:, ds[0] < 3 * err_ds[0]] = np.nan
    return ds[0], ds[1], ds[2], err_ds[0], err_ds[1], err_ds[2]

def load_freq(ar_name):
    load_archive = psrchive.Archive_load(ar_name)
    I = load_archive.get_Integration(0)
    band_inverted = load_archive.get_bandwidth() < 0
    x = np.array([I.get_centre_frequency(i) for i in range(I.get_nchan())])
    x *= I.get_doppler_factor()
    if band_inverted: x = x[::-1]
    return x


def load_QU(x, Q, U, err_Q, err_U, check_NaN=True):
    L = np.sqrt(U**2+Q**2)
    err_L = np.sqrt((err_Q * Q / L)**2 + (err_U * U / L)**2)
    QL = Q / L
    UL = U / L
    err_QL = np.sqrt((err_Q / L)**2 + (err_L * Q / L**2)**2)
    err_UL = np.sqrt((err_U / L)**2 + (err_L * U / L**2)**2)
    QU = np.hstack([QL, UL])
    err_QU = np.hstack([err_QL, err_UL])
    if check_NaN:
      idx = np.where(~np.isnan(QU))[0]
      x = x[idx[idx<x.size]]
      QU = QU[idx]
      err_QU = err_QU[idx]
    return x, QU, err_QU


def cos(x, *p):
    if len(p) == 3:
      ph0, RM, exp = p
    elif len(p) == 2:
      exp = 2.
      ph0, RM = p
    return np.cos(RM*cc.c.value**2/(x*1e6)**exp *2 + ph0)

def sin(x, *p):
    if len(p) == 3:
      ph0, RM, exp = p                                    
    elif len(p) == 2:
      exp = 2.
      ph0, RM = p
    return np.sin(RM*cc.c.value**2/(x*1e6)**exp *2 + ph0)

def RM_fit(x, *p):
    if len(p) == 3:
      ph0, RM, exp = p                                    
    elif len(p) == 2:
      exp = 2.
      ph0, RM = p
    y = np.exp(1j*(RM*cc.c.value**2/(x*1e6)**exp *2 + ph0))
    return np.hstack([y.real, y.imag])


def multiple_plots():
  ar_list = glob('/data/FRB121102/analyses/C-band/*_puppi_*.DM2')
  ar_list = [val for val in ar_list if path.splitext(path.basename(val))[0] in ar_pars.keys()]

  for i,ar in enumerate(ar_list):
    fig, (ax1, ax2) = plt.subplots(2, sharey=True, sharex=True, figsize=(20,4))
 
    _, Q, U, _, err_Q, err_U = load_UQ_freq(ar)

    xQU = load_freq(ar)
    xQU, QU, err_QU = load_QU(xQU, Q, U, err_Q, err_U)

    ax1.errorbar(xQU, QU[:QU.size/2], yerr=err_QU[:QU.size/2], fmt='r.')
    ax2.errorbar(xQU, QU[QU.size/2:], yerr=err_QU[QU.size/2:], fmt='r.')

    fit_exp = False
    p0 = [0.,102500.,2.]
    if fit_exp:
      par, pcov = curve_fit(RM_fit, xQU, QU, p0=p0, sigma=err_QU, bounds=[(-np.inf,-np.inf,1),(np.inf,np.inf,3)])
    else:
      p0 = [0.,102500.]
      par, pcov = curve_fit(RM_fit, xQU, QU, p0=p0, sigma=err_QU)
    err_par = np.sqrt(np.diag(pcov))
    model = RM_fit(xQU, *par)
    chi2 = np.sum(((model-QU)/err_QU)**2)
    red_chi2 = np.mean(((model-QU)/err_QU)**2)

    x = np.linspace(xQU.min(), xQU.max(), 1e4)
    ax1.plot(x, cos(x, *par), 'k-')
    ax2.plot(x, sin(x, *par), 'k-')

    
    print "\nBurst {} - {}".format(i, ar)
    print "RM: {:.0f} +- {:.0f}, PA: {:.3f} +- {:.3f}".format(par[1], err_par[1], np.mod(par[0], np.pi), err_par[0])
    if len(p0) == 3: 
      print "exp: {} +- {} ".format(par[2], err_par[2])
      print par[2] - 3*err_par[2], par[2], par[2] + 3*err_par[2]
    print "red chi2:", red_chi2
    

    #print par[2] - 3*err_par[2], par[2], par[2] + 3*err_par[2]
    #print "{}: {:.0f}, {:.0f}, {:.3f}, {:.3f}".format(path.basename(ar)[:2],par[1], err_par[1], np.mod(par[0], np.pi), err_par[0])
    #print red_chi2

    ax1.set_xlim([4140, 4860])
    ax1.set_ylim([-2,2])
    ax1.set_ylabel("Q/L")
    ax2.set_ylabel("U/L")
    ax2.set_xlabel("Frequency (MHz)")
    ax1.yaxis.set_major_locator(MultipleLocator(1))
    ax2.yaxis.set_major_locator(MultipleLocator(1))

    fig.tight_layout()
    plt.show()

    
    fig, ax1 = plt.subplots(1, sharey=True, sharex=True, figsize=(8,8))
    _, Q, U, _, err_Q, err_U = load_UQ_freq(ar, RM=par[1])
    xQU = load_freq(ar)
    idx = np.where(~np.logical_or(np.isnan(Q), np.isnan(U)))[0]
    xQU = xQU[idx]
    U = U[idx]
    Q = Q[idx]
    err_U = err_U[idx]
    err_Q = err_Q[idx]
    ang = np.mod(np.rad2deg(np.arctan2(Q, U) / 2.), 180)
    ang = np.mod(ang + np.median(ang), 180)
    err_ang = np.rad2deg(np.sqrt((err_Q*U)**2+(err_U*Q)**2)/2./(Q**2+U**2))
    ax1.errorbar(xQU, ang, yerr=err_ang, ecolor='k', marker='None', capsize=0, linestyle='None')
    par, pcov = np.polyfit(xQU, ang, 1, w=1/err_ang, cov=True)
    err_par = np.sqrt(np.diag(pcov))
    print "Slope PA: {:.3f} +- {:.3f}".format(par[0], err_par[0])
    lin = np.poly1d(par)
    ax1.plot(xQU, lin(xQU), 'k-')
    ax1.set_xlim([4140, 4860])
    ax1.set_ylabel("PA (deg)")
    ax1.set_xlabel("Frequency (MHz)")
    fig.tight_layout()
    plt.show()

  return



def fit(idx):
  def RM_fit_global(x, *p):
    exp = 2.
    #ph0, RM = np.reshape(p, [2,-1])  #PA+RM
    ph0 = np.array([p[0]],np.newaxis); RM = np.array(p[1:])  #RM
    y = np.exp(1j*(RM[:,np.newaxis]*cc.c.value**2/(x*1e6)**exp *2 + ph0[:,np.newaxis]))
    y = np.dstack([y.real, y.imag])
    return y.flatten()[idx]
  return RM_fit_global

def fit2(idx):
  #None
  def RM_fit_global(x, *p):
    exp = 2.
    ph0, RM = p
    y = np.exp(1j*(RM*cc.c.value**2/(x*1e6)**exp *2 + ph0))
    y = np.dstack([y.real, y.imag])
    return y.flatten()[idx]
  return RM_fit_global

def global_fit():
  ar_list = glob('*_puppi_5774*.DM2')
  ar_list = [val for val in ar_list if path.splitext(path.basename(val))[0] in ar_pars.keys()]
  fig, (ax1, ax2) = plt.subplots(2, sharey=True, sharex=True, figsize=[183*mm_to_in,183/3.*mm_to_in])

  xQU = load_freq(ar_list[0])
  xsize = xQU.size
  n_ar = len(ar_list)
  QU = np.zeros([n_ar, xsize, 2])
  err_QU = np.zeros([n_ar, xsize, 2])
  for i,ar in enumerate(ar_list):
    _, Q, U, _, err_Q, err_U = load_UQ_freq(ar)
    _ , QU_t, err_QU_t = load_QU(xQU, Q, U, err_Q, err_U, check_NaN=False)
    ax1.plot(xQU, QU_t[:QU_t.size/2], '.', zorder=2, ms=5)
    ax2.plot(xQU, QU_t[QU_t.size/2:], '.', zorder=2, ms=5)
    QU[i, :, :] = QU_t.reshape([2, -1]).T
    err_QU[i, :, :] = err_QU_t.reshape([2, -1]).T

  y = QU.flatten()
  yerr = err_QU.flatten()
  x = np.tile(xQU, (QU.shape[0],1))
  idx = np.where(~np.isnan(y))[0]
  y = y[idx]
  yerr = yerr[idx]

  #p0 = [0.]*n_ar + [102500.]*n_ar  #PA+RM
  #p0 = [1.92] + [103500.]*n_ar  #RM
  p0 = [1.92, 102500.]  #None
  par, pcov = curve_fit(fit2(idx), x, y, p0=p0, sigma=yerr)#, maxfev=10000)
  err_par = np.sqrt(np.diag(pcov))
  model = fit(idx)(xQU, *par)
  chi2 = np.sum(((model-y)/yerr)**2)
  red_chi2 = np.mean(((model-y)/yerr)**2)

  #PA + RM
  #par = np.array(par).reshape([2,-1])
  #par[0] = np.rad2deg(np.mod(par[0], np.pi))
  #err_par = np.array(err_par).reshape([2,-1])
  #print "PA:\n", np.array_str(np.vstack([par[0],err_par[0]]).T, precision=2, suppress_small=True)
  #print "RM:\n", np.array_str(np.vstack([par[1],err_par[1]]).T, precision=0, suppress_small=True)

  #RM
  #print "PA: {:.2f} +- {:.2f}".format(np.rad2deg(np.mod(par[0], np.pi)), np.rad2deg(err_par[0])) 
  #print "RM:\n", np.array_str(np.vstack([par[1:],err_par[1:]]).T, precision=0, suppress_small=True)

  #None
  print "PA: {:.2f} +- {:.2f}".format(np.rad2deg(np.mod(par[0], np.pi)), np.rad2deg(err_par[0])) 
  print "RM: {:.0f} +- {:.0f}".format(par[1], err_par[1])


  print "chi2red:", red_chi2



  """


  QU = QL.flatten()
  UL = UL.flatten()
  err_QL = err_QL.flatten()
  err_UL = err_UL.flatten()

  x = np.linspace(4500.78125-400, 4500.78125+400, 1e4)
  p0 = [0.,102500]
  Q_par, pcov = curve_fit(cos, xQ, QL, p0=p0, sigma=err_QL)
  err_Q_par = np.sqrt(np.diag(pcov))
  U_par, pcov = curve_fit(sin, xU, UL, p0=p0, sigma=err_UL)
  err_U_par = np.sqrt(np.diag(pcov))
  ax1.plot(x, cos(x, *Q_par), '-', zorder=1, color='grey')
  ax2.plot(x, sin(x, *U_par), '-', zorder=1, color='grey')

  print "Q - RM: {:.0f} +- {:.0f}, PA: {:.3f} +- {:.3f}".format(Q_par[1], err_Q_par[1], Q_par[0], err_Q_par[0])
  print "U - RM: {:.0f} +- {:.0f}, PA: {:.3f} +- {:.3f}".format(U_par[1], err_U_par[1], U_par[0], err_U_par[0])

  for i in range(QL.shape[0]):
    ax3.errorbar(xQ[i], cos(xQ[i], *Q_par) - QL[i], yerr=err_QL[i], fmt='None', ecolor='r', capsize=0, zorder=2)
    ax3.errorbar(xU[i], sin(xU[i], *U_par) - UL[i], yerr=err_UL[i], fmt='None', ecolor='b', capsize=0, xorder=2)
  ax3.axhline(0, color='grey', zorder=1)


  #ax1.plot(x, cos(x, 0), 'k-')
  #ax2.plot(x, sin(x, 0), 'k-')

  ax1.set_xlim([4140, 4860])
  ax1.set_ylim([-1.99,1.99])
  ax1.set_ylabel("Q/L")
  ax2.set_ylabel("U/L")
  ax3.set_ylabel("res.")
  ax3.set_xlabel("Frequency (MHz)")
  ax1.yaxis.set_major_locator(MultipleLocator(1))
  ax2.yaxis.set_major_locator(MultipleLocator(1))

  #ax1.set_xlim([4200,4500])


  fig.subplots_adjust(hspace=0.0, left=0.06,right=.99,bottom=.15,top=.98)

  plt.savefig('UQ.pdf', format='pdf', dpi=300)

  """

  plt.show()
  return

  


def average(val, sig):
  val = np.array(val, dtype=float)
  sig = np.array(sig, dtype=float)
  avg = np.sum(val*sig**-2) / np.sum(sig**-2)
  std = np.sqrt(1. / np.sum(sig**-2))
  return avg, std


def RMevolution():
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
  mjd = 57747
  idx = [i for i,n in enumerate(ar_prop.keys()) if n.find(str(mjd))>-1]
  y_mjd, y_err_mjd = average(y[idx], y_err[idx], returned=True)
  ax.errorbar(0, y_mjd, yerr=y_err_mjd, fmt='r+', zorder=10)
  mjd = 57748
  idx = [i for i,n in enumerate(ar_prop.keys()) if n.find(str(mjd))>-1]
  y_mjd, y_err_mjd = average(y[idx], y_err[idx], returned=True)
  ax.errorbar(1, y_mjd, yerr=y_err_mjd, fmt='r+', zorder=10)
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


  fig.subplots_adjust(left=0.15,right=.98,bottom=.1,top=.98)
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


def LI_RM():
  fig, ax = plt.subplots(1, figsize=[89*mm_to_in,120*mm_to_in])
  RM_range = np.linspace(52500, 152500, 1e3)
  LI = np.load('/data/FRB121102/analyses/C-band/LI_RM.npy')
  d = LI.shape[0]
  c = ['b'] * 6 + ['r'] * 4 + ['g']
  for i,n in enumerate(LI):
    #ax.fill_between(RM_range, 0, n-np.median(n)+d, facecolor='w', edgecolor=c[i])
    ax.plot(RM_range, n-np.median(n)+d, color=c[i])
    d -= 0.3

  RM_avg = RM_range[np.argmax(LI, axis=1)].mean()
  print 'Average RM: ',RM_avg

  RM_avg = RM_range[np.argmax(LI[:6], axis=1)].mean()

  ax.axvline(RM_avg, c='grey',lw=2, zorder=0)

  ax.set_xlabel('RM (rad$\,$m$^{-2}$)')
  ax.set_ylabel('L/I')
  ax.set_ylim([d+0.2, LI.shape[0]+1-0.1])
  ax.set_xlim([RM_avg-15000, RM_avg+15000])
  ax.tick_params(axis='y', labelleft='off', left='off', right='off')
  ax.ticklabel_format(useOffset=False)

  fig.subplots_adjust(left=0.1,right=.98,bottom=.1,top=.98)
  plt.savefig('LI_RM.pdf', format='pdf', dpi=300)
  plt.show()
  return


def Faraday_spectrum():
  ar_list = glob('/data/FRB121102/analyses/C-band/*_puppi_*.calibP')
  ar_list = [val for val in ar_list if path.splitext(path.basename(val))[0] in ar_pars.keys()]

  fig, ax = plt.subplots(1, figsize=[89*mm_to_in,89*mm_to_in])

  d = len(ar_list)
  RM_range = np.linspace(52500, 152500, 1e3)
  LI = np.zeros([d, RM_range.size])
  err_LI = np.zeros_like(LI)
  for j,ar_name in enumerate(ar_list):
    print "Archive processing n.",j

    def load_LI(archive, rm):
      I, Q, U, err_I, err_Q, err_U = load_UQ_time(ar_name, RM=rm)
      L = np.sqrt(U**2 + Q**2)
      err_L = np.sqrt((err_Q*Q)**2 + (err_U*U)**2) / L
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

    ax.fill_between(RM_range, 0, LI[j]-np.median(LI[j])+d, facecolor='w', edgecolor='k')
    d -= 0.5

  np.save('LI_RM_onlyL', LI)

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

  
  ax.plot([55,225], 1.1927**2*np.array([102677.45,102677.45]), 'g-', lw=3)

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
  ar_list = glob('/data/FRB121102/analyses/C-band/*_puppi_*.DM2')
  ar_list = [val for val in ar_list if path.splitext(path.basename(val))[0] in ar_pars.keys()]

  """
  fig, (ax1, ax2, ax3) = plt.subplots(3, figsize=[89*mm_to_in,89*mm_to_in])
  x = np.linspace(4500.78125-400+1.5625/2, 4500.78125+400-1.5625/2, 512)
  xt_flag = 0
  obs_end = []
  c = ['blue', 'orange', 'green', 'red', 'purple', 'brown', 'pink', 'gray', 'olive', 'cyan', 'k']
  for i,ar_name in enumerate(ar_list):
    RM = ar_prop[path.splitext(path.basename(ar_name))[0]][1]
    I, Q, U, err_I, err_Q, err_U = load_UQ_freq(ar_name, RM=RM)

    L = np.sqrt(U**2+Q**2)
    err_L = np.sqrt((err_Q*Q.T)**2 + (err_U*U.T)**2).T / L
    LI = L/I
    err_LI = np.sqrt((err_L/I)**2+(err_I*L.T/I.T**2).T**2)
    ax1.plot(x, LI, c[i])

    ang = np.rad2deg(np.arctan2(Q, U) / 2.)
    err_ang = np.rad2deg(np.sqrt((err_Q*U)**2+(err_U*Q)**2)/2./(Q**2+U**2))
    ax2.errorbar(x, ang, yerr=err_ang, ecolor=c[i], marker='None', capsize=0, linestyle='None')

    I, Q, U, err_I, err_Q, err_U = load_UQ_time(ar_name, RM=RM)
    ang = np.rad2deg(np.arctan2(Q, U) / 2.)
    err_ang = np.rad2deg(np.sqrt((err_Q*U)**2+(err_U*Q)**2)/2./(Q**2+U**2))
    xt = np.arange(ang.size) + xt_flag
    #ax3.errorbar(xt, ang, yerr=err_ang, ecolor='k', fmt='None', capsize=0, lw=3, elinewidth=3)
    ax3.errorbar(xt, ang, yerr=err_ang, marker='None', ecolor=c[i], capsize=0, linestyle='None')

    xt_flag += ang.size
    if (i==5) or (i==9): obs_end.append(xt_flag)

  ax3.axvline(obs_end[0], color='r')
  ax3.axvline(obs_end[1], color='r')

  ax1.set_xlim([4500.78125-400+1.5625/2, 4500.78125+400-1.5625/2])
  #ax.set_ylim([0,2])
  ax1.set_ylabel("L/I")
  ax1.set_xlabel("Frequency (MHz)")
  ax2.set_ylabel("PA (deg)")
  ax2.set_xlabel("Frequency (MHz)")
  ax3.set_ylabel("PA (deg)")
  ax3.set_xlabel("Time (arb.)")
  ax3.tick_params(axis='x', labelbottom='off', bottom='off')

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


  """
  RM_PA[2, RM_PA[2] < np.pi/4] += np.pi / 2.
  RM_PA[2, RM_PA[0] < 102555] += np.pi / 2.
  #RM_PA[2, RM_PA[2] < 60] += np.pi / 2.
  RM_PA[2, RM_PA[0] > 102800] -= np.pi / 2.
  """

  ax1.errorbar(RM_PA[0], np.rad2deg(RM_PA[2]), xerr=RM_PA[1], yerr=np.rad2deg(RM_PA[3]), fmt='ko')
  ax1.set_ylabel("PA (deg)")
  ax1.set_xlabel("RM (rad/m2)")
  plt.savefig('PA_RM.pdf', format='pdf', dpi=300)
  plt.show()



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
  }

  def load_ar_table(ar_name, ar_list):
    load_archive = psrchive.Archive_load(ar_name)
    load_archive.dedisperse()
    load_archive.pscrunch()
    load_archive.remove_baseline()
    ds = load_archive.get_data().squeeze()
    w = load_archive.get_weights().squeeze()
    ds = ds * w[:, np.newaxis]
    prof = ds.sum(axis=0)
    err_lim = ar_list[path.splitext(path.basename(ar_name))[0]]
    err_prof = prof[err_lim[0]:err_lim[1]].copy()
    mean = err_prof.mean()
    err_prof -= mean
    prof -= mean
    std = err_prof.std()
    err_prof /= std
    prof /= std
    return prof


  ar_list = glob('/data/FRB121102/analyses/C-band/*_puppi_*.DM2')
  ar_list = [val for val in ar_list if path.splitext(path.basename(val))[0] in ar_pars_all.keys()]

  prof = np.zeros([len(ar_list), 4096])

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
    print n.max() * 25./8./(2*800*10.24)**.5

  print ''
  print 'F'
  for n in prof:
    print n.sum() * 25./8./(2*800*10.24)**.5


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
  ang = np.mod(np.arctan2(Q, U)/2., 2*np.pi) / np.pi 
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
  #multiple_plots() 
  #RMevolution()
  #LI_plot()
  #Faraday_spectrum()
  #LI_RM()
  #DM_RM_all()
  #PA_f()
  #fill_table()
  #IQUV_plot()
  global_fit()








