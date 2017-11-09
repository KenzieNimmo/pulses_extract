import psrchive
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.optimize import brent, brentq, least_squares
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
import numpy.ma as ma

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


def load_ar(ar, use_psrchive_baseline=False):
  load_archive = psrchive.Archive_load(ar)
  load_archive.tscrunch()
  load_archive.convert_state('Stokes')
  load_archive.dedisperse()
  if use_psrchive_baseline: load_archive.remove_baseline()
  ds = ma.asarray(load_archive.get_data().squeeze())
  w = load_archive.get_weights().flatten()
  idx = np.where(w==0)[0]
  ds = np.multiply(ds, w[np.newaxis,:,np.newaxis])
  ds[:,idx,:] = ma.masked
  if load_archive.get_bandwidth() < 0: ds = ds[:,::-1,:]  #Inverted band
  lim = ar_pars[path.basename(ar).split('.')[0]]
  err_ds = ds.copy()
  err_ds[lim[0] - ds.shape[-1] / 20 : lim[1] + ds.shape[-1] / 20] = ma.masked
  mean = err_ds.mean(axis=2)
  if not use_psrchive_baseline: ds -= mean[:,:,np.newaxis]
  err_ds -= mean[:,:,np.newaxis]
  err_ds = err_ds.std(axis=2)
  return ds, err_ds


def load_Q_U(ar, use_psrchive_baseline=True, rms_level=None):
  ds, err_ds = load_ar(ar)
  lim = ar_pars[path.basename(ar).split('.')[0]]
  ds = ds[:, :, lim[0]:lim[1]].sum(axis=2)
  err_ds *= np.sqrt(lim[1] - lim[0])
  if rms_level is not None:
    L = np.sqrt(ds[1]**2+ds[2]**2)
    err_L = np.sqrt((err_ds[1] * ds[1] / L)**2 + (err_ds[2] * ds[2] / L)**2)
    ds[:, L < rms_level * err_L] = ma.masked
  return ds[0], ds[1], ds[2], err_ds[0], err_ds[1], err_ds[2]


def load_freq(ar, fscrunch=None):
  load_archive = psrchive.Archive_load(ar)
  if fscrunch is not None: load_archive.fscrunch(fscrunch)
  I = load_archive.get_Integration(0)
  x = np.array([I.get_centre_frequency(i) for i in range(I.get_nchan())])
  x *= I.get_doppler_factor()
  if load_archive.get_bandwidth() < 0: x = x[::-1]
  return x


def get_dRM(ar):
  load_archive = psrchive.Archive_load(ar)
  f_min = load_archive.get_centre_frequency() - abs(load_archive.get_bandwidth()) / 2.
  return 2. * np.pi * (f_min * 1e6)**2 / cc.c.value**2 / 16.


def get_complex(arrA, arrB):
    return arrA + 1j * arrB


def RM_model(x, *p):
  exp = 2.
  ph0, RM = p
  y = np.exp(2j*RM_model_approx(x, *p))
  return y


def get_fit_goodness(RM, x, QU):
  p = [0., RM]
  QU_model = RM_model(x, *p)
  goodness = np.sum(QU / QU_model)
  return -np.absolute(goodness)


def get_fscrunch(fscrunch, arr, s_arr=None):
  if s_arr is None:
    scrunch = np.reshape(arr, [-1, fscrunch]).mean(axis=1)
    return scrunch
  else:
    mean, err = get_weighted_average(arr, s_arr, scrunch=fscrunch)
    return mean, err


def get_weighted_average(arr, s_arr, scrunch=None):
  w = s_arr**-2
  if scrunch is not None:
    w = np.reshape(w, [-1, fscrunch])
    arr = np.reshape(arr, [-1, fscrunch])
  mean = np.sum(w * arr, axis=-1) / np.sum(w, axis=-1)
  err = np.sum(w * (arr - mean)**2, axis=-1) / np.sum(w, axis=-1)
  return mean, err


def wrap_ang(ang):
  return np.mod(ang + np.pi/2., np.pi) - np.pi/2.


def RM_model_approx(x, *p):
  return p[1] * cc.c.value**2 / (x*1e6)**2 + p[0]


def get_RM_PA(ar, fscrunch=None, rms_level=None, plot_diagnostic=False):
  #Load quantities
  I, Q, U, err_I, err_Q, err_U = load_Q_U(ar)
  QU = get_complex(Q, U)
  err_QU = get_complex(err_Q, err_U)
  err_QU = np.absolute(err_QU) / np.sqrt(2.)
  x = load_freq(ar)

  #Brute-force RM
  RM_steps = np.arange(85000, 110000, get_dRM(ar))

  fit_goodness = np.zeros_like(RM_steps)
  for i,RM in enumerate(RM_steps):
    fit_goodness[i] = get_fit_goodness(RM, x, QU)

  if plot_diagnostic:
    #Diagnostic plot
    plt.plot(RM_steps, fit_goodness, 'k.')
    plt.xlabel('RM (rad/m2)')
    plt.ylabel('Fit goodness')
    plt.show()

  #Refined RM
  best_RM = fit_goodness.argmin()
  bracket = RM_steps[best_RM-1 : best_RM+2]
  RM = brent(get_fit_goodness, args=(x,QU), brack=bracket)
  print 'Best RM must be in the range [{:.0f}, {:.0f}]'.format(bracket[0], bracket[2])

  #PA
  p = [0., RM]
  QU_model = RM_model(x, *p)
  PA = wrap_ang(np.angle(np.sum(QU / QU_model)) / 2.)

  #
  p = [PA, RM]
  QU_model = RM_model(x, *p)
  defar = QU / QU_model
  err_defar = err_QU
  if fscrunch is not None:
    #Scrunch in frequency
    defar, err_defar = get_fscrunch(fscrunch, defar, err_defar)
    x = get_fscrunch(fscrunch, x)
  if rms_level is not None:
    #Thresholding
    L = np.absolute(QU)
    err_L = np.sqrt((err_Q * Q / L)**2 + (err_U * U / L)**2)
    if fscrunch is not None: L, err_L = get_fscrunch(fscrunch, L, err_L)
    idx = np.where(L > rms_level * err_L)[0]
    defar = defar[idx]
    err_defar = err_defar[idx]
    x = x[idx]
    I = I[idx]
    err_I = err_I[idx]
    QU = QU[idx]
    err_QU = err_QU[idx]
  
  #Estimate deviations and errors
  PA_f = wrap_ang(np.angle(defar) / 2.)  #Wrap on -pi/2,pi/2
  err_PA_f = err_defar / 2. / np.absolute(defar)
  p0 = [0., 0.]  #PA, RM
  par, pcov = curve_fit(RM_model_approx, x, PA_f, p0=p0, sigma=err_PA_f)
  err_par = np.sqrt(np.diag(pcov))
  model = RM_model_approx(x, *par)
  red_chi2 = np.mean(((model - PA_f) / err_PA_f)**2)
  print 'red_chi2 = {:.2f}'.format(red_chi2)

  if plot_diagnostic:
    #Diagnostic plot
    x_p = np.linspace(x.min(), x.max(), 1e4)
    plt.errorbar(x, np.rad2deg(PA_f), yerr=np.rad2deg(err_PA_f), fmt='ro')
    plt.plot(x_p, np.rad2deg(RM_model_approx(x_p, *par)), 'k-')
    plt.xlabel('Frequency (MHz)')
    plt.ylabel('PA (deg)')
    plt.show()

  RM += par[1]
  err_RM = err_par[1]
  PA += par[0]
  err_PA = err_par[0]

  #Polarization fracion
  LI, err_LI = get_pol_fract(PA, RM, x, QU, err_QU, I, err_I, ar)
  print "Pol. fraction = {:.2f} +- {:.2f}".format(LI, err_LI)
  return RM, err_RM, np.rad2deg(PA), np.rad2deg(err_PA)


def get_pol_fract_OLD(PA, RM, x, QU, err_QU, I, err_I, ar):
  p = [PA, RM]
  QU_model = RM_model(x, *p)
  defar = QU / QU_model
  err_defar = err_QU / QU_model

  load_archive = psrchive.Archive_load(ar)
  chan_width = abs(load_archive.get_bandwidth() / load_archive.get_nchan())
  pol_f = np.sinc(RM * cc.c.value**2 * (chan_width*1e6) / (x * 1e6)**3)
  defar /= pol_f
  err_defar /= pol_f

  LI_f =  np.absolute(defar) / I
  err_LI_f = np.sqrt((err_defar.real * defar.real / I / np.absolute(defar))**2 + (err_defar.imag * defar.imag / I / np.absolute(defar))**2 + (err_I * np.absolute(defar) / I**2)**2)

  LI, err_LI = get_weighted_average(LI_f, err_LI_f)
  return LI, err_LI


def get_pol_fract(PA, RM, x, QU, err_QU, I, err_I, ar):
  p = [PA, RM]
  QU_model = RM_model(x, *p)
  defar = QU / QU_model
  err_defar = err_QU
  #I_t, err_I_t = get_weighted_average(I, err_I)
  I_t = I.sum()
  err_I_t = np.sqrt(np.sum(err_I**2))

  load_archive = psrchive.Archive_load(ar)
  chan_width = abs(load_archive.get_bandwidth() / load_archive.get_nchan())
  pol_f = np.sinc(RM * cc.c.value**2 * (chan_width*1e6) / (x * 1e6)**3)
  defar /= pol_f
  err_defar /= pol_f

  #P, err_P = get_weighted_average(defar, err_defar)
  P = defar.sum()
  err_P = np.sqrt(np.sum(err_defar**2))
  LI = np.absolute(P) / I_t
  #err_LI = np.sqrt((err_P.real * P.real / I_t / np.absolute(P))**2 + (err_P.imag * P.imag / I_t / np.absolute(P))**2 + (err_I_t * np.absolute(P) / I_t**2)**2)
  err_LI = np.sqrt((err_P / I_t)**2 + (err_I_t * np.absolute(P) / I_t**2)**2)
  return LI, err_LI


def main():
  ar_list = glob('*.clean')
  ar_list = [val for val in ar_list if path.basename(val).split('.')[0] in ar_pars.keys()]

  for i,ar in enumerate(ar_list):
    print '\nFitting for', ar
    RM, err_RM, PA, err_PA = get_RM_PA(ar, rms_level=5, fscrunch=8)
    print "RM = {:.0f}+-{:.0f}, PA = {:.3f}+-{:.3f}".format(RM, err_RM, PA, err_PA)

  return



if __name__ == '__main__':
  main()



