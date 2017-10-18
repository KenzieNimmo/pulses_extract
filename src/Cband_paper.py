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


mpl.rcParams['font.size'] = 7
mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['axes.linewidth'] = 1
mpl.rcParams['ps.useafm'] = True
mpl.rcParams['pdf.use14corefonts'] = True
mpl.rcParams['text.usetex'] = True
mm_to_in = 0.0393701


ar_pars = {
  '02_puppi_57747_C0531+33_0558_1183.dm559.72': [2030,2190],
  '03_puppi_57747_C0531+33_0558_1202.dm559.72': [2195,2258],
  '06_puppi_57747_C0531+33_0558_3687.dm559.72': [1383,1424],
  '07_puppi_57747_C0531+33_0558_3688.dm559.72': [1760,1886],
  '08_puppi_57747_C0531+33_0558_3689.dm559.72': [2112,2400],
  '01_puppi_57747_C0531+33_0558_5.dm559.72':    [1672,1866],
  '15_puppi_57748_C0531+33_0594_1269.dm559.72': [1750,1917],
  '12_puppi_57748_C0531+33_0594_48.dm559.72':   [1200,1292],
  '13_puppi_57748_C0531+33_0594_49.dm559.72':   [2060,2238],
  '14_puppi_57748_C0531+33_0594_50.dm559.72':   [1790,1840],
  '16_puppi_57772_C0531+33_0007_2695.dm559.72': [870,1050],

  'pulse_81845749.calibP': [680,730],
  'pulse_82398320.calibP': [225,270],
  'pulse_81822537.calibP': [680,730],
  'pulse_82398320.calibP.RM': [225,270],
}



ar_prop = {
  '01_puppi_57747_C0531+33_0558_5.dm559.72':    [57747.1232313932, 102368,  8, 3.745, 0.069, 102376,  8, 3.680, 0.076, 567.5],
  '02_puppi_57747_C0531+33_0558_1183.dm559.72': [57747.1308531025, 102405, 40, 3.582, 0.374, 102383, 31, 3.764, 0.281, 572.0],
  '03_puppi_57747_C0531+33_0558_1202.dm559.72': [57747.1399375004, 102373, 18, 3.747, 0.162, 102329, 15, 4.166, 0.137, 563.5],
  '04_puppi_57747_C0531+33_0558_25437.dm559.72':[57747.1452406758, 102350, 100,    0,     0, 102350, 100,    0,     0, 567.0],
  '05_puppi_57747_C0531+33_0558_3683.dm559.72': [57747.1481344927, 102350, 100,    0,     0, 102350, 100,    0,     0, 564.0],
  '06_puppi_57747_C0531+33_0558_3687.dm559.72': [57747.1539561231, 102369, 32, 3.695, 0.293, 102356, 29, 3.852, 0.265, 559.7],
  '07_puppi_57747_C0531+33_0558_3688.dm559.72': [57747.1540104922, 102278, 20, 4.636, 0.201, 102316, 33, 4.263, 0.328, 560.2],
  '08_puppi_57747_C0531+33_0558_3689.dm559.72': [57747.1594948018, 102358, 23, 3.796, 0.193, 102307, 20, 4.217, 0.171, 580.0],
  '09_puppi_57747_C0531+33_0558_3690.dm559.72': [57747.1600418535, 102350, 100,    0,     0, 102350, 100,    0,     0, 565.0],
  '10_puppi_57747_C0531+33_0558_12568.dm559.72':[57747.1696345353, 102350, 100,    0,     0, 102350, 100,    0,     0, 561.0],
  '11_puppi_57748_C0531+33_0594_2.dm559.72':    [57748.1193274896, 102200, 100,    0,     0, 102200, 100,    0,     0, 560.2],
  '12_puppi_57748_C0531+33_0594_48.dm559.72':   [57748.1472089520, 102093, 36, 4.661, 0.355, 102168, 38, 3.947, 0.372, 563.0],
  '13_puppi_57748_C0531+33_0594_49.dm559.72':   [57748.1488993263, 102187, 15, 3.720, 0.122, 102205, 14, 3.579, 0.115, 562.0],
  '14_puppi_57748_C0531+33_0594_50.dm559.72':   [57748.1512921576, 102134, 17, 4.197, 0.159, 102176, 19, 3.796, 0.169, 562.0],
  '15_puppi_57748_C0531+33_0594_1269.dm559.72': [57748.1693814308, 102295, 24, 2.559, 0.240, 102225, 38, 3.206, 0.377, 560.2],
  '16_puppi_57772_C0531+33_0007_2695.dm559.72': [57772.1236149924, 102669, 10,-2.311, 0.094, 102633, 10,-1.979, 0.093, 565.0],
}


def load_ar(ar_name):
    load_archive = psrchive.Archive_load(ar_name)
    load_archive.remove_baseline()
    load_archive.dedisperse()
    ds = load_archive.get_data().squeeze()
    ds = ds[:,::-1,:]
    lim = ar_pars[path.splitext(path.basename(ar_name))[0]]
    w = load_archive.get_weights().flatten()[::-1]
    idx = np.where(w==0)[0]
    ds = np.multiply(w, np.rollaxis(ds,2))
    ds = np.rollaxis(np.rollaxis(ds, 1, 0), 2, 1)
    ds[:,idx,:] = np.nan
    if lim[0] > ds.shape[-1]/2: err_lim = [0, lim[0]-ds.shape[-1]/20]
    else: err_lim = [lim[1]+ds.shape[-1]/20, ds.shape[-1]]
    err_ds = ds[:,:,err_lim[0]:err_lim[1]].copy()
    mean = np.nanmean(np.reshape(err_ds, [err_ds.shape[0], err_ds.shape[1]*err_ds.shape[2]]), axis=1)
    err_ds -= mean[:, np.newaxis, np.newaxis]
    err_ds = np.nanstd(np.reshape(err_ds, [err_ds.shape[0], err_ds.shape[1]*err_ds.shape[2]]), axis=1)
    return ds, err_ds


def load_UQ(ar_name):
    ds, err_ds = load_ar(ar_name)
    lim = ar_pars[path.splitext(path.basename(ar_name))[0]]
    ds = np.nansum(ds[:, :, lim[0]:lim[1]], axis=2)
    err_ds *= np.sqrt(lim[1] - lim[0])
    ds[:, ds[0] < 5 * err_ds[0]] = np.nan
    return ds[0], ds[1], ds[2], err_ds[0], err_ds[1], err_ds[2]


def load_UQ_time(ar_name):
    ds, err_ds = load_ar(ar_name)
    lim = ar_pars[path.splitext(path.basename(ar_name))[0]]
    err_ds *= np.sqrt(ds.shape[1])
    ds = np.nansum(ds[:, :, lim[0]:lim[1]], axis=1)
    ds[:, ds[0] < 5 * err_ds[0]] = np.nan
    return ds[0], ds[1], ds[2], err_ds[0], err_ds[1], err_ds[2]


def cos(x, *p):
    ph0, RM = p                                      
    return np.cos(RM*3e8**2/(x*1e6)**2 * 2 + ph0)

def sin(x, *p):
    ph0, RM = p                                      
    return np.sin(RM*3e8**2/(x*1e6)**2 * 2 + ph0)

def RM_fit(x, *p):
  ph0, RM = p
  #ang = np.mod(2.*(RM*3e8**2/(x*1e6)**2+ph0), 2*np.pi)
  #ang[ang > np.pi] = np.pi - ang[ang > np.pi]
  #return ang
  return 2.*(RM*3e8**2/(x*1e6)**2+ph0)

def multiple_plots():
  ar_list = glob('/data/FRB121102/analyses/C-band/*_puppi_*.calibP')

  for i,ar in enumerate(ar_list):
    if not ar_pars.has_key(path.splitext(path.basename(ar))[0]): continue
    fig, (ax1, ax2) = plt.subplots(2, sharey=True, sharex=True, figsize=(20,4))
 
    _, Q, U, _, err_Q, err_U = load_UQ(ar)

    L = np.sqrt(U**2+Q**2)
    err_L = np.sqrt((err_Q*Q.T)**2 + (err_U*U.T)**2).T / L
    QL = Q / L
    UL = U / L
    err_QL = np.sqrt(err_Q**2+(err_L*Q/L).T**2).T/L
    err_UL = np.sqrt(err_U**2+(err_L*U/L).T**2).T/L

    xUQ = np.linspace(4500.78125-400, 4500.78125+400, L.shape[0])
    ax1.errorbar(xUQ, QL, yerr=err_QL, fmt='r.')
    ax2.errorbar(xUQ, UL, yerr=err_UL, fmt='r.')

    QU = Q + 1j*U
    err_QU = err_Q + 1j*err_U
    idx = np.where(~np.isnan(QU))[0]
    xUQ = xUQ[idx]
    QU = QU[idx]

    y = np.angle(QU)

    p0 = [-1000.,100000.]
    par, pcov = curve_fit(RM_fit, xUQ, y, p0=p0)#, maxfev = 10000)#, sigma=err_QU)#, bounds=[[-2*np.pi,-np.inf],[2*np.pi,np.inf]])
    err_par = np.sqrt(np.diag(pcov))

    x = np.linspace(4500.78125-400, 4500.78125+400, 1e4)
    ax1.plot(x, cos(x, *par), 'k-')
    ax2.plot(x, sin(x, *par), 'k-')

    print "\nBurst {} - {}".format(i, ar)
    print "RM: {:.0f} +- {:.0f}, PA: {:.3f} +- {:.3f}".format(par[1], err_par[1], par[0], err_par[0])

    ax1.set_xlim([4140, 4860])
    ax1.set_ylim([-2,2])
    ax1.set_ylabel("Q/L")
    ax2.set_ylabel("U/L")
    ax2.set_xlabel("Frequency (MHz)")
    ax1.yaxis.set_major_locator(MultipleLocator(1))
    ax2.yaxis.set_major_locator(MultipleLocator(1))

    fig.tight_layout()
    plt.show()

  return


def single_plot():
  ar_list = glob('/data/FRB121102/analyses/C-band/*_puppi_57748*.calibP')
  fig, (ax1, ax2) = plt.subplots(2, sharey=True, sharex=True, figsize=(20,4))

  for i,ar in enumerate(ar_list):
    if not ar_pars.has_key(path.basename(ar)): continue

    I, Q, U, err_I, err_Q, err_U = load_UQ(ar)

    xUQ = np.linspace(4500.78125-400, 4500.78125+400, U.size)

    L = np.sqrt(U**2+Q**2)
    err_L = np.sqrt((err_Q*Q.T)**2 + (err_U*U.T)**2).T / L
    QL = Q / L
    UL = U / L
    err_QL = np.sqrt(err_Q**2+(err_L*Q/L).T**2).T/L
    err_UL = np.sqrt(err_U**2+(err_L*U/L).T**2).T/L

    ang = np.arccos(QL)
    QL = QL*np.cos(ar_prop[path.basename(ar)][3]) + (((np.mod(ang, 2*np.pi) < np.pi/2) | (np.mod(ang, 2*np.pi) > 3*np.pi/2)) * 2 - 1) * np.sqrt(1-QL**2)*np.sin(ar_prop[path.basename(ar)][3])


    idx = np.where(~np.isnan(L))
    #idx = np.where(L>3*err_L)[0]
    xUQ = xUQ[idx]
    QL = QL[idx]
    UL = UL[idx]
    err_QL = err_QL[idx]
    err_UL = err_UL[idx]


    plot(ax1, ax2, ar, xUQ, QL, UL, err_QL, err_UL)

  fig.tight_layout()
  plt.savefig('UQ.png', format='png')
  plt.show()

  return



def plot(ax1, ax2, ar_name, xUQ, yQ, yU, err_yQ, err_yU, fit=False):
  x = np.linspace(4500.78125-400, 4500.78125+400, xUQ.size)

  ax1.errorbar(xUQ, yQ, yerr=err_yQ, fmt='.', label=ar_name)
  ax2.errorbar(xUQ, yU, yerr=err_yU, fmt='.', label=ar_name)

  if fit:
    p0 = [0.,102500]
    Q_par, _ = curve_fit(cos, xUQ, yQ, p0=p0, sigma=err_yQ)
    U_par, _ = curve_fit(sin, xUQ, yU, p0=p0, sigma=err_yU)
    ax1.plot(x, cos(x, *Q_par), 'k-')
    ax2.plot(x, sin(x, *U_par), 'k-')

    print Q_par
    print U_par

    #ax1.plot(x, cos(x, 0), 'k-')
    #ax2.plot(x, sin(x, 0), 'k-')

  #ax1.set_xlim([4140, 4860])
  ax1.set_ylim([-2,2])
  ax1.set_ylabel("Q/L")
  ax2.set_ylabel("U/L")
  ax2.set_xlabel("Frequency (MHz)")
  ax1.yaxis.set_major_locator(MultipleLocator(1))
  ax2.yaxis.set_major_locator(MultipleLocator(1))

  return


def global_fit():
  ar_list = glob('/data/FRB121102/analyses/C-band/[0-9][0-9]_puppi_57748*.calibP')
  fig, (ax1, ax2, ax3) = plt.subplots(3, sharey=True, sharex=True, figsize=[183*mm_to_in,183/3.*mm_to_in])

  n_ar = len(ar_pars.keys())
  U = np.zeros([n_ar,512])
  Q = np.zeros([n_ar,512])
  err_U = np.zeros(n_ar)
  err_Q = np.zeros(n_ar)
  i = 0
  for ar in ar_list:
    if not ar_pars.has_key(path.basename(ar)): continue
    _, Q[i], U[i], _, err_Q[i], err_U[i] = load_UQ(ar)
    i += 1
  
  L = np.sqrt(U**2+Q**2)
  err_L = np.sqrt((err_Q*Q.T)**2 + (err_U*U.T)**2).T / L
  QL = Q / L
  UL = U / L
  err_QL = np.sqrt(err_Q**2+(err_L*Q/L).T**2).T/L
  err_UL = np.sqrt(err_U**2+(err_L*U/L).T**2).T/L

  #QL = np.nansum(QL*err_QL**-2, axis=0)/np.nansum(err_QL**-2, axis=0)
  #UL = np.nansum(UL*err_UL**-2, axis=0)/np.nansum(err_UL**-2, axis=0)
  #err_QL = np.nansum(err_QL**-2, axis=0)
  #err_UL = np.nansum(err_UL**-2, axis=0)
  #xUQ = np.linspace(4500.78125-400, 4500.78125+400, U.size)

  QL = QL.flatten()
  UL = UL.flatten()
  err_QL = err_QL.flatten()
  err_UL = err_UL.flatten()

  xUQ = np.tile(np.linspace(4500.78125-400, 4500.78125+400, L.shape[1]), L.shape[1])
  idx = np.where(~np.isnan(QL))[0]
  xQ = xUQ[idx]
  QL = QL[idx]
  err_QL = err_QL[idx]
  idx = np.where(~np.isnan(UL))[0]
  xU = xUQ[idx]
  UL = UL[idx]
  err_UL = err_UL[idx]


  #idx = np.where(QL>3*err_QL)[0]
  #xQ = xUQ[idx]
  #QL = QL[idx]
  #idx = np.where(UL>3*err_UL)[0]
  #xU = xUQ[idx]
  #UL = UL[idx]

  for i in range(QL.shape[0]):
    #ax1.errorbar(xQ[i], QL[i], yerr=err_QL[i], fmt='.', capsize=0, zorder=2)
    #ax2.errorbar(xU[i], UL[i], yerr=err_UL[i], fmt='.', capsize=0, zorder=2)
    ax1.plot(xQ[i], QL[i], '.', zorder=2, ms=5)
    ax2.plot(xU[i], UL[i], '.', zorder=2, ms=5)

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

  plt.savefig('UQ.pdf', format='pdf')
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
    y[i], y_err[i] = average([n[1],n[5]], [n[2],n[6]])

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
  y_mjd, y_err_mjd = average(y[idx], y_err[idx])
  ax.errorbar(0, y_mjd, yerr=y_err_mjd, fmt='None', zorder=10, lw=5, capsize=0, ecolor='r')

  idx = np.where(x.astype(int) == 1)[0]
  y_mjd, y_err_mjd = average(y[idx], y_err[idx])
  ax.errorbar(1, y_mjd, yerr=y_err_mjd, fmt='None', zorder=10, lw=5, capsize=0, ecolor='r')

  
  #Inset
  axins = inset_axes(ax, width="70%", height="60%", loc=4)
  ax.add_patch(patches.Rectangle((-0.5, 102250), 1, 200, fill=False )) 
  idx = np.where(x.astype(int) == 0)[0]
  x = x[idx] * 24*60
  y = y[idx]
  y_err = y_err[idx]
  y_mjd, y_err_mjd = average(y, y_err)
  axins.errorbar(x, y, yerr=y_err, fmt='None', ecolor='k', capsize=0, lw=3)
  axins.errorbar(np.mean(x), y_mjd, yerr=y_err_mjd, fmt='None', zorder=10, lw=5, capsize=0, ecolor='r')
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
  plt.savefig('RMevolution.pdf', format='pdf')
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
  plt.savefig('RMevolution_zoom.png', format='png')
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
  plt.savefig('LI_f.pdf', format='pdf')
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
      archive = psrchive.Archive_load(ar_name)
      archive.dedisperse()
      archive.set_rotation_measure(rm)
      archive.defaraday()
      archive.fscrunch()
      archive.remove_baseline()
      [I,Q,U,V] = archive.get_data().squeeze()
      lim = ar_pars[path.basename(ar_name)]
      if lim[0] > I.size/2: err_lim = [0, lim[0]-I.size/20]
      else: err_lim = [lim[1]+I.size/20, I.size-1]
      err_I = np.std(I[err_lim])
      err_Q = np.std(Q[err_lim])
      err_U = np.std(U[err_lim])
      L = np.sqrt(U**2 + Q**2)
      #err_L = np.sqrt((err_Q*Q)**2 + (err_U*U)**2) / L
      err_L = np.std(L[err_lim])
      L = L[lim[0]:lim[1]].sum()
      I = I[lim[0]:lim[1]].sum()
      #err_L = np.sqrt(np.sum(err_L[lim[0]:lim[1]]**2))
      err_L *= np.sqrt(lim[1]-lim[0])
      err_I *= np.sqrt(lim[1]-lim[0])
      LI = L / I
      err_LI = np.sqrt((err_L/I)**2+(err_I*L/I**2)**2)
      return LI, err_LI
    for i,rm in enumerate(RM_range):
      if (i*len(RM_range))%len(RM_range) == 0: print "Finished: {:.2f}%".format(float(i)/len(RM_range)*100)
      LI[j,i], err_LI[j,i] = load_LI(ar_name, rm)

    ax.fill_between(RM_range, 0, LI[j]-np.median(LI[j])+d, facecolor='w', edgecolor='k')
    d -= 0.5

  np.save('LI_RM', LI)

  #ax.plot(RM_range, LI, 'k.')
  #ax.errorbar(RM_range, LI, yerr=err_LI, fmt='k.', ecolor='k', capsize=0)
  ax.set_xlabel('RM (rad$\,$m$^{-2}$)')
  ax.set_ylabel('L/I')
  ax.set_ylim([d, len(ar_list)+1.1])
  ax.tick_params(axis='y', labelleft='off', left='off', right='off')
  ax.set_xlim([RM_range[0], RM_range[-1]])
  ax.ticklabel_format(useOffset=False)

  fig.subplots_adjust(left=0.1,right=.98,bottom=.1,top=.98)
  plt.savefig('LI_RM.pdf', format='pdf')
  #plt.show()
  
  return



def DM_RM_all():
  fig, ax = plt.subplots(1, figsize=[89*mm_to_in,89*mm_to_in])
  DM_RM = np.loadtxt("ATNF.dat", usecols=[2,3], delimiter=';', skiprows=1).T
  PSR_list = np.loadtxt("ATNF.dat", usecols=[1], delimiter=';', skiprows=1, dtype='str')
  ax.plot(DM_RM[0], np.abs(DM_RM[1]), 'ko')

  magnetars = ['J1550-5418', 'J1809-1943', 'J1745-2900', '1622-4950']
  sorter = np.argsort(PSR_list)
  idx = sorter[np.searchsorted(PSR_list, magnetars, sorter=sorter)]
  for i in idx:
    ax.plot(DM_RM[0,i], np.abs(DM_RM[1,i]), 'ro', zorder=10)

  
  ax.plot([55,225], 1.1927**2*np.array([102300,102300]), 'g-', lw=3)


  ax.annotate('J1746-2849', xy=(1456,10104), ha='left', va='baseline')
  ax.annotate('J1746-2856', xy=(1168,13253), ha='left', va='baseline')
  ax.annotate('J1745-2900', xy=(1778,66080), ha='left', va='baseline')
  ax.annotate('Sgr A*', xy=(0,66080), ha='left', va='baseline')
  ax.annotate('J1745-2900', xy=(1778,66080), ha='left', va='baseline')


  ax.set_xscale("log", nonposx='clip')
  ax.set_yscale("log", nonposx='clip')
  ax.set_xlim([2,2300])
  ax.set_ylim([-.5,2e5])
  ax.set_xlabel('DM (pc$\,$cm$^{-3}$')

  ax.set_ylabel('RM (rad$\,$m$^{-2}$')
  fig.subplots_adjust(left=0.15,right=.98,bottom=.15,top=.98)
  plt.savefig('DM_RM_all.pdf', format='pdf')
  plt.show()
  return


def PA_f():
  ar_list = glob('/data/FRB121102/analyses/C-band/*_puppi_*.RM')
  ar_list = [val for val in ar_list if path.splitext(path.basename(val))[0] in ar_pars.keys()]

  fig, (ax1, ax2, ax3) = plt.subplots(3, figsize=[89*mm_to_in,89*mm_to_in])
  x = np.linspace(4500.78125-400, 4500.78125+400, 512)
  xt_flag = 0

  c = ['blue', 'orange', 'green', 'red', 'purple', 'brown', 'pink', 'gray', 'olive', 'cyan', 'k']
  for i,ar_name in enumerate(ar_list):
    I, Q, U, err_I, err_Q, err_U = load_UQ(ar_name)

    L = np.sqrt(U**2+Q**2)
    err_L = np.sqrt((err_Q*Q.T)**2 + (err_U*U.T)**2).T / L
    LI = L/I
    err_LI = np.sqrt((err_L/I)**2+(err_I*L.T/I.T**2).T**2)
    ax1.plot(x, LI, c[i])

    ang = np.rad2deg(np.arctan2(Q, U) / 2.)
    err_ang = np.rad2deg(np.sqrt((err_Q*U)**2+(err_U*Q)**2)/2./(Q**2+U**2))
    ax2.errorbar(x, ang, yerr=err_ang, ecolor=c[i], fmt='None')

    I, Q, U, err_I, err_Q, err_U = load_UQ_time(ar_name)
    ang = np.rad2deg(np.arctan2(Q, U) / 2.)
    err_ang = np.rad2deg(np.sqrt((err_Q*U)**2+(err_U*Q)**2)/2./(Q**2+U**2))
    xt = np.arange(ang.size) + xt_flag
    ax3.errorbar(xt, ang, yerr=err_ang, ecolor=c[i], fmt='None', capsize=0)
    xt_flag += ang.size

  ax1.set_xlim([4500.78125-400, 4500.78125+400])
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
  plt.savefig('PA.pdf', format='pdf')
  plt.show()


if __name__ == '__main__':
  #multiple_plots() 
  #single_plot()
  #RMevolution()
  LI_plot()
  #Faraday_spectrum()
  #LI_RM()
  #DM_RM_all()
  #PA_f()



