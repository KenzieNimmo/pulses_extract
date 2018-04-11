import psrchive
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from matplotlib.ticker import MultipleLocator
from glob import glob
from os import path

mpl.rcParams['font.size'] = 7
mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['axes.linewidth'] = 1
mpl.rcParams['ps.useafm'] = True
mpl.rcParams['pdf.use14corefonts'] = True
mpl.rcParams['text.usetex'] = True
mm_to_in = 0.0393701


ar_pars = {
  '02_puppi_57747_C0531+33_0558_1183.dm559.72.calibP': [2030,2190],
  '03_puppi_57747_C0531+33_0558_1202.dm559.72.calibP': [2195,2258],
  '06_puppi_57747_C0531+33_0558_3687.dm559.72.calibP': [1383,1424],
  '07_puppi_57747_C0531+33_0558_3688.dm559.72.calibP': [1760,1886],
  '08_puppi_57747_C0531+33_0558_3689.dm559.72.calibP': [2112,2400],
  '01_puppi_57747_C0531+33_0558_5.dm559.72.calibP':    [1672,1866],
  '15_puppi_57748_C0531+33_0594_1269.dm559.72.calibP': [1750,1917],
  '12_puppi_57748_C0531+33_0594_48.dm559.72.calibP':   [1200,1292],
  '13_puppi_57748_C0531+33_0594_49.dm559.72.calibP':   [2060,2238],
  '14_puppi_57748_C0531+33_0594_50.dm559.72.calibP':   [1790,1840],
  '16_puppi_57772_C0531+33_0007_2695.dm559.72.calibP': [870,1050],
}


ar_prop = {
  '02_puppi_57747_C0531+33_0558_1183.dm559.72.calibP': [57747.1308531025, 102405, 40, 3.582, 0.374, 102383, 31, 3.764, 0.281],
  '03_puppi_57747_C0531+33_0558_1202.dm559.72.calibP': [57747.1399375004, 102373, 18, 3.747, 0.162, 102329, 15, 4.166, 0.137],
  '06_puppi_57747_C0531+33_0558_3687.dm559.72.calibP': [57747.1539561231, 102369, 32, 3.695, 0.293, 102356, 29, 3.852, 0.265],
  '07_puppi_57747_C0531+33_0558_3688.dm559.72.calibP': [57747.1540104922, 102278, 20, 4.636, 0.201, 102316, 33, 4.263, 0.328],
  '08_puppi_57747_C0531+33_0558_3689.dm559.72.calibP': [57747.1594948018, 102358, 23, 3.796, 0.193, 102307, 20, 4.217, 0.171],
  '01_puppi_57747_C0531+33_0558_5.dm559.72.calibP':    [57747.1232313932, 102368, 8, 3.745, 0.069, 102376, 8, 3.680, 0.076],
  '15_puppi_57748_C0531+33_0594_1269.dm559.72.calibP': [57748.1693814308, 102295, 24, 2.559, 0.240, 102225, 38, 3.206, 0.377],
  '12_puppi_57748_C0531+33_0594_48.dm559.72.calibP':   [57748.1472089520, 102093, 36, 4.661, 0.355, 102168, 38, 3.947, 0.372],
  '13_puppi_57748_C0531+33_0594_49.dm559.72.calibP':   [57748.1488993263, 102187, 15, 3.720, 0.122, 102205, 14, 3.579, 0.115],
  '14_puppi_57748_C0531+33_0594_50.dm559.72.calibP':   [57748.1512921576, 102134, 17, 4.197, 0.159, 102176, 19, 3.796, 0.169],
  '16_puppi_57772_C0531+33_0007_2695.dm559.72.calibP': [57772.1236149924, 102669, 10, -2.311, 0.094, 102633, 10, -1.979, 0.093],
}


def load_UQ(ar_name):
    load_archive = psrchive.Archive_load(ar_name)
    load_archive.remove_baseline()
    load_archive.dedisperse()
    ds = load_archive.get_data().squeeze()
    ds = ds[:,::-1,:]
    lim = ar_pars[path.basename(ar_name)]
    w = load_archive.get_weights().flatten()[::-1]
    idx = np.where(w==0)[0]
    ds = np.multiply(w, np.rollaxis(ds,2))
    ds = np.rollaxis(np.rollaxis(ds, 1, 0), 2, 1)
    if lim[0] > ds.shape[-1]/2: err_lim = [0, lim[0]-ds.shape[-1]/20]
    else: err_lim = [lim[1]+ds.shape[-1]/20, ds.shape[-1]]
    err_I, err_Q, err_U = load_err(ds, err_lim, lim[1]-lim[0])
    ds = ds[:, :, lim[0]:lim[1]]
    ds = ds.sum(axis=2)
    ds[:,idx] = np.nan
    ds[:, ds[0] < 5 * err_I] = np.nan
    return ds[0], ds[1], ds[2], err_I, err_Q, err_U

def load_err(ds, lim, bins):
  bl = ds[:,:,lim[0]:lim[1]].copy()
  new_size = (lim[1] - lim[0]) / bins * bins
  bl = bl[:,:,:new_size]
  bl = bl.reshape([bl.shape[0], bl.shape[1], new_size/bins, bins]).sum(axis=-1)
  return np.nanstd(bl[0]), np.nanstd(bl[1]), np.nanstd(bl[2])


def cos(x, *p):
    ph0, RM = p                                      
    return np.cos(RM*3e8**2/(x*1e6)**2 * 2 + ph0)

def sin(x, *p):
    ph0, RM = p                                      
    return np.sin(RM*3e8**2/(x*1e6)**2 * 2 + ph0)


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



def multiple_plots():
  ar_list = glob('/data/FRB121102/analyses/C-band/*_puppi_*.calibP')

  for i,ar in enumerate(ar_list):
    if not ar_pars.has_key(path.basename(ar)): continue
    fig, (ax1, ax2) = plt.subplots(2, sharey=True, sharex=True, figsize=(20,4))
 
    _, Q, U, _, err_Q, err_U = load_UQ(ar)

    L = np.sqrt(U**2+Q**2)
    err_L = np.sqrt((err_Q*Q.T)**2 + (err_U*U.T)**2).T / L
    QL = Q / L
    UL = U / L
    err_QL = np.sqrt(err_Q**2+(err_L*Q/L).T**2).T/L
    err_UL = np.sqrt(err_U**2+(err_L*U/L).T**2).T/L

    xUQ = np.linspace(4500.78125-400, 4500.78125+400, L.shape[0])
    idx = np.where(~np.isnan(QL))[0]
    xQ = xUQ[idx]
    QL = QL[idx]
    err_QL = err_QL[idx]
    idx = np.where(~np.isnan(UL))[0]
    xU = xUQ[idx]
    UL = UL[idx]
    err_UL = err_UL[idx]

    ax1.errorbar(xQ, QL, yerr=err_QL, fmt='r.')
    ax2.errorbar(xU, UL, yerr=err_UL, fmt='r.')

    x = np.linspace(4500.78125-400, 4500.78125+400, 1e4)
    p0 = [0.,102500]
    Q_par, pcov = curve_fit(cos, xQ, QL, p0=p0, sigma=err_QL)#, bounds=[[-2*np.pi,-np.inf],[2*np.pi,np.inf]])
    err_Q_par = np.sqrt(np.diag(pcov))
    U_par, pcov = curve_fit(sin, xU, UL, p0=p0, sigma=err_UL)#, bounds=[[-2*np.pi,-np.inf],[2*np.pi,np.inf]])
    err_U_par = np.sqrt(np.diag(pcov))
    ax1.plot(x, cos(x, *Q_par), 'k-')
    ax2.plot(x, sin(x, *U_par), 'k-')

    print "\nBurst {} - {}".format(i, ar)
    print "Q - RM: {:.0f} +- {:.0f}, PA: {:.3f} +- {:.3f}".format(Q_par[1], err_Q_par[1], Q_par[0], err_Q_par[0])
    print "U - RM: {:.0f} +- {:.0f}, PA: {:.3f} +- {:.3f}".format(U_par[1], err_U_par[1], U_par[0], err_U_par[0])

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


def RMevolution():
  y = np.zeros(len(ar_prop.keys()))
  y_err = np.zeros_like(y)
  for i,n in enumerate(ar_prop.values()):
    y[i], y_err[i] = np.average([n[1],n[5]], weights=[n[2],n[6]], returned=True)

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
  ax.errorbar(x, y, yerr=y_err, fmt='k.')
  ax.set_xlim([-4.9,30])
  ax.set_xlabel('Days')
  ax.set_ylabel('RM (rad$\,$m$^{-2}$)')
  ax.ticklabel_format(useOffset=False)

  fig.subplots_adjust(left=0.15,right=.98,bottom=.1,top=.98)
  plt.savefig('RMevolution.png', format='png')
  plt.show()


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


  return


def LI_plot():
  ar_list = glob('/data/FRB121102/analyses/C-band/[0-9][0-9]_puppi_57748*.calibP')
  fig, ax = plt.subplots(1, figsize=[89*mm_to_in,89*mm_to_in])

  n_ar = len(ar_pars.keys())
  I = np.zeros([n_ar,512])
  U = np.zeros([n_ar,512])
  Q = np.zeros([n_ar,512])
  err_I = np.zeros(n_ar)
  err_U = np.zeros(n_ar)
  err_Q = np.zeros(n_ar)
  i = 0
  for ar in ar_list:
    if not ar_pars.has_key(path.basename(ar)): continue
    I[i], Q[i], U[i], err_I[i], err_Q[i], err_U[i] = load_UQ(ar)
    i += 1

  L = np.sqrt(U**2+Q**2)
  err_L = np.sqrt((err_Q*Q.T)**2 + (err_U*U.T)**2).T / L

  LI = L/I
  err_LI = np.sqrt((err_L/I)**2+(err_I*L.T/I.T**2).T**2) 

  x = np.linspace(4500.78125-400, 4500.78125+400, L.shape[1])
  for i,n in enumerate(LI):
    ax.errorbar(x, n, yerr=err_LI[i], fmt='-', capsize=0, zorder=2)

  x = np.tile(np.linspace(4500.78125-400, 4500.78125+400, LI.shape[1]), LI.shape[1])
  LI = LI.flatten()
  err_LI = err_LI.flatten()
  idx = np.where(~np.isnan(LI))[0]
  x = x[idx]
  LI = LI[idx]
  err_LI = err_LI[idx]

  lin = np.poly1d(np.polyfit(x, LI, 1, w=1/err_LI))
  squ = np.poly1d(np.polyfit(x, LI, 2, w=1/err_LI))

  x_f = np.linspace(4500.78125-500, 4500.78125+500, 1e5)
  ax.plot(x_f, lin(x_f), color='k', zorder=3)
  #ax.plot(x_f, squ(x_f), color='k', zorder=3)

  ax.set_xlim([4140, 4860])
  #ax.set_ylim([-1.99,1.99])
  ax.set_ylabel("L/I")
  ax.set_xlabel("Frequency (MHz)")
  #ax.yaxis.set_major_locator(MultipleLocator(1))
  #ax.yaxis.set_major_locator(MultipleLocator(1))
  

  fig.subplots_adjust(left=.11,right=.98,bottom=.1,top=.98)
  plt.savefig('LI_f.pdf', format='pdf')
  plt.show()


  return


def Faraday_spectrum():
  ar_name = '/data/FRB121102/analyses/C-band/08_puppi_57747_C0531+33_0558_3689.dm559.72.calibP'

  RM_range = np.linspace(50000, 150000, 1e2)

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
    else: err_lim = [lim[1]+I.size/20, I.size]
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

  LI = np.zeros_like(RM_range)
  err_LI = np.zeros_like(RM_range)
  for i,rm in enumerate(RM_range):
    if (i*len(RM_range))%len(RM_range) == 0: print "Finished: {:.2f}%".format(float(i)/len(RM_range)*100)
    LI[i], err_LI[i] = load_LI(ar_name, rm)

  fig, ax = plt.subplots(1, figsize=[89*mm_to_in,89*mm_to_in])
  
  #ax.plot(RM_range, LI, 'k.')
  ax.errorbar(RM_range, LI, yerr=err_LI, fmt='k.', ecolor='k', capsize=0)
  ax.set_xlabel('RM (rad$\,$m$^{-2}$)')
  ax.set_ylabel('L/I')
  ax.set_xlim([RM_range[0], RM_range[-1]])
  ax.ticklabel_format(useOffset=False)

  fig.subplots_adjust(left=0.13,right=.95,bottom=.1,top=.98)
  plt.savefig('LI_RM.pdf', format='pdf')
  plt.show()
  
  return


if __name__ == '__main__':
  #multiple_plots() 
  Faraday_spectrum()
  #single_plot()
  #RMevolution()
  #LI_plot()


