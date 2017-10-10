import psrchive
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from matplotlib.ticker import MultipleLocator
from glob import glob
from os import path


ar_pars = {
  'puppi_57747_C0531+33_0558_1183.dm559.72.calibP': [2030,2190],
  'puppi_57747_C0531+33_0558_1202.dm559.72.calibP': [2195,2258],
  'puppi_57747_C0531+33_0558_3687.dm559.72.calibP': [1383,1424],
  'puppi_57747_C0531+33_0558_3688.dm559.72.calibP': [1760,1886],
  'puppi_57747_C0531+33_0558_3689.dm559.72.calibP': [2112,2400],
  'puppi_57747_C0531+33_0558_5.dm559.72.calibP':    [1672,1866],
  'puppi_57748_C0531+33_0594_1269.dm559.72.calibP': [1750,1917],
  'puppi_57748_C0531+33_0594_48.dm559.72.calibP':   [1200,1292],
  'puppi_57748_C0531+33_0594_49.dm559.72.calibP':   [2060,2238],
  'puppi_57748_C0531+33_0594_50.dm559.72.calibP':   [1790,1840],
  'puppi_57772_C0531+33_0007_2695.dm559.72.calibP': [870,1050],
}


ar_prop = {
  'puppi_57747_C0531+33_0558_1183.dm559.72.calibP': [57747.1308531025, 102405, 40, 3.582, 0.374, 102383, 31, 3.764, 0.281],
  'puppi_57747_C0531+33_0558_1202.dm559.72.calibP': [57747.1399375004, 102373, 18, 3.747, 0.162, 102329, 15, 4.166, 0.137],
  'puppi_57747_C0531+33_0558_3687.dm559.72.calibP': [57747.1539561231, 102369, 32, 3.695, 0.293, 102356, 29, 3.852, 0.265],
  'puppi_57747_C0531+33_0558_3688.dm559.72.calibP': [57747.1540104922, 102278, 20, 4.636, 0.201, 102316, 33, 4.263, 0.328],
  'puppi_57747_C0531+33_0558_3689.dm559.72.calibP': [57747.1594948018, 102358, 23, 3.796, 0.193, 102307, 20, 4.217, 0.171],
  'puppi_57747_C0531+33_0558_5.dm559.72.calibP':    [57747.1232313932, 102368, 8, 3.745, 0.069, 102376, 8, 3.680, 0.076],
  'puppi_57748_C0531+33_0594_1269.dm559.72.calibP': [57748.1693814308, 102295, 24, 2.559, 0.240, 102225, 38, 3.206, 0.377],
  'puppi_57748_C0531+33_0594_48.dm559.72.calibP':   [57748.1472089520, 102093, 36, 4.661, 0.355, 102168, 38, 3.947, 0.372],
  'puppi_57748_C0531+33_0594_49.dm559.72.calibP':   [57748.1488993263, 102187, 15, 3.720, 0.122, 102205, 14, 3.579, 0.115],
  'puppi_57748_C0531+33_0594_50.dm559.72.calibP':   [57748.1512921576, 102134, 17, 4.197, 0.159, 102176, 19, 3.796, 0.169],
  'puppi_57772_C0531+33_0007_2695.dm559.72.calibP': [57772.1236149924, 102669, 10, -2.311, 0.094, 102633, 10, -1.979, 0.093],
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
    ds[:, ds[0] < 3 * err_I] = np.nan
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
  ar_list = glob('/data/FRB121102/analyses/C-band/RM/UQ/puppi_5777*.calibP')
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
  ar_list = glob('/data/FRB121102/analyses/C-band/RM/UQ/puppi_*.calibP')

  for i,ar in enumerate(ar_list):
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

  ax1.set_xlim([4140, 4860])
  ax1.set_ylim([-2,2])
  ax1.set_ylabel("Q/L")
  ax2.set_ylabel("U/L")
  ax2.set_xlabel("Frequency (MHz)")
  ax1.yaxis.set_major_locator(MultipleLocator(1))
  ax2.yaxis.set_major_locator(MultipleLocator(1))

  return


def global_fit():
  ar_list = glob('/data/FRB121102/analyses/C-band/RM/UQ/puppi_57748*.calibP')
  fig, (ax1, ax2) = plt.subplots(2, sharey=True, sharex=True, figsize=(20,4))

  n_ar = len(ar_list)
  U = np.zeros([n_ar,512])
  Q = np.zeros([n_ar,512])
  err_U = np.zeros(n_ar)
  err_Q = np.zeros(n_ar)
  for i,ar in enumerate(ar_list):
    _, Q[i], U[i], _, err_Q[i], err_U[i] = load_UQ(ar)
  
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

  ax1.errorbar(xQ, QL, yerr=err_QL, fmt='r.')
  ax2.errorbar(xU, UL, yerr=err_UL, fmt='r.')


  x = np.linspace(4500.78125-400, 4500.78125+400, 1e4)
  p0 = [0.,102500]
  Q_par, pcov = curve_fit(cos, xQ, QL, p0=p0, sigma=err_QL)
  err_Q_par = np.sqrt(np.diag(pcov))
  U_par, pcov = curve_fit(sin, xU, UL, p0=p0, sigma=err_UL)
  err_U_par = np.sqrt(np.diag(pcov))
  ax1.plot(x, cos(x, *Q_par), 'k-')
  ax2.plot(x, sin(x, *U_par), 'k-')

  print "Q - RM: {:.0f} +- {:.0f}, PA: {:.3f} +- {:.3f}".format(Q_par[1], err_Q_par[1], Q_par[0], err_Q_par[0])
  print "U - RM: {:.0f} +- {:.0f}, PA: {:.3f} +- {:.3f}".format(U_par[1], err_U_par[1], U_par[0], err_U_par[0])

  #ax1.plot(x, cos(x, 0), 'k-')
  #ax2.plot(x, sin(x, 0), 'k-')

  ax1.set_xlim([4140, 4860])
  ax1.set_ylim([-2,2])
  ax1.set_ylabel("Q/L")
  ax2.set_ylabel("U/L")
  ax2.set_xlabel("Frequency (MHz)")
  ax1.yaxis.set_major_locator(MultipleLocator(1))
  ax2.yaxis.set_major_locator(MultipleLocator(1))


  fig.tight_layout()
  plt.savefig('UQ_fit.png', format='png')
  plt.show()
  return



if __name__ == '__main__':
  multiple_plots() 
  #global_fit()
  #single_plot()
