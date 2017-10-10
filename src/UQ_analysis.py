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
  'burst16.ar': [870,1050]
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


def main():
  ar_list = glob('/data/FRB121102/analyses/C-band/RM/UQ/puppi_57747*.calibP')
  fig, (ax1, ax2) = plt.subplots(2, sharey=True, sharex=True, figsize=(20,4))

  for i,ar in enumerate(ar_list):
    I, Q, U, err_I, err_Q, err_U = load_UQ(ar)

    xUQ = np.linspace(4500.78125-400, 4500.78125+400, U.size)

    L = np.sqrt(U**2+Q**2)
    err_L = np.sqrt(err_Q**2*Q**2 + err_U**2*U**2) / L

    idx = np.where(~np.isnan(L))
    #idx = np.where(L>3*err_L)[0]
    xUQ = xUQ[idx]
    Q = Q[idx]
    U = U[idx]
    L = L[idx]
    err_L = err_L[idx]

    err_QL = np.sqrt(err_Q**2+(err_L*Q/L)**2)/L
    err_UL = np.sqrt(err_U**2+(err_L*U/L)**2)/L


    plot(ax1, ax2, ar, Q, U, err_Q, err_U)

  fig.tight_layout()
  plt.savefig('UQ.png', format='png')
  plt.show()

  return


def plot(ax1, ax2, ar_name, xUQ, yQ, yU, err_yQ, err_yU, fit=False):
  x = np.linspace(4500.78125-400, 4500.78125+400, U.size)

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
  Q_par, _ = curve_fit(cos, xQ, QL, p0=p0, sigma=err_QL)
  U_par, _ = curve_fit(sin, xU, UL, p0=p0, sigma=err_UL)
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


  fig.tight_layout()
  plt.savefig('UQ_fit.png', format='png')
  plt.show()
  return



if __name__ == '__main__':
  #main() 
  global_fit()

