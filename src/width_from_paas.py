import numpy as np
from glob import glob

period = 41.94304  #ms

def von_mises(x, p):
  #To open .m files from paas
  #mu, k and A are the three parameters in the fil
  #x is defined in phase
  mu, k, A = p
  vm = np.exp(k*np.cos((x-mu)*2.*np.pi))/2./np.pi/np.i0(k)
  return vm / vm.max() * A

def von_mises_approx(x, p):
  #Approximate von Mises function with normal distribution
  #x is defined in phase
  mu, k, A = p
  s2 = 1./k
  vm = 1./(2*s2*np.pi)**0.5*np.exp(-((x-mu)*2*np.pi)**2/2/s2) 
  return vm / vm.max() * A

def curve(x, m_name):
  m = np.loadtxt(m_name)
  if len(m.shape) == 1:  y = von_mises_approx(x, m)
  else:
    y = np.zeros_like(x)
    for mi in m:
      y += von_mises_approx(x, mi)
  return y

x = np.linspace(0,1,1e6)

ar_list = glob('puppi_*.m')
for ar in ar_list:
  y = curve(x, ar)
  HM = x[y >= y.max()/2.]
  FWHM = (HM.max() - HM.min()) * period
  print "{:.30}: w = {:.3f}".format(ar, FWHM)


