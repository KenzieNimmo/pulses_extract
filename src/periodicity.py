import numpy as np
import matplotlib.pyplot as plt
import os
import argparse
from matplotlib.ticker import MaxNLocator


def parser():
  # Command-line options
  parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                  description="Brute-force code to find periodicity.")
  parser.add_argument('-Pmin', help="Minimum period to search (ms).", default=0.1, type=float)
  parser.add_argument('-Pmax', help="Maximum period to search (ms).", default=100., type=float)
  parser.add_argument('-bin_res', help="Number of bins in the profile.", default=10., type=float)
  parser.add_argument('-step_res', help="Time resolution of the steps.", default=1e5, type=float)
  return parser.parse_args()


#PSR J0139+33, p=1247.95ms
Time_RRAT_B = np.array([367.72415, 918.07129, 1005.4286, 1064.0923, 1065.3408, 2629.0305, 2897.3308, 3059.5774, 3126.9622])
Time_RRAT_A = np.array([3.3460729, 215.49582, 324.06693, 433.88669, 868.18079, 1091.5649, 1138.9696, 1433.505, 1535.8436, 1918.9585, 2077.4385, 2078.6909, 2260.8911])

#puppi_57614_C0531+33_0803
Time_Miraculus = np.array([317.576187644, 655.057103252, 1970.12925811, 2200.80836418, 2366.69169011, 2546.00205032, 2617.06354011, 2639.36052864, 2709.56064772, 2781.85758264, 3035.79115318, 3046.17221511, 3095.19773518, 3779.10300011, 4026.48944264, 4408.66385311, 4686.96856554, 4729.77584832, 4921.753436, 5026.17730918, 6032.31547218, 6123.56886064, 6780.10195564])

#puppi_57747_C0531+33_0558
Time_Cband = np.array([515.07245744, 1173.61328249, 1958.50752531, 2416.69374149, 2666.70009955, 3169.70825108, 3174.41030384, 3648.23245769, 3695.51994455, 4524.31952084])

#Breakthrough
Time_BTL = np.array([16.252, 263.415, 285.463, 323.381, 344.799, 356.065, 580.683, 597.667, 691.892, 704.128, 769.934, 841.012, 993.339, 1036.519, 1142.513, 1257.52, 1454.681, 1789.513])


def main(args):
  # Defines the two observations to use
  obs_A = Time_Miraculus.copy()
  obs_B = Time_BTL.copy()
  obs_C = Time_Cband.copy()

  #Convert times in ms and place the first in the middle of phace
  obs_A *= 1000
  obs_A -= obs_A.min()
  #obs_A += args.Pmin / 2.
  obs_B *= 1000
  obs_B -= obs_B.min()
  #obs_B += args.Pmin / 2.
  obs_C *= 1000
  obs_C -= obs_C.min()
  #obs_C += args.Pmin / 2.

  # Get periodicities from obs A
  period_A = TOAs(obs_A, min_period=args.Pmin, max_period=args.Pmax, bin_res=args.bin_res, step_res=args.step_res, lim=obs_A.size/3.)

  print "Periodicities detected in the first dataset:", len(period_A)

  period_B = []
  for i,pA in enumerate(period_A):
    p = TOAs(obs_B, min_period=pA-pA/args.step_res, max_period=pA+pA/args.step_res, bin_res=args.bin_res, step_res=args.step_res*10, lim=obs_B.size/3.)
    pB = get_maxTOA(p, obs_B, bin_res=10., n_max=1)
    period_B.extend(pB)
  print "Periodicities detected in the second dataset:", len(period_B)

  period_C = []
  for i,pB in enumerate(period_B):
    p = TOAs(obs_C, min_period=pB-pB/args.step_res, max_period=pB+pB/args.step_res, bin_res=args.bin_res, step_res=args.step_res*10, lim=obs_C.size/3.)
    pC = get_maxTOA(p, obs_C, bin_res=10., n_max=1)
    period_C.extend(pC)
  print "Periodicities detected in the third dataset:", len(period_C)

  period = period_C
  obs = obs_C

  if len(period) == 0: 
    print "No periodicity detected."
    return
  elif len(period) > 10:  
    print "Detected {} periods, only the most significat 10 are plotted".format(len(period_C))
    period = get_maxTOA(period, obs, bin_res=args.bin_res, n_max=10)
  plot_period(period, obs, bin_res=args.bin_res)
  
  return


def get_profile(p, bin_res, Time):
  profile = np.zeros(int(np.ceil(bin_res)), dtype=np.int8)  #ceil?
  t_bin = (Time % p) / p * bin_res
  for n in t_bin.astype(int):
    profile[n:n+1] += 1
  return profile


def get_maxTOA(period, Time, bin_res=10., n_max=10):
  max_list = []
  for p in period:
    max_list.append(max(get_profile(p, bin_res, Time)))
  max_list = np.array(max_list)

  idx = np.argsort(-max_list)
  period = np.array(period)[idx]
  period = period[:n_max]
  idx = np.argsort(period)
  period = period[idx]

  return period


def TOAs(Time, min_period=0.1, max_period=1e2, bin_res=1., step_res=1e5, lim=1):
  p = min_period
  period = []
  while p < max_period:
    profile = get_profile(p, bin_res, Time)
    if profile[profile>1].size > 0:
      if profile.max() >= lim:
        period.append(p)
    p += (p / step_res)
  return np.array(period)


def plot_period(period, Time, bin_res=10.):
  profile = []
  for p in period:
    profile.append(get_profile(p, bin_res, Time))

  #Plot the result
  if len(period) == 1:
    fig, ax = plt.subplots(1, figsize=(10,10))
    ax.bar(np.arange(profile[0].size), profile[0], width=1., color='k', linewidth=0)
    ax.annotate("p = {:.2f} ms".format(period[0]), (.8,.9), xycoords='axes fraction')
    ax.tick_params(axis='both', left='off', right='off', labelleft='on', top='off', bottom='off', labelbottom='off')
    ax.yaxis.set_major_locator(MaxNLocator(nbins=5, integer=True))

  else:
    fig, ax_arr = plt.subplots(len(period), 1, figsize=(10,20))
    for i,ax in enumerate(ax_arr):
      ax.bar(np.arange(profile[i].size), profile[i], width=1., color='k', linewidth=0)
      ax.annotate("p = {:.2f} ms".format(period[i]), (.8,.9), xycoords='axes fraction', va='top')
      ax.tick_params(axis='both', left='off', right='off', labelleft='on', labelright='on', top='off', bottom='off', labelbottom='off')
      ax.set_ylim([0, profile[i].max()+.9])
      ax.yaxis.set_major_locator(MaxNLocator(nbins=5, integer=True))

  ax.set_xlabel('Pulse phase')
  fig.tight_layout()
  fig.subplots_adjust(wspace=0, hspace=0)
  plt.show()

  return


def plot_profile(period, Time, bin_res=10.):
  profile = get_profile(period, bin_res, Time)
  fig, ax = plt.subplots(1, figsize=(10,10))
  ax.bar(np.arange(profile[0].size), profile[0], width=1., color='k', linewidth=0)
  ax.annotate("p = {:.2f} ms".format(period[0]), (.8,.9), xycoords='axes fraction')
  ax.tick_params(axis='both', left='off', right='off', labelleft='on', top='off', bottom='off', labelbottom='off')
  ax.yaxis.set_major_locator(MaxNLocator(nbins=5, integer=True))
  ax.set_xlabel('Pulse phase')
  fig.tight_layout()
  fig.subplots_adjust(wspace=0, hspace=0)
  plt.show()
  return



if __name__ == '__main__':
  args = parser()
  main(args)

