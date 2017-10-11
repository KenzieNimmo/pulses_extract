import numpy as np
import matplotlib.pyplot as plt
from glob import glob
#import psrchive
import os
import argparse


def parser():
  # Command-line options
  parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                  description="Brute-force code to find periodicity.")
  parser.add_argument('-Pmin', help="Minimum period to search (ms).", default=0., type=float)
  parser.add_argument('-Pmax', help="Maximum period to search (ms).", default=100., type=float)
  parser.add_argument('-bin_num', help="Number of bins in the profile.", default=1., type=float)
  parser.add_argument('-step_res', help="Time resolution of the step (ms).", default=1e5, type=float)
  return parser.parse_args()



def main():
  data, puls_n = load_matrix()
  out = np.zeros([data.shape[1], data.shape[1]])
  for i in range(data.shape[1]):
    for j,n in enumerate(data):
      roll_n = -np.mod(np.int(i * puls_n[j]), data.shape[1])
      out[i] += np.roll(n, roll_n)
    #out[i,:-i] = 0
  
  return


def prepare_matrix():
  ar_list = glob('/path/to*.sp')  # Archives with 10.24us resolution
  data = []
  for ar in ar_list:
    load_archive = psrchive.Archive_load(ar)
    load_archive.tscrunch()
    load_archive.fscrunch()
    load_archive.pscrunch()
    load_archive.remove_baseline()
    w = load_archive.get_weights()
    archive = load_archive.get_data().squeeze() * w
    data.append(archive)
  data = np.array(data)
  #puls_n = #extract pulse number from filename
  puls_n -= puls_n.min()
  data = np.vstack([puls_n, data.T]).T
  np.save(data_file, data)
  return


def load_matrix():  
  data = np.load(data_file)  #Matrix with one pulse per row
  puls_n = data[:,0]
  #puls_n /= puls_n.max()
  data = data[:,1:]
  #data = np.reshape(data, [data.shape[0], data.shape[1]/10, 10]).sum(axis=2) #Reduce time resolution to 0.1024 ms
  #data = np.roll(data, -np.argmax(data[0])+30, axis=1)  #Bring the first pulse 0.3072 ms from the origin
  data = np.roll(data, -np.argmax(data[0])+data.shape[1]/2, axis=1)  #Bring the first pulse at the center
  return data, puls_n



def TOAs(max_period = 1e2):
  #puppi_57614_C0531+33_0803
  #Time = np.array([318.310236, 655.795077, 1970.861998, 2201.543721, 2367.42443, 2546.742641, 2617.79628, 2640.094577, 2710.297313, 2782.591631, 3036.52651, 3046.904955, 3095.933092, 3779.83574, 4027.2234909999997, 4409.396593, 4687.717007, 4730.516439, 4921.753436, 5026.912665999999, 6033.050829, 6124.302909, 6780.836004000001])

  #puppi_57747_C0531+33_0558
  #Time = np.array([515.1714509999999, 1173.711585, 1958.604964, 2416.7920440000003, 2666.797711, 3169.8055170000002, 3174.507397, 3648.334561, 3695.617556, 4524.416614])
  
  #Breakthrough
  #Time = np.array([16.2505, 285.436, 323.352, 344.768, 356.033, 597.612, 691.83, 769.864, 1036.42, 1142.41, 1454.55, 1904.922, 2453.071, 3326.39])
  
  #PSR J0139+33
  Time = np.array([367.72409, 918.07233, 1005.4276, 1064.0929, 1065.3397, 1469.8207, 2161.0476, 2897.3311, 3059.5791])

  
  
  Time *= 1000
  Time -= Time.min()

  t_res = max_period / 1000.
  
  Time += max_period / 2.
  puls_n = np.floor(Time / max_period)
  puls_n /= puls_n.max()
  
  size_x = np.floor(max_period / t_res / puls_n[1]).astype(int)
  size_y = np.ceil(max_period / t_res).astype(int)
  data = np.zeros((Time.size, size_y), dtype=np.int8)
  t_bin = ((Time % max_period) / t_res).astype(int)
  interval = np.ceil(2. / t_res).astype(int)
  for i,n in enumerate(t_bin):
    data[i, n-interval:n+interval] = 1
  out = np.zeros((size_x, size_y), dtype=np.int8)

  for i in range(size_x):
    if i%100==0: print "%.2f" % (float(i)/size_x*100)
    for j,n in enumerate(data):
      roll_n = int(-i * puls_n[j])
      out[i] += np.roll(n, roll_n)
  
  f, (ax1, ax2) = plt.subplots(1, 2, sharey=True, figsize=(10,10))
  overlap = np.clip(out-1, 0, np.inf).sum(axis=1)
  ax2.plot(overlap, np.linspace(max_period,0,out.shape[0]), 'k')
  ax1.imshow(out,origin="lower",aspect='auto',interpolation='nearest',cmap='Greys', extent=[0, max_period, max_period, 0])
  ax1.set_xlabel('Phase (ms)')
  ax2.set_xlabel('Power (arb.)')
  ax1.set_ylabel('Period (ms)')
  #plt.savefig('periodMax_{:.0f}ms.png'.format(max_period), papertype = 'a4', orientation = 'portrait', format = 'png')
  plt.show()
  return



def TOAs2(min_period=0., max_period=1e2, puls_width=1.):
  #PSR J0139+33
  Time = np.array([367.72409, 918.07233, 1005.4276, 1064.0929, 1065.3397, 2161.0476, 2897.3311, 3059.5791])  #p=1247.95ms
  
  #Convert times in ms and place the first in the middle of phace
  Time *= 1000
  Time -= Time.min()
  Time += max_period / 2.
  
  #Create the matrix with TOAs
  length = max_period / puls_width
  data = np.zeros([Time.size, int(np.ceil(length))], dtype=np.int8)  #ceil?
  
  t_bin = (Time % max_period) / puls_width
  for i,n in enumerate(t_bin.astype(int)):
    data[i, n:n+1] = 1

  puls_n = np.floor(Time / max_period)
  puls_rel = puls_n/puls_n.max()

  #Identify peaks
  out = []
  period = []
  n_iter = int((np.ceil((max_period - min_period) / puls_width) * puls_n.max())) #int((np.ceil((max_period - min_period) / puls_width)) * (Time.max() / max_period))
  print "  0 % completed",
  for i in range(n_iter):
    progress = int((float(i) / n_iter * 100))
    if progress % 1 == 0: print " \r{:2}".format(progress),
    p = max_period - i / puls_n.max() * puls_width  #p = max_period - i * puls_width / Time.max() * max_period
    fold_n = data.shape[1] - int(i/puls_n.max())   
    profile = np.zeros(data.shape[1])  
    for j,n in enumerate(data):
      roll_n = int(np.round(-i * puls_rel[j]))
      profile += np.roll(n, roll_n)
    profile.resize([int(np.ceil(float(profile.size) / fold_n)), fold_n])
    profile = profile.sum(axis=0)
    #profile.resize(bins, profile.size/bins+1)
    #profile = profile.sum(axis=1)
    profile.resize(profile.size/4+1, 4)
    profile = profile.sum(axis=1)
    if profile[profile>1].size > 0:
      if profile.max() > Time.size / 2:
        period.append(p)
        out.append(profile)
  
  print "\n"

  if len(period) == 0: 
    print "No periodicity detected."
    exit()

  max_list = np.array([max(n) for n in out])
  idx = np.argwhere(max_list == np.amax(max_list))[0]
  print "Best periods: ", [period[i] for i in idx]

  if len(period) > 30: 
    print "Too many periods!", len(period)
    exit()

  elif len(period) == 1:
    fig, ax = plt.subplots(1, figsize=(10,10))
    ax.bar(np.arange(out[i].size), out[i], width=1., color='k', linewidth=0)
    ax.annotate("p = {:.2f} ms".format(period[i]), (.8,.9), xycoords='axes fraction')
    ax.tick_params(axis='both', left='off', right='off', labelleft='on', top='off', bottom='off', labelbottom='off')
    fig.subplots_adjust(wspace=0, hspace=0)
    plt.show()

  #Plot the result
  else:
    fig, ax_arr = plt.subplots(len(period), 1, figsize=(10,10))
    for i,ax in enumerate(ax_arr):
      ax.bar(np.arange(out[i].size), out[i], width=1., color='k', linewidth=0)
      ax.annotate("p = {:.2f} ms".format(period[i]), (.8,.9), xycoords='axes fraction')
      ax.tick_params(axis='both', left='off', right='off', labelleft='on', top='off', bottom='off', labelbottom='off')
      ax.set_ylim([0, out[i].max()+.5])
    fig.subplots_adjust(wspace=0, hspace=0)
    plt.show()

  return



def TOAs3(min_period=0., max_period=1e2, bin_res=1., step_res=1e5):
  #PSR J0139+33
  #Time = np.array([367.72409, 918.07233, 1005.4276, 1064.0929, 1065.3397, 2161.0476, 2897.3311, 3059.5791])  #p=1247.95ms

  #puppi_57614_C0531+33_0803
  #Time = np.array([318.310236, 655.795077, 1970.861998, 2201.543721, 2367.42443, 2546.742641, 2617.79628, 2640.094577, 2710.297313, 2782.591631, 3036.52651, 3046.904955, 3095.933092, 3779.83574, 4027.2234909999997, 4409.396593, 4687.717007, 4730.516439, 4921.753436, 5026.912665999999, 6033.050829, 6124.302909, 6780.836004000001])

  #puppi_57747_C0531+33_0558
  Time = np.array([515.1714509999999, 1173.711585, 1958.604964, 2416.7920440000003, 2666.797711, 3169.8055170000002, 3174.507397, 3648.334561, 3695.617556, 4524.416614])
  
  #Breakthrough
  #Time = np.array([16.2505, 285.436, 323.352, 344.768, 356.033, 597.612, 691.83, 769.864, 1036.42, 1142.41, 1454.55, 1904.922, 2453.071, 3326.39])
  
  #Convert times in ms and place the first in the middle of phace
  Time *= 1000
  Time -= Time.min()
  #Time += max_period / 2.

  def calculate_profile(p, bin_res, Time):
    profile = np.zeros(int(np.ceil(bin_res)), dtype=np.int8)  #ceil?
    t_bin = (Time % p) / p * bin_res
    for n in t_bin.astype(int):
      profile[n:n+1] += 1
    return profile


  p = min_period
    
  period = []
  out = []
  while p < max_period:
    profile = calculate_profile(p, bin_res, Time)
    if profile[profile>1].size > 0:
      if profile.max() > Time.size / 2:
        period.append(p)
        out.append(profile)
    p += (p / step_res)

  if len(period) == 0: 
    print "No periodicity detected."
    exit()

  max_list = np.array([max(n) for n in out])
  if len(period) > 10: 
    print "Detected {} periods, only the most significat 10 are plotted".format(len(period))
    idx = np.argsort(-max_list)
    period = np.array(period)[idx]
    out = np.array(out)[idx]
    period = period[:10]
    out = out[:10]
    idx = np.argsort(period)
    period = period[idx]
    out = out[idx]

  max_list = np.array([max(n) for n in out])
  idx = np.argwhere(max_list == np.amax(max_list))
  print "Best periods: ", period

  if len(period) == 1:
    fig, ax = plt.subplots(1, figsize=(10,10))
    ax.bar(np.arange(out[i].size), out[i], width=1., color='k', linewidth=0)
    ax.annotate("p = {:.2f} ms".format(period[i]), (.8,.9), xycoords='axes fraction')
    ax.tick_params(axis='both', left='off', right='off', labelleft='on', top='off', bottom='off', labelbottom='off')
    fig.subplots_adjust(wspace=0, hspace=0)
    plt.show()

  #Plot the result
  else:
    fig, ax_arr = plt.subplots(len(period), 1, figsize=(10,20))
    for i,ax in enumerate(ax_arr):
      ax.bar(np.arange(out[i].size), out[i], width=1., color='k', linewidth=0)
      ax.annotate("p = {:.2f} ms".format(period[i]), (.8,.9), xycoords='axes fraction', va='top')
      ax.tick_params(axis='both', left='off', right='off', labelleft='on', labelright='on', top='off', bottom='off', labelbottom='off')
      ax.set_ylim([0, out[i].max()+.9])
    fig.subplots_adjust(wspace=0, hspace=0)
    plt.show()



  return



def plot(arr, puls_width=10., f=10, max_period=1.8967887942982455e3):
          Time = np.array(arr)
          #Convert times in ms and place the first in the middle of phace
          Time *= 1000
          Time -= Time.min()
          Time += max_period / 2.
          length = max_period / puls_width
          data = np.zeros([Time.size, int(np.ceil(length))], dtype=np.int8)  #ceil?
          t_bin = (Time % max_period) / puls_width
          for i,n in enumerate(t_bin.astype(int)):
              data[i, n:n+1] = 1
          puls_n = (np.floor(Time / max_period) / f).astype(int)
          t = np.zeros((puls_n.max()+5, data.shape[1]))
          for i,j in enumerate(puls_n):
            t[j] = data[i]
          plt.imshow(t, aspect='auto', interpolation='nearest', cmap='Greys', origin='lower', extent=[0, max_period, 0, max_period*puls_n.max()*f/1000.])
          plt.xlabel('Phase (ms)')
          plt.ylabel('Time (s)')
          plt.title('period {:.5} ms'.format(max_period))
          plt.show()


def plot(arr, length=20., f=20., max_period=1.8967887942982455e3):
          Time = np.array(arr)
          #Convert times in ms and place the first in the middle of phace
          Time *= 1000
          Time -= Time.min()
          Time += max_period / 2.
          puls_width = max_period / length
          data = np.zeros([Time.size, int(np.ceil(length))], dtype=np.int8)  #ceil?
          t_bin = (Time % max_period) / puls_width
          for i,n in enumerate(t_bin.astype(int)):
              data[i, n:n+1] = 1
          puls_n = np.floor(Time / max_period)
          puls_n = (puls_n / puls_n.max() * f).astype(int)
          t = np.zeros((puls_n.max()+1, data.shape[1]))
          for i,j in enumerate(puls_n):
            t[j] += data[i]
          t = np.clip(t, 0, 1)
          plt.imshow(t, aspect='auto', interpolation='nearest', cmap='Greys', origin='lower', extent=[0, max_period, 0, max_period*puls_n.max()*f/1000.])
          plt.xlabel('Phase (ms)')
          plt.ylabel('Time (s)')
          plt.title('period {:.5} ms'.format(max_period))
          plt.show()


def plot_hist(arr, length=20., max_period=1.8967887942982455e3):
          Time = np.array(arr)
          #Convert times in ms and place the first in the middle of phace
          Time *= 1000
          Time -= Time.min()
          Time += max_period / 2.
          puls_width = max_period / length
          data = np.zeros([Time.size, int(np.ceil(length))], dtype=np.int8)  #ceil?
          t_bin = (Time % max_period) / puls_width
          for i,n in enumerate(t_bin.astype(int)):
              data[i, n:n+1] = 1
          plt.bar(np.linspace(0., max_period, length), data.sum(axis=0), width=puls_width, color='k', linewidth=0)
          plt.xlabel('Phase (ms)')
          plt.ylabel('N')
          plt.title('period {:.5} ms'.format(max_period))
          #plt.xlim([0, max_period])
          plt.show()



if __name__ == '__main__':
  args = parser()
  try:
    os.system('setterm -cursor off')
    TOAs3(min_period=args.Pmin, max_period=args.Pmax, bin_res=args.bin_num, step_res=args.step_res)
  finally:
    os.system('setterm -cursor on')

