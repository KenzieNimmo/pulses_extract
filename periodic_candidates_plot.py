import sys
import re
import glob
import subprocess
import argparse
import os

import numpy as np

import sifting

#From ACCEL_sift.py

# glob for ACCEL files
globaccel = "*ACCEL_*0"
# glob for .inf files
globinf = "*DM*.inf"
# In how many DMs must a candidate be detected to be considered "good"
min_num_DMs = 2
# Lowest DM to consider as a "real" pulsar
low_DM_cutoff = 2.0
# Ignore candidates with a sigma (from incoherent power summation) less than this
sifting.sigma_threshold = 4.0
# Ignore candidates with a coherent power less than this
sifting.c_pow_threshold = 100.0
# How close a candidate has to be to another candidate to                
# consider it the same candidate (in Fourier bins)
sifting.r_err = 1.1
# Shortest period candidates to consider (s)
sifting.short_period = 0.0005
# Longest period candidates to consider (s)
sifting.long_period = 15.0
# Ignore any candidates where at least one harmonic does exceed this power
sifting.harm_pow_cutoff = 8.0

def load_cands(folder='.'):
  # Try to read the .inf files first, as _if_ they are present, all of
  # them should be there.  (if no candidates are found by accelsearch
  # we get no ACCEL files...
  inffiles = glob.glob(folder + '/TEMP/' + globinf)
  candfiles = glob.glob(folder + '/TEMP/' + globaccel)
  # Check to see if this is from a short search
  if len(re.findall("_[0-9][0-9][0-9]M_" , inffiles[0])):
      dmstrs = [x.split("DM")[-1].split("_")[0] for x in candfiles]
  else:
      dmstrs = [x.split("DM")[-1].split(".inf")[0] for x in inffiles]
  dms = map(float, dmstrs)
  dms.sort()
  dmstrs = ["%.2f"%x for x in dms]

  # Read in all the candidates
  cands = sifting.read_candidates(candfiles)

  # Remove candidates that are duplicated in other ACCEL files
  if len(cands):
      cands = sifting.remove_duplicate_candidates(cands)

  # Remove candidates with DM problems
  if len(cands):
      cands = sifting.remove_DM_problems(cands, min_num_DMs, dmstrs, low_DM_cutoff)

  # Remove candidates that are harmonically related to each other
  # Note:  this includes only a small set of harmonics
  if len(cands):
      cands = sifting.remove_harmonics(cands)

  # Read dm and period of best candidates
  if len(cands):
      cands.sort(sifting.cmp_sigma)

      n_cands = 10  #Number of candidates to return
      dm = np.zeros(n_cands)
      p = np.zeros(n_cands)
      for i, goodcand in enumerate(cands):
        dm[i] = goodcand.DM
        p[i] = goodcand.p
        if i == n_cands-1: break
      
      return dm, p

  print "No periodic candidates found!"
  return None, None


def parser():
    # Command-line options
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                    description="The program plots candidates from the periodicity search.")
    parser.add_argument('-folder', help="Directory to store the output.", default='.')
    parser.add_argument('-fits', help="Full name (with path) of .fits file.", default='*.fits')
    return parser.parse_args()
  
  
  
if __name__ == '__main__':
  args = parser()
  basename = os.path.splitext(os.path.basename(args.fits))[0]
  dm_list, p_list = load_cands(args.folder)
  if not isinstance(dm,np.ndarray): exit()
  
  for idx,[dm, p] in enumerate(zip(dm_list, p_list)):
    if subprocess.call(['prepfold', '-nsub', '64', '-p', str(p), '-dm', str(dm), '-o', '{}/{}_periodic_cand_{}'.format(args.folder, basename, idx), args.fits]):
     print "Error in prepfold!"


