import subprocess
from glob import glob
import os
import argparse
import multiprocessing as mp

import psrchive
import pyfits
import numpy as np
import pandas as pd
from presto import psr_utils




def parser():
    # Command-line options
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                    description="The program create a sinlge-pulse psrchive file.")
    parser.add_argument('db_name', help="Name of the HDF5 database containing single pulses.", default='.')
    parser.add_argument('-fits_file', help="Name of fits file.", default='*')
    parser.add_argument('-obsPATH', help="Path of the observation output folder.", default='.')
    parser.add_argument('-par_file', help="Name of the parameter file for the puls.", default='*.par')
    parser.add_argument('-profile_bins', help="Number of bins within the profile.", default=4096, type=int)
    return parser.parse_args()
  
def burst_nsub(puls, profile_bins, fits_file=False):
  if not fits_file: fits_file = glob('/psr_archive/hessels/hessels/AO-FRB/pipeline_products/puppi_57614_C0531+33_0803/pulses/{}/*.fits'.format(puls.name))[0]
  with pyfits.open(fits_file) as fits:
    header = fits['SUBINT'].header + fits['Primary'].header
    #MJD seconds of the beginning of the fits file
    mjds_chop = header['STT_SMJD'] + header['STT_OFFS'] + header['NSUBOFFS'] * header['NSBLK'] * header['TBIN']
    #Time resolution of the fits file
    time_resolution = header['TBIN']
    
  burst_start = (puls.SMJD - mjds_chop) % 84000  #PRESTO reference with respect beginning of raw chopped fits file

  period = time_resolution * profile_bins  #Folding period for the archive
  nsub = int(burst_start / period)  #Subintegration number in the archive
  return nsub  



def dspsr(puls, par_file, fits_file, profile_bins=4096, parallel=False):
  nsub = burst_nsub(puls, profile_bins, fits_file)  #Subintegration number containing the pulse
  archive_name = os.path.splitext(fits_file)[0]
  
  #Create unfolded archive for timing
  if parallel:
    subprocess.call(['dspsr', '-t', str(mp.cpu_count()), '-K', '-b', str(profile_bins), '-s', '-A', '-E', par_file, '-O', archive_name+'_p1', fits_file])
  else:
    subprocess.call(['dspsr', '-N',  '-A', '-O', archive_name+'_timing', fits_file])
      
  #Fold the fits file to create the archives at two different phases
  if not (os.path.isfile(archive_name+'_p1.ar') or os.path.isfile(archive_name+'_p2.ar')):
    if parallel:
      subprocess.call(['dspsr', '-t', str(mp.cpu_count()), '-K', '-b', str(profile_bins), '-s', '-A', '-E', par_file, '-O', archive_name+'_p1', fits_file])
      subprocess.call(['dspsr', '-t', str(mp.cpu_count()), '-K', '-b', str(profile_bins), '-s', '-A', '-E', par_file, '-O', archive_name+'_p2', fits_file])      
    else:
      subprocess.call(['dspsr', '-K', '-b', str(profile_bins), '-s', '-A', '-E', par_file, '-O', archive_name+'_p1', fits_file])
      subprocess.call(['dspsr', '-K', '-b', str(profile_bins), '-s', '-A', '-E', par_file, '-O', archive_name+'_p2', fits_file])

    #Select the signle pulse at the two phases
    subprocess.call(['pam', '-m', '-x', '"{} {}"'.format(nsub, nsub), archive_name+'_p1.ar'])
    subprocess.call(['pam', '-m', '-x', '"{} {}"'.format(nsub-1, nsub-1), archive_name+'_p2.ar'])

  #Clean the archives
  if not (os.path.isfile(archive_name+'_p1.paz') or os.path.isfile(archive_name+'_p2.paz')):
    subprocess.call(['paz', '-e', 'paz', '-r', archive_name+'_p1.ar', archive_name+'_p2.ar'])
  
  #Create compressed archives
  if not (os.path.isfile(archive_name+'_p1.FTp') or os.path.isfile(archive_name+'_p2.FTp')):
    subprocess.call(['pam', '-e', 'FTp', '-FTp', archive_name+'_p1.paz', archive_name+'_p2.paz'])
  
  #Select the right phase
  if not (os.path.isfile(archive_name+'_p1.downsamp') or os.path.isfile(archive_name+'_p2.downsamp')):
    #Downsample at the closest factor
    downfact = puls.Downfact
    while profile_bins % downfact != 0:
      downfact -= 1
    subprocess.call(['pam', '-e', 'downsamp', '-b', str(downfact), archive_name+'_p1.FTp', archive_name+'_p2.FTp'])
    #Find maximum peak
    ar_p1 = psrchive.Archive_load(archive_name+'_p1.downsamp')
    ar_p2 = psrchive.Archive_load(archive_name+'_p2.downsamp')
    ar_p1.remove_baseline()
    ar_p2.remove_baseline()
    ar_p1 = ar_p1.get_data().squeeze()
    ar_p2 = ar_p2.get_data().squeeze()
    #Remove wrong-phase archive
    if ar_p1.max() > ar_p2.max(): 
      os.remove(archive_name+'_p2.ar')
      os.remove(archive_name+'_p2.paz')
      os.remove(archive_name+'_p2.FTp')
      os.remove(archive_name+'_p2.downsamp')
    else:
      os.remove(archive_name+'_p1.ar')
      os.remove(archive_name+'_p1.paz')
      os.remove(archive_name+'_p1.FTp')
      os.remove(archive_name+'_p1.downsamp')

  return

  

  
if __name__ == '__main__':
  args = parser()
  
  pulses = pd.read_hdf(args.db_name,'pulses')
  pulses = pulses[pulses.Pulse <= 1]
  
  if '*' in args.par_file: par_file = glob(args.par_file)[0]
  else: par_file = args.par_file

  fits_path = os.path.join(args.obsPATH, '{}', args.fits_file)
    
  for idx_p, puls in pulses.iterrows():
    fits_file = fits_path.format(idx_p) 
    if '*' in fits_file: fits_file = glob(fits_file)[0]
    dspsr(puls, par_file, fits_file, profile_bins=args.profile_bins)



#TO_TEST
burst_start = (puls.SMJD - mjds_chop) % 84000

offset = period / 2.

fits_duration = 3.0408704
starting_phase = (burst_start - offset) / fits_duration

DM_delay = psr_utils.delay_from_DM(560, 1380.78125-400) - psr_utils.delay_from_DM(560, 1380.78125+400)
duration = DM_delay + offset

dspsr -E 0531+33.par -A -O test -p starting_phase -T duration puppi_57614_C0531+33_0803_33.fits
#Check if it is fine for timing



