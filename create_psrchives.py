import subprocess
from glob import glob
import os
import argparse
import multiprocessing as mp
import shutil

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

def read_fits(fits_file):
  with pyfits.open(fits_file) as fits:
    header = fits['SUBINT'].header + fits['Primary'].header
    #MJD seconds of the beginning of the fits file
    mjds_chop = header['STT_SMJD'] + header['STT_OFFS'] + header['NSUBOFFS'] * header['NSBLK'] * header['TBIN']
    #Time resolution of the fits file
    time_resolution = header['TBIN']
    freq_c = header['OBSFREQ']
    bandwidth = abs(header['OBSBW'])
    
  return {'time_resolution': time_resolution, 'freq_c': freq_c, 'bandwidth': bandwidth}

  
def dspsr(puls, par_file, fits_file, profile_bins=4096, parallel=False):
  if not os.path.isfile(archive_name + '.ar'):
    archive_name = os.path.splitext(fits_file)[0]
    readfile = read_fits(fits_file)
    period = readfile['time_resolution'] * profile_bins
    
    #Fold the fits file to create single-pulse archives
    subprocess.call(['dspsr', '-D', str(puls.DM), '-K', '-b', str(profile_bins), '-s', '-E', par_file, fits_file])
    
    #Lists of archive names and starting times (s)
    archive_list = np.array(glob('pulse_*.ar'))
    archive_time_list = np.array([psrchive.Archive_load(ar).start_time().get_secs() + psrchive.Archive_load(ar).start_time().get_fracsec() for ar in archive_list])
    idx_sorted = np.argsort(archive_list)
    archive_list = archive_list[idx_sorted]
    archive_time_list = archive_time_list[idx_sorted]
    
    #Find archive where dispersed pulse would start
    start_dispersed_puls = puls.SMJD - archive_time_list
    idx_puls = np.where( (start_dispersed_puls > 0) & (start_dispersed_puls < period))[0]
    
    #Check that puls is centered
    phase = start_dispersed_puls[idx_puls] / period
    if abs(phase - 0.5) > 0.5: 
      for ar in archive_list: os.remove(ar)
      subprocess.call(['dspsr', '-S', '0.5', '-D', str(puls.DM), '-K', '-b', str(profile_bins), '-s', '-E', par_file, '-O', archive_name, fits_file])
      archive_list = np.array(glob('pulse_*.ar'))
      archive_time_list = np.array([psrchive.Archive_load(ar).start_time().get_secs() + psrchive.Archive_load(ar).start_time().get_fracsec() for ar in archive_list])
      idx_sorted = np.argsort(archive_list)
      archive_list = archive_list[idx_sorted]
      archive_time_list = archive_time_list[idx_sorted]
      start_dispersed_puls = puls.SMJD - archive_time_list
      idx_puls = np.where( (start_dispersed_puls > 0) & (start_dispersed_puls < period))[0]
    
    #Delay between maximum and central frequencies
    DM_delay = psr_utils.delay_from_DM(puls.DM, readfile['freq_c']) - psr_utils.delay_from_DM(puls.DM, readfile['freq_c'] + readfile['bandwidth'] / 2. )
    
    n_puls = int(DM_delay / period)
    idx_puls += n_puls
    
    shutil.copyfile(archive_list[idx_puls], archive_name + '.ar')
    for ar in archive_list: os.remove(ar)
    
  #Clean the archive
  if not os.path.isfile(archive_name + '.paz'):
    subprocess.call(['paz', '-e', 'paz', '-r', archive_name + '.ar'])
  
  #Create compressed archive
  if not os.path.isfile(archive_name + '.FTp'):
    subprocess.call(['pam', '-e', 'FTp', '-FTp', archive_name + '.paz'])
  
  #Create downsampled archive at the closest factor
  if not os.path.isfile(archive_name + '.downsamp'):
    downfact = puls.Downfact
    while profile_bins % downfact != 0:
      downfact -= 1
    subprocess.call(['pam', '-e', 'downsamp', '-b', str(downfact), archive_name + '.FTp'])
    
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


