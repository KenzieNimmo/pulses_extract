import argparse
import os
from glob import glob

import pyfits

from extract_psrfits_subints import extract_subints_from_observation
from create_psrchives import dspsr



RAW_DIR="/psr_archive/hessels/hessels/AO-FRB/raw_data"



def parser():
    # Command-line options
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                    description="The program extracts chunks of data from fits file into psarchives around a specified time list.")
    parser.add_argument('obsID', help="ID of the observation to process.")
    parser.add_argument('time_list', help="List of times to extract, either in seconds from beginning of observation or MJD.", nargs='+', type=float)
    parser.add_argument('-output_folder', help="Path of the folder to store the bursts.", default='.')
    parser.add_argument('-time_window', help="Time window to extract around the pulse.", default=0.04194304, type=float)
    return parser.parse_args()
 


def main():
  args = parser()
  
  if max(args.time_list) < 36000:
    print "Seconds from the start of observation inserted."
    t_sec = args.time_list
    t_mjd, mjd = secs2mjd(t_sec, args.obsID)
    print "Equivalent MJDs: {}".format(mjd + t_mjd / 24. / 3600.)
  else:
    print "MJD values inserted."
    t_mjd = args.time_list
    t_sec = mjd2sec(t_mjd, args.obsID)
    print "Equivalent seconds from the start of observation: {}".format(t_sec)
    
  #Output directory
  out_dir = args.output_folder
  if not os.path.exists(out_dir):
    os.makedirs(out_dir)
  
  for i in len(args.time_list):
    t_folder = os.path.join(out_dir,'%.2f'%t)
    os.makedirs(t_folder)
    extract_archive(t_sec[i], t_mjd[i], args.obsID, out_dir=t_folder, width=args.time_window)
  
  return


  
def sec2mjd(t_sec, obsID):
  obs_start, mjd = start_of_obs(obsID)
  return [t + obs_start for t in t_sec], mjd

def mjd2sec(t_mjd, obsID):
  obs_start, mjd = start_of_obs(obsID)
  return [t - obs_start for t in t_mjd] 
  
def start_of_obs(obsID):
  raw_files = os.path.join(RAW_DIR, obsID)
  raw_fits_files = glob(raw_files+'*')
  file_starts = []
  for obs in raw_fits_files:
    with pyfits.open(obs,memmap=True) as fits:
      header = fits['SUBINT'].header + fits['PRIMARY'].header
      file_starts.append(header['STT_SMJD'] + header['STT_OFFS'] + header['NSUBOFFS'] * header['NSBLK'] * header['TBIN'])
  mjd = header['STT_IMJD']
  return min(file_starts), mjd
  
  
  
def extract_archive(t_sec, t_mjd, obsID, out_dir='.', width=0.04194304):
  #Create fits file
  raw_files = os.path.join(RAW_DIR, obsID)
  extract_subints_from_observation(raw_files, out_dir, [t_sec,], -2, 8, pulseID='%.2f'%t_sec)
  
  #Create psrchive
  fits_file = os.path.join(out_dir,'%.2f.fits'%t_sec)
  dspsr(fits_file, SMJD=t_mjd, width=width)
  
  return



if __name__ == '__main__':
  main()
