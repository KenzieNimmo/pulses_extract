import argparse
import glob
import os
import sys

import pandas as pd
import numpy as np
from scipy import special
import pyfits

from src import C_Funct
import auto_waterfaller
from extract_psrfits_subints import extract_subints_from_observation


DM_low = 461.           #Lowest de-dispersed value
DM_high = 660.          #Highest de-dispersed value
SNR_peak_min = 8.       #Minumum peak SNR value
SNR_min = 6.            #Minumum SNR value
Downfact_max = 100      #Maximum Downfact value



def main():
  def parser():
    # Command-line options
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                    description="The program groupes events from single_pulse_search.py in pulses.")
    parser.add_argument('-fits', help="Filename of .fits file.", default='*.fits')
    parser.add_argument('-idL', help="Basename of .singlepulse files.", default='*')
    parser.add_argument('-folder', help="Path of the folder containig the .singlepulse files.", default='.')
    parser.add_argument('-store_events', help="Store events in a SinglePulses.hdf5 database.", action='store_true')
    parser.add_argument('-events_dt', help="Duration in sec within two events are related to the same pulse.", default=20e-3,
                        type=int)
    parser.add_argument('-events_dDM', help="Number of DM steps within two events are related to the same pulse.", 
                        default=5, type=float)
    parser.add_argument('-DM_step', help="Value of the DM step between timeseries.", 
                        default=1., type=float)
    parser.add_argument('-events_database', help="Load events from this database.", default='')
    parser.add_argument('-pulses_database', help="Load pulses from this database.", default='')
    parser.add_argument('-store_dir', help="Path of the folder to store the output.", default='.')
    parser.add_argument('-plot_pulses', help="Save plots of detected pulses.", action='store_true')
    parser.add_argument('-extract_raw', help="Extract raw data specified in this path around detected pulses.", default='')
    parser.add_argument('-pulses_checked', help="Path of a text file containig a list of pulse identifiers to label as RFI.", default='')
    parser.add_argument('-plot_statistics', help="Produce plots with statistics of the pulses.", action='store_true')
    return parser.parse_args()
  args = parser()
  
  if args.pulses_database: pulses = pd.read_hdf(args.pulses_database,'pulses')
  else: 
    header = fits_header(args.fits)
    pulses = pulses_database(args, header)
    store = pd.HDFStore('{}/SinglePulses.hdf5'.format(args.store_dir), 'a')
    store.append('pulses',pulses)
    #store.append('pulses_bu',pulses) #Create a back up table in the database
    store.close()
    pulses.to_csv('{}/pulses_list.txt'.format(args.store_dir), sep='\t', cols=['Pulse',], header=['Rank',], index_label='#PulseID')

  if args.pulses_checked: 
    pulses_checked(pulses, args.pulses_checked)
    store = pd.HDFStore('{}/SinglePulses.hdf5'.format(args.store_dir), 'r+')
    store.remove('pulses')
    store.append('pulses',pulses)
    store.close()
    
  if args.plot_pulses: 
    real_pulses = pulses[pulses.Pulse < 2]
    auto_waterfaller.main(args.fits, np.array(real_pulses.Time), np.array(real_pulses.DM), np.array(real_pulses.Sigma), \
                                             np.array(real_pulses.Duration), top_freq=real_pulses.top_Freq.iloc[0], directory=args.store_dir)
  
  if args.extract_raw: 
    real_pulses = pulses[pulses.Pulse < 2]
    extract_subints_from_observation(args.extract_raw, args.store_dir+'/fits', np.array(real_pulses.Time), -2, 8)
  
  if args.plot_statistics: plot_statistics(pulses[pulses.Pulse == 0])

  return


def events_database(args, header):
  #Create events database
  sp_files = glob.glob("{}/{}*.singlepulse".format(args.folder, args.idL))
  events = pd.concat(pd.read_csv(f, delim_whitespace=True, dtype=np.float64) for f in sp_files)
  events.reset_index(drop=True, inplace=True)
  events.columns = ['a','DM','Sigma','Time','Sample','Downfact','b']
  events = events.ix[:,['DM','Sigma','Time','Sample','Downfact']]
  events.index.name = 'idx'
  events['Pulse'] = 0
  events.Pulse = events.Pulse.astype(np.int32)
  events.Downfact = events.Downfact.astype(np.int16)
  events.Sample = events.Sample.astype(np.int32)
  events.sort(['DM','Time'],inplace=True)

  if args.store_events:
    store = pd.HDFStore('{}/SinglePulses.hdf5'.format(args.store_dir), 'w')
    store.append('events',events,data_columns=['Pulse','SAP','BEAM','DM','Time'])
    store.close()
  
  #Remove last 10s of data
  obs_length = header['NSBLK'] * header['NAXIS2'] * header['TBIN']
  pulses = pulses[pulses.Time < obs_length-10.]

  C_Funct.Get_Group(events.DM.values, events.Sigma.values, events.Time.values, events.Pulse.values, 
                    args.events_dDM, args.events_dt, args.DM_step)
  
  return events[events.Pulse >= 0]
  
  
def pulses_database(args, header, events=None):
  #Create pulses database
  if args.events_database: events = pd.read_hdf(args.events_database,'events')
  elif not isinstance(events, pd.DataFrame): events = events_database(args, header)
  gb = events.groupby('Pulse',sort=False)
  pulses = events.loc[gb.Sigma.idxmax()]
  pulses.index = pulses.Pulse
  pulses.index.name = None
  pulses = pulses.loc[:,['DM','Sigma','Time','Sample','Downfact']]
  pulses.index.name = 'idx'
  pulses['IMJD'] = header['STT_IMJD']
  pulses['SMJD'] = header['STT_SMJD'] + pulses.Time
  pulses['Duration'] = pulses.Downfact * header['TBIN']
  pulses['top_Freq'] = header['OBSFREQ'] + abs(header['OBSBW']) / 2.
  pulses['Pulse'] = 0
  pulses.Pulse = pulses.Pulse.astype(np.int8)
  pulses['dDM'] = (gb.DM.max() - gb.DM.min()) / 2.
  pulses.dDM=pulses.dDM.astype(np.float32)
  pulses['dTime'] = (gb.Time.max() - gb.Time.min()) / 2.
  pulses.dTime=pulses.dTime.astype(np.float32)
  pulses['N_events'] = gb.DM.count()
  pulses.N_events = pulses.N_events.astype(np.int16)

  pulses = pulses[pulses.N_events > 5]
  pulses = pulses[pulses.Sigma >= SNR_peak_min]
  pulses = pulses[pulses.Downfact <= Downfact_max]
  
  RFIexcision(events, pulses)
  
  pulses.sort('Sigma', ascending=False, inplace=True)
  return pulses[pulses.Pulse == 0]
  

def RFIexcision(events, pulses):
  events = events[events.Pulse.isin(pulses.index)]
  events.sort('DM',inplace=True)
  gb = events.groupby('Pulse')
  pulses.sort_index(inplace=True)
  
  #Remove flat SNR pulses. Minimum ratio to have weakest pulses with SNR = 8
  pulses.Pulse[pulses.Sigma / gb.Sigma.min() <= SNR_peak_min / SNR_min] += 1
  
  #Remove flat duration pulses. Minimum ratio to have weakest pulses with SNR = 8 (from Eq.6.21 of Pulsar Handbook)
  pulses.Pulse[gb.Downfact.max() / pulses.Downfact < (SNR_peak_min / SNR_min)**2] += 1
  
  #Remove pulses peaking near the DM edges
  DM_frac = (DM_high - DM_low) * 0.2  #Remove 5% of DM range from each edge
  pulses.Pulse[(pulses.DM < DM_low+DM_frac) | (pulses.DM > DM_high-DM_frac)] += 1
  
  #Remove pulses intersecting half the maximum SNR different than 2 or 4 times
  def crosses(sig):
    diff = sig - (sig.max() + sig.min()) / 2.
    count = np.count_nonzero(np.diff(np.sign(diff)))
    return (count != 2) & (count != 4) & (count != 6)
  pulses.Pulse += gb.apply(lambda x: crosses(x.Sigma)).astype(np.int8)
  
  #Remove weaker pulses within a temporal window
  def simultaneous(p):                            
    puls = pulses.Pulse[np.abs(pulses.Time-p.Time) < 0.02]
    if puls.shape[0] == 1: return 0
    if p.name == puls.index[0]: return 0
    else: return 1
  pulses.Pulse += pulses.apply(lambda x: simultaneous(x), axis=1)
  
  return


def fits_header(filename):
  with pyfits.open(filename,memmap=True) as fits:
    header = fits['SUBINT'].header + fits['PRIMARY'].header
  return header
  
  
def pulses_checked(pulses, filename):
  RFI_list = np.genfromtxt(filename, dtype=int).T

  #print "Folowing pulses will be marked as RFI: ", RFI_list[0,RFI_list[1]==0]
  #sys.stdout.write("Proceed? [y/n]")
  #choice = raw_input().lower()
  #if not (choice == 'y') | (choice == 'yes'): 
    #print "Aborting..."
    #return
  
  pulses.Pulse.loc[RFI_list[0]] = RFI_list[1]
  return pulses


def plot_statistics(pulses):
  auto_waterfaller.histogram(pulses.DM)
  auto_waterfaller.histogram(pulses.Sigma)

  

if __name__ == '__main__':
  main()
