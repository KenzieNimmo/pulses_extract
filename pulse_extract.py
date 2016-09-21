import argparse
import glob
import os

import pandas as pd
import numpy as np
from scipy import special
import pyfits

from src import C_Funct
import auto_waterfaller

def main():
  def parser():
    # Command-line options
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                    description="The program groupes events from single_pulse_search.py in pulses.")
    parser.add_argument('fits', help="Filename of .fits file.", default='')
    parser.add_argument('-idL', help="Basename of .singlepulse files.", default='')
    parser.add_argument('-folder', help="Path of the folder containig the .singlepulse files.", default='')
    parser.add_argument('-store_events', help="Store events in a SinglePulses.hdf5 database.",action='store_true')
    parser.add_argument('-remove_singlepulse', help="Remove .singlepulse files at the end of the processing.",action='store_true')
    parser.add_argument('-events_dt', help="Duration in sec within two events are related to the same pulse.", default=20e-3,
                        type=int)
    parser.add_argument('-events_dDM', help="Number of DM steps within two events are related to the same pulse.", 
                        default=5, type=float)
    parser.add_argument('-DM_step', help="Value of the DM step between timeseries.", 
                        default=1., type=float)
    parser.add_argument('-events_database', help="Path of the eventual SinglePulses.hdf5 database to load events.")
    parser.add_argument('-pulses_database', help="Path of the eventual SinglePulses.hdf5 database to load pulses.")
    parser.add_argument('-store_dir', help="Path of the folder to store the SinglePulses.hdf5 database.", default='')
    return parser.parse_args()
  args = parser()
  if args.folder != '': args.folder += '/'
  if args.store_dir != '': args.store_dir += '/'

  header = fits_header(args.fits)
  
  if isinstance(args.pulses_database, str): pulses = pd.read_hdf(args.pulses_database + '/SinglePulses.hdf5','pulses')
  else: 
    pulses = pulses_database(args, header)
    store = pd.HDFStore('{}{}/SinglePulses.hdf5'.format(args.folder, args.store_dir), 'a')
    store.append('pulses',pulses)
    store.close()
  
  auto_waterfaller.main(args.fits, np.array(pulses.Time), np.array(pulses.DM), directory=args.store_dir)
  
  return


  
def RFIexcision(events, pulses):
  events = events[events.Pulse.isin(pulses.index)]
  events.DM = events.DM.astype(np.float64)
  events.Sigma = events.Sigma.astype(np.float64)
  events.sort('DM',inplace=True)
  gb = events.groupby('Pulse')
  pulses.sort_index(inplace=True)
  
  #Remove flat SNR pulses. Minimum ratio to have weakest pulses with SNR = 8
  pulses.Pulse[pulses.Sigma / gb.Sigma.min() <= 8./6.] += 1
  
  #Remove flat duration pulses. Minimum ratio to have weakest pulses with SNR = 8 (from Eq.6.21 of Pulsar Handbook)
  pulses.Pulse[gb.Downfact.max() / pulses.Downfact < 1.8] += 1
  
  #Remove pulses peaking near the DM edges
  pulses.Pulse[(pulses.DM < 471.) | (pulses.DM > 651.)] += 1
  
  def crosses(sig):
    diff = sig - (sig.max() + sig.min()) / 2.
    count = np.count_nonzero(np.diff(np.sign(diff)))
    return (count != 2) & (count != 4)
  
  pulses.Pulse += gb.apply(lambda x: crosses(x.Sigma)).astype(np.int8)
    
  return


def events_database(args):
  #Create events database
  sp_files = glob.glob("{}{}*.singlepulse".format(args.folder, args.idL))
  events = pd.concat(pd.read_csv(f, delim_whitespace=True, dtype=np.float32) for f in sp_files)
  events.reset_index(drop=True, inplace=True)
  events.columns = ['a','DM','Sigma','Time','Sample','Downfact','b']
  events = events.ix[:,['DM','Sigma','Time','Sample','Downfact']]
  events.index.name = 'idx'
  events['Pulse'] = 0
  events.Pulse = events.Pulse.astype(np.int32)
  events.sort(['DM','Time'],inplace=True)
  C_Funct.Get_Group(events.DM.values, events.Sigma.values, events.Time.values, events.Pulse.values, 
                    args.events_dDM, args.events_dt, args.DM_step)
  if args.store_events:
    store = pd.HDFStore('{}{}SinglePulses.hdf5'.format(args.folder, args.store_dir), 'a')
    store.append('events',events,data_columns=['Pulse','SAP','BEAM','DM','Time'])
    store.close()
  if args.remove_singlepulse:
    for f in sp_files:
      #os.remove(f)
      print "Files would have been removed!!!"
      pass
  return events[events.Pulse >= 0]
  
  
def pulses_database(args, header, events=None):
  #Create pulses database
  if isinstance(args.events_database, str): events = pd.read_hdf(args.events_database + '/SinglePulses.hdf5','events')
  elif not isinstance(events, pd.DataFrame): events = events_database(args)
  gb = events.groupby('Pulse',sort=False)
  pulses = events.loc[gb.Sigma.idxmax()]  
  pulses.index = pulses.Pulse
  pulses.index.name = None
  pulses = pulses.loc[:,['DM','Sigma','Time','Sample','Downfact']]
  pulses.index.name = 'idx'
  pulses['MJD'] = header['STT_IMJD'] + (header['STT_SMJD'] + pulses.Time) / 86400.
  pulses['Duration'] = pulses.Downfact * header['TBIN']
  pulses['Pulse'] = 0
  pulses.Pulse = pulses.Pulse.astype(np.int8)
  pulses['dDM'] = (gb.DM.max() - gb.DM.min()) / 2.
  pulses.dDM=pulses.dDM.astype(np.float32)
  pulses['dTime'] = (gb.Time.max() - gb.Time.min()) / 2.
  pulses.dTime=pulses.dTime.astype(np.float32)
  pulses['N_events'] = gb.DM.count()
  pulses.N_events = pulses.N_events.astype(np.int16)

  pulses = pulses[pulses.N_events > 5]
  pulses = pulses[pulses.Sigma >= 8.]
  pulses = pulses[pulses.Downfact <= 100]
  
  RFIexcision(events, pulses)  

  pulses.sort('Sigma', ascending=False, inplace=True)
  return pulses[pulses.Pulse == 0]
  
   
  
def fits_header(filename):
  with pyfits.open(filename,memmap=True) as fits:
    header = fits['SUBINT'].header + fits['PRIMARY'].header
  return header
  

if __name__ == '__main__':
  main()
