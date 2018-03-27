#####################################
#	FRB130628 arecibo search pipeline
#	
#	Kelly Gourdji 2017
#
#	Must be run in specific pointing's directory
#   labeled e.g. 4bit-p2030.20160702.FRB130628_0
#   pulses_extract must be git cloned to this directory.
#		python setup.py build_ext --inplace
#   Must create a subbanded_data subdirectory with soft links to fits files
#####################################

import os
import subprocess
from glob import glob #use for wildcards in subprocess?
import sys
import time
import pandas as pd
import numpy as np

def execute(command,working_dir=None):
	print command
	p = subprocess.Popen(command, shell=True, executable='/bin/bash',cwd=working_dir)
	p.wait() #careful. necessary in order to allow this child process to finish before script continues. look at communicate.
def make_prepsubband(infile,downsamp,lodm,dmstep,numdms,maskfile,base,beam,subband,working_dir):
	execute("prepsubband -nsub 120 -noscales -nooffsets -noweights -nobary -downsamp %d -lodm %f\
			-dmstep %f -numdms %d -zerodm -mask %s -o\
			 %s_b%ds%d_ZERO %s"\
			%(downsamp,lodm,dmstep,numdms,maskfile,base,beam,subband,infile), working_dir)


general_dir = "/data/gourdji/FRB130628_pipeline"#remove test later #Where everything pipeline related is stored 
script_dir = os.getcwd() + "/pulses_extract/src"
fits_dir = os.getcwd() + "/subbanded_data"

###test variables ###
#base = "4bit-p2030.20160702.FRB130628_0"
#beam = 0
#subband = 0
#cwd = os.getcwd()
#infile = cwd + "/" + glob("%s.b%ds%d*.fits"%(base,beam,subband))[0]
########################

### ONE SUBBAND VERSION ###
base_path = os.getcwd()
base = os.path.basename(base_path)#"4bit-p2030.20160702.FRB130628_1"
#pointing = int(base[-1])
beams = range(7)
subbands = range(2)
dmstep = 1.00
numdms = 96
downsamp = 3
lodm = 0.
dsubDM = 96.
calls = 7
calls = range(calls)

for beam in beams:
	for subband in subbands:
		print "NOW PROCESSING SUBBAND %d of BEAM %d"%(subband,beam)
		execute("mkdir %s_b%ds%d"%(base,beam,subband))
		outdir = '%s/%s_b%ds%d'%(base_path,base,beam,subband)		
		execute("mkdir %s/obs_data"%outdir)
		execute("mkdir %s/pulses"%outdir)
		execute("mkdir %s/TEMP"%outdir)
		#infile = glob("%s/%s.b%ds%d*.fits"%(fits_dir,base,beam,subband))[0]
		infile = glob("%s/*b%ds%d*.fits"%(fits_dir,beam,subband))[0]
		execute("rfifind -time 2.0 -psrfits -noscales -nooffsets -noweights -o %s_b%ds%d %s"%(base,beam,subband,infile), working_dir="%s/TEMP"%outdir)
		maskfile = glob("%s/TEMP/%s_b%ds%d*_rfifind.mask"%(outdir,base,beam,subband))[0]
		lodm = 0.
		if subband == 0:
			dmstep = 0.50
			numdms = 50
			downsamp = 2
			dsubDM = 25.
			calls = 19
			calls = range(calls)
			for call in calls:
				make_prepsubband(infile,downsamp,lodm,dmstep,numdms,maskfile,base,beam,subband,working_dir="%s/TEMP"%outdir)
				lodm +=  dsubDM
			execute("ls %s_b%ds%d_ZERO*.dat | xargs -n 1 single_pulse_search.py --noplot -m 100 -t 5.0 -b"%(base,beam,subband), working_dir="%s/TEMP"%outdir)
			dmstep = 1.00
			numdms = 50
			downsamp = 3
			dsubDM = 50.
			calls = 2
			calls = range(calls)
			for call in calls:
				make_prepsubband(infile,downsamp,lodm,dmstep,numdms,maskfile,base,beam,subband,working_dir="%s/TEMP"%outdir)
				lodm +=  dsubDM
			execute("ls %s_b%ds%d_ZERO_DM47[5-9]*.dat | xargs -n 1 single_pulse_search.py --noplot -m 70 -t 5.0 -b"%(base,beam,subband), working_dir="%s/TEMP"%outdir)
			execute("ls %s_b%ds%d_ZERO_DM4[8-9][0-9]*.dat | xargs -n 1 single_pulse_search.py --noplot -m 70 -t 5.0 -b"%(base,beam,subband), working_dir="%s/TEMP"%outdir)
			execute("ls %s_b%ds%d_ZERO_DM5[0-9][0-9]*.dat | xargs -n 1 single_pulse_search.py --noplot -m 70 -t 5.0 -b"%(base,beam,subband), working_dir="%s/TEMP"%outdir)

		else:
			dmstep = 0.30
			numdms = 50
			downsamp = 2
			dsubDM = 15.
			calls = 21
			calls = range(calls)
			for call in calls:
				make_prepsubband(infile,downsamp,lodm,dmstep,numdms,maskfile,base,beam,subband,working_dir="%s/TEMP"%outdir)
				lodm +=  dsubDM
			execute("ls %s_b%ds%d_ZERO*.dat | xargs -n 1 single_pulse_search.py --noplot -m 100 -t 5.0 -b"%(base,beam,subband), working_dir="%s/TEMP"%outdir)
			dmstep = 0.50
			numdms = 50
			downsamp = 3
			dsubDM = 25.
			calls = 9
			calls = range(calls)
			for call in calls:
				make_prepsubband(infile,downsamp,lodm,dmstep,numdms,maskfile,base,beam,subband,working_dir="%s/TEMP"%outdir)
				lodm +=  dsubDM
			execute("ls %s_b%ds%d_ZERO_DM31[5-9]*.dat | xargs -n 1 single_pulse_search.py --noplot -m 70 -t 5.0 -b"%(base,beam,subband), working_dir="%s/TEMP"%outdir)	
			execute("ls %s_b%ds%d_ZERO_DM3[2-9][0-9]*.dat | xargs -n 1 single_pulse_search.py --noplot -m 70 -t 5.0 -b"%(base,beam,subband), working_dir="%s/TEMP"%outdir)	
			execute("ls %s_b%ds%d_ZERO_DM4[0-9][0-9]*.dat | xargs -n 1 single_pulse_search.py --noplot -m 70 -t 5.0 -b"%(base,beam,subband), working_dir="%s/TEMP"%outdir)	
			execute("ls %s_b%ds%d_ZERO_DM5[0-3][0-9]*.dat | xargs -n 1 single_pulse_search.py --noplot -m 70 -t 5.0 -b"%(base,beam,subband), working_dir="%s/TEMP"%outdir)	
			dmstep = 1.00
			numdms = 50
			downsamp = 6
			dsubDM = 50.
			calls = 1
			make_prepsubband(infile,downsamp,lodm,dmstep,numdms,maskfile,base,beam,subband,working_dir="%s/TEMP"%outdir)
			execute("ls %s_b%ds%d_ZERO_DM5[4-8][0-9].dat | xargs -n 1 single_pulse_search.py --noplot -m 30 -t 5.0 -b"%(base,beam,subband), working_dir="%s/TEMP"%outdir)

		#execute("ls %s_b%ds%d_ZERO*.dat | xargs -n 1 single_pulse_search.py --noplot -m 150 -t 5.0 -b"%(base,beam,subband), working_dir="%s/TEMP"%outdir)
		execute("single_pulse_search.py -t 10 %s_b%ds%d_ZERO*singlepulse"%(base,beam,subband), working_dir="%s/TEMP"%outdir)
		
		#execute("mv %s_b%ds%d_ZERO* %s_b%ds%d_TEST_proc"%(base,beam,subband,base,beam,subband))
		#execute("mv %s_b%ds%d_rfifind.* %s_b%ds%d_TEST_proc"%(base,beam,subband,base,beam,subband))
		#this was used to perform actions on files already moved out of the "test" directory.
		#if not os.path.isfile("%s/%s_b%ds%d_TEST_proc/%s_b%ds%d_ZERO_SNR5_singlepulse.ps"\
														#%(path,base,beam,subband,base,beam,subband)):
			#execute("mv %s/%s_b%ds%d_TEST_proc/%s_b%ds%d_ZERO_singlepulse.ps\
			 #%s/%s_b%ds%d_TEST_proc/%s_b%ds%d_ZERO_SNR5_singlepulse.ps"%(path,base,beam,subband,base,beam,subband,path,base,beam,subband,base,beam,subband))
		#execute("single_pulse_search.py -t 10 %s/%s_b%ds%d_TEST_proc/*.singlepulse"%(path,base,beam,subband))
		#subprocess.Popen("single_pulse_search.py -t 10 *.singlepulse", shell=True, cwd='%s_b%ds%d_TEST_proc'%(base,beam,subband))


##### STEP 2: ONCE SINGLEPULSE FILES ARE CREATED FOR EACH DM #####
		execute("python %s/pulses_extract.py -db_name %s_b%ds%d_proc.hdf5 -fits %s\
		 		-store_events -idL %s_b%ds%d_ZERO_DM -store_dir %s \
					-beam_num %d -group_num %d -plot_statistics -parameters_id FRB130628_Alfa_s%d > /dev/null"\
				%(script_dir,base,beam,subband,infile,base,beam,subband,outdir,beam,subband,subband), working_dir="%s/TEMP"%outdir)
		execute("cp %s/TEMP/%s_b%ds%d_ZERO_singlepulse.ps %s/obs_data"%(outdir,base,beam,subband,outdir))
		execute("cp %s/TEMP/%s_b%ds%d_ZERO_DM470.00.dat %s/obs_data"%(outdir,base,beam,subband,outdir))
		execute("cp %s/TEMP/%s_b%ds%d_ZERO_DM470.00.inf %s/obs_data"%(outdir,base,beam,subband,outdir))
		execute("cp %s/%s_b%ds%d_proc.hdf5 %s "%(outdir,base,beam,subband,base_path))
		#execute("rm -rf %s/TEMP"%outdir)
		sys.exit() #REMOVE THIS. This is so that only one file is processed.

execute("mkdir pulses")
execute("python %s/pulses_extract.py -beam_comparison *_proc.hdf5"%script_dir)
execute("python %s/pulses_extract.py -fits %s -pulses_database -store_dir %s/pulses\
   -plot_pulses -plot_statistics -parameters_id FRB130628_Alfa_s0"%(script_dir,fits_dir,base)) #doesn't matter which subband param ID to use.
execute("python %s/beams_plot.py > /dev/null"%script_dir)
#Use this for debugging so can print messages within pulses_extract.py 
#remove > dev/null since otherwise won't print output on command line
#process = subprocess.Popen("python %s/pulses_extract.py -db_name %s_b%ds%d_SinglePulses.hdf5 -fits %s\
#		-store_events -idL %s_b%ds%d_ZERO_DM -store_dir pulses \
#			-plot_pulses -plot_statistics -parameters_id FRB130628_Alfa_s%d"\ 
#		%(script_dir,base,beam,subband,infile,base,beam,subband,subband), shell=True, stdout=subprocess.PIPE)
#stdout = process.communicate()[0]
#print "STDOUT: %s"%stdout
		
"""
#concatenate all singlepulse files to create one master file
#"find -maxdepth 1 -name "%s_b%ds%d_ZERO_DM*" | xargs -n 1 "%(base,beam,subband)
execute("cat %s_b%ds%d_ZERO_DM*.singlepulse |\
 head -n 1 > %s_b%ds%d.singlepulse"%(base,beam,subband))
execute("cat %s_b%ds%d_ZERO_DM*.singlepulse | grep -v \#\ >> %s_b%ds%d.singlepulse"\
	%(base,beam,subband,base,beam,subband))

path = "/psr_archive/hessels/hessels/AO-FRB/P3055/New_data/FRB130628_pipeline/test"
base = "4bit-p2030.20160702.FRB130628_0"
beam = 0
subband = 0
#create events database
events = pd.read_csv("p2030.20160702.FRB130628_0_b0s0.singlepulse", dtype=np.float64,\
						delim_whitespace=True, usecols=[1,2,3,4,5])
events.rename(columns={'(s)':'Sample', 'Sample':'Downfact'}, inplace=True)
events.index.name = 'idx' #gives index column a name
events['Pulse'] = 0 #Place holder for grouping events (Daniele's script in C)
events.Pulse = events.Pulse.astype(np.int32)
events.Downfact = events.Downfact.astype(np.int16)
events.Sample = events.Sample.astype(np.int32)
events.sort(['DM','Time'],inplace=True) #sort_values later version of pandas... see which one to use.
C_Funct.Get_Group(events.DM.values, events.Sigma.values, events.Time.values, events.Pulse.values, \
					5., 20e-3, 1.) #label events by group ('Pulses' value). Use these in groupby to actually do the grouping.
events = events[events.Pulse >= 0] 
store = pd.HDFStore('%s_b%ds%d_SinglePulses.hdf5'%(base,beam,subband))
store.append('events',events)
store.close()

#create pulses database, append to hdfstore
gb = events.groupby('Pulse',sort=False) #actually group the events that are in the same group ('Pulse' value)
pulses = events.loc[gb.Sigma.idxmax()] #set the metadata of each group to be those corresponding to the event with highest Sigma.
pulses.index = pulses.Pulse #let the Pulse val be the new index (makes sense since the old indices are all out of order after the grouping step)
#store.append('pulses',pulses)


#removing last 10s of data?
#C_funct.Get_Group?
#def events_database(args, header):
  #Create events database
  #sp_files = glob.glob(os.path.join(args.folder,'{}*.singlepulse'.format(args.idL)))
  #events = pd.concat(pd.read_csv(f, delim_whitespace=True, dtype=np.float64) for f in sp_files)
  #events.reset_index(drop=True, inplace=True)
  #events.columns = ['DM','Sigma','Time','Sample','Downfact','a','b']
  #events = events.ix[:,['DM','Sigma','Time','Sample','Downfact']] #nec?
  #events.index.name = 'idx' #necessary? doesn't come back.
  #events['Pulse'] = 0
  #events.Pulse = events.Pulse.astype(np.int32)
  #events.Downfact = events.Downfact.astype(np.int16)
  #events.Sample = events.Sample.astype(np.int32)
  #events.sort_values(['DM','Time'],inplace=True)
   #Remove last 10s of data
  #obs_length = header['NSBLK'] * header['NAXIS2'] * header['TBIN']
  #events = events[events.Time < obs_length-10.]

  #C_Funct.Get_Group(events.DM.values, events.Sigma.values, events.Time.values, events.Pulse.values, 
                    #args.events_dDM, args.events_dt, args.DM_step) #what does this do?

  #events = events[events.Pulse >= 0]
"""
 




