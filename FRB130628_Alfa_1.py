#####################################
#	FRB130628 arecibo search pipeline
#	
#	Kelly Gourdji 2017
#
#	Must be run in same dir as pulses_extract
#####################################

import os
import subprocess
from glob import glob #use for wildcards in subprocess?
import sys
import time
import pandas as pd
import numpy as np
#sys.path.insert(0, "/psr_archive/hessels/hessels/AO-FRB/pipeline_products/pulses_extract")
#from src import C_Funct

def execute(command):
	print command
	subprocess.call(command, shell=True)

def make_prepsubband(infile,downsamp,lodm,dmstep,numdms,maskfile,base,beam,subband):
	execute("prepsubband -nsub 120 -noscales -nooffsets -downsamp %d -lodm %f\
			-dmstep %f -numdms %d -zerodm -mask %s -o\
			 %s_b%ds%d_ZERO %s"\
			%(downsamp,lodm,dmstep,numdms,maskfile,base,beam,subband,infile))

cwd = os.getcwd()
script_dir = cwd + "/pulses_extract/src"
### test variables ###
#base = "4bit-p2030.20160702.FRB130628_0"
#beam = 0
#subband = 0
#infile = cwd + "/" + glob("%s.b%ds%d*.fits"%(base,beam,subband))[0]
########################

### ONE SUBBAND VERSION ###
path = "/data/gourdji/FRB130628_pipeline/test" #come back to this
base = "4bit-p2030.20160702.FRB130628_0"
#beams = range(7)
#subbands = range(2)
beams = [0]
subbands = [1]
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
		infile = cwd + "/" + glob("%s.b%ds%d*.fits"%(base,beam,subband))[0]
		execute("rfifind -time 2.0 -psrfits -noscales -nooffsets -o %s_b%ds%d %s"%(base,beam,subband,infile))
		maskfile = glob("%s_b%ds%d*_rfifind.mask"%(base,beam,subband))[0]
		lodm = 0.
		if subband == 0:
			dmstep = 0.50
			numdms = 50
			downsamp = 2
			dsubDM = 25.
			calls = 19
			calls = range(calls)
			for call in calls:
				make_prepsubband(infile,downsamp,lodm,dmstep,numdms,maskfile,base,beam,subband)
				lodm +=  dsubDM
			dmstep = 1.00
			numdms = 50
			downsamp = 3
			dsubDM = 50.
			calls = 2
			calls = range(calls)
			for call in calls:
				make_prepsubband(infile,downsamp,lodm,dmstep,numdms,maskfile,base,beam,subband)
				lodm +=  dsubDM	

		else:
			dmstep = 0.30
			numdms = 50
			downsamp = 2
			dsubDM = 15.
			calls = 21
			calls = range(calls)
			for call in calls:
				make_prepsubband(infile,downsamp,lodm,dmstep,numdms,maskfile,base,beam,subband)
				lodm +=  dsubDM
			dmstep = 0.50
			numdms = 50
			downsamp = 3
			dsubDM = 25.
			calls = 9
			calls = range(calls)
			for call in calls:
				make_prepsubband(infile,downsamp,lodm,dmstep,numdms,maskfile,base,beam,subband)
				lodm +=  dsubDM	
			dmstep = 1.00
			numdms = 50
			downsamp = 6
			dsubDM = 50.
			calls = 1
			make_prepsubband(infile,downsamp,lodm,dmstep,numdms,maskfile,base,beam,subband)

		execute("ls %s_b%ds%d_ZERO*.dat | xargs -n 1 single_pulse_search.py --noplot -m 70 -t 5.0"%(base,beam,subband))
		execute("single_pulse_search.py -t 10 %s_b%ds%d_ZERO*singlepulse"%(base,beam,subband))
		execute("mkdir %s_b%ds%d_TEST_proc"%(base,beam,subband))
		execute("mv %s_b%ds%d_ZERO* %s_b%ds%d_TEST_proc"%(base,beam,subband,base,beam,subband))
		execute("mv %s_b%ds%d_rfifind.* %s_b%ds%d_TEST_proc"%(base,beam,subband,base,beam,subband))
		#this was used to perform actions on files already moved out of the "test" directory.
		#if not os.path.isfile("%s/%s_b%ds%d_TEST_proc/%s_b%ds%d_ZERO_SNR5_singlepulse.ps"\
														#%(path,base,beam,subband,base,beam,subband)):
			#execute("mv %s/%s_b%ds%d_TEST_proc/%s_b%ds%d_ZERO_singlepulse.ps\
			 #%s/%s_b%ds%d_TEST_proc/%s_b%ds%d_ZERO_SNR5_singlepulse.ps"%(path,base,beam,subband,base,beam,subband,path,base,beam,subband,base,beam,subband))
		#execute("single_pulse_search.py -t 10 %s/%s_b%ds%d_TEST_proc/*.singlepulse"%(path,base,beam,subband))
		#subprocess.Popen("single_pulse_search.py -t 10 *.singlepulse", shell=True, cwd='%s_b%ds%d_TEST_proc'%(base,beam,subband))


##### STEP 2: ONCE SINGLEPULSE FILES ARE CREATED FOR EACH DM #####
#PUT BACK INDENTS (2) (should be inside subband loop)
		os.chdir("%s/%s_b%ds%d_TEST_proc"%(cwd,base,beam,subband)) #hard code path later
		#execute("mkdir obs_data")
		execute("mkdir pulses")
		#execute("mkdir periodic_cands")
		#execute("mkdir TEMP")
		#execute("cd TEMP")
		execute("python %s/pulses_extract.py -db_name %s_b%ds%d_TEST_proc.hdf5 -fits %s\
		 		-store_events -idL %s_b%ds%d_ZERO_DM -store_dir pulses \
					-plot_pulses -plot_statistics -parameters_id FRB130628_Alfa_s%d > /dev/null"\
				%(script_dir,base,beam,subband,infile,base,beam,subband,subband))
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
 ###combined subbands version###
 
"""
#base = argv[1] #base name of fits file, up until beam number "b".
#base = argv[1] #up until channel "s".
base = "4bit-p2030.20160702.FRB130628_0"
beam = 0

#combine subbands into single file
#execute("combine_mocks -o %s_b%d %s.b%ds%d*000.fits %s.b%ds1*000.fits"\ #* stands for 0 or more of the preceding character (not a wildcard)
		#%(base,beam,subband,base,beam,subband,base,beam,subband))
#RFI find
#rfifind -time 2.0 -psrfits -noscales -nooffsets -o p2030.20160702.FRB130628_0_b0_0001 4bit-p2030.20160702.FRB130628_0_b0_0001.fits
# Create a Dedispersion plan. Only needs to be run once; can be applied to all obs.
#with open("%s_b%d_DDplan.txt"%(base, beam), "w+") as output:
#	subprocess.call("DDplan.py -d 600 -f 1375.43 -b 322.73 -n 960 -t 0.00006548\
#	 -s 120 -r 0.5 -o %s_b%d_DDplan %s_b%d_*.fits"%(base,beam,subband,base,beam,subband), shell=True, stdout=output)
#unsure about -r option. ask Jason.

dmstep = 0.50
numdms = 96
downsamp = 3
lodm = 0.
dsubDM = 48.
calls = 13
calls = range(calls)

for call in calls:
	execute("prepsubband -nsub 120 -noscales -nooffsets -downsamp %d -lodm %f\
			-dmstep %f -numdms %d -zerodm -mask *.mask -o %s_b%d_ZERO %s_b%d*.fits"\
			%(downsamp,lodm,dmstep,numdms,base,beam,subband,base,beam,subband))
	lodm +=  dsubDM

execute("mkdir %s_b%d_proc"%(base,beam,subband))
execute("mv %s_b%d_ZERO* %s_b%d_proc"%(base,beam,subband,base,beam,subband))

#one subband version
base = "4bit-p2030.20160702.FRB130628_0"
beam = 0

for call in calls:
	execute("prepsubband -nsub 120 -noscales -nooffsets -downsamp %d -lodm %f\
			-dmstep %f -numdms %d -zerodm -mask 4bit-p2030.20160702.FRB130628_0.b0s%dg0_rfifind.mask -o\
			 %s_b%ds%d_ZERO 4bit-p2030.20160702.FRB130628_0.b0s%dg0.00000.fits"\
			%(downsamp,lodm,dmstep,numdms,base,beam,subband))
	lodm +=  dsubDM

execute("mkdir %s_b%ds%d_TEST_proc"%(base,beam,subband))
execute("mv %s_b%ds%d_ZERO* %s_b%ds%d_TEST_proc"%(base,beam,subband,base,beam,subband))
"""




