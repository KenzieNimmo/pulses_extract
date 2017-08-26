import pandas as pd
from matplotlib import pyplot as plt
import numpy as np
from matplotlib.collections import EllipseCollection
from matplotlib.collections import CircleCollection
import glob
import os

def get_markers(array):
	marker_size = array.Sigma
	marker_size[(marker_size>=10) & (marker_size<15)] = 10
	marker_size[(marker_size>=15) & (marker_size<40)] = 40
	marker_size[marker_size>=40] = 200
	return marker_size

def make_plot(subband):
	databases = glob.glob("*s%d*proc.hdf5"%subband)
	fig, ax = plt.subplots(nrows=3,ncols=3, figsize=(16, 9))#, facecolor='w', edgecolor='k')

	ax = ax.ravel()
	for i,file in enumerate(databases):
		events = pd.read_hdf(file, 'events')
		#Events = events[(events.Sigma > 10)]# & (events.Sigma < 26)]
		events = events[events.Sigma>=6]
		marker_size = 5
		ax[i].axhline(470,c='c',lw=1) #DM of FRB
		ax[i].scatter(events.Time, events.DM, s=marker_size, marker='o', facecolors='none', edgecolors ='k', alpha =0.2, label='Events')
		try:
			pulses = pd.read_hdf(file, 'pulses')
			pulses = pulses[pulses.Sigma>=6]
			#all grouped pulses
			marker_size = get_markers(pulses)
			ax[i].scatter(pulses.Time, pulses.DM, s=marker_size, marker='o', facecolors='none', edgecolors ='r', label='All pulses')
			#RFI
			#RFI = pulses[pulses.Pulse > 9]
			#if len(RFI)>0:
				#marker_size = get_markers(RFI)
				#ax[i].scatter(RFI.Time, RFI.DM, s=marker_size, marker='o', facecolors='none', edgecolors ='y', label='RFI pulses')
			#Remaining pulse candidates
			pulse_cands = pulses[pulses.Pulse == -1]
			if len(pulse_cands)>0:
				marker_size = get_markers(pulse_cands)
				ax[i].scatter(pulse_cands.Time, pulse_cands.DM, s=marker_size, marker='o', facecolors='none', edgecolors ='g', label='pulse candidates')
		except KeyError: 
			pass
		ax[i].set_xlim([-1, 600])
		ax[i].set_ylim([-1, 575])
		ax[i].set_xlabel('Time (s)')
		ax[i].set_ylabel(r'DM (pc cm$^{-3}$)')
		#ax[i].legend(fontsize=8)
		try:
			num_cands = len(pulse_cands)
		except UnboundLocalError:
			num_cands = 0

		ax[i].set_title('BEAM %d: %d pulse candidates'%(i,num_cands), fontsize=12)
	base_path = os.getcwd()
	base = os.path.basename(base_path)
	fig.delaxes(ax[7])
	fig.delaxes(ax[8])
	#ax[8].axis('off')
	red_artist = plt.Line2D((0,1),(0,0),marker='o', ls='None', mec='red', mfc='none')
	green_artist = plt.Line2D((0,1),(0,0), marker='o', ls='None', mec='green', mfc='none')
	#yellow_artist = plt.Line2D((0,1),(0,0), marker='o', ls='None', mec='#D4AC0D', mfc='none')
	black_artist = plt.Line2D((0,1),(0,0), marker='o', ls='None', mec='black', mfc='none')
	#small_marker = plt.Line2D((0,1),(0,0),marker='o', ls='None', markersize = 10, mec='blue', mfc='none')
	#medium_marker = plt.Line2D((0,1),(0,0), marker='o', markersize = 40, ls='None', mec='blue', mfc='none')
	#large_marker = plt.Line2D((0,1),(0,0), marker='o', markersize = 200, ls='None', mec='blue', mfc='none')
	DM_line = plt.Line2D((0,1),(0,0), lw=2,color='c')
	#ax[6].scatter(300,400,s=10, marker='o', facecolors='none', edgecolors='none', label="10 <= S/N < 15")
	#ax[6].scatter(310,400,s=40, marker='o', facecolors='none', edgecolors='none')
	#ax[6].scatter(320,400,s=200, marker='o', facecolors='none', edgecolors='none')
	ax[6].legend([green_artist, red_artist, black_artist, DM_line],['Pulse candidates',\
			'Grouped Events', 'Events', r'FRB DM = 470 pc cm$^{-3}$'], bbox_to_anchor=(1.05, 1), loc=2, ncol=2,  borderaxespad=0., frameon=False) # '10 <= S/N < 15','15 <= S/N < 40', 'S/N >= 40'],\
	fig.tight_layout(w_pad = 2, h_pad = 2)
	plt.subplots_adjust(hspace = 0.5, wspace=.3)
	plt.savefig('overview_%s_subband_%d.png'%(base,subband),bbox_inches='tight', dpi=200)
	fig.clf()
	plt.clf()
	return

make_plot(0)
make_plot(1)

#offsets = list(zip(Events.Time, Events.DM))

#fig, ax = plt.subplots(1,2,figsize=(16,6), sharex=True,sharey=True)
#ax[0].scatter(Events.Time, Events.DM, s=marker_size, marker='o', facecolors='none', edgecolors ='k')
#ax[1].add_collection(EllipseCollection(widths=marker_size, heights=marker_size,angles=0,units='xy', facecolors='none',edgecolors='k',offsets=offsets, transOffset=ax[1].transData))
#ax[1].axis('equal')

#plt.scatter(300,300, s=25, marker='o', facecolors='none', edgecolors ='r')
#plt.scatter(400,300, s=10, marker='o', facecolors='none', edgecolors ='r')
#ax[1].add_collection(EllipseCollection(widths=25,heights=25,angles=0,units='xy', facecolors='none',edgecolors='r',offsets=(300,300), transOffset=ax[1].transData))
#ax[1].add_collection(EllipseCollection(widths=10,heights=10,angles=0,units='xy', facecolors='none',edgecolors='r',offsets=(400,300), transOffset=ax[1].transData))
