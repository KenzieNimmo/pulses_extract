import os
import sys
import Image
from PIL import ImageTk
import pandas as pd
import numpy as np
import pyfits
import h5py
import Tkinter
import argparse

#DM_min = 530
#DM_max = 570
#Sigma_min = 10
#duration_max = 8 #duration is in seconds in hd file
#Need to add these to command line options for select cands
def select_cands(filename, DM_min=None, DM_max=None, Sigma_min=None,duration_max=None, sort=None):
	cands = pd.read_hdf(filename, 'pulses')
	if DM_min:
		cands = cands[cands.DM > DM_min]
	if DM_max:
		cands = cands[cands.DM < DM_max]
	if Sigma_min:
		cands = cands[cands.Sigma > Sigma_min]
	if duration_max:
		cands = cands[cands.Duration < duration_max]

	new_filename = os.path.splitext(filename)[0] + '_pulses_unranked.txt'
	if sort:
		cands.sort(sort, ascending=False, inplace=True)
		return cands.to_csv(new_filename, sep='\t', cols=['Pulse',], header=['Rank',], index_label='#PulseID')
	else:
		#default sorting (by pulse ID)
		return cands.sort_index().to_csv(new_filename, sep='\t', cols=['Pulse',], header=['Rank',], index_label='#PulseID')

def next_image(event):
	print("Hit the enter key to view the next figure.")
	event.widget.quit() # this will cause mainloop to unblock.

def next_image_quiet(event):
	event.widget.quit()

def viewer(txtfile, path_to_pulses, obs, view_only_mode=False,ranks_to_view=None):
	#put txtfile columns into numpy arrays:
	with open(txtfile, 'r') as tf:
		pulse_IDs = np.array([]).reshape(0,1)
		ranks = np.array([]).reshape(0,1)
		next(tf)
		for row in tf:
			cols = np.array(row.split())
			pulse_IDs = np.vstack([pulse_IDs, int(cols[0])])
			ranks = np.vstack([ranks, cols[1]])

	cands = np.loadtxt(txtfile, usecols=(0,), dtype='int')
	if view_only_mode:
		root = Tkinter.Tk()
		root.bind('<Return>', next_image_quiet)
		for i, cand in enumerate(cands):
			cand_rank = ranks[i]
			cand_rank = cand_rank.astype(int)
			if cand_rank in ranks_to_view:
				path = path_to_pulses + str(cand)
				for file in os.listdir(path):
					if file.endswith('.png'):
						img = Image.open(os.path.join(path, file))
						img = img.resize((800, 500), Image.ANTIALIAS)
						root.geometry('%dx%d' % (img.size[0],img.size[1]))
						tkpi = ImageTk.PhotoImage(img)
						label_image = Tkinter.Label(root, image=tkpi)
						label_image.place(x=0,y=0,width=img.size[0],height=img.size[1])
						root.title(file)				
						root.mainloop()
					else:
						pass
			else:
				pass
	else:
		for i, cand in enumerate(cands):
			path = path_to_pulses + str(cand)
			root = Tkinter.Tk()
			#root.geometry("800x500")
			root.bind('<Return>', next_image)
			for file in os.listdir(path):
				if file.startswith('FRB'): #avoid looking at diagnostic plot
					img = Image.open(os.path.join(path, file))
					img = img.resize((800, 500), Image.ANTIALIAS)
					root.geometry('%dx%d' % (img.size[0],img.size[1]))
					tkpi = ImageTk.PhotoImage(img)
					label_image = Tkinter.Label(root, image=tkpi)
					label_image.place(x=0,y=0,width=img.size[0],height=img.size[1])
					root.title(file)				
					root.mainloop()
				else:
					pass
			print "Now ranking candidate %d from %s"%(cand,obs)
			root.destroy()
			input = raw_input("Provide a ranking for this candidate:")
			while input.isdigit() == False:
				input = raw_input("Oups. Please enter an integer:")
			ranks[i] = input

		pulse_IDs=pulse_IDs.flatten()
		pulse_IDs = pulse_IDs.astype(int)
		ranks = ranks.flatten()
		ranks = ranks.astype(int)
		path, txtfilename = os.path.split(txtfile)
		ranked_file = path + "/" + obs + '_pulses.txt'
		np.savetxt(ranked_file, np.column_stack([pulse_IDs, ranks]), delimiter='\t', header='PulseID\tRank',fmt="%s")
	
if __name__ == '__main__':
	#filename = 'puppi_57607_C0531+33_0603/pulses/puppi_57607_C0531+33_0603.hdf5'
	#select_cands(filename, DM_min=530, DM_max=570, Sigma_min=10, sort='Sigma')
	def parser():
		parser = argparse.ArgumentParser(description='View or rank an observation. Hit enter to view the next image.')
		parser.add_argument('-view_only_mode', help="View candidates without ranking them", action='store_true')
		parser.add_argument('-rank', help="View only candidates with these/this rank(s). E.g. If you want to view \
		 					candidates of ranks 1 and 2 type: -rank 12. Default: Rank 0 only", default=[0], type=list)
		parser.add_argument('observation', help="Name of observation to view/rank.")
		parser.add_argument('-create_cand_list', help="Create a list of candidates to rank", action='store_true')
		parser.add_argument('-not_dop263', help="If you want to run script outside of dop263 (script must be in same dir as observation dir).", action='store_true')

		return parser.parse_args()
	args = parser()

	obs = args.observation
	if args.not_dop263:
		obs_path = ''
	else:
		obs_path = '/psr_archive/hessels/hessels/AO-FRB/pipeline_products/'
	path_to_pulses = obs_path + obs + '/pulses/'
	filename = path_to_pulses + obs + '.hdf5' #for select_cands

	if args.create_cand_list:
		select_cands(filename)

	if args.view_only_mode: 
		txtfile = path_to_pulses + obs + '_pulses.txt'
		ranks_to_view = args.rank
		ranks_to_view = map(int, ranks_to_view)
		viewer(txtfile, path_to_pulses, obs, view_only_mode=True, ranks_to_view=ranks_to_view)

	else:
		txtfile = path_to_pulses + obs + '_pulses_unranked.txt'
		viewer(txtfile, path_to_pulses, obs)



	#txtfile = path_to_pulses + obs + '_pulses_unranked.txt'
	#obs = sys.argv[1]

	#path, txtfilename = os.path.split(txtfile)
	#select_cands(filename)
	#viewer(txtfile, path_to_pulses, obs)

	#viewer('puppi_57638_C0531+33_1218/pulses/puppi_57638_C0531+33_1218_pulses.txt','puppi_57638_C0531+33_1218/pulses/')
	#ranked_file = os.path.splitext(txtfile)[0] + '_ranked.txt'



