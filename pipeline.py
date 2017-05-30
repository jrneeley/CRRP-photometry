#!/usr/bin/env python

import sys
import glob
import daophot_setup
import shutil
import daophot
import re
import allstar
import daomatch
import daomaster
import coordinates
import optical
import numpy as np

target_name = sys.argv[1]
channel = sys.argv[2]

file_list = glob.glob('all/'+channel+'*[0-9].fits')
dn_list = glob.glob('all/'+channel+'*_dn.fits')

# Identify relevant paths
dao_folder = '/Users/jrneeley/Daophot/'
optical_folder = '/Users/jrneeley/CRRP/OpticalCatalogs/'

# Set daophot.opt file to appropriate channel
daophot_setup.set_opt_files(channel)


# Convert images to counts
if (len(dn_list) == 0):
	print "Converting images to DN..."
	for image in file_list:
		daophot_setup.spitzer_flux2dn(image, "")
	print "Done."
else:
	print "Files already converted to counts."

dn_list = glob.glob('all/'+channel+'*dn.fits')

# Define fields (Note: you need to create lists ahead of time for maps)
num_fields, fields = daophot_setup.find_fields(file_list, channel)
print str(num_fields)+' fields defined with '+str(len(fields[0]))+' images each.'


ask_skip = raw_input("Do you want to skip any steps? [y/n]: ")
if ask_skip == 'n':
	start = 0
if ask_skip == 'y':
	print "Enter 0 to begin with intial DAOPHOT photometry."
	print "Enter 1 to begin with master PSF."
	print "Enter 2 to begin with ALLSTAR."
	print "Enter 3 to begin with DAOMATCH (between Spitzer images)."
	print "Enter 4 to begin with DAOMASTER (between Spitzer images)."
	print "Enter 5 to begin with calculating optical field boundaries."
	start = input("Your choice: ")
if (start <= 0):
## Run DAOPHOT
	print "Starting DAOPHOT"
	for image in dn_list:
		daophot.init_phot(dao_folder, target_name, image)

	print "Initial aperture photometry complete."

if (start <= 1):
## Find PSF on first frame only (manually)
	master_frame = raw_input("Identify master frame: ")
#	daophot.find_psf(master_frame)
## Copy this PSF to each epoch

	master_file = re.sub(".fits",".psf",master_frame)
	shutil.copy(master_file, 'master.psf')

	shutil.copy(master_frame, re.sub("all","master",master_frame))
	shutil.copy(master_file, re.sub("all","master",master_file))

	for file in dn_list:
		psfname = re.sub(".fits",".psf", file)
		shutil.copy('master.psf', psfname)
if (start <= 2):
## Run ALLSTAR
	print "Starting ALLSTAR..."
	for file in dn_list:
		allstar.allstar_init(dao_folder, target_name, file)
	print "ALLSTAR complete."

if (start <= 3):
## Run DAOMATCH
## You may want to check the .mch files after this step
	daomatch.daomatch_init(dao_folder, channel, target_name, fields, num_fields)
	print "DAOMATCH complete."

## Run DAOMASTER
if (start <= 4):
	for x in range(0,num_fields):
#		print "Working on field "+ str(x+1)
		daomaster.daomaster_init(dao_folder, channel+"_field"+str(x+1)+".mch")

## Find appropriate window in source catalog
if (start <=5):

	ids, xcat, ycat, ra, dec = optical.read_optical_catalog(optical_folder, target_name)

	xmin = np.zeros(num_fields)
	xmax = np.zeros(num_fields)
	ymin = np.zeros(num_fields)
	ymax = np.zeros(num_fields)
	f = np.zeros(num_fields, dtype='S8')

	for ind in range(0,num_fields):
		print "Calculating field " + str(ind+1)+ " boundaries..."
		x1, x2, y1, y2 = coordinates.find_coord_window(fields[ind])
		print x1, x2, y1, y2
		c1, c2, c3, c4 = coordinates.radec2pix(target_name, x1, x2, y1, y2, xcat, ycat, ra, dec)
		print "Xmin, Xmax, Ymin, Ymax for optical catalog:"
		print c1, c2, c3, c4
		xmin[ind] = c1
		xmax[ind] = c2
		ymin[ind] = c3
		ymax[ind] = c4
		f[ind] = 'Field'+str(ind+1)

# Save boundary window for each field into a text file (e.g. I1-catalog-cuts.txt)
	data_save = np.array(zip(f, xmin, xmax, ymin, ymax), dtype=[('c1', 'S8'),
		('c2', float), ('c3', float), ('c4', float), ('c5', float)])
	np.savetxt(channel+'-catalog-cuts.txt', data_save, comments='', fmt='%s %0.3f %0.3f %0.3f %0.3f')
