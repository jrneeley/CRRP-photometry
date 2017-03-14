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


target_name = sys.argv[1]
channel = sys.argv[2]

file_list = glob.glob('all/'+channel+'*[0-9].fits')
dn_list = glob.glob('all/'+channel+'*_dn.fits')

daophot_setup.set_opt_files(channel)


# Convert images to counts
if (len(dn_list) == 0):
	for image in file_list:
		daophot_setup.spitzer_flux2dn(image, "")
	print "All files converted to counts."
else:
	print "Files already converted to counts."

dn_list = glob.glob('all/'+channel+'*dn.fits')

# Define fields (Note: you need to create lists ahead of time for maps)
num_fields, fields = daophot_setup.find_fields(dn_list, channel)
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
		daophot.init_phot(target_name, image)

	print "Initial aperture photometry complete."

if (start <= 1):
## Find PSF on first frame only (manually)
	master_frame = raw_input("Identify master frame: ")
#	daophot.find_psf(master_frame)
## Copy this PSF to each epoch
	shutil.copy(master_frame, 'master_'+channel+'.fits')
	master_file = re.sub(".fits",".psf",master_frame)
	master_list = re.sub(".fits",".lst",master_frame)
	master_psf='master_'+channel+'.psf'
	shutil.copy(master_file, master_psf)
	shutil.copy(master_list, 'master_'+channel+'.lst')

	for file in dn_list:
		psfname = re.sub(".fits",".psf", file)
		shutil.copy(master_psf, psfname)
if (start <= 2):
## Run ALLSTAR
	for file in dn_list:
		allstar.allstar_init(target_name, file)
	print "ALLSTAR complete."

if (start <= 3):
## Run DAOMATCH
## You may want to check the .mch files after this step
	daomatch.daomatch_init(channel, target_name, fields, num_fields)
	print "DAOMATCH complete."

## Run DAOMASTER
if (start <= 4):
	for x in range(0,num_fields):
		print "Working on field "+ str(x+1)
		daomaster.daomaster_init(channel+"_field"+str(x+1)+".mch")

## Find appropriate window in source catalog
if (start <=5):

	ids, xcat, ycat, v_mags, ra, dec = optical.read_optical_catalog(target_name)

	for ind in range(0,num_fields):
		print "Calculating field " + str(ind+1)+ " boundaries..."
		x1, x2, y1, y2 = coordinates.find_coord_window(fields[ind])
		print x1, x2, y1, y2
		c1, c2, c3, c4 = coordinates.radec2pix(target_name, x1, x2, y1, y2, xcat, ycat, ra, dec)
		print "Xmin, Xmax, Ymin, Ymax for optical catalog:"
		print c1, c2, c3, c4

	#print "Calculating off field boundaries..."
	#x1, x2, y1, y2 = coordinates.find_coord_window(off_list)
	#c1, c2, c3, c4 = coordinates.radec2pix(target_name, x1, x2, y1, y2, xcat, ycat, ra, dec)
	#print "Xmin, Xmax, Ymin, Ymax for optical catalog:"
	#print c1, c2, c3, c4
