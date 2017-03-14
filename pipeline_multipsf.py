#!/usr/bin/env python

import sys
import glob
import daophot_setup
import shutil
import daophot
import re
import allstar
import daomatch


target_name = sys.argv[1]
channel = sys.argv[2]
start = 3

file_list = glob.glob('all/'+channel+'*[0-9].fits')
dn_list = glob.glob('all/'+channel+'*_dn.fits')

daophot_setup.set_opt_files(channel)

if (start <= 0):
#Convert images to counts
	if (len(dn_list) == 0):
		for image in file_list:
			daophot_setup.spitzer_flux2dn(image, "")
		print "All files converted to counts."
	else:
		print "Files already converted to counts."

	dn_list = glob.glob('all/'+channel+'*dn.fits')

if (start <= 1):
## Run DAOPHOT
	print "Starting DAOPHOT"
	for image in dn_list:
		daophot.init_phot(image)
		daophot.find_psf(image)
	print "Initial aperture photometry and PSF complete."

if (start <= 2):
## Run ALLSTAR
	for file in dn_list:
		allstar.allstar_init(file)
	print "ALLSTAR complete."

if (start <= 3):
## Run DAOMATCH
	daomatch.daomatch_init(channel)
	print "DAOMATCH complete."
