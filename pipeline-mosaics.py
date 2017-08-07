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
from astropy.io import fits

target = sys.argv[1]
channel = sys.argv[2]
exptime = sys.argv[3]

file_list = glob.glob('mosaics/'+channel+'*.fits')
dn_list = glob.glob('mosaics/'+channel+'*_dn.fits')

# Identify relevant paths
dao_dir = '/usr/local/phot/'
optical_dir = '/Volumes/Annie/CRRP/OpticalCatalogs/'
opt_dir = '/Volumes/Annie/CRRP/OPTfiles/'

# Set daophot.opt file to appropriate channel
daophot_setup.set_opt_files_mosaic(opt_dir, channel, exptime)

# Convert images to counts
if (len(dn_list) == 0):
	print "Converting images to DN..."
	for image in file_list:
		daophot_setup.spitzer_flux2dn(image, "")
	print "Done."
else:
	print "Files already converted to counts."

dn_list = glob.glob('mosaics/'+channel+'*_dn.fits')

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
		daophot.init_phot(dao_dir, target, image)

	print "Initial aperture photometry complete."

if (start <= 1):

	hdulist = fits.open(dn_list[0])
	data = hdulist[0].data
	hdulist.close()
	nx = len(data[0])
	ny = len(data)
	newdata = data
	if channel == 'I1':
		newdata[:,int(nx/3):] = np.nan
	if channel == 'I2':
		newdata[:,0:int(nx*2/3)] = np.nan
	outfile = channel+'-half-mosaic.fits'
	hdu = fits.PrimaryHDU(newdata)
	hdu.writeto(outfile, clobber=True)


## Find PSF on first frame only (manually)
#	master_frame = raw_input("Identify master frame: ")
	ready = raw_input("Generate a PSF and copy to the first file. Type continue when ready: ")
#	daophot.find_psf(master_frame)

## Copy this PSF to each epoch
	master_file = re.sub(".fits",".psf", dn_list[0])

	for dn in dn_list:
		if dn == dn_list[0]:
			continue
		else:
			psfname = re.sub(".fits",".psf", dn)
			shutil.copy('master.psf', psfname)
if (start <= 2):
## Run ALLSTAR
	print "Starting ALLSTAR..."
	for dn in dn_list:
		allstar.allstar_mosaic(dao_dir, target, dn)
	print "ALLSTAR complete."

if (start <= 3):
## Run DAOMATCH
## You may want to check the .mch files after this step
	daomatch.daomatch_mosaic(dao_dir, channel, target, dn_list)
	print "DAOMATCH complete."


if (start <= 4):
## Run DAOMASTER
	daomaster.daomaster_mosaic(dao_dir, channel+"_mosaic.mch")

## Find appropriate window in source catalog
if (start <=5):

	dtype1 = np.dtype([('x', float), ('y', float)])
	data = np.loadtxt(channel+'_mosaic.mag', dtype=dtype1, skiprows=3, usecols=(1,2))
	xmin = np.min(data['x'])
	xmax = np.max(data['x'])
	ymin = np.min(data['y'])
	ymax = np.max(data['y'])

	ids, xcat, ycat, ra, dec = optical.read_optical_fnl(optical_dir, target)

	print "Calculating field boundaries..."
	x1, x2, y1, y2 = coordinates.find_coord_window_mosaic(dn_list[0], xmin, xmax, ymin, ymax)
	print x1, x2, y1, y2
	c1, c2, c3, c4 = coordinates.radec2pix(target, x1, x2, y1, y2, xcat, ycat, ra, dec)
	print "Xmin, Xmax, Ymin, Ymax for optical catalog:"
	print c1, c2, c3, c4
	xmin = [c1]
	xmax = [c2]
	ymin = [c3]
	ymax = [c4]
	f = ['Field1']

# Save boundary window for each field into a text file (e.g. I1-catalog-cuts.txt)
	data_save = np.array(zip(f, xmin, xmax, ymin, ymax), dtype=[('c1', 'S8'),
		('c2', float), ('c3', float), ('c4', float), ('c5', float)])
	np.savetxt(channel+'-mosaic-catalog-cuts.txt', data_save, comments='', fmt='%s %0.3f %0.3f %0.3f %0.3f')
