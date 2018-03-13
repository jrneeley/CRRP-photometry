#!/usr/bin/env python

import sys
import glob
import daophot_setup
import shutil
import daophot
import re
import allstar
import daomatch_crrp
import daomaster_crrp
#import coordinates
#import optical
import numpy as np
import os
#import deep_phot
#import pexpect

target = sys.argv[1]
channel = sys.argv[2]
exptime = sys.argv[3]

if (channel == '1' or channel == '3p6um' or channel == 'I1'): channel = 'I1'
elif (channel == '2' or channel == '4p5um' or channel == 'I2'): channel = 'I2'
else:
	print 'invalid channel'
	sys.exit()

file_list = glob.glob('mosaics/'+channel+'*[0-9].fits')
dn_list = glob.glob('mosaics/'+channel+'*_dn.fits')

# Identify relevant paths
dao_dir = '/usr/local/phot/'
optical_dir = '/Volumes/Annie/CRRP/OpticalCatalogs/'
opt_dir = '/Volumes/Annie/CRRP/OPTfiles/'
working_dir = os.getcwd()

deep_mosaic_fits = target+'_'+channel+'_deep.fits'

print '\nWorking on cluster {} in the {} band....'.format(target, channel)
print 'The deep mosaic is {} \n'.format(deep_mosaic_fits)

#dao_dir, optical_dir, opt_dir = daophot_setup.folder_setup()

# Copy appropriate opt files to current directory
daophot_setup.set_opt_files_mosaic(opt_dir, channel, exptime, warm=1)

# Convert BCDs to counts
if (len(dn_list) == 0):
	print "Converting images to DN..."
	for image in file_list:
		daophot_setup.spitzer_flux2dn(image, "")
	print "Done.\n"
else:
	print "Files already converted to counts.\n"

dn_list = glob.glob('mosaics/'+channel+'*dn.fits')


ask_skip = raw_input("Do you want to skip any steps? [y/n]: ")
if ask_skip == 'n':
	start = 0
if ask_skip == 'y':
	print "Enter 0 to begin with intial DAOPHOT photometry."
	print "Enter 1 to begin with copying the master PSF."
	print "Enter 2 to begin with initial ALLSTAR."
	print "Enter 3 to begin with defining the MIR star list."
	start = input("Your choice: ")
if (start <= 0):
## Run DAOPHOT
	print "\nStarting DAOPHOT..."
	for image in dn_list:
		daophot.init_phot(dao_dir, target, image, mosaics=1)

	print "Initial aperture photometry complete. \n"

if (start <= 1):

	mosaic_dns = glob.glob('mosaics/'+channel+'*dn.fits')
	mosaic_dn = mosaic_dns[0]

	psf_stars = np.loadtxt('mosaics/'+channel+'-psf.reg', dtype=np.dtype([('x', float), ('y', float)]))
	coo_file = re.sub('.fits', '.coo', mosaic_dn)
	lst_file = re.sub('.fits', '.lst', mosaic_dn)
	f = open(coo_file, 'r')
	f2 = open(lst_file, 'w')
	for ind in range(3):
		line = f.readline()
		f2.write(line)
	f.close()
	f2.close()
	dtype = np.dtype([('id', int), ('x', float), ('y', float), ('c1', float),
    	('c2', float), ('c3', float), ('c4', float)])
	all_stars = np.loadtxt(coo_file, dtype=dtype, skiprows=3)
	f_handle = open(lst_file, 'a')
	for x, y in zip(psf_stars['x'], psf_stars['y']):
		dist = np.sqrt((all_stars['x'] - x)**2 + (all_stars['y'] - y)**2)
		line = all_stars[dist == np.min(dist)]
		np.savetxt(f_handle, line, fmt='%8i %8.2f %8.2f %8.3f %8.3f %8.3f %8.3f')
	f_handle.close()

	master_psf = raw_input("Identify master PSF: ")

## Copy the master PSF to each epoch
	for file in dn_list:
		psfname = re.sub(".fits",".psf", file)
		shutil.copy(master_psf, psfname)

if (start <= 2):
## Run ALLSTAR on individual BCDs
	print "\nStarting ALLSTAR..."
	for file in dn_list:
		allstar.allstar_init(dao_dir, target, file, mosaics=1)
	print "ALLSTAR complete.\n"


if (start <=3):
	# check what directory we are in
	if os.getcwd() == working_dir: os.chdir(working_dir+'/DeepMosaic')
## Run DAOMATCH
## You may want to check the .mch files after this step
	print '\nRunning DAOMATCH/DAOAMASTER between deep mosaic and individual epoch mosaics...'
	daomatch_crrp.deep_mosaic(dao_dir, deep_mosaic_fits, target, dn_list, channel=channel, mosaics=1)
	print "DAOMATCH complete."

## Run DAOMASTER
	mch_file = channel+'_mosaics.mch'
	daomaster_crrp.deep_mosaic(dao_dir, mch_file, mosaics=1)
	print 'DAOMASTER complete.\n'


print 'Finished. Its time to run ALLFRAME!!'
