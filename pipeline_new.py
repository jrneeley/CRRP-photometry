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
import os
import deep_phot
import pexpect

target = sys.argv[1]
channel = sys.argv[2]
exptime = sys.argv[3]

if (channel == '1' or channel == '3p6um' or channel == 'I1'): channel = 'I1'
elif (channel == '2' or channel == '4p5um' or channel == 'I2'): channel = 'I2'
else:
	print 'invalid channel'
	sys.exit()

file_list = glob.glob('data/'+channel+'*[0-9].fits')
dn_list = glob.glob('data/'+channel+'*_dn.fits')

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
daophot_setup.set_opt_files(opt_dir, channel, exptime, warm=1)

# Convert BCDs to counts
if (len(dn_list) == 0):
	print "Converting images to DN..."
	for image in file_list:
		daophot_setup.spitzer_flux2dn(image, "")
	print "Done.\n"
else:
	print "Files already converted to counts.\n"

dn_list = glob.glob('data/'+channel+'*dn.fits')


ask_skip = raw_input("Do you want to skip any steps? [y/n]: ")
if ask_skip == 'n':
	start = 0
if ask_skip == 'y':
	print "Enter 0 to begin with intial DAOPHOT photometry."
	print "Enter 1 to begin with copying the master PSF."
	print "Enter 2 to begin with initial ALLSTAR."
	print "Enter 3 to begin with defining the MIR star list."
	print "Enter 4 to begin with DAOMATCH/DAOMASTER (between Spitzer images)."
	print "Enter 5 to begin with calculating optical catalog boundaries."
	print "Enter 6 to begin with matching the optical and MIR catalogs"
	start = input("Your choice: ")
if (start <= 0):
## Run DAOPHOT
	print "\nStarting DAOPHOT..."
	for image in dn_list:
		daophot.init_phot(dao_dir, target, image)

	print "Initial aperture photometry complete. \n"

if (start <= 1):
## Find PSF on first frame only (manually)
	master_psf = raw_input("Identify master PSF: ")

## Copy the master PSF to each epoch
	for file in dn_list:
		psfname = re.sub(".fits",".psf", file)
		shutil.copy(master_psf, psfname)

if (start <= 2):
## Run ALLSTAR on individual BCDs
	print "\nStarting ALLSTAR..."
	for file in dn_list:
		allstar.allstar_init(dao_dir, target, file)
	print "ALLSTAR complete.\n"

if (start <=3):
##### Define star list on deep mosaic
	deep_phot.mosaic_phot(dao_dir, opt_dir, deep_mosaic_fits, channel, exptime)

if (start <=4):

	# check what directory we are in
	if os.getcwd() == working_dir: os.chdir(working_dir+'/DeepMosaic')
## Run DAOMATCH
## You may want to check the .mch files after this step
	print '\nRunning DAOMATCH/DAOAMASTER between deep mosaic and individual BCDs...'
	daomatch.deep_mosaic(dao_dir, deep_mosaic_fits, target, dn_list)
	print "DAOMATCH complete."

## Run DAOMASTER
	mch_file = re.sub('.fits', '_dn.mch', deep_mosaic_fits)
	daomaster.deep_mosaic(dao_dir, mch_file)
	print 'DAOMASTER complete.\n'
## Find appropriate window in source catalog
if (start <=5):

	curr_dir = os.getcwd()
	if curr_dir == working_dir: os.chdir(working_dir+'/DeepMosaic')
#	print 'Matching with optical catalog...'
	ids, catalog_x, catalog_y, catalog_ra, catalog_dec = optical.read_optical_fnl(optical_dir, target)

	als_file = re.sub('.fits', '_dn.als', deep_mosaic_fits)
	dtype1 = np.dtype([('x', float), ('y', float)])
	data = np.loadtxt(als_file, dtype=dtype1, skiprows=3, usecols=(1,2))
	xmin = np.min(data['x'])
	xmax = np.max(data['x'])
	ymin = np.min(data['y'])
	ymax = np.max(data['y'])

#	f = 'Deep'
	print "Calculating optical boundaries..."

	ra1, ra2, dec1, dec2 = coordinates.find_coord_window_mosaic(deep_mosaic_fits, xmin, xmax, ymin, ymax)
	print ra1, ra2, dec1, dec2
	min_x, min_y = coordinates.radec2catalogpix(ra1, dec1, catalog_x, catalog_y, catalog_ra, catalog_dec)
	max_x, max_y = coordinates.radec2catalogpix(ra2, dec2, catalog_x, catalog_y, catalog_ra, catalog_dec)
#		c1, c2, c3, c4 = coordinates.radec2pix(target, x1, x2, y1, y2, xcat, ycat, ra, dec)
	print "Xmin, Xmax, Ymin, Ymax for optical catalog:"
	print min_x, max_x, min_y, max_y

	xmin = [min_x]
	xmax = [max_x]
	ymin = [min_y]
	ymax = [max_y]
	f = ['Deep']

# Save boundary window for each field into a text file (e.g. I1-catalog-cuts.txt)
	data_save = np.array(zip(f, xmin, xmax, ymin, ymax), dtype=[('c1', 'S8'),
		('c2', float), ('c3', float), ('c4', float), ('c5', float)])
	np.savetxt(channel+'-deep-cuts.txt', data_save, comments='', fmt='%s %0.3f %0.3f %0.3f %0.3f')

if (start <= 6):

# DAOMATCH between optical catalog and MIR catalog
	mir_data = re.sub('.fits', '_dn.als', deep_mosaic_fits)
	optical_data = 'optical:'+target+'-I.mag'
	limits = '{:.2f}, {:.2f}, {:.2f}, {:.2f}'.format(min_x, max_x, min_y, max_y)
	daomatch = pexpect.spawn(dao_dir+'daomatch')
#    daomatch.logfile = sys.stdout
	daomatch.expect("Master input file")
	daomatch.sendline(optical_data+'*')
	daomatch.expect('Ymin, Ymax')
	daomatch.sendline(limits)
	daomatch.expect("Output file")
	daomatch.sendline('op-'+channel+'-deep.mch')
	check = daomatch.expect(["Next input file", "Write this transformation"])
	if check == 0: daomatch.sendline(mir_data+'!') # / forces scale to be 1
	if check == 1: daomatch.sendline('y')
	daomatch.expect("Next input file")
	daomatch.sendline("")
	daomatch.expect("Good bye")
	daomatch.close()

	daomaster = pexpect.spawn(dao_dir+'daomaster')
	daomaster.logfile = sys.stdout

	daomaster.expect("File with list of input files")
	daomaster.sendline('op-'+channel+'-deep.mch')
	daomaster.expect("Minimum number")
	daomaster.sendline("1,1,1")
	daomaster.expect("Maximum sigma")
	daomaster.sendline("10")
	daomaster.expect("Your choice")
	daomaster.sendline("20")
	daomaster.expect("Critical match-up radius")
	daomaster.sendline("-5")
	daomaster.expect("New match-up radius")
	daomaster.sendline("4")
	daomaster.expect("New match-up radius")
	daomaster.sendline("3")
	daomaster.expect("New match-up radius")
	daomaster.sendline("2")
	daomaster.expect("New match-up radius")
	daomaster.sendline("2")
	daomaster.expect("New match-up radius")
	daomaster.sendline("2")
	daomaster.expect("New match-up radius")
	daomaster.sendline("1")
	daomaster.expect("New match-up radius")
	daomaster.sendline("1")
	daomaster.expect("New match-up radius")
	daomaster.sendline("1")
	daomaster.expect("New match-up radius")
	daomaster.sendline("1")
	daomaster.expect("New match-up radius")
	daomaster.sendline("1")
	daomaster.expect("New match-up radius")
	daomaster.sendline("0")
	daomaster.expect("Assign new star IDs")
	daomaster.sendline("n")
	daomaster.expect("A file with mean magnitudes")
	daomaster.sendline("n")
	daomaster.expect("A file with corrected magnitudes")
	daomaster.sendline("n")
	daomaster.expect("A file with raw magnitudes")
	daomaster.sendline("n")
	daomaster.expect("A file with the new transformations")
	daomaster.sendline("y")
	daomaster.expect("Output file name")
	daomaster.sendline("")
	daomaster.expect("New output file name")
	daomaster.sendline("")
	daomaster.expect("A file with the transfer table")
	daomaster.sendline("e")
	daomaster.close(force=True)

	combine = pexpect.spawn(dao_dir+'combine')
	combine.expect('master transformations')
	combine.sendline('op-'+channel+'-deep.mch')
	combine.expect('Expand')
	combine.sendline('n')
	combine.expect('EXIT')
	combine.sendline('')
	combine.close()

print 'Finished. Its time to run ALLFRAME!!'
