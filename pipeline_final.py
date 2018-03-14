#!/usr/bin/env python

import sys
sys.path.insert(0,'/home/jill/python/CRRP-photometry/')
import glob
import daophot_setup
import shutil
import re
import dao
import numpy as np
import os

target = sys.argv[1]
channel = sys.argv[2]
exptime = sys.argv[3]

dao_dir = '/apps/daophot32/'
optical_dir = '/mnt/data/public/jill/CRRP/OpticalCatalogs/'
opt_dir = '/mnt/data/public/jill/CRRP/OPTfiles/'
current_dir = os.getcwd()
data_dir = '/mnt/data/public/jill/CRRP/'+target+'/mosaics/'

deep_mosaic_fits = target+'_'+channel+'_deep.fits'


print 'Working on cluster {} in the {} band....'.format(target, channel)
print 'The deep mosaic is {} \n'.format(deep_mosaic_fits)

# copy appropriate daophot options files to the current directory
os.chdir(data_dir)
daophot_setup.set_opt_files(opt_dir, channel, exptime, warm=1, mosaic=1)

# get list of images
file_list = glob.glob(channel+'*[0-9].fits')
dn_list = glob.glob(channel+'*_dn.fits')

if (len(dn_list) == 0):
    print "Converting images to DN..."
    for image in file_list:
        daophot_setup.spitzer_flux2dn(image, "")
    print "Done.\n"
else:
    print "Files already converted to counts.\n"

dn_list = glob.glob(channel+'*dn.fits')
print dn_list


## Run DAOPHOT
print "\nStarting DAOPHOT..."
for image in dn_list:
    dao.daophot(dao_dir, image)

print "Initial aperture photometry complete. \n"

# locate PSF stars used for the deep mosaic
first_frame = dn_list[0]
psf_stars = np.loadtxt(channel+'-psf.reg', dtype=np.dtype([('x', float), ('y', float)]))
coo_file = re.sub('.fits', '.coo', first_frame)
lst_file = re.sub('.fits', '.lst', first_frame)
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

daophot.find_psf(dao_dir, first_frame)
repeat = raw_input('Do you want to remove any PSF stars? [y/n]: ')
if repeat == 'y':
    ready = raw_input('Delete from .lst file and type continue when ready: ')
    daophot.find_psf(dao_dir, first_frame)
if repeat == 'n':
    print 'PSF ok...continuing...\n'

# copy master PSF to every other frame
frame1_psf = re.sub('.fits', '.psf', first_frame)
master_psf = channel+'-master.psf'
shutil.move(frame1_psf, master_psf)

## Copy the master PSF to each epoch
for file in dn_list:
    psfname = re.sub(".fits",".psf", file)
    shutil.copy(master_psf, psfname)

## Run ALLSTAR on individual BCDs
print "Starting ALLSTAR..."
for img in dn_list:
    dao.allstar(dao_dir, img)
print "ALLSTAR complete."


# copy deep mosaic into current directory
#os.chdir(data_dir)
deep = '../DeepMosaic/'+target+'_'+channel+'_deep_dn'
deep_new = channel+'-deep'
shutil.copy(deep+'.fits', deep_new+'.fits')
shutil.copy(deep+'.psf', deep_new+'.psf')
shutil.copy(deep+'.alf', deep_new+'.alf')
als_list = glob.glob(channel+'*.als')
als_list.insert(0, deep_new+'.alf')

print '\nRunning DAOMATCH/DAOAMASTER between deep mosaic and individual epoch mosaics...'
mch_file = channel+'-mosaics.mch'
dao.daomatch(dao_dir, als_list, mch_file, force_scale_rot=1)
print "DAOMATCH complete."
dao.daomaster(dao_dir, mch_file)
print 'DAOMASTER complete.'

print 'Finished. Its time to run ALLFRAME!!'
