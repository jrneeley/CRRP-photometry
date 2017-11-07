#!/usr/bin/env python

import calibration
import sys
import optical
import daophot_setup
import matplotlib.pyplot as mp
from astropy.io import fits
import shutil
import glob
import coordinates
import read_dao
import lightcurves
import numpy as np
import os
import re

target_name = sys.argv[1]
channel = sys.argv[2]
dao_dir = '/usr/local/phot/'
optical_dir = '/Volumes/Annie/CRRP/OpticalCatalogs/'
opt_dir = '/Volumes/Annie/CRRP/OPTfiles/'
working_dir = os.getcwd()

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
# measure the PSF
#daophot.find_psf(dao_dir, mosaic_dn)
#repeat = raw_input('Do you want to remove any PSF stars? [y/n]: ')
#if repeat == 'y':
#    ready = raw_input('Delete from .lst file and type continue when ready: ')
#    daophot.find_psf(dao_dir, mosaic_dn)
#if repeat == 'n':
#    print 'PSF ok...continuing...\n'
