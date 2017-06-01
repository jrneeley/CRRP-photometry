#!/usr/bin/env python

import calibration
import sys
import optical
import daophot_setup
import matplotlib.pyplot as mp
from astropy.io import fits
from matplotlib.colors import LogNorm
import shutil
import glob
import coordinates
import read_dao
import lightcurves
import numpy as np
target_name = sys.argv[1]
#channel = sys.argv[2]

#optical_folder = raw_input("Enter path to optical catalog: ")
optical_folder = '/Users/jrneeley/CRRP/OpticalCatalogs/'

dtype1 = np.dtype([('id', int), ('period', float)])
data = np.loadtxt(target_name+'-clement.txt', dtype=dtype1, usecols=(0,3))

for ind, lcv in enumerate(data['id']):
    lcv_file = optical_folder+target_name+'lcvs/'+target_name+'V'+str(lcv)+'.lcv'
    try:
        U, B, V, R, I = lightcurves.read_optical_lcv(lcv_file)
#    lightcurves.find_period(V[0], V[1], V[2])
#    lightcurves.plot_raw_optical_lcv(U, B, V, R, I)
        lightcurves.plot_phased_optical_lcv(U, B, V, R, I, data['period'][ind], 'V'+str(lcv))
    except:
        print 'Star not found.'
#calibration.find_stars_in_cat(optical_folder, target_name, channel)
#calibration.find_zp(channel)
#hdulist = fits.open('mosaic_dn.fits', mode='update')
#prihdr = hdulist[0].header
#scidata = hdulist[0].data
#newdata = scidata[250:762,250:762]
#outfile = 'slice_dn.fits'
#hdu = fits.PrimaryHDU(newdata)
#hdu.writeto(outfile, clobber=True)
#infile = 'alf-ap-phot.dat'
#calibration.find_zp_single_frame(infile)
