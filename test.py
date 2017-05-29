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

target_name = sys.argv[1]
channel = sys.argv[2]

#optical_folder = raw_input("Enter path to optical catalog: ")
optical_folder = '/Users/jrneeley/CRRP/OpticalCatalogs/'


#calibration.find_stars_in_cat(optical_folder, target_name, channel)


calibration.find_zp(channel)

#hdulist = fits.open('mosaic_dn.fits', mode='update')
#prihdr = hdulist[0].header
#scidata = hdulist[0].data
#newdata = scidata[250:762,250:762]

#outfile = 'slice_dn.fits'
#hdu = fits.PrimaryHDU(newdata)
#hdu.writeto(outfile, clobber=True)
#infile = 'alf-ap-phot.dat'
#calibration.find_zp_single_frame(infile)
