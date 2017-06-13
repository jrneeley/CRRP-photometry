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

optical_folder = '/Users/jrneeley/CRRP/OpticalCatalogs/'

datasets = optical.compile_datasets(optical_folder, target_name)
print '\n\nDatasets: '
print datasets

#legend =
## Find periods of optical light curves
dtype1 = np.dtype([('id', int), ('period', float)])
data = np.loadtxt(target_name+'-clement.txt', dtype=dtype1, usecols=(0,3))


for ind, lcv in enumerate(data['id']):
    lcv_file = optical_folder+target_name+'lcvs/'+target_name+'V'+str(lcv)+'.lcv'
    try:
        U, B, V, R, I = lightcurves.read_optical_lcv(lcv_file, old=1)

        period = lightcurves.find_period(V[0], V[1], V[2], data['period'][ind])
        print lcv, period, 1/period, 1/data['period'][ind]

#    lightcurves.plot_raw_optical_lcv(U, B, V, R, I)
        lightcurves.plot_phased_optical_lcv(U, B, V, R, I, period, 'V'+str(lcv), datasets)
    except:
        print lcv_file
        e = sys.exc_info()[0]
        print str(e)
    mp.show()
    sys.exit()
#######

## Find periods of IRAC lightcurves
lcvs = glob.glob('lcvs/matches/*.lcv')
for lcv in lcvs:
    dtype1 = np.dtype([('filter', 'S2'), ('aor', int), ('mjd', float),
        ('mag', float), ('err', float)])
    data = np.loadtxt(lcv, dtype=dtype1, usecols=(0,1,3,6,7))

    filters = np.unique(data['filter'])
    for filt in filters:
        mag_all = data['mag'][data['filter'] == filt]
        err_all = data['err'][data['filter'] == filt]
        mjd_all = data['mjd'][data['filter'] == filt]
        aor_all = data['aor'][data['filter'] == filt]
        mag = mag_all[~np.isnan(mag_all)]
        err = err_all[~np.isnan(mag_all)]
        mjd = mjd_all[~np.isnan(mag_all)]
        period = lightcurves.find_period(mag, err, mjd)
        print lcv, period
###########




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
