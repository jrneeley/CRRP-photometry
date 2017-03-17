#!/usr/bin/env python

import read_dao
import matplotlib.pyplot as mp
import numpy as np
import optical
from astropy.io import ascii
import math

def find_cal_star_coords(target, channel):
    optical_folder = '/Users/Jill/CRRP/OpticalCatalogs/'
    #read in lst file of calibration stars
    psf_stars, psf_x, psf_y = read_dao.read_lst(channel+'_cal.lst')

    cat_ids, x, y, v, ra, dec = optical.read_optical_catalog(optical_folder, target)
    for star in psf_stars:
        index_match = np.argwhere(cat_ids == star)
        cal_ra = ra[index_match]
        cal_dec = dec[index_match]
        print star, cal_ra, cal_dec

def find_cal_stars(target, channel):
    #identify calibration stars

    #read in lst file of calibration stars
    psf_stars, psf_x, psf_y = read_dao.read_lst(channel+'_cal.lst')
    # read in raw file with PSF photometry in all frames
    star_ids, mags = read_dao.read_raw('optical_alf.raw')

    # For each calibration star, find all of the PSF photometry
    bad_index = []
    for star in psf_stars:
        try:
            index_match = star_ids.index(str(star))

        except ValueError:
            bad_index.append(np.argwhere(psf_stars == star))

    psf_stars = np.delete(psf_stars, bad_index)
    print 'Using '+str(len(psf_stars))+' calibration stars.'


    return psf_stars
def find_zp(psf_stars, channel):

    # read in raw file with PSF photometry in all frames
    star_ids, mags = read_dao.read_raw('optical_alf.raw')

    zp = []
    zp_er = []
    avg_mag = []
    for star in psf_stars:
        # Find PSF photometry
        index_match = star_ids.index(str(star))
        mag_match = mags[index_match]
        psf_mag = mag_match[5:len(mag_match)-2:2]
        psf_err = mag_match[6:len(mag_match)-2:2]
        f_num = np.arange(len(psf_mag)/2)

        # Find aperture photometry
        ap_phot = 'lcvs/'+channel+'_'+str(star)+'.txt'
        data = ascii.read(ap_phot, delimiter=' ')
        aor = data['col1']
        frame = data['col2']
        ap_mag = data['col3']
        ap_err = data['col4']

        ap_mag = np.array(ap_mag, dtype='f')
        ap_err = np.array(ap_err, dtype='f')
        I1_psf_mag = np.array(psf_mag[0:len(psf_mag)/2], dtype='f')
        I1_psf_err = np.array(psf_err[0:len(psf_err)/2], dtype='f')
        I2_psf_mag = np.array(psf_mag[len(psf_mag)/2:], dtype='f')
        I2_psf_err = np.array(psf_err[len(psf_err)/2:], dtype='f')
        if channel == 'I1':
            residual = ap_mag - I1_psf_mag
            no_match = np.argwhere(abs(residual) > 10.0)
            residual[no_match] = float('NaN')
        #    res_err = math.sqrt(ap_err**2 + I1_psf_err**2)
        if channel == 'I2':
            residual = ap_mag - I2_psf_mag
            no_match = np.argwhere(abs(residual) > 10.0)
            residual[no_match] = float('NaN')
        #    res_err = math.sqrt(ap_err**2 + I2_psf_err**2)
        zp.append(np.nanmean(residual))
        zp_er.append(np.nanstd(residual))
        avg_mag.append(np.nanmean(ap_mag))



    #    hist, bins = np.histogram(residual)
    #    width = 0.7 * (bins[1] - bins[0])
    #    center = (bins[:-1] + bins[1:]) / 2
    #    mp.bar(center, hist, align='center')
    #    mp.bar(bins, hist)
    #    mp.show()
    keep_zp = []
    keep_zp_er = []
    keep_avg_mag = []

    if channel == 'I1':
        sat_limit = 10.1
    if channel == 'I2':
        sat_limit = 9.9
    for ind in range(len(avg_mag)):
        if avg_mag[ind] > sat_limit:
            keep_zp.append(zp[ind])
            keep_zp_er.append(zp_er[ind])
            keep_avg_mag.append(avg_mag[ind])

    print np.mean(keep_zp), np.std(keep_zp_er)
    mp.plot(avg_mag, zp, 'ro', keep_avg_mag, keep_zp, 'bo')
#    mp.plot(keep_avg_mag, keep_zp, 'bo')
    mp.show()
