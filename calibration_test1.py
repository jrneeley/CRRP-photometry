#!/usr/bin/env python

import read_dao
import matplotlib.pyplot as mp
import numpy as np
import optical
from astropy.io import ascii
import math

def find_cal_star_coords(optical_folder, target, channel):
#    optical_folder = '/Users/Jill/CRRP/OpticalCatalogs/'
    #read in lst file of calibration stars
    psf_stars, psf_x, psf_y = read_dao.read_lst(channel+'_cal.lst')
#    psf_stars, psf_mags = read_dao.read_alf(channel+'.alf')

    args = np.argsort(psf_mags)
    psf_mags = psf_mags[args]
    psf_stars = psf_stars[args]
    cal_ra = []
    cal_dec = []
    stars = []
    cat_ids, x, y, v, ra, dec = optical.read_optical_catalog(optical_folder, target)
    for star in psf_stars:
        index_match = np.argwhere(cat_ids == star)
        if len(index_match):
            stars.append(star)
            cal_ra.append(ra[index_match])
            cal_dec.append(dec[index_match])

    ascii.write([stars, cal_ra, cal_dec], channel+'_cal_stars.txt')

def find_cal_stars(target, channel):
    #identify calibration stars

    #read in lst file of calibration stars
    psf_stars, psf_x, psf_y = read_dao.read_lst('master_'+channel+'.lst')
#    psf_stars, psf_m = read_dao.read_alf(channel+'.alf')
    # read in raw file with PSF photometry in all frames
    star_ids, mags = read_dao.read_raw('optical-'+channel+'-alf.raw')

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
    #star_ids, mags = read_dao.read_raw('optical_alf.raw')
    star_ids, mags = read_dao.read_raw('optical-'+channel+'-alf.raw')

    zp = []
    zp_er = []
    avg_mag = []
    avg_mag_er = []
    skipped = 0
    for star in psf_stars:
        # Find PSF photometry
        index_match = star_ids.index(str(star))
        mag_match = mags[index_match]
        psf_mag = np.array(mag_match[5:len(mag_match)-2:2], dtype='f')
        psf_err = np.array(mag_match[6:len(mag_match)-2:2], dtype='f')
        f_num = np.arange(len(psf_mag)/2)

        # Find aperture photometry
        ap_phot = 'lcvs/psf_stars/'+channel+'_'+str(star)+'.txt'
        data = ascii.read(ap_phot, delimiter=' ')
        aor = data['col1']
        frame = data['col2']
        try:
            ap_mag = np.array(data['col3'], dtype='f')
            ap_err = np.array(data['col4'], dtype='f')
        except:
            skipped += 1
            continue
        # Remove some bad aperture photometry mags
        median_mag = np.nanmedian(ap_mag)
        diff = abs(ap_mag - median_mag)
        bad_phot = np.argwhere(diff > 3*np.nanstd(ap_mag))
        ap_mag[bad_phot] = float('NaN')

        if channel == 'I1':
        #    psf_mag = np.array(psf_mag[0:len(psf_mag)/2], dtype='f')
        #    psf_err = np.array(psf_err[0:len(psf_err)/2], dtype='f')
            psf_mag[psf_mag > 90] = float('NaN')
        #    res_err = math.sqrt(ap_err**2 + I1_psf_err**2)
        if channel == 'I2':
        #    psf_mag = np.array(psf_mag[len(psf_mag)/2:], dtype='f')
        #    psf_err = np.array(psf_err[len(psf_err)/2:], dtype='f')
            psf_mag[psf_mag > 90] = float('NaN')

        residual = ap_mag - psf_mag


    #    if np.nanstd(ap_mag) < 0.1:
    #        mp.plot(f_num, residual, 'ro')
    #        mp.show()
        zp.append(np.nanmean(residual))
        zp_er.append(np.nanstd(residual))
        avg_mag.append(np.nanmean(ap_mag))
        avg_mag_er.append(np.nanstd(ap_mag))

    #    good_values = residual[~np.isnan(residual)]
    #    bin_num = np.sqrt(len(good_values))
    #    print bin_num
    #    n, bins, patches = mp.hist(good_values, bins=bin_num)
    #    mp.show()


    print 'Skipped '+str(skipped)+' stars.'



    #mp.plot(avg_mag, avg_mag_er, 'ro')
    #mp.show()

    keep_zp = []
    keep_zp_er = []
    keep_avg_mag = []

    if channel == 'I1':
        sat_limit = 11.1
        err_limit = 0.1
    if channel == 'I2':
        sat_limit = 10.4
        err_limit = 0.1
    for ind in range(len(avg_mag)):
        #if avg_mag[ind] > sat_limit:
        if avg_mag_er[ind] < err_limit:
            keep_zp.append(zp[ind])
            keep_zp_er.append(zp_er[ind])
            keep_avg_mag.append(avg_mag[ind])
#    keep_zp = zp[avg_mag_er < err_limit]
#    keep_zp_er = zp[avg_mag_er < err_limit]
#    keep_avg_mag = zp[avg_mag_er < err_limit]

    linex = [sat_limit,20]
    liney = [np.nanmedian(keep_zp), np.nanmedian(keep_zp)]
    print np.nanmean(keep_zp), np.nanstd(keep_zp_er)
    mp.plot(avg_mag, zp, 'ro', keep_avg_mag, keep_zp, 'bo')
    mp.plot(linex, liney, 'k-')
#    mp.plot(keep_avg_mag, keep_zp, 'bo')
    mp.show()

    cutoff = np.nanmedian(keep_zp)

    final_zp = []
    final_avg_mag = []

    for ind in range(len(keep_zp)):
        if keep_zp[ind] > cutoff:
            final_zp.append(keep_zp[ind])
            final_avg_mag.append(keep_avg_mag[ind])


    mp.plot(final_avg_mag, final_zp, 'ro')
    mp.show()
