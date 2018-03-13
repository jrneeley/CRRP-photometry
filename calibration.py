#!/usr/bin/env python

import read_dao
import matplotlib.pyplot as mp
import numpy as np
import optical
from astropy.io import ascii
import math
import coordinates
import StringIO
import glob
import sys
import re
from astropy.stats import sigma_clip
import daophot_setup
import pexpect
import os



def find_stars_in_cat(optical_dir, target, channel, data_dir=''):

    cat_ids, x, y, ra, dec = optical.read_optical_fnl(optical_dir, target)

    reg_file = open(data_dir+channel+'.reg').read().replace(':', ' ')
    dtype1 = np.dtype([('ra_h', int), ('ra_m', int), ('ra_s', float),
        ('dec_d', int), ('dec_m', int), ('dec_s', float)])
    data = np.loadtxt(StringIO.StringIO(reg_file), dtype=dtype1)
    ra_h = data['ra_h']
    ra_m = data['ra_m']
    ra_s = data['ra_s']
    dec_d = data['dec_d']
    dec_m = data['dec_m']
    dec_s = data['dec_s']

    cal_ra, cal_dec = coordinates.hms2deg(ra_h, ra_m, ra_s, dec_d, dec_m, dec_s)

    alf_list = glob.glob(data_dir+'data/'+channel+'*.alf')

    phot_data = np.zeros(len(cal_ra), dtype=[('id', 'S8'), ('ra', float), ('dec', float),
        ('neigh', int), ('aor', int, len(alf_list)), ('f_num', int, len(alf_list)),
        ('x', float, len(alf_list)), ('y', float, len(alf_list)),
        ('psf_mag', float, len(alf_list)), ('psf_err', float, len(alf_list))])
    neighbors = np.zeros(len(cal_ra), dtype=object)

# Find id numbers of selected stars in optical catalog - and all stars within 6 arcsec
    for obj in range(0,len(cal_ra)):

        dist = coordinates.radial_dist(cal_ra[obj], cal_dec[obj], ra, dec)

        cat_match = np.argmin(dist)
        nei_match = np.argwhere(dist < 6.0)
        num_neighbors = len(nei_match)

        phot_data['id'][obj] = cat_ids[cat_match]
        phot_data['ra'][obj] = ra[cat_match]
        phot_data['dec'][obj] = dec[cat_match]
        phot_data['neigh'][obj] = num_neighbors
        neighbors[obj] = cat_ids[nei_match]
# Find the selected stars in the individual alf files
    for ind in range(0,len(alf_list)):

        alf_id, x, y, alf_mag, alf_err = read_dao.read_alf(alf_list[ind])

        for ind2 in range(0,len(cal_ra)):
            alf_match = np.argwhere(alf_id == int(phot_data['id'][ind2]))
            alf_nei_match = np.argwhere(alf_id == neighbors[ind2])
            trash1, aor_num, frame_num, trash2 = alf_list[ind].split('_')
            phot_data['aor'][ind2,ind] = int(aor_num)
            phot_data['f_num'][ind2, ind] = int(frame_num)

            if len(alf_match):
                phot_data['x'][ind2, ind] = x[alf_match]
                phot_data['y'][ind2, ind] = y[alf_match]
                phot_data['psf_mag'][ind2, ind] = alf_mag[alf_match]
                phot_data['psf_err'][ind2, ind] = alf_err[alf_match]
            else:
                phot_data['x'][ind2,ind] = float('NaN')
                phot_data['y'][ind2,ind] = float('NaN')
                phot_data['psf_mag'][ind2,ind] = float('NaN')
                phot_data['psf_err'][ind2,ind] = float('NaN')

    print 'Writing files...'
#    print phot_data['x'][1]
    for ind in range(0,len(cal_ra)):
        if np.all(np.isnan(phot_data['psf_mag'][ind])):
            continue

        np.savetxt(data_dir+'cal_stars/'+channel+'_'+phot_data['id'][ind]+'.coo',
            np.c_[phot_data['aor'][ind], phot_data['f_num'][ind],
            phot_data['x'][ind],phot_data['y'][ind],
            phot_data['psf_mag'][ind],
            phot_data['psf_err'][ind]],
            header=str(phot_data['id'][ind])+' '+str(phot_data['ra'][ind])+' '+str(phot_data['dec'][ind])+' '+str(phot_data['neigh'][ind])+'\n',
            comments='', fmt='%8i %2i %7.3f %7.3f %6.3f %6.4f')




def find_zp_old(channel, verbose=0, data_dir=''):

    phot_list = glob.glob(data_dir+'cal_stars/'+channel+'*.phot')

    n_stars = len(phot_list)
    zp = np.zeros(n_stars)
    zp_er = np.zeros(n_stars)
    avg_mag = np.zeros(n_stars)
    avg_mag_er = np.zeros(n_stars)
    id_num = np.zeros(n_stars, dtype=int)
    ra = np.zeros(n_stars)
    dec = np.zeros(n_stars)
    num_neighbors = np.zeros(n_stars, dtype=int)
    for ind, phot in enumerate(phot_list):


        f = open(phot, 'r')
        header = f.readline()
        f.close()
        head1, head2, head3, head4 = header.split()

        id_num[ind] = int(head1)
        ra[ind] = float(head2)
        dec[ind] = float(head3)
        num_neighbors[ind] = int(head4)

        dtype1 = np.dtype([('psf_mag',float), ('psf_err', float),
            ('ap_mag', float), ('ap_err', float)])
        data = np.loadtxt(phot, dtype=dtype1, usecols=(4,5,6,7), skiprows=1)

        if np.all(np.isnan(data['ap_mag'])) :
            zp[ind] = float('NaN')
            zp_er[ind] = float('NaN')
            avg_mag[ind] = float('NaN')
            avg_mag_er[ind] = float('NaN')
            continue
        median_ap = np.nanmedian(data['ap_mag'])
        mean_ap = np.nanmean(data['ap_mag'])
        median_psf = np.nanmedian(data['psf_mag'])
        mean_psf = np.nanmean(data['psf_mag'])
        std_ap = np.nanstd(data['ap_mag'])
        std_psf = np.nanstd(data['psf_mag'])

        filtered_ap = np.array(data['ap_mag'])
        filtered_psf = np.array(data['psf_mag'])
        filtered_ap_er = np.array(data['ap_err'])
        filtered_psf_er = np.array(data['psf_err'])
        for i in range(5):
            outliers_ap = np.argwhere(abs(filtered_ap - mean_ap) > 3*std_ap)
            filtered_ap[outliers_ap] = float('NaN')
            filtered_ap_er[outliers_ap] = float('NaN')
            mean_ap = np.nanmean(filtered_ap)
            std_ap = np.nanstd(filtered_ap)
            outliers_psf = np.argwhere(abs(filtered_psf - mean_psf) > 3*std_psf)
            filtered_psf[outliers_psf] = float('NaN')
            filtered_psf_er[outliers_psf] = float('NaN')
            mean_psf = np.nanmean(filtered_psf)
            std_psf = np.nanstd(filtered_psf)

        # Show aperture photometry and PSF photometry for all stars
        if verbose == 1:
            f, axarr = mp.subplots(2, sharex=True)
            axarr[0].plot(data['ap_mag'], 'ro')
            axarr[0].plot(filtered_ap, 'bo')
            axarr[1].plot(data['psf_mag'], 'ro')
            axarr[1].plot(filtered_psf, 'bo')
            axarr[1].set_xlabel('Observation number')
            axarr[1].set_ylabel('PSF mag')
            axarr[0].set_ylabel('Aperture mag')
            axarr[0].set_ylim(np.nanmax(data['ap_mag'])+0.05, np.nanmin(data['ap_mag'])-0.05)
            axarr[1].set_ylim(np.nanmax(data['psf_mag'])+0.05, np.nanmin(data['psf_mag'])-0.05)
            mp.show()

        residual = filtered_ap - filtered_psf
        res_err = np.sqrt( filtered_ap_er**2 + filtered_psf_er**2)
        # Show residuals for each bcd for all stars
    #    mp.errorbar(np.arange(len(residual)), residual, yerr=res_err , fmt='o')
    #    mp.show()
        zp[ind] = np.nanmean(residual)
        zp_er[ind] = np.nanstd(residual)
        avg_mag[ind] = mean_ap
        avg_mag_er[ind] = std_ap



    mp.errorbar(avg_mag, zp, yerr=zp_er, fmt='o', color='r')
    mp.xlabel('Aperture mag')
    mp.ylabel('Aperture - PSF mag')
    mp.savefig(data_dir+channel+'-raw-zp.eps', format='eps')
#    mp.show()
    mp.gcf().clear()

    f, ax = mp.subplots(2)
    ax[0].plot(num_neighbors, zp, 'ro')
    ax[0].set_xlabel('Number of neighbors')
    ax[0].set_ylabel('ZP')
    ax[1].plot(zp_er, zp, 'ro')
    ax[1].set_xlabel('ZP error')
    ax[1].set_ylabel('ZP')
    mp.savefig(data_dir+channel+'-quality-check.eps', format='eps')
    mp.show()
    mp.gcf().clear()
#    mp.show()
    neighbor_threshold = input('Maximum number of neighbors?: ')
    error_threshold = input('Maximum standard deviation?:')

    good_values = np.argwhere((num_neighbors < neighbor_threshold) &
        (zp_er < error_threshold)).flatten()
    good_zp = zp[good_values]
    good_zp_er = zp_er[good_values]
    good_avg_mag = avg_mag[good_values]

    remo=[0]
    while len(remo):
        mean_zp = np.nanmean(good_zp)
        std_zp = np.nanstd(good_zp)
        remo = (abs(good_zp - mean_zp) > 2*std_zp).nonzero()[0]
        good_zp = np.delete(good_zp, remo)
        good_zp_er = np.delete(good_zp_er, remo)
        good_avg_mag = np.delete(good_avg_mag, remo)
    print str(len(good_zp))+' final calibration stars.'
    print 'Mean, median zero point and standard deviation:'
    print np.nanmean(good_zp), np.nanmedian(good_zp), np.nanstd(good_zp)
    #mp.plot(avg_mag, zp, 'ro', good_avg_mag, good_zp, 'bo')
    mp.errorbar(avg_mag, zp, yerr=zp_er, fmt='o', color='r')
    mp.errorbar(good_avg_mag, good_zp, yerr=good_zp_er, fmt='o', color='b')
    mp.axhline(np.nanmean(good_zp))
#    liney = [np.nanmean(good_zp), np.nanmean(good_zp)]
#    mp.plot(linex, liney, 'k-')
    mp.ylabel('Aperture - PSF')
    mp.xlabel('Aperture mag')
#    mp.show()
    mp.savefig(data_dir+channel+'-final-zp.eps', format='eps')
    mp.gcf().clear()
#    mp.plot(avg_mag, zp_er, 'ro', good_avg_mag, good_zp_er, 'bo')
#    mp.show()
def find_zp_single_frame_old(infile):


    dtype1 = np.dtype([('id', int), ('x', float), ('y', float), ('ap_mag', float),
        ('ap_er', float), ('psf_mag', float), ('psf_er', float)])
    data = np.loadtxt(infile, dtype = dtype1)

    bad_phot = (data['ap_mag'] < 0).nonzero()[0]
    bad_phot2 = (data['ap_mag'] > 30).nonzero()[0]
    bad_phot = np.append(bad_phot, bad_phot2)
    if len(bad_phot):
        data['ap_mag'][bad_phot] = float('NaN')
        data['ap_er'][bad_phot] = float('NaN')

    zp = data['ap_mag'] - data['psf_mag']
    zp_er = np.sqrt( data['ap_er']**2 + data['psf_er']**2 )

    mean_zp = np.nanmean(zp)
    std_zp = np.nanstd(zp)

    remo = [0]
    while len(remo):
        remo = (abs(zp - mean_zp) > 2*std_zp).nonzero()[0]
        zp[remo] = float('NaN')
        zp_er[remo] = float('NaN')

    print np.nanmean(zp), np.nanstd(zp)
    mp.plot(data['ap_mag'], zp, 'ro')
    mp.show()

def find_zp_bcds(channel, data_dir='', mosaics=0):

    if mosaics == 0:
        fits_list = glob.glob(data_dir+'data/'+channel+'*[0-9].fits')
        if channel == 'I1': zmag = 17.12
        if channel == 'I2': zmag = 16.65
    if mosaics == 1:
        fits_list = glob.glob(data_dir+'mosaics/'+channel+'*[0-9].fits')
        if channel == 'I1': zmag = 18.67
        if channel == 'I2': zmag = 18.19

    diff1 = np.array([], dtype=float)
    diff2 = np.array([], dtype=float)
    names = np.array([], dtype=int)
    ap_mags = np.array([], dtype=float)
    obs_num = np.array([], dtype=float)
    err = np.array([], dtype=float)

    for img in fits_list:
        coo_file = re.sub('.fits', '.coo', img)
        ap_file = re.sub('.fits', '.ap', img)

        ids, ap_phot, ap_err = read_dao.read_ap(ap_file)
        ap_phot = ap_phot - 25.0 + zmag
        data = read_dao.read_coo_new(coo_file)

        diff1 = np.append(diff1, ap_phot - data['mag1'])
        diff2 = np.append(diff2, ap_phot - data['mag2'])
        err = np.append(err, np.sqrt(ap_err**2 + data['err']**2))
        obs_num = np.append(obs_num, np.arange(len(diff1)))
        ap_mags = np.append(ap_mags, ap_phot)
        names = np.append(names, ids)

    # delete residuals where the aperture photometry was missing
    ap_mags = np.delete(ap_mags, np.argwhere(diff1 > 50))
    obs_num = np.delete(obs_num, np.argwhere(diff1 > 50))
    diff2 = np.delete(diff2, np.argwhere(diff1 > 50))
    err = np.delete(err, np.argwhere(diff1 > 50))
    names = np.delete(names, np.argwhere(diff1 > 50))
    diff1 = np.delete(diff1, np.argwhere(diff1 > 50))

    # delete saturated stars
    obs_num = np.delete(obs_num, np.argwhere(ap_mags < 10))
    diff2 = np.delete(diff2, np.argwhere(ap_mags < 10))
    err = np.delete(err, np.argwhere(ap_mags < 10))
    names = np.delete(names, np.argwhere(ap_mags < 10))
    diff1 = np.delete(diff1, np.argwhere(ap_mags < 10))
    ap_mags = np.delete(ap_mags, np.argwhere(ap_mags < 10))

    # delete faint stars
    obs_num = np.delete(obs_num, np.argwhere(ap_mags > 14))
    diff2 = np.delete(diff2, np.argwhere(ap_mags > 14))
    err = np.delete(err, np.argwhere(ap_mags > 14))
    names = np.delete(names, np.argwhere(ap_mags > 14))
    diff1 = np.delete(diff1, np.argwhere(ap_mags > 14))
    ap_mags = np.delete(ap_mags, np.argwhere(ap_mags > 14))

    # sigma clip full sample
    filtered_diff1 = sigma_clip(diff1, sigma=3, iters=5)
    filtered_diff2 = sigma_clip(diff2, sigma=3, iters=5)


    mp.plot(ap_mags, diff2, 'ro')
    mp.plot(ap_mags, filtered_diff2, 'o')
    mp.show()
    print np.mean(filtered_diff1), np.std(filtered_diff1)
    print np.mean(filtered_diff2), np.std(filtered_diff2)

    psf_stars = np.unique(names)
    print 'There are {} psf stars.'.format(len(psf_stars))

    avg_zp1 = np.zeros(len(psf_stars))
    avg_zp2 = np.zeros(len(psf_stars))
    std_zp1 = np.zeros(len(psf_stars))
    std_zp2 = np.zeros(len(psf_stars))
    avg_mag = np.zeros(len(psf_stars))
    for ind, star in enumerate(psf_stars):
        x = obs_num[names == star]
        y3 = ap_mags[names == star]
        y1 = diff1[names == star]
        y2 = diff2[names == star]
        psf = y3 - y2
        filtered_y3 = sigma_clip(y3, sigma=2, iters=5)
        if mosaics == 0:
            mp.plot(filtered_y3, y1, 'ro')
            mp.plot(filtered_y3, y2, 'bo')
        if mosaics == 1:
            mp.plot(y3, y1, 'ro')
            mp.plot(y3, y2, 'bo')
        mp.xlabel('Aperture phot mag')
        mp.ylabel('Residual')
        mp.title(star)
        mp.show()
        ap_range = np.max(filtered_y3) - np.min(filtered_y3)
        psf_range = np.max(filtered_y3 - y2) - np.min(filtered_y3 - y2)
        print ap_range, psf_range
        print np.std(filtered_y3), np.std(filtered_y3 - y2)
        avg_zp2[ind] = np.mean(y2[y3 == filtered_y3])
        avg_zp1[ind] = np.mean(y1[y3 == filtered_y3])
        std_zp1[ind] = np.std(y1[y3 == filtered_y3])
        std_zp2[ind] = np.std(y2[y3 == filtered_y3])
        avg_mag[ind] = np.mean(y3[y3 == filtered_y3])


    filtered_avzp1 = sigma_clip(avg_zp1, sigma=3, iters=5)
    filtered_avzp2 = sigma_clip(avg_zp2, sigma=3, iters=5)
    mp.errorbar(avg_mag, filtered_avzp1, yerr=std_zp1, fmt='ro')
    mp.errorbar(avg_mag, filtered_avzp2, yerr=std_zp2, fmt='o')
    mp.xlabel('Average mag')
    mp.ylabel('Zero point')
    mp.show()

    print np.mean(filtered_avzp1), np.std(filtered_avzp1)
    print np.mean(filtered_avzp2), np.std(filtered_avzp2)

    return np.mean(filtered_avzp2), np.std(filtered_avzp2)

def find_zp_deep_mosaic(target, channel, data_dir=''):

    if channel == 'I1': zmag = 18.67
    if channel == 'I2': zmag = 18.19
    deep_mosaic = data_dir+'DeepMosaic/'+target+'_'+channel+'_deep.fits'


    coo_file = re.sub('.fits', '.coo', deep_mosaic)
    ap_file = re.sub('.fits', '.ap', deep_mosaic)

    ids, ap_mags, ap_err = read_dao.read_ap(ap_file)
    ap_mags = ap_mags - 25.0 + zmag
    data = read_dao.read_coo_new(coo_file)

    diff1 = ap_mags - data['mag1']
    diff2 = ap_mags - data['mag2']
    err = np.sqrt(ap_err**2 + data['err']**2)

    # delete residuals where the aperture photometry was missing
    ap_mags = np.delete(ap_mags, np.argwhere(diff1 > 50))
    diff2 = np.delete(diff2, np.argwhere(diff1 > 50))
    err = np.delete(err, np.argwhere(diff1 > 50))
    diff1 = np.delete(diff1, np.argwhere(diff1 > 50))

    # delete saturated stars
    diff2 = np.delete(diff2, np.argwhere(ap_mags < 10))
    err = np.delete(err, np.argwhere(ap_mags < 10))
    diff1 = np.delete(diff1, np.argwhere(ap_mags < 10))
    ap_mags = np.delete(ap_mags, np.argwhere(ap_mags < 10))

    # delete faint stars
    diff2 = np.delete(diff2, np.argwhere(ap_mags > 14))
    err = np.delete(err, np.argwhere(ap_mags > 14))
    diff1 = np.delete(diff1, np.argwhere(ap_mags > 14))
    ap_mags = np.delete(ap_mags, np.argwhere(ap_mags > 14))

    # sigma clip full sample
    filtered_diff1 = sigma_clip(diff1, sigma=3, iters=5)
    filtered_diff2 = sigma_clip(diff2, sigma=3, iters=5)

    zp = np.mean(filtered_diff2)
    zp_er = np.std(filtered_diff2)

    #mp.plot(ap_mags, diff2, 'ro')
    fig, ax = mp.subplots(1)
    ax.plot(ap_mags, filtered_diff2, 'o')
    ax.axhline(zp, color='k')
    ax.axhline(zp+zp_er, color='k', alpha=0.5)
    ax.axhline(zp-zp_er, color='k', alpha=0.5)
    mp.show()
    print np.mean(filtered_diff1), np.std(filtered_diff1)
    print np.mean(filtered_diff2), np.std(filtered_diff2)

    psf_stars = ids
    print 'There are {} psf stars.'.format(len(psf_stars))

    return zp, zp_er

def apply_calibration_bcds(target, channel, zp, data_dir='', mosaics=0):

    if mosaics == 0:
        alf_list = glob.glob(data_dir+'data/'+channel+'*.alf')
    if mosaics == 1:
        alf_list = glob.glob(data_dir+'mosaics/'+channel+'*.alf')

    alc_params = [0.98397385,-8.6387206e-05,3.6310403e-05,-5.6284359e-07,1.6023687e-06,1.3574380e-06]

    for alf in alf_list:

        new_file = re.sub('.alf', '.cal', alf)
        f = open(alf, 'r')
        header1 = f.readline()
        header2 = f.readline()
        header3 = f.readline()
        f.close()


        dtype1 = np.dtype([('c1',int), ('c2', float), ('c3', float), ('c4', float),
            ('c5', float), ('c6', float), ('c7', 'S4'), ('c8', float), ('c9', float)])
        data = np.loadtxt(alf, dtype=dtype1, skiprows=3)

        flux = 10**(data['c4'] / -2.5)
        x = np.round(data['c2'])
        y = np.round(data['c3'])
        if mosaics == 0:
            alc = alc_params[0] + alc_params[1]*(x-127) + alc_params[2]*(y-127) + alc_params[3]*(x-127)*(y-127) + alc_params[4]*(x-127)**2 + alc_params[5]*(y-127)**2
            cflux = flux/alc
        if mosaics == 1:
            cflux = flux

        cmag = -2.5*np.log10(cflux)
        data['c4'] = cmag + zp

        f = open(new_file, 'w')
        f.write(header1)
        f.write(header2)
        f.write(header3)
        f.close()

        f = open(new_file, 'a')
        np.savetxt(f, data, fmt='%7i %8.3f %8.3f %8.3f %8.4f %8.2f %8s %8.2f %8.3f')
        f.close()

def apply_calibration_deep_mosaic(target, channel, zp, data_dir=''):

    deep_alf = data_dir+'DeepMosaic/'+target+'_'+channel+'_deep_dn.alf'


    new_file = re.sub('.alf', '.cal', deep_alf)
    f = open(deep_alf, 'r')
    header1 = f.readline()
    header2 = f.readline()
    header3 = f.readline()
    f.close()


    dtype1 = np.dtype([('c1',int), ('c2', float), ('c3', float), ('c4', float),
        ('c5', float), ('c6', float), ('c7', 'S4'), ('c8', float), ('c9', float)])
    data = np.loadtxt(deep_alf, dtype=dtype1, skiprows=3)

    x = np.round(data['c2'])
    y = np.round(data['c3'])

    data['c4'] = data['c4'] + zp

    f = open(new_file, 'w')
    f.write(header1)
    f.write(header2)
    f.write(header3)
    f.close()

    f = open(new_file, 'a')
    np.savetxt(f, data, fmt='%7i %8.3f %8.3f %8.3f %8.4f %8.2f %8s %8.2f %8.3f')
    f.close()

def find_lst_stars(target, channel, data_dir='', mosaics=0):

    lst_file = data_dir+'DeepMosaic/'+target+'_'+channel+'_deep_dn.lst'
    alf_file = data_dir+'DeepMosaic/'+target+'_'+channel+'_deep_dn.alf'

    if mosaics == 0:
        alf_all_list = glob.glob(data_dir+'data/'+channel+'*.alf')
    if mosaics == 1:
        alf_all_list = glob.glob(data_dir+'mosaics/'+channel+'*.alf')

    id_init, x_init, y_init = read_dao.read_lst(lst_file)
    id_alf, x_new, y_new, psf_mag, psf_err = read_dao.read_alf(alf_file)

    n_psf_stars = len(x_init)
    psf_stars = np.zeros(n_psf_stars, dtype=object)
    psf_star_mags = np.zeros(n_psf_stars, dtype=float)
    psf_star_errs = np.zeros(n_psf_stars, dtype=float)
    psf_star_cmags = np.zeros(n_psf_stars, dtype=float)
    psf_star_x = np.zeros(n_psf_stars, dtype=float)
    psf_star_y = np.zeros(n_psf_stars, dtype=float)
    neigh_dists = np.zeros(n_psf_stars, dtype=object)
    psf_star_id = np.zeros(n_psf_stars, dtype=int)
    for idx, (x, y) in enumerate(zip(x_init, y_init)):

        dist = np.sqrt((x_new - x)**2 + (y_new - y)**2)
        match_x = x_new[dist == np.min(dist)]
        match_y = y_new[dist == np.min(dist)]
        match_id = id_alf[dist == np.min(dist)]
        match_mag = psf_mag[dist == np.min(dist)]
        match_err = psf_err[dist == np.min(dist)]
        # find all stars within 10 (mosiac) pixels of the calibrators
        neighbor_ids = id_alf[dist < 10]
        neighbor_mags = psf_mag[dist < 10]
        neighbor_dists = dist[dist < 10]
        cmag = match_mag
        if len(neighbor_ids) > 0 :
            neighbor_ids = np.delete(neighbor_ids, np.argwhere(neighbor_ids == match_id))
            neighbor_mags = np.delete(neighbor_mags, np.argwhere(neighbor_mags == match_mag))
            neighbor_dists = np.delete(neighbor_dists, np.argwhere(neighbor_dists == np.min(neighbor_dists)))
            if len(neighbor_ids) > 0:
                flux = 10**(match_mag/ -2.5)
                for n, d in zip(neighbor_mags, neighbor_dists):
                    fluxratio = get_flux_ratio(channel, d/2.0)
                    flux_neigh = 10**(n/ -2.5)
                    flux += flux_neigh*fluxratio
                cmag = -2.5*np.log10(flux)
        psf_stars[idx] = np.append(match_id, neighbor_ids)
        psf_star_mags[idx] = match_mag
        psf_star_errs[idx] = match_err
        psf_star_x[idx] = match_x
        psf_star_y[idx] = match_y
        psf_star_cmags[idx] = cmag
        psf_star_id[idx] = match_id
        neigh_dists[idx] = neighbor_dists
        #print match_id, neighbor_ids

    # we only need this step if it wasn't already done.
    if mosaics == 1:
        # make coo file for deep mosaic
        coo_mosaic = re.sub('_dn.lst', '.coo', lst_file)
        f_handle = open(coo_mosaic, 'w')
        f = open(lst_file, 'r')
        for ii in range(3):
            line = f.readline()
            f_handle.write(line)
        f.close()
        f_handle.close()
        f_handle = open(coo_mosaic, 'a')
        data_save = np.array(zip(psf_star_id, psf_star_x, psf_star_y, psf_star_mags,
            psf_star_errs, psf_star_cmags, np.repeat(-9.999, n_psf_stars)), dtype=[('id', int), ('x', float), ('y', float),
            ('mag', float), ('err', float), ('mag2', float), ('oth2', float)])
        np.savetxt(f_handle, data_save, fmt='%8i %8.2f %8.2f %8.3f %8.3f %8.3f %8.3f')
        f_handle.close()

    # collect multi-epoch allframe photometry
    for alf_file in alf_all_list:
        coo_file = re.sub('_dn.alf', '.coo', alf_file)
        f_handle = open(coo_file, 'w')
        dn_coo = re.sub('.coo', '_dn.coo', coo_file)
        f = open(dn_coo, 'r')
        for ii in range(3):
            line = f.readline()
            f_handle.write(line)
        f.close()
        f_handle.close()

        id_epoch, x_epoch, y_epoch, psf_epoch, err_epoch = read_dao.read_alf(alf_file)

        for star in range(n_psf_stars):
    #        print psf_stars[star][0], psf_stars[star][1:], neigh_dists[star]
            f_handle = open(coo_file, 'a')

            this_x = x_epoch[id_epoch == psf_stars[star][0]]
            this_y = y_epoch[id_epoch == psf_stars[star][0]]
            this_mag = psf_epoch[id_epoch == psf_stars[star][0]]
            this_err = err_epoch[id_epoch == psf_stars[star][0]]
            if len(this_mag) > 0:
                flux = 10**(this_mag[0]/ -2.5)
            #    ds = np.array(neigh_dists[star], dtype=float)

                for neigh, d in zip(psf_stars[star][1:], neigh_dists[star]):
                    neigh_mag = psf_epoch[id_epoch == neigh]
                    if len(neigh_mag) > 0:
                    #    print 'Adding {} to {}...'.format(neigh_mag, this_mag)
                        fluxratio = get_flux_ratio(channel, d/2.0)
                        flux_neigh = 10**(neigh_mag[0]/ -2.5)
                        flux += flux_neigh*fluxratio
                corr_mag = -2.5*np.log10(flux)
        #    print this_mag, corr_mag

                data_save = np.array(zip([psf_stars[star][0]], this_x, this_y, this_mag,
                    this_err, [corr_mag], [-9.999]), dtype=[('id', int), ('x', float), ('y', float),
                    ('mag', float), ('err', float), ('mag2', float), ('oth2', float)])
                np.savetxt(f_handle, data_save, fmt='%8i %8.2f %8.2f %8.3f %8.3f %8.3f %8.3f')
            f_handle.close()


def do_ap_phot(target, channel, exptime, data_dir='', dao_dir='', mosaics=0, opt_dir='/Volumes/Annie/CRRP/OPTfiles/'):

    # list original [in flux units] fits files
    if mosaics == 0:
        fits_list = glob.glob(data_dir+'data/'+channel+'*[0-9].fits')
    if mosaics == 1:
        fits_list = glob.glob(data_dir+'mosaics/'+channel+'*[0-9].fits')
    deep_mosaic = target+'_'+channel+'_deep.fits'

    # copy appropriate OPT files for invidual BCDs

    if mosaics == 0:
        daophot_setup.set_opt_files(opt_dir, channel, exptime, warm=1)
    if mosaics == 1:
        daophot_setup.set_opt_files_mosaic(opt_dir, channel, exptime, warm=1)
    # start DAOPHOT
    daophot = pexpect.spawn(dao_dir+'daophot')
    daophot.logfile = sys.stdout
    for image in fits_list:
        if mosaics == 0:
            img = re.sub(data_dir+'data/', target+':', image)
        if mosaics == 1:
            img = re.sub(data_dir+'mosaics/', target+'m:', image)
        daophot.expect("Command:")
        daophot.sendline("at " + img)
        daophot.expect("Command:")
        daophot.sendline("phot")
        daophot.expect("File with aperture radii")
        daophot.sendline("")
        daophot.expect("PHO>")
        daophot.sendline("")
        daophot.expect("Input position file")
        daophot.sendline('')
        daophot.expect("Output file")
        daophot.sendline('')
    daophot.expect('Command:')
    daophot.sendline('exit')
    daophot.close(force=True)

    # only need to do this if is isn't already done
    if mosaics == 1:
        os.chdir(data_dir+'/DeepMosaic')
        print 'Changed directory to DeepMosaic\n'
        # Now do aperture photometry on the deep mosaic
        daophot_setup.set_opt_files_mosaic(opt_dir, channel, exptime, warm=1)
        lst_file = re.sub('.fits', '.coo', deep_mosaic)
        ap_file = re.sub('.fits', '.ap', deep_mosaic)
        # start DAOPHOT
        daophot = pexpect.spawn(dao_dir+'daophot')
        daophot.logfile = sys.stdout
        img = re.sub(data_dir+'data/', target+':', image)
        daophot.expect("Command:")
        daophot.sendline("at " + deep_mosaic)
        daophot.expect("Command:")
        daophot.sendline("phot")
        daophot.expect("File with aperture radii")
        daophot.sendline("")
        daophot.expect("PHO>")
        daophot.sendline("")
        daophot.expect("Input position file")
        daophot.sendline(lst_file)
        daophot.expect("Output file")
        daophot.sendline(ap_file)
        daophot.expect('Command:')
        daophot.sendline('exit')
        daophot.close(force=True)
        os.chdir(data_dir)


def get_flux_ratio(channel, pixel_distance):
    # pixel distance is in Native IRAC pixels
    pix_from_centroid = np.arange(51)
    pix_from_centroid = pix_from_centroid/5.0
    if channel == 'I1':
        flux_ratio = np.array([1.0000, 0.9994, 0.9980, 0.9955, 0.9911, 0.9838,
            0.9719, 0.9537, 0.9270, 0.8902, 0.8396, 0.7707, 0.6834, 0.5833,
            0.4802, 0.3837, 0.3007, 0.2340, 0.1813, 0.1391, 0.1056, 0.0797,
            0.0602, 0.0463, 0.0365, 0.0296, 0.0248, 0.0215, 0.0191, 0.0173,
            0.0158, 0.0146, 0.0135, 0.0125, 0.0117, 0.0109, 0.0102, 0.0096,
            0.0091, 0.0086, 0.0081, 0.0077, 0.0073, 0.0069, 0.0066, 0.0062,
            0.0059, 0.0057, 0.0054 ,0.0052 ,0.0051])
    if channel == 'I2':
        flux_ratio = np.array([1.0000, 0.9990, 0.9967, 0.9921, 0.9842, 0.9720,
            0.9547, 0.9320, 0.9031, 0.8650, 0.8115, 0.7366, 0.6417, 0.5358,
            0.4288, 0.3316, 0.2537, 0.1969, 0.1553, 0.1218, 0.0934, 0.0701,
            0.0519, 0.0388, 0.0300, 0.0243, 0.0206, 0.0181, 0.0163, 0.0150,
            0.0139, 0.0129, 0.0120, 0.0111, 0.0103, 0.0095, 0.0089, 0.0084,
            0.0079, 0.0075, 0.0072, 0.0070, 0.0067, 0.0065, 0.0063, 0.0060,
            0.0058, 0.0056, 0.0053, 0.0050, 0.0048])

    idx = np.searchsorted(pix_from_centroid, pixel_distance)

    return flux_ratio[idx]
