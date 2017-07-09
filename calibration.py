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


def find_stars_in_cat(optical_folder, target, channel, folder=''):

    cat_ids, x, y, ra, dec = optical.read_optical_fnl(optical_folder, target)

    reg_file = open(folder+channel+'.reg').read().replace(':', ' ')
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

    alf_list = glob.glob(folder+'all/'+channel+'*.alf')

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

        np.savetxt(folder+'cal_stars/'+channel+'_'+phot_data['id'][ind]+'.coo',
            np.c_[phot_data['aor'][ind], phot_data['f_num'][ind],
            phot_data['x'][ind],phot_data['y'][ind],
            phot_data['psf_mag'][ind],
            phot_data['psf_err'][ind]],
            header=str(phot_data['id'][ind])+' '+str(phot_data['ra'][ind])+' '+str(phot_data['dec'][ind])+' '+str(phot_data['neigh'][ind])+'\n',
            comments='', fmt='%8i %2i %7.3f %7.3f %6.3f %6.4f')




def find_zp(channel, verbose=0, folder=''):

    phot_list = glob.glob(folder+'cal_stars/'+channel+'*.phot')

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
    mp.savefig(folder+channel+'-raw-zp.eps', format='eps')
#    mp.show()
    mp.gcf().clear()

    f, ax = mp.subplots(2)
    ax[0].plot(num_neighbors, zp, 'ro')
    ax[0].set_xlabel('Number of neighbors')
    ax[0].set_ylabel('ZP')
    ax[1].plot(zp_er, zp, 'ro')
    ax[1].set_xlabel('ZP error')
    ax[1].set_ylabel('ZP')
    mp.savefig(folder+channel+'-quality-check.eps', format='eps')
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
    mp.savefig(folder+channel+'-final-zp.eps', format='eps')
    mp.gcf().clear()
#    mp.plot(avg_mag, zp_er, 'ro', good_avg_mag, good_zp_er, 'bo')
#    mp.show()
def find_zp_single_frame(infile):


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

def apply_calibration(channel, zp, folder=''):

    alf_list = glob.glob(folder+'all/'+channel+'*.alf')

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

        data['c4'] = data['c4'] + zp

        f = open(new_file, 'w')
        f.write(header1)
        f.write(header2)
        f.write(header3)
        f.close()

        f = open(new_file, 'a')
        np.savetxt(f, data, fmt='%7i %8.3f %8.3f %8.3f %8.4f %8.2f %8s %8.2f %8.3f')
        f.close()
