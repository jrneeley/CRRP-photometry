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

def find_stars_in_cat(optical_folder, target, channel):

    cat_ids, x, y, ra, dec = optical.read_optical_catalog(optical_folder, target)

    dtype1 = np.dtype([('ra_h', 'i2'), ('ra_m', 'i2'), ('ra_s', 'f6'), ('dec_d', 'i3'), ('dec_m', 'i2'), ('dec_s', 'f5')])
    data = np.loadtxt(channel+'.reg', dtype=dtype1)
    ra_h = data['ra_h']
    ra_m = data['ra_m']
    ra_s = data['ra_s']
    dec_d = data['dec_d']
    dec_m = data['dec_m']
    dec_s = data['dec_s']

    cal_ra, cal_dec = coordinates.hms2deg(ra_h, ra_m, ra_s, dec_d, dec_m, dec_s)


    ra1 = np.radians(ra)
    dec1 = np.radians(dec)
    id_match = []
    ra_match = []
    dec_match = []
    num_neighbors = []

    for obj in range(0,len(cal_ra)):
        ra2 = np.radians(cal_ra[obj])
        dec2 = np.radians(cal_dec[obj])
        x1 = np.cos(dec1)*np.cos(ra1)
        y1 = np.cos(dec1)*np.sin(ra1)
        z1 = np.sin(dec1)
        x2 = np.cos(dec2)*np.cos(ra2)
        y2 = np.cos(dec2)*np.sin(ra2)
        z2 = np.sin(dec2)

        dist = np.sqrt( (x2 - x1)**2 + (y2 - y1)**2 + (z2 - z1)**2)
        dist = np.degrees(dist)*3600.

        neighbors = dist[dist < 3.6]
        #if len(neighbors) == 1:
        ind_match = np.argmin(dist)

        id_match.append(cat_ids[ind_match])
        ra_match.append(ra[ind_match])
        dec_match.append(dec[ind_match])
        num_neighbors.append(len(neighbors))


    np.savetxt(channel+'_cal_star.txt', np.c_[id_match, ra_match, dec_match, num_neighbors], fmt='%i %f %f %i')

def find_stars_in_cat2(optical_folder, target, channel):

    cat_ids, x, y, ra, dec = optical.read_optical_catalog(optical_folder, target)

    reg_file = open(channel+'.reg').read().replace(':', ' ')
    dtype1 = np.dtype([('ra_h', 'i2'), ('ra_m', 'i2'), ('ra_s', 'f6'), ('dec_d', 'i3'), ('dec_m', 'i2'), ('dec_s', 'f5')])
    data = np.loadtxt(StringIO.StringIO(reg_file), dtype=dtype1)
    ra_h = data['ra_h']
    ra_m = data['ra_m']
    ra_s = data['ra_s']
    dec_d = data['dec_d']
    dec_m = data['dec_m']
    dec_s = data['dec_s']

    cal_ra, cal_dec = coordinates.hms2deg(ra_h, ra_m, ra_s, dec_d, dec_m, dec_s)

    alf_list = glob.glob('all/'+channel+'*.alf')

    phot_data = np.zeros(len(cal_ra), dtype=[('id', 'S8'), ('ra', float), ('dec', float),
        ('neigh', int), ('aor', 'i8', len(alf_list)), ('f_num', 'i2', len(alf_list)),
        ('x', float, len(alf_list)), ('y', float, len(alf_list)),
        ('psf_mag', float, len(alf_list)), ('psf_err', float, len(alf_list))])

    for obj in range(0,len(cal_ra)):

        dist = coordinates.radial_dist(cal_ra[obj], cal_dec[obj], ra, dec)
        neighbors = dist[dist < 3.6]
        num_neighbors = len(neighbors)
        cat_match = np.argmin(dist)

        phot_data['id'][obj] = cat_ids[cat_match]
        phot_data['ra'][obj] = ra[cat_match]
        phot_data['dec'][obj] = dec[cat_match]
        phot_data['neigh'][obj] = len(neighbors)


    for ind in range(0,len(alf_list)):

        alf_id, x, y, alf_mag, alf_err = read_dao.read_alf(alf_list[ind])
        for ind2 in range(0,len(cal_ra)):
            alf_match = np.argwhere(alf_id == int(phot_data['id'][ind2]))
            trash1, aor_num, frame_num, trash2 = alf_list[ind].split('_')
            phot_data['aor'][ind2,ind] = int(aor_num)
            phot_data['f_num'][ind2, ind] = int(frame_num)

            if len(alf_match):
                phot_data['x'][ind2,ind] = x[alf_match]
                phot_data['y'][ind2,ind] = y[alf_match]
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

        np.savetxt('cal_stars/'+channel+'_'+phot_data['id'][ind]+'.coo',
            np.c_[phot_data['aor'][ind], phot_data['f_num'][ind],
            phot_data['x'][ind],phot_data['y'][ind], phot_data['psf_mag'][ind],
            phot_data['psf_err'][ind]],
            header=str(phot_data['id'][ind])+' '+str(phot_data['ra'][ind])+' '+str(phot_data['dec'][ind])+' '+str(phot_data['neigh'][ind]),
            comments='', fmt='%8i %2i %7.3f %7.3f %6.3f %6.4f')

def find_cal_star_coords(optical_folder, target, channel):
#    optical_folder = '/Users/Jill/CRRP/OpticalCatalogs/'
    #read in lst file of calibration stars
    psf_stars, psf_mags, psf_errs = read_dao.read_alf(channel+'.alf')

    #args = np.argsort(psf_mags)
    #psf_mags = psf_mags[args]
    #psf_stars = psf_stars[args]
    cal_ra = []
    cal_dec = []
    stars = []
    num_neighbors = []
    cat_ids, x, y, ra, dec = optical.read_optical_catalog(optical_folder, target)

    for star in psf_stars:
        index_match = np.argwhere(cat_ids == star)
        if len(index_match):
            stars.append(star)
            cal_ra.append(ra[index_match])
            cal_dec.append(dec[index_match])

            dist = coordinates.radial_dist(ra[index_match], dec[index_match], ra, dec)
            neighbors = dist[dist < 3.6]
            num_neighbors.append(len(neighbors))

    np.savetxt(channel+'_cal_stars.txt', np.c_[stars, cal_ra, cal_dec, num_neighbors], fmt='%i %f %f %i')


def find_cal_stars(target, channel):
    #identify calibration stars

    #read in lst file of calibration stars
#    psf_stars, psf_x, psf_y = read_dao.read_lst('master_'+channel+'.lst')
#    psf_stars, psf_m, psf_e = read_dao.read_alf(channel+'.alf')
    dtype1 = np.dtype([('id', 'i8'), ('neighbors', 'i2')])
    data = np.loadtxt(channel+'_cal_star.txt', usecols=(0,3), dtype=dtype1)
    psf_stars = data['id']
    num_neighbors = data['neighbors']

    # read in raw file with PSF photometry in all frames
    star_ids, mags = read_dao.read_raw('optical-'+channel+'-alf.raw')

    # For each calibration star, find all of the PSF photometry
    bad_index = []
    for star in psf_stars:
        try:
            index_match = star_ids.index(str(star))

        except ValueError:
            bad_index.append(np.argwhere(psf_stars == star))
    if len(bad_index):
        psf_stars = np.delete(psf_stars, bad_index)
        num_neighbors = np.delete(num_neighbors, bad_index)
    print 'Using '+str(len(psf_stars))+' calibration stars.'
    print psf_stars
    print num_neighbors

    return psf_stars, num_neighbors


def find_zp(psf_stars, num_neighbors, channel):

    # read in raw file with PSF photometry in all frames
    #star_ids, mags = read_dao.read_raw('optical_alf.raw')
    star_ids, mags = read_dao.read_raw('optical-'+channel+'-alf.raw')

    zp = []
    zp_er = []
    avg_mag = []
    avg_mag_er = []
    num_nei = []
    skipped = 0

    for star in psf_stars:
        # Find PSF photometry
        index_match = star_ids.index(str(star))
        mag_match = mags[index_match]
        psf_mag = np.array(mag_match[5:len(mag_match)-2:2], dtype='f')
        psf_err = np.array(mag_match[6:len(mag_match)-2:2], dtype='f')
        f_num = np.arange(len(psf_mag)/2)
        if star < 5000000:
            try:
                # Find aperture photometry
                ap_phot = 'lcvs/hand_selected/with_corrections/'+channel+'_'+str(star)+'.txt'
                dtype1 = np.dtype([('mag', 'f7'), ('err', 'f6')])
                data = np.loadtxt(ap_phot, dtype=dtype1, usecols=(2,3))
#        data = ascii.read(ap_phot, delimiter=' ')
#        aor = data['col1']
#        frame = data['col2']

                ap_mag=data['mag']
                ap_err=data['err']
            #ap_mag = np.array(data['col3'], dtype='f')
            #ap_err = np.array(data['col4'], dtype='f')
            except:
                skipped += 1
                print star
                continue
        else:
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
        num_nei.append(num_neighbors[np.where(psf_stars == star)])
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
        #if avg_mag_er[ind] < err_limit:
        if num_nei[ind] < 2:
            keep_zp.append(zp[ind])
            keep_zp_er.append(zp_er[ind])
            keep_avg_mag.append(avg_mag[ind])
#    keep_zp = zp[avg_mag_er < err_limit]
#    keep_zp_er = zp[avg_mag_er < err_limit]
#    keep_avg_mag = zp[avg_mag_er < err_limit]

    linex = [sat_limit,20]
    liney = [np.nanmedian(keep_zp), np.nanmedian(keep_zp)]
    print len(keep_zp)
    print np.nanmean(keep_zp), np.nanstd(keep_zp_er)
    mp.plot(avg_mag, zp, 'ro', keep_avg_mag, keep_zp, 'bo')
    mp.plot(linex, liney, 'k-')
#    mp.plot(keep_avg_mag, keep_zp, 'bo')
    mp.show()
    mp.plot(avg_mag, zp_er, 'ro', keep_avg_mag, keep_zp_er, 'bo')

#    mp.plot(keep_avg_mag, keep_zp, 'bo')
    mp.show()
#    cutoff = np.nanmedian(keep_zp)

#    final_zp = []
#    final_avg_mag = []

#    for ind in range(len(keep_zp)):
#        if keep_zp[ind] > cutoff:
#            final_zp.append(keep_zp[ind])
#            final_avg_mag.append(keep_avg_mag[ind])


#    mp.plot(final_avg_mag, final_zp, 'ro')
#    mp.show()

def find_zp2(channel):

    phot_list = glob.glob('cal_stars/'+channel+'*.phot')

    n_stars = len(phot_list)
    zp = np.array(n_stars)
    zp_er = np.array(n_stars)
    avg_mag = np.array(n_stars)
    avg_mag_er = np.array(n_stars)
    id_num = np.array(n_stars, dtype=np.str_)
    ra = np.array(n_stars)
    dec = np.array(n_stars)
    num_neighbors = np.array(n_stars, dtype=int)
    for phot in phot_list:

        current_ind = phot_list.index(phot)
        print current_ind
        f = open(phot, 'r')
        header = f.readline()
        f.close()
        head1, head2, head3, head4 = header.split()
        id_num[current_ind] = str(head1)
        ra[current_ind] = float(head2)
        dec[current_ind] = float(head3)
        num_neighbors[current_ind] = int(head4)

        dtype1 = np.dtype([('psf_mag',float), ('psf_err', float), ('ap_mag', float), ('ap_err', float)])
        data = np.loadtxt(phot, dtype=dtype1, usecols=(4,5,6,7), skiprows=1)

        residual = data['ap_mag'] - data['psf_mag']


        median_ap = np.nanmedian(data['ap_mag'])
        mean_ap = np.nanmean(data['ap_mag'])
        median_psf = np.nanmedian(data['psf_mag'])
        mean_psf = np.nanmean(data['psf_mag'])
        std_ap = np.nanstd(data['ap_mag'])
        std_psf = np.nanstd(data['psf_mag'])


        zp[current_ind] = np.nanmean(residual)
        zp_er[current_ind] = np.nanstd(residual)
        avg_mag[current_ind] = mean_ap
        avg_mag_er[current_ind] = std_ap

    print num_neighbors
    print np.argwhere(num_neighbors < 2)

    good_zp = zp[num_neighbors < 2]
    mp.plot(avg_mag, zp, 'ro')
    mp.show()
    mp.plot(avg_mag, zp_er, 'ro')
    mp.show()
