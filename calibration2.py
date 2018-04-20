#!/usr/bin/env python


import matplotlib.pyplot as mp
import glob
import config
import os
import numpy as np
from astropy.io import fits

from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from matplotlib.colors import LogNorm

import aper
import dao
from astropy.stats import sigma_clip
import re



def do_ap_phot(target, channel, exptime):

    print 'Doing aperture photometry on deep mosaic.'
    data_dir = config.top_dir+target
    os.chdir(data_dir+'/newmosaics')
    print 'Changed directory to {}'.format(data_dir+'/newmosaics')

    hdu = fits.open(channel+'-deep.fits')[0]
    image = hdu.data.astype(float)
    wcs = WCS(hdu.header)

    coords = np.loadtxt(channel+'-calibrators.reg', dtype=[('ra', 'S15'), ('dec', 'S15')])

    sky_coords = SkyCoord(coords['ra'], coords['dec'], frame='fk5', unit=(u.hour, u.deg))

    fig = mp.figure(figsize=(20,15))
    ax = fig.add_subplot(111, projection=wcs)
    mp.imshow(hdu.data, origin='lower', cmap='gray_r', norm=LogNorm())
    ax.scatter(sky_coords.ra, sky_coords.dec, s=100, edgecolor='cyan', alpha=0.5, facecolor='none', transform=ax.get_transform('fk5'))
    mp.show()

    pixels = sky_coords.to_pixel(wcs)
    xpix = pixels[0]
    ypix = pixels[1]

    alf_data = dao.read_alf(channel+'-deep.alf')

    cal_data = np.zeros(len(xpix), dtype=[('id', int), ('x', float), ('y', float), ('mag', float), ('err', float)])
    ii = 0
    for x, y in zip(xpix, ypix):

        dist = np.sqrt((alf_data['x'] - x)**2 + (alf_data['y'] - y)**2)
        match = np.argmin(dist)
        cal_data['id'][ii] = alf_data['id'][match]
        cal_data['x'][ii] = alf_data['x'][match]
        cal_data['y'][ii] = alf_data['y'][match]
        cal_data['mag'][ii] = alf_data['mag'][match]
        cal_data['err'][ii] = alf_data['err'][match]
        ii += 1

    if channel == 'I1':
        f_vega = 280.9
        gain = 3.7
        ap_corr = 1.1257 # aperture correction for 3,3,7 pixel aperture
        fluxconv = 0.1257*4
    if channel == 'I2':
        f_vega = 179.7
        gain = 3.71
        ap_corr = 1.1226
        fluxconv = 0.1447*4
    zmag = 2.5*np.log10(f_vega*exptime/(2.3504e-5*ap_corr*0.6*0.6*fluxconv))

    deep_mag, e_deep_mag, flux, e_flux, sky, e_sky, badflag, outstr  = \
        aper.aper(image, xpix, ypix, phpadu=gain, apr=6, skyrad=[6,14], zeropoint=zmag)

    deep_mag = deep_mag.flatten()
    e_deep_mag = e_deep_mag.flatten()

    #Save results to file
    dt = np.dtype([('id', int), ('x', float), ('y', float), ('psf_mag', float),
        ('psf_err', float), ('ap_mag', float), ('ap_err', float)])
    data_save = np.array(zip(cal_data['id'], cal_data['x'], cal_data['y'], \
        cal_data['mag'], cal_data['err'], deep_mag, e_deep_mag), dtype=dt)
    np.savetxt(channel+'-calibrators.phot', data_save, fmt='%10i %9.2f %9.2f '+2*'%6.3f %5.3f ')


def calculate_zeropoint(targets):

    data_dir = config.top_dir
    os.chdir(data_dir)
    print 'Changed directory to {}'.format(data_dir)

    fig = mp.figure(figsize=(15,6))

    channels = ['I1', 'I2']
    for channel in channels:
        if channel == 'I1': ax = fig.add_subplot(121)
        if channel == 'I2': ax = fig.add_subplot(122)
        all_residuals = np.array([])
        all_errs = np.array([])
        all_mags = np.array([])

        for target in targets:

            dt = np.dtype([('id', int), ('x', float), ('y', float), ('psf_mag', float),
                ('psf_err', float), ('ap_mag', float), ('ap_err', float)])
            cluster_data = np.loadtxt(target+'/newmosaics/'+channel+'-calibrators.phot', dtype=dt)


            residual = cluster_data['ap_mag'] - cluster_data['psf_mag']
            e_residual = np.sqrt(cluster_data['ap_err']**2 + cluster_data['psf_err']**2)

            filtered = sigma_clip(residual, sigma=3, iters=5)
            all_residuals = np.append(all_residuals, residual[~filtered.mask])
            all_errs = np.append(all_errs, e_residual[~filtered.mask])
            all_mags = np.append(all_mags, cluster_data['ap_mag'][~filtered.mask])

            ax.errorbar(cluster_data['ap_mag'][~filtered.mask], residual[~filtered.mask], \
                yerr=e_residual[~filtered.mask], fmt='o', label=target)

        zp = np.mean(all_residuals)
        zp_err = np.std(all_residuals)

        ax.axhline(zp, color='k')
        ax.axhspan(zp-zp_err, zp+zp_err, color='k', alpha=0.2)
        ax.set_ylabel('ap mag - psf_mag')


        if channel == 'I1':
            ax.legend(loc=0)
            ax.set_xlabel('[3.6] mag')
            I1zp = zp
            I1zp_err = zp_err
            print '[3.6] calibration zp is {:.3f} +/- {:.3f}'.format(zp, zp_err)
        if channel == 'I2':
            ax.set_xlabel('[4.5] mag')
            I2zp = zp
            I2zp_err = zp_err
            print '[4.5] calibration zp is {:.3f} +/- {:.3f}'.format(zp, zp_err)
    mp.savefig('calibration.pdf', format='pdf')

    return I1zp, I1zp_err, I2zp, I2zp_err

def apply_calibration(targets, I1zp, I2zp):

    channels = ['I1', 'I2']
    zps = [I1zp, I2zp]

    for target in targets:

        for ii, channel in enumerate(channels):
            data_dir = config.top_dir+target
            alf_list = glob.glob(data_dir+'/newmosaics/'+channel+'_*.alf')
            deep_alf = data_dir+'/newmosaics/'+channel+'-deep.alf'

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

            data['c4'] = data['c4'] + zps[ii]

            f = open(new_file, 'w')
            f.write(header1)
            f.write(header2)
            f.write(header3)
            f.close()

            f = open(new_file, 'a')
            np.savetxt(f, data, fmt='%7i %8.3f %8.3f %8.3f %8.4f %8.2f %8s %8.2f %8.3f')
            f.close()

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

                cflux = flux

                cmag = -2.5*np.log10(cflux)
                data['c4'] = cmag + zps[ii]

                f = open(new_file, 'w')
                f.write(header1)
                f.write(header2)
                f.write(header3)
                f.close()

                f = open(new_file, 'a')
                np.savetxt(f, data, fmt='%7i %8.3f %8.3f %8.3f %8.4f %8.2f %8s %8.2f %8.3f')
                f.close()

def apply_calibration_deep_mosaic(target, channel, zp):

    data_dir = config.top_dir+target
    deep_alf = data_dir+'/newmosaics/'+channel+'-deep.alf'


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
