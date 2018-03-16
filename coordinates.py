#!/usr/bin/env python

import sys
#import glob
import numpy as np
#from astropy.io import ascii
from astropy.wcs import WCS
import optical


def find_coord_window(img_list):

# initialize boundary conditions
    bxmin = 9999
    bxmax = -9999
    bymin = 9999
    bymax = -9999

    for image in img_list:
        w = WCS(image)
        x1, y1 = w.wcs_pix2world(0, 0, 0)
        x2, y2 = w.wcs_pix2world(0, 255, 0)
        x3, y3 = w.wcs_pix2world(255, 0, 0)
        x4, y4 = w.wcs_pix2world(255, 255, 0)

        xmin = min(x1, x2, x3, x4)
        xmax = max(x1, x2, x3, x4)
        ymin = min(y1, y2, y3, y4)
        ymax = max(y1, y2, y3, y4)

        if (xmin < bxmin):
            bxmin = min(x1, x2, x3, x4)
        if (xmax > bxmax):
            bxmax = max(x1, x2, x3, x4)
        if (ymin < bymin):
            bymin = min(y1, y2, y3, y4)
        if (ymax > bymax):
            bymax = max(y1, y2, y3, y4)

    return(bxmin, bxmax, bymin, bymax)

def find_coord_window_mosaic(img, xmin, xmax, ymin, ymax):

    w = WCS(img)
    x1, y1 = w.wcs_pix2world(xmin, ymin, 0)
    x2, y2 = w.wcs_pix2world(xmin, ymax, 0)
    x3, y3 = w.wcs_pix2world(xmax, ymin, 0)
    x4, y4 = w.wcs_pix2world(xmax, ymax, 0)

    bxmin = min(x1, x2, x3, x4)
    bxmax = max(x1, x2, x3, x4)
    bymin = min(y1, y2, y3, y4)
    bymax = max(y1, y2, y3, y4)

    return(bxmin, bxmax, bymin, bymax)

# Convert celestial coordinates to pixels using optical catalog
def radec2catalogpix(ra_in, dec_in, catalog_x, catalog_y, catalog_ra, catalog_dec):

    ra_diff = np.abs(catalog_ra - ra_in)
    dec_diff = np.abs(catalog_dec - dec_in)
    index_x = np.argmin(ra_diff)
    index_y = np.argmin(dec_diff)

    target_x = catalog_x[index_x]

    target_y = catalog_y[index_y]

    return(target_x, target_y)

# Convert coordinates in sexagesimal to decimal units
def hms2deg(ra_h, ra_m, ra_s, dec_d, dec_m, dec_s):

    if type(dec_d) != np.ndarray:
        if dec_d < 0:
            dec = dec_d - dec_m/60. - dec_s/3600.
        if dec_d >= 0 :
            dec = dec_d + dec_m/60. + dec_s/3600.
        ra = 15*(ra_h+ra_m/60.+ra_s/3600.)
    else:
        dec = np.zeros(len(dec_d))
        for ind, d in enumerate(dec_d):
            if dec_d[ind] < 0:
            #    dec.append(dec_d[ind] - dec_m[ind]/60. - dec_s[ind]/3600.)
                dec[ind] = dec_d[ind] - dec_m[ind]/60. - dec_s[ind]/3600.
            if dec_d[ind] >= 0:
            #    dec.append(dec_d[ind] + dec_m[ind]/60. + dec_s[ind]/3600.)
                dec[ind] = dec_d[ind] + dec_m[ind]/60. + dec_s[ind]/3600.
#    if len(ra_h) > 1:
        ra = np.zeros(len(ra_h))
        for ind, r in enumerate(ra_h):
            #ra.append(15.*(ra_h[ind]+ra_m[ind]/60.+ra_s[ind]/3600.))
            ra[ind] = 15.*(ra_h[ind]+ra_m[ind]/60.+ra_s[ind]/3600.)

    return ra, dec

def radec_string2deg(ra, dec):

    if type(ra) != np.ndarray:

        ra_sep = ra.split(':')
        dec_sep = dec.split(':')
        ra_new, dec_new = hms2deg(float(ra_sep[0]), float(ra_sep[1]),
            float(ra_sep[2]), float(dec_sep[0]), float(dec_sep[1]), float(dec_sep[2]))
    else:
        num_stars = len(ra)

        ra_new = np.zeros(num_stars)
        dec_new = np.zeros(num_stars)
        for ind, star in enumerate(ra):

            ra_sep = ra[ind].split(':')
            dec_sep = dec[ind].split(':')
            #print ra_sep, dec_sep
            ra_deg, dec_deg = hms2deg(float(ra_sep[0]), float(ra_sep[1]),
                float(ra_sep[2]), float(dec_sep[0]), float(dec_sep[1]), float(dec_sep[2]))
            ra_new[ind] = ra_deg
            dec_new[ind] = dec_deg

    return ra_new, dec_new

# Finds radial distance between coordinates in arcsec
# second RA/DEC should be a scalar, first can be scalar or array
def radial_dist(ra1, dec1, ra2, dec2):

    ra1 = np.radians(ra1)
    dec1 = np.radians(dec1)
    ra2 = np.radians(ra2)
    dec2 = np.radians(dec2)


    x1 = np.cos(dec1)*np.cos(ra1)
    y1 = np.cos(dec1)*np.sin(ra1)
    z1 = np.sin(dec1)
    x2 = np.cos(dec2)*np.cos(ra2)
    y2 = np.cos(dec2)*np.sin(ra2)
    z2 = np.sin(dec2)

    dist = np.sqrt( (x2 - x1)**2 + (y2 - y1)**2 + (z2 - z1)**2)
    dist = np.degrees(dist)*3600.

    return dist

def find_deep_mos_coords(deep_mosaic, bcd):


    pixcord = np.array([[0,0], [0,255], [255, 0], [255, 255]], np.float_)
    w1 = WCS(bcd)
    world = w1.wcs_pix2world(pixcord, 0)
#    print world

    w2 = WCS(deep_mosaic)
    pixnew = w2.wcs_world2pix(world, 0)
#    print pixnew

    xmin = np.min(pixnew[:,0])
    xmax = np.max(pixnew[:,0])
    ymin = np.min(pixnew[:,1])
    ymax = np.max(pixnew[:,1])

    return (xmin, xmax, ymin, ymax)
