#!/usr/bin/env python

import sys
#import glob
import numpy as np
#from astropy.io import ascii
from astropy.wcs import WCS
import optical


def find_coord_window(img_list):

# initialize boundary conditions
    bxmin = 999
    bxmax = 0
    bymin = 999
    bymax = -999

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
            bxmin = xmin
        if (xmax > bxmax):
            bxmax = xmax
        if (ymin < bymin):
            bymin = ymin
        if (ymax > bymax):
            bymax = ymax

    return(bxmin, bxmax, bymin, bymax)

# Convert celestial coordinates to pixels using optical catalog
def radec2pix(target, xmin, xmax, ymin, ymax, x, y, ra, dec):

    index_xmin = np.searchsorted(ra, xmin)
    index_xmax = np.searchsorted(ra, xmax)
    dec2 = np.sort(dec)
    y2 = y[np.argsort(dec)]
    index_ymin = np.searchsorted(dec2, ymin)
    index_ymax = np.searchsorted(dec2, ymax)

    cat_xmin = x[index_xmin]
    cat_xmax = x[index_xmax]
    cat_ymin = y2[index_ymin]
    cat_ymax = y2[index_ymax]

    return(cat_xmin, cat_xmax, cat_ymin, cat_ymax)

# Convert coordinates in sexagesimal to decimal units
def hms2deg(ra_h, ra_m, ra_s, dec_d, dec_m, dec_s):

    if len(dec_d) == 1 :
        if dec_d < 0:
            dec = dec_d - dec_m/60. - dec_s/3600.
        if dec_d >= 0 :
            dec = dec_d + dec_m/60. + dec_s/3600.
    if len(dec_d) > 1:
        dec = []
        for ind in range(0,len(ra_h)):
            if dec_d[ind] < 0:
                dec.append(dec_d[ind] - dec_m[ind]/60. - dec_s[ind]/3600.)
            if dec_d[ind] > 0:
                dec.append(dec_d[ind] + dec_m[ind]/60. + dec_s[ind]/3600.)
    if len(ra_h) > 1:
        ra = []
        for ind in range(0,len(ra_h)):
            ra.append(15.*(ra_h[ind]+ra_m[ind]/60.+ra_s[ind]/3600.))
    if len(ra_h) == 1:
        ra = 15*(ra_h+ra_m/60.+ra_s/3600.)

    return ra, dec

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
