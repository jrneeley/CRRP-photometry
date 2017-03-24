#!/usr/bin/env python

import numpy as np
from astropy.io import ascii
import coordinates


def read_optical_catalog(optical_folder, target):

    catalog=optical_folder+target+'.dat'

    print "Reading optical catalog for "+target+"..."

    dtype1 = np.dtype([('id', 'I5'), ('x', 'f8'), ('y', 'f8'), ('ra_h', 'i2'), ('ra_m', 'i2'), ('ra_s', 'f5'), ('dec_d', 'i3'), ('dec_m', 'i2'), ('dec_s', 'f4')])
    data = np.loadtxt(catalog, dtype=dtype1, skiprows=1, usecols=[0,1,2,22,23,24,25,26,27])

    id_num = data['id']
    x = data['x']
    y = data['y']
    ra_h = data['ra_h']
    ra_m = data['ra_m']
    ra_s = data['ra_s']
    dec_d = data['dec_d']
    dec_m = data['dec_m']
    dec_s = data['dec_s']

    ra, dec = coordinates.hms2deg(ra_h, ra_m, ra_s, dec_d, dec_m, dec_s)

    print "Finished reading optical catalog."
    return(id_num, x, y, ra, dec)

def read_optical_catalog_old(target):

    catalog='test_old.dat'

    print "Reading optical catalog for "+target+"..."
    data = ascii.read(catalog, delimiter=' ', data_start=3)

    id_num = np.array(data["col1"])
    x = np.array(data["col2"])
    y = np.array(data["col3"])
    u_mags = np.array(data["col4"])
    b_mags = np.array(data["col6"])
    v_mags = np.array(data["col8"])
    r_mags = np.array(data["col10"])
    i_mags = np.array(data["col11"])
    ra_h = np.array(data["col26"])
    ra_m = np.array(data["col27"])
    ra_s = np.array(data["col28"])
    dec_d = np.array(data["col29"])
    dec_m = np.array(data["col30"])
    dec_s = np.array(data["col31"])

    ra, dec = coordinats.hms2deg(ra_h, ra_m, ra_s, dec_d, dec_m, dec_s)

    print "Finished reading optical catalog."
    return(id_num, x, y, v_mags, ra, dec)
