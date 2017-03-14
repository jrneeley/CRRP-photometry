#!/usr/bin/env python

import numpy as np
from astropy.io import ascii
import coordinates


def read_optical_catalog(target):

    catalog='/Users/jrneeley/CRRP/OpticalCatalogs/'+target+'.dat'

    print "Reading optical catalog for "+target+"..."
    data = ascii.read(catalog, delimiter=' ', data_start=1)

    id_num = np.array(data["col1"])
    x = np.array(data["col2"])
    y = np.array(data["col3"])
    u_mags = np.array(data["col4"])
    b_mags = np.array(data["col7"])
    v_mags = np.array(data["col10"])
    r_mags = np.array(data["col13"])
    i_mags = np.array(data["col16"])
    ra_h = np.array(data["col23"])
    ra_m = np.array(data["col24"])
    ra_s = np.array(data["col25"])
    dec_d = np.array(data["col26"])
    dec_m = np.array(data["col27"])
    dec_s = np.array(data["col28"])

    ra, dec = coordinates.hms2deg(ra_h, ra_m, ra_s, dec_d, dec_m, dec_s)

    print "Finished reading optical catalog."
    return(id_num, x, y, v_mags, ra, dec)

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
