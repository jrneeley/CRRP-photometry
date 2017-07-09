#!/usr/bin/env python

import numpy as np
from astropy.io import ascii
import coordinates
import shutil
import lightcurves
import glob
from astropy.time import Time
import plotting_utilities

def read_optical_catalog(optical_folder, target):

    catalog=optical_folder+target+'.dat'

    print "Reading optical catalog for "+target+"..."

    dtype1 = np.dtype([('id', int), ('x', float), ('y', float), ('ra_h', int),
        ('ra_m', int), ('ra_s', float), ('dec_d', int), ('dec_m', int),
        ('dec_s', float)])
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

def read_optical_fnl(optical_folder, target):

    catalog=optical_folder+target+'.fnl'

    print "Reading optical catalog for "+target+"..."

    dtype1 = np.dtype([('id', int), ('x', float), ('y', float), ('ra_h', int),
        ('ra_m', int), ('ra_s', float), ('dec_d', int), ('dec_m', int),
        ('dec_s', float)])
    data = np.loadtxt(catalog, dtype=dtype1, skiprows=3, usecols=[0,1,2,25,26,27,28,29,30])

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

def read_fnl_for_cmd(optical_folder, target, center_ra, center_dec):

    catalog=optical_folder+target+'.fnl'

    print "Reading optical catalog for "+target+"..."

    dtype1 = np.dtype([('id', int), ('V', float), ('B', float), ('I', float),
        ('R', float), ('U', float), ('ra_h', int), ('ra_m', int), ('ra_s', float),
        ('dec_d', int), ('dec_m', int), ('dec_s', float)])
    data = np.loadtxt(catalog, dtype=dtype1, skiprows=3, usecols=[0,3,5,7,9,11,25,26,27,28,29,30])

    ra_h = data['ra_h']
    ra_m = data['ra_m']
    ra_s = data['ra_s']
    dec_d = data['dec_d']
    dec_m = data['dec_m']
    dec_s = data['dec_s']

    ra, dec = coordinates.hms2deg(ra_h, ra_m, ra_s, dec_d, dec_m, dec_s)

    dist = coordinates.radial_dist(ra, dec, center_ra, center_dec)

    print "Finished reading optical catalog."
    return(data, dist)

def find_variables_fnl(optical_folder, target):

    catalog=optical_folder+target+'.fnl'

    print "Reading optical catalog for "+target+"..."

    f = open(catalog, 'r')
    lines = f.readlines()
    master_id=[]
    variable_id=[]

    for line in lines:
        temp = line.split()
        if len(temp) == 32:
            master_id.append(temp[0])
            variable_id.append(temp[-1])

    data_save = np.array(zip(variable_id, master_id), dtype=[('c1', 'S10'), ('c2', int)])
    np.savetxt('PeterIDs.txt', data_save, comments='', fmt='%10s %8i')



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

def compile_datasets(target, old=0, folder=''):

    lcvs = glob.glob(folder+'lcvs/optical/*.lcv')
#    all_datasets = np.zeros(1, dtype='S30')
#    all_jds = np.zeros(1, dtype=float)
    for lcv in lcvs:

        U, B, V, R, I = lightcurves.read_optical_lcv(lcv, old=old)
        if lcv == lcvs[0]:
            all_datasets = np.array(U[3], dtype='S30')
            all_jds = np.array(U[2], dtype=float)
        else:
            all_datasets = np.append(all_datasets, U[3])
            all_jds = np.append(all_jds, U[2])
        all_datasets = np.append(all_datasets, B[3])
        all_datasets = np.append(all_datasets, V[3])
        all_datasets = np.append(all_datasets, R[3])
        all_datasets = np.append(all_datasets, I[3])
        all_jds = np.append(all_jds, B[2])
        all_jds = np.append(all_jds, V[2])
        all_jds = np.append(all_jds, R[2])
        all_jds = np.append(all_jds, I[2])
    all_jds = all_jds + 2400000.5
    datasets_prefix = np.zeros(len(all_datasets), dtype='S30')
    for ind, string in enumerate(all_datasets):
        datasets_prefix[ind] = string.split(':')[0]
    unique, counts = np.unique(datasets_prefix, return_counts=True)

    dataset_names = unique[np.argsort(counts)[::-1]]
    dataset_counts = counts[np.argsort(counts)[::-1]]
    print '\n\nDatasets:\n'
    for ind, dataset in enumerate(dataset_names):
        jds = all_jds[datasets_prefix == dataset]
        t_min = Time(jds.min(), format='jd')
        t_max = Time(jds.max(), format='jd')
        t_min.out_subfmt = 'date'
        t_max.out_subfmt = 'date'
        print '%10s %6i %s %s %s' % (dataset, dataset_counts[ind], t_min.iso, t_max.iso, plotting_utilities.get_color(ind))

    return unique[np.argsort(counts)[::-1]]
