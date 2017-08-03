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

def read_fnl_w_radial_dist(optical_folder, target, center_ra, center_dec):

    catalog=optical_folder+target+'.fnl'

    print "Reading optical catalog for "+target+"..."

    dtype1 = np.dtype([('id', int), ('x', float), ('y', float), ('V', float),
        ('Ver', float), ('B', float), ('Ber', float), ('I', float), ('Ier', float),
        ('R', float), ('Rer', float), ('U', float), ('Uer', float), ('nV', int),
        ('nB', int), ('nI', int), ('nR', int), ('nU', int), ('chi', float),
        ('sharp', float), ('var1', float), ('var2', float), ('var3', float),
        ('var4', float), ('var5', float), ('ra_h', int), ('ra_m', int),
        ('ra_s', float), ('dec_d', int), ('dec_m', int), ('dec_s', float)])
    data = np.loadtxt(catalog, dtype=dtype1, skiprows=3, usecols=range(31))

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

def read_fnl(optical_folder, target):

    catalog=optical_folder+target+'.fnl'

    print "Reading optical catalog for "+target+"..."

    dtype1 = np.dtype([('id', int), ('x', float), ('y', float), ('V', float),
        ('Ver', float), ('B', float), ('Ber', float), ('I', float), ('Ier', float),
        ('R', float), ('Rer', float), ('U', float), ('Uer', float), ('nV', int),
        ('nB', int), ('nI', int), ('nR', int), ('nU', int), ('chi', float),
        ('sharp', float), ('var1', float), ('var2', float), ('var3', float),
        ('var4', float), ('var5', float), ('ra_h', int), ('ra_m', int),
        ('ra_s', float), ('dec_d', int), ('dec_m', int), ('dec_s', float)])
    data = np.loadtxt(catalog, dtype=dtype1, skiprows=3, usecols=range(31))

    print "Finished reading optical catalog."
    return data

def find_variables_fnl(optical_folder, target, center_ra, center_dec, folder=''):

    catalog=optical_folder+target+'.fnl'

    print "Reading optical catalog for "+target+"..."

    f = open(catalog, 'r')
    lines = f.readlines()
    master_id=[]
    variable_id=[]

    for line in lines:
        temp = line.split()
        if len(temp) == 32:
            master_id.append(int(temp[0]))
            variable_id.append(temp[-1])

    data, dist = read_fnl_w_radial_dist(optical_folder, target, center_ra, center_dec)
    mags = np.zeros(len(master_id))
    errs = np.zeros(len(master_id))
    rad_dist = np.zeros(len(master_id))
    var1 = np.zeros(len(master_id))
    var2 = np.zeros(len(master_id))
    var3 = np.zeros(len(master_id))
    var4 = np.zeros(len(master_id))
    var5 = np.zeros(len(master_id))
    ra = np.zeros(len(master_id), dtype='S13')
    dec = np.zeros(len(master_id), dtype='S13')
    for ind, star in enumerate(master_id):

        mags[ind] = data['V'][data['id'] == star]
        errs[ind] = data['Ver'][data['id'] == star]
        rad_dist[ind] = dist[data['id'] == star]
        var1[ind] = data['var1'][data['id'] == star]
        var2[ind] = data['var2'][data['id'] == star]
        var3[ind] = data['var3'][data['id'] == star]
        var4[ind] = data['var4'][data['id'] == star]
        var5[ind] = data['var5'][data['id'] == star]
        ra_h = data['ra_h'][data['id'] == star]
        ra_m = data['ra_m'][data['id'] == star]
        ra_s = data['ra_s'][data['id'] == star]
        dec_d = data['dec_d'][data['id'] == star]
        dec_m = data['dec_m'][data['id'] == star]
        dec_s = data['dec_s'][data['id'] == star]
    #    ra_d, dec_de = coordinates.hms2deg(ra_h, ra_m, ra_s, dec_d, dec_m, dec_s)
    #    ra_deg[ind] = ra_d[0]
    #    dec_deg[ind] = dec_de[0]
        ra[ind] = str(ra_h[0])+':'+str(ra_m[0])+':'+str(ra_s[0])
        dec[ind] = str(dec_d[0])+':'+str(dec_m[0])+':'+str(dec_s[0])
    data_save = np.array(zip(variable_id, master_id, mags, errs, var1, var2, var3,
        var4, var5, ra, dec, rad_dist),
        dtype=[('c1', 'S10'), ('c2', int), ('c3', float), ('c4', float), ('c5', float),
        ('c6', float), ('c7', float), ('c8', float), ('c9', float), ('c10', 'S13'),
        ('c11', 'S13'), ('c12', float)])
    np.savetxt(folder+'PeterIDs.txt', data_save, comments='',
        fmt='%10s %8i %7.3f %6.3f %7.3f %7.3f %7.3f %6.1f %7.4f %13s %13s %10.3f')



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

def compile_datasets(target, old=0, returnColors=True):

    lcvs = glob.glob(target+'lcvs/optical/*.lcv')
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
    colors = np.zeros(len(dataset_names), dtype='S25')
    print '\n\nDatasets:\n'
    for ind, dataset in enumerate(dataset_names):
        jds = all_jds[datasets_prefix == dataset]
        t_min = Time(jds.min(), format='jd')
        t_max = Time(jds.max(), format='jd')
        t_min.out_subfmt = 'date'
        t_max.out_subfmt = 'date'
        colors[ind] = plotting_utilities.get_color(ind)
        print '%10s %6i %s %s %s' % (dataset, dataset_counts[ind], t_min.iso, t_max.iso, plotting_utilities.get_color(ind))

    if returnColors == True:
        return dataset_names, colors
    if returnColors == False:
        return dataset_names
