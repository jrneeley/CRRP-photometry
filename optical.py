#!/usr/bin/env python

import numpy as np
from astropy.io import ascii
import coordinates
import shutil
import lightcurves
import glob
from astropy.time import Time
import plotting_utilities
import dao
import config

def read_optical_catalog(optical_dir, target):

    catalog=optical_dir+target+'.dat'

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

def read_optical_fnl(target):

    catalog = config.optical_dir+target+'.fnl'

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

def read_fnl_w_radial_dist(target, center_ra, center_dec):

    catalog = config.optical_dir+target+'.fnl'

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

def read_fnl(target):

    catalog = config.optical_dir+target+'.fnl'

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

def read_nir_fnl(target):

    catalog = config.optical_dir+target+'ir.fnl'

    print "Reading NIR catalog for "+target+"..."

    dtype1 = np.dtype([('id', int), ('x', float), ('y', float), ('J', float),
        ('Jer', float), ('H', float), ('Her', float), ('K', float), ('Ker', float),
        ('nJ', int), ('nH', int), ('nK', int), ('chi', float),
        ('sharp', float), ('var1', float), ('var2', float), ('var3', float),
        ('var4', float), ('var5', float), ('ra_h', int), ('ra_m', int),
        ('ra_s', float), ('dec_d', int), ('dec_m', int), ('dec_s', float)])
    data = np.loadtxt(catalog, dtype=dtype1, skiprows=3, usecols=range(25))

    print "Finished reading NIR catalog."
    return data

def find_variables_fnl(target, center_ra, center_dec):
    data_dir = config.top_dir+target
    catalog = config.optical_dir+target+'.vary'

    f = open(catalog, 'r')
    lines = f.readlines()
    master_id=[]
    variable_id=[]

    for line in lines:
        temp = line.split()
        if len(temp) == 32:
            master_id.append(int(temp[0]))
            variable_id.append(temp[-1])
    f.close()
    print 'Found {} candidate variables.'.format(len(variable_id))

    data, dist = read_fnl_w_radial_dist(target, center_ra, center_dec)
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
        match = data['id'] == star

        mags[ind] = data['V'][match]
        errs[ind] = data['Ver'][match]
        rad_dist[ind] = dist[match]
        var1[ind] = data['var1'][match]
        var2[ind] = data['var2'][match]
        var3[ind] = data['var3'][match]
        var4[ind] = data['var4'][match]
        var5[ind] = data['var5'][match]
        ra_h = data['ra_h'][match]
        ra_m = data['ra_m'][match]
        ra_s = data['ra_s'][match]
        dec_d = data['dec_d'][match]
        dec_m = data['dec_m'][match]
        dec_s = data['dec_s'][match]
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
    np.savetxt(data_dir+'/candidate-vars.txt', data_save, comments='',
        fmt='%10s %8i %7.3f %6.3f %7.3f %7.3f %7.3f %6.1f %7.4f %13s %13s %10.3f')

def find_variables(target, center_ra, center_dec):
    data_dir = config.top_dir+target
    catalog = config.optical_dir+target+'.vary'

    f = open(catalog, 'r')
    f2 = open('candidate-variables.txt', 'w')
    lines = f.readlines()
    master_id=[]
    variable_id=[]

    for line in lines:
        temp = line.split()
        if len(temp) == 32:
            master_id.append(int(temp[0]))
            variable_id.append(temp[-1])
            f2.write(line)
    f.close()
    print 'Found {} candidate variables.'.format(len(variable_id))

def compile_datasets(target, old=0, returnColors=True):

    lcvs = glob.glob(target+'lcvs/optical/*.lcv')
#    all_datasets = np.zeros(1, dtype='S30')
#    all_jds = np.zeros(1, dtype=float)
    num_images = 0
    for lcv in lcvs:

        U, B, V, R, I = lightcurves.read_optical_lcv(lcv, old=old)
        if lcv == lcvs[0]:
            all_datasets = np.array(U[3], dtype='S35')
            all_jds = np.array(U[2], dtype=float)
            num_images += len(U[0])+len(B[0])+len(V[0])+len(R[0])+len(I[0])
        else:
            all_datasets = np.append(all_datasets, U[3])
            all_jds = np.append(all_jds, U[2])
            num_images += len(U[0])+len(B[0])+len(V[0])+len(R[0])+len(I[0])
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
    print 'Total number of images: {} '.format(num_images)
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
        print '%3i %10s %6i %s %s %s' % (ind+1, dataset, dataset_counts[ind], t_min.iso, t_max.iso, plotting_utilities.get_color(ind))

    if returnColors == True:
        return dataset_names, colors
    if returnColors == False:
        return dataset_names
