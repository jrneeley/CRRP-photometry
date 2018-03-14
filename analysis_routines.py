import numpy as np
import matplotlib.pyplot as mp
import variables
import glob
import read_dao
import re
import os
import optical
import coordinates
from astropy.io import fits
#import progressbar
from time import sleep
import config




def merge_catalogs(target, optical_dir, working_dir, cluster_coord=None):

    if cluster_coord == None: opt_data = optical.read_fnl(optical_dir, target)
    else:
        opt_data, rad_dist = optical.read_fnl_w_radial_dist(optical_dir, target,
            cluster_coord[0], cluster_coord_[1])

    dtype1 = np.dtype([('id', 'S8'), ('3.6', float), ('3.6err', float), ('4.5', float),
    ('4.5err', float), ('n3.6', int), ('n4.5', int), ('var1', float), ('var2', float)])
    print 'Reading MIR catalog for '+target+'...'
    mir_data = np.loadtxt(working_dir+'/mir-catalog.txt', dtype=dtype1)

    # change nan values to 9.999 and 99.999
    mir_data['3.6'][mir_data['n3.6'] == 0] = 99.999
    mir_data['3.6err'][mir_data['n3.6'] == 0] = 9.9999
    mir_data['4.5'][mir_data['n4.5'] == 0] = 99.999
    mir_data['4.5err'][mir_data['n4.5'] == 0] = 9.9999

    # add header
    head = ' 5 FILTERS:                 V             B             I             R             U           [3.6]         [4.5]           n    n    n    n    n    n    n    chi  sharp |---------- variability ----------|--- RA  (2000)  Dec ----'


    dtype_comb = np.dtype([('id', int), ('x', float), ('y', float), ('V', float),
        ('Ver', float), ('B', float), ('Ber', float), ('I', float), ('Ier', float),
        ('R', float), ('Rer', float), ('U', float), ('Uer', float), ('3.6', float),
        ('3.6er', float), ('4.5', float), ('4.5er', float), ('nV', int),
        ('nB', int), ('nI', int), ('nR', int), ('nU', int), ('n3.6', int),
        ('n4.5', int), ('chi', float), ('sharp', float), ('var1', float),
        ('var2', float), ('var3', float), ('var4', float), ('var5', float),
        ('ra_h', int), ('ra_m', int), ('ra_s', float),
        ('dec_d', int), ('dec_m', int), ('dec_s', float)])
    dtype_comb2 = np.dtype([('id', int), ('x', float), ('y', float), ('V', float),
        ('Ver', float), ('B', float), ('Ber', float), ('I', float), ('Ier', float),
        ('R', float), ('Rer', float), ('U', float), ('Uer', float), ('3.6', float),
        ('3.6er', float), ('4.5', float), ('4.5er', float), ('nV', int),
        ('nB', int), ('nI', int), ('nR', int), ('nU', int), ('n3.6', int),
        ('n4.5', int), ('chi', float), ('sharp', float), ('var1', float),
        ('var2', float), ('var3', float), ('var4', float), ('var5', float),
        ('ra_h', int), ('ra_m', int), ('ra_s', float),
        ('dec_d', int), ('dec_m', int), ('dec_s', float), ('rad_dist', float)])

    if cluster_coord == None:
        data_save = np.array(zip(opt_data['id'], opt_data['x'], opt_data['y'],
            opt_data['V'], opt_data['Ver'], opt_data['B'], opt_data['Ber'],
            opt_data['I'], opt_data['Ier'], opt_data['R'], opt_data['Rer'],
            opt_data['U'], opt_data['Uer'], mir_data['3.6'], mir_data['3.6err'],
            mir_data['4.5'], mir_data['4.5err'], opt_data['nV'], opt_data['nB'],
            opt_data['nI'], opt_data['nR'], opt_data['nU'], mir_data['n3.6'],
            mir_data['n4.5'], opt_data['chi'], opt_data['sharp'], opt_data['var1'],
            opt_data['var2'], opt_data['var3'], opt_data['var4'], opt_data['var5'],
            opt_data['ra_h'], opt_data['ra_m'], opt_data['ra_s'], opt_data['dec_d'],
            opt_data['dec_m'], opt_data['dec_s']), dtype=dtype_comb)
    else:
        data_save = np.array(zip(opt_data['id'], opt_data['x'], opt_data['y'],
            opt_data['V'], opt_data['Ver'], opt_data['B'], opt_data['Ber'],
            opt_data['I'], opt_data['Ier'], opt_data['R'], opt_data['Rer'],
            opt_data['U'], opt_data['Uer'], mir_data['3.6'], mir_data['3.6err'],
            mir_data['4.5'], mir_data['4.5err'], opt_data['nV'], opt_data['nB'],
            opt_data['nI'], opt_data['nR'], opt_data['nU'], mir_data['n3.6'],
            mir_data['n4.5'], opt_data['chi'], opt_data['sharp'], opt_data['var1'],
            opt_data['var2'], opt_data['var3'], opt_data['var4'], opt_data['var5'],
            opt_data['ra_h'], opt_data['ra_m'], opt_data['ra_s'], opt_data['dec_d'],
            opt_data['dec_m'], opt_data['dec_s'], rad_dist), dtype=dtype_comb2)

    np.savetxt(working_dir+'merged-catalog.txt', data_save,
        fmt='%8i %8.2f %8.2f %6.3f %6.4f %6.3f %6.4f %6.3f %6.4f %6.3f %6.4f %6.3f %6.4f %6.3f %6.4f %6.3f %6.4f %4i %4i %4i %4i %4i %4i %4i %6.3f %6.3f %6.3f %6.3f %6.3f %6.1f %6.3f %3i %02i %05.2f %+03i %02i %04.1f',
        header=head)

def merge_opt_deep_catalogs(target, cluster_coord=None):

    optical_dir = config.optical_dir
    data_dir = config.top_dir+target

    if cluster_coord == None: opt_data = optical.read_fnl(target)
    else:
        opt_data, rad_dist = optical.read_fnl_w_radial_dist(target,
            cluster_coord[0], cluster_coord[1])
    include_nir = 0
    nir_file = optical_dir+target+'ir.fnl'
    if os.path.isfile(nir_file):
        nir_data = optical.read_nir_fnl(target)
        include_nir = 1
    dtype1 = np.dtype([('id', int), ('3.6', float), ('3.6err', float)])
    dtype2 = np.dtype([('id', int), ('4.5', float), ('4.5err', float)])
    print 'Reading MIR catalog for '+target+'...'
    data3p6 = np.loadtxt(data_dir+'/DeepMosaic/'+target+'_I1_deep_dn.cal', dtype=dtype1, usecols=(0,3,4), skiprows=3)
    data4p5 = np.loadtxt(data_dir+'/DeepMosaic/'+target+'_I2_deep_dn.cal', dtype=dtype2, usecols=(0,3,4), skiprows=3)

    mir_data = np.zeros(len(opt_data['id']), dtype=([('3.6', float), ('3.6err', float),
        ('4.5', float), ('4.5err', float), ('n3.6', int), ('n4.5', int)]))

    mir_data['n3.6'][:] = 1
    mir_data['n4.5'][:] = 1
    for ind, star in enumerate(opt_data['id']):
        match3p6 = np.argwhere(data3p6['id'] == star)
        match4p5 = np.argwhere(data4p5['id'] == star)
        if len(match3p6) > 0:
            mir_data['3.6'][ind] = data3p6['3.6'][data3p6['id'] == star]
            mir_data['3.6err'][ind] = data3p6['3.6err'][data3p6['id'] == star]
        else:
            mir_data['3.6'][ind] = 99.999
            mir_data['3.6err'][ind] = 9.9999
        if len(match4p5) > 0:
            mir_data['4.5'][ind] = data4p5['4.5'][data4p5['id'] == star]
            mir_data['4.5err'][ind] = data4p5['4.5err'][data4p5['id'] == star]
        else:
            mir_data['4.5'][ind] = 99.999
            mir_data['4.5err'][ind] = 9.9999

    if include_nir == 0:
        # add header
        head = ' 7 FILTERS:                 U             B             V             R             I           [3.6]         [4.5]           n    n    n    n    n    n    n    chi  sharp |---------- variability ----------|--- RA  (2000)  Dec ----'

        dtype_comb = np.dtype([('id', int), ('x', float), ('y', float), ('U', float),
            ('Uer', float), ('B', float), ('Ber', float), ('V', float), ('Ver', float),
            ('R', float), ('Rer', float), ('I', float), ('Ier', float), ('3.6', float),
            ('3.6er', float), ('4.5', float), ('4.5er', float), ('nU', int),
            ('nB', int), ('nV', int), ('nR', int), ('nI', int), ('n3.6', int),
            ('n4.5', int), ('chi', float), ('sharp', float), ('var1', float),
            ('var2', float), ('var3', float), ('var4', float), ('var5', float),
            ('ra_h', int), ('ra_m', int), ('ra_s', float),
            ('dec_d', int), ('dec_m', int), ('dec_s', float), ('rad_dist', float)])

        data_save = np.array(zip(opt_data['id'], opt_data['x'], opt_data['y'],
            opt_data['U'], opt_data['Uer'], opt_data['B'], opt_data['Ber'],
            opt_data['V'], opt_data['Ver'], opt_data['R'], opt_data['Rer'],
            opt_data['I'], opt_data['Ier'], mir_data['3.6'], mir_data['3.6err'],
            mir_data['4.5'], mir_data['4.5err'], opt_data['nU'], opt_data['nB'],
            opt_data['nV'], opt_data['nR'], opt_data['nI'], mir_data['n3.6'],
            mir_data['n4.5'], opt_data['chi'], opt_data['sharp'], opt_data['var1'],
            opt_data['var2'], opt_data['var3'], opt_data['var4'], opt_data['var5'],
            opt_data['ra_h'], opt_data['ra_m'], opt_data['ra_s'], opt_data['dec_d'],
            opt_data['dec_m'], opt_data['dec_s'], rad_dist), dtype=dtype_comb)

        np.savetxt(data_dir+'/merged-deep-catalog.txt', data_save,
            fmt='%8i %8.2f %8.2f %6.3f %6.4f %6.3f %6.4f %6.3f %6.4f %6.3f %6.4f %6.3f %6.4f %6.3f %6.4f %6.3f %6.4f %4i %4i %4i %4i %4i %4i %4i %6.3f %6.3f %6.3f %6.3f %6.3f %6.1f %6.3f %3i %02i %05.2f %+03i %02i %04.1f %0.2f',
            header=head)
    if include_nir == 1:

        head = '10 FILTERS:                 U             B             V             R             I             J             H             K           [3.6]         [4.5]           n    n    n    n    n    n    n    chi  sharp |---------- variability ----------|--- RA  (2000)  Dec ----'

        dtype_comb = np.dtype([('id', int), ('x', float), ('y', float), ('U', float),
            ('Uer', float), ('B', float), ('Ber', float), ('V', float), ('Ver', float),
            ('R', float), ('Rer', float), ('I', float), ('Ier', float), ('J', float),
            ('Jer', float), ('H', float), ('Her', float), ('K', float), ('Ker', float),
            ('3.6', float), ('3.6er', float), ('4.5', float), ('4.5er', float),
            ('nU', int), ('nB', int), ('nV', int), ('nR', int), ('nI', int),
            ('nJ', int), ('nH', int), ('nK', int), ('n3.6', int), ('n4.5', int),
            ('chi', float), ('sharp', float), ('var1', float),
            ('var2', float), ('var3', float), ('var4', float), ('var5', float),
            ('ra_h', int), ('ra_m', int), ('ra_s', float),
            ('dec_d', int), ('dec_m', int), ('dec_s', float), ('rad_dist', float)])


        data_save = np.array(zip(opt_data['id'], opt_data['x'], opt_data['y'],
            opt_data['U'], opt_data['Uer'], opt_data['B'], opt_data['Ber'],
            opt_data['V'], opt_data['Ver'], opt_data['R'], opt_data['Rer'],
            opt_data['I'], opt_data['Ier'], nir_data['J'], nir_data['Jer'],
            nir_data['H'], nir_data['Her'], nir_data['K'], nir_data['Ker'],
            mir_data['3.6'], mir_data['3.6err'], mir_data['4.5'], mir_data['4.5err'],
            opt_data['nU'], opt_data['nB'], opt_data['nV'], opt_data['nR'], opt_data['nI'],
            nir_data['nJ'], nir_data['nH'], nir_data['nK'], mir_data['n3.6'],
            mir_data['n4.5'], opt_data['chi'], opt_data['sharp'], opt_data['var1'],
            opt_data['var2'], opt_data['var3'], opt_data['var4'], opt_data['var5'],
            opt_data['ra_h'], opt_data['ra_m'], opt_data['ra_s'], opt_data['dec_d'],
            opt_data['dec_m'], opt_data['dec_s'], rad_dist), dtype=dtype_comb)

        np.savetxt(data_dir+'/merged-deep-catalog.txt', data_save,
            fmt='%8i %8.2f %8.2f %6.3f %6.4f %6.3f %6.4f %6.3f %6.4f %6.3f %6.4f %6.3f %6.4f %6.3f %6.4f %6.3f %6.4f %6.3f %6.4f %6.3f %6.4f %6.3f %6.4f %4i %4i %4i %4i %4i %4i %4i %4i %4i %4i %6.3f %6.3f %6.3f %6.3f %6.3f %6.1f %6.3f %3i %02i %05.2f %+03i %02i %04.1f %0.2f',
            header=head)


def read_merged_catalog(target):#, center_ra, center_dec):

    data_dir = config.top_dir+target
    dtype_comb = np.dtype([('id', int), ('x', float), ('y', float), ('U', float),
        ('Uer', float), ('B', float), ('Ber', float), ('V', float), ('Ver', float),
        ('R', float), ('Rer', float), ('I', float), ('Ier', float), ('3.6', float),
        ('3.6er', float), ('4.5', float), ('4.5er', float), ('nU', int),
        ('nB', int), ('nV', int), ('nR', int), ('nI', int), ('n3.6', int),
        ('n4.5', int), ('chi', float), ('sharp', float), ('var1', float),
        ('var2', float), ('var3', float), ('var4', float), ('var5', float),
        ('ra_h', int), ('ra_m', int), ('ra_s', float),
        ('dec_d', int), ('dec_m', int), ('dec_s', float), ('rad_dist', float)])

    data = np.loadtxt(data_dir+'/merged-deep-catalog.txt', dtype=dtype_comb)

    return data

def read_merged_catalog2(target):#, center_ra, center_dec):

    data_dir = config.top_dir+target
    dtype_comb = np.dtype([('id', int), ('x', float), ('y', float), ('U', float),
        ('Uer', float), ('B', float), ('Ber', float), ('V', float), ('Ver', float),
        ('R', float), ('Rer', float), ('I', float), ('Ier', float), ('J', float),
        ('Jer', float), ('H', float), ('Her', float), ('K', float), ('Ker', float),
        ('3.6', float), ('3.6er', float), ('4.5', float), ('4.5er', float),
        ('nU', int), ('nB', int), ('nV', int), ('nR', int), ('nI', int),
        ('nJ', int), ('nH', int), ('nK', int), ('n3.6', int), ('n4.5', int),
        ('chi', float), ('sharp', float), ('var1', float),
        ('var2', float), ('var3', float), ('var4', float), ('var5', float),
        ('ra_h', int), ('ra_m', int), ('ra_s', float),
        ('dec_d', int), ('dec_m', int), ('dec_s', float), ('rad_dist', float)])

    data = np.loadtxt(data_dir+'/merged-deep-catalog.txt', dtype=dtype_comb)

    return data
