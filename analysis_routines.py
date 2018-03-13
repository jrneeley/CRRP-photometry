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


def make_mir_catalog(channels, target, dao_ids, data_dir=''):

    img_folder = data_dir+'data/'
    all_files = np.array([], dtype='S30')
    #for channel in channels:
    temp = glob.glob(img_folder+'*.cal')

    if len(temp) == 0:
        all_files = np.append(all_files, glob.glob(img_folder+'*.alf'))
    else:
        all_files = np.append(all_files, temp)

    master_data = np.zeros(len(dao_ids), dtype=[('filter', 'S2', len(all_files)),
        ('id', 'S10'), ('mjd', float, len(all_files)),
        ('psf_mag', float, len(all_files)), ('psf_err', float, len(all_files))])

    prev_num = 0
    for num_channel, channel in enumerate(channels):
        file_list = glob.glob(img_folder+channel+'*.cal')
        if len(file_list) == 0:
            file_list = glob.glob(img_folder+channel+'*.alf')
        phot_data = np.zeros(len(dao_ids), dtype=[('filter', 'S2', len(file_list)),
            ('id', 'S10'), ('mjd', float, len(file_list)),
            ('psf_mag', float, len(file_list)), ('psf_err', float, len(file_list))])

        phot_data['id'] = dao_ids
#        bar = progressbar.ProgressBar(maxval=len(file_list), \
#            widgets=[progressbar.Bar('=', '[',']'), ' ', progressbar.Percentage()])
        bar.start()
        for ind in range(0,len(file_list)):

            alf_id, x, y, alf_mag, alf_err = read_dao.read_alf(file_list[ind])
            fits_file = re.sub(".cal",".fits", file_list[ind])
            if fits_file == file_list[ind]:
                fits_file = re.sub(".alf", ".fits", file_list[ind])
            hdulist = fits.open(fits_file)
            prihdr = hdulist[0].header
            mjd = prihdr['mjd_obs']

            for ind2 in range(0,len(dao_ids)):
                alf_match = np.argwhere(alf_id == dao_ids[ind2])
                phot_data['filter'][ind2, ind] = channel
                phot_data['mjd'][ind2,ind] = mjd
                if len(alf_match):
                    phot_data['psf_mag'][ind2, ind] = alf_mag[alf_match]
                    phot_data['psf_err'][ind2, ind] = alf_err[alf_match]
                else:
                    phot_data['psf_mag'][ind2,ind] = float('NaN')
                    phot_data['psf_err'][ind2,ind] = float('NaN')
            bar.update(ind+1)
            sleep(0.1)
        first_index = prev_num
        second_index = first_index + len(file_list)
        master_data['filter'][:,first_index:second_index] = phot_data['filter']
        master_data['id'] = phot_data['id']
        master_data['mjd'][:,first_index:second_index] = phot_data['mjd']
        master_data['psf_mag'][:,first_index:second_index] = phot_data['psf_mag']
        master_data['psf_err'][:,first_index:second_index] = phot_data['psf_err']
        prev_num = len(file_list)


    mean_mag = np.zeros((len(dao_ids), len(channels)))
    mean_err = np.zeros((len(dao_ids), len(channels)))
    n_obs = np.zeros((len(dao_ids), len(channels)))
    kvar = np.zeros(len(dao_ids))
    ivar = np.zeros(len(dao_ids))
    for ind, star in enumerate(dao_ids):

        all_mags = master_data['psf_mag'][ind]
        all_errs = master_data['psf_err'][ind]
        all_filts = master_data['filter'][ind]
        all_times = master_data['mjd'][ind]

        for num_ch, channel in enumerate(channels):
            mags_in_ch = all_mags[all_filts == channel]
            errs_in_ch = all_errs[all_filts == channel]
            times_in_ch = all_times[all_filts == channel]
            good_mags = mags_in_ch[~np.isnan(mags_in_ch)]
            good_errs = errs_in_ch[~np.isnan(mags_in_ch)]
            good_times = times_in_ch[~np.isnan(mags_in_ch)]
            N = len(good_mags)
            n_obs[ind, num_ch] = N
            if N > 0:
                mean_mag[ind, num_ch], mean_err[ind, num_ch] = variables.robust_weighted_mean(good_mags, good_errs)
        #    kvar[ind], ivar[ind] = variables.welch_stetson_indices(good_mags, good_errs, good_times)
            else:
                mean_mag[ind] = np.nan
                mean_err[ind] = np.nan

        mags_no_nan = all_mags[~np.isnan(all_mags)]
        errs_no_nan = all_errs[~np.isnan(all_mags)]
        times_no_nan = all_times[~np.isnan(all_mags)]
        if len(mags_no_nan) > 0:
        #    kvar[ind], ivar[ind] = variables.welch_stetson_indices(mags_no_nan, errs_no_nan, times_no_nan)
            kvar[ind] = 0.0
            ivar[ind] = 0.0
        else:
            kvar[ind] = np.nan
            ivar[ind] = np.nan


    data_save = np.array(zip(dao_ids, mean_mag[:,0], mean_err[:,0],
        mean_mag[:,1], mean_err[:,1], n_obs[:,0], n_obs[:,1], kvar, ivar),
        dtype=[('c1', 'S8'), ('c2', float), ('c3', float), ('c4', float),
        ('c5', float), ('c6', int), ('c7', int), ('c8', float), ('c9', float)])
    np.savetxt(data_dir+'mir-catalog.txt', data_save,
        fmt='%8s %6.3f %5.3f %6.3f %5.3f %4i %4i %5.3f %7.2f')


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

def merge_opt_deep_catalogs(target, optical_dir, working_dir, cluster_coord=None):

    if cluster_coord == None: opt_data = optical.read_fnl(optical_dir, target)
    else:
        opt_data, rad_dist = optical.read_fnl_w_radial_dist(optical_dir, target,
            cluster_coord[0], cluster_coord[1])
    include_nir = 0
    nir_file = optical_dir+target+'ir.fnl'
    if os.path.isfile(nir_file):
        nir_data = optical.read_nir_fnl(optical_dir, target)
        include_nir = 1
    dtype1 = np.dtype([('id', int), ('3.6', float), ('3.6err', float)])
    dtype2 = np.dtype([('id', int), ('4.5', float), ('4.5err', float)])
    print 'Reading MIR catalog for '+target+'...'
    data3p6 = np.loadtxt(working_dir+'/DeepMosaic/'+target+'_I1_deep_dn.cal', dtype=dtype1, usecols=(0,3,4), skiprows=3)
    data4p5 = np.loadtxt(working_dir+'/DeepMosaic/'+target+'_I2_deep_dn.cal', dtype=dtype2, usecols=(0,3,4), skiprows=3)

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
        head = ' 7 FILTERS:                 V             B             I             R             U           [3.6]         [4.5]           n    n    n    n    n    n    n    chi  sharp |---------- variability ----------|--- RA  (2000)  Dec ----'

        dtype_comb = np.dtype([('id', int), ('x', float), ('y', float), ('V', float),
            ('Ver', float), ('B', float), ('Ber', float), ('I', float), ('Ier', float),
            ('R', float), ('Rer', float), ('U', float), ('Uer', float), ('3.6', float),
            ('3.6er', float), ('4.5', float), ('4.5er', float), ('nV', int),
            ('nB', int), ('nI', int), ('nR', int), ('nU', int), ('n3.6', int),
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

        np.savetxt(working_dir+'merged-deep-catalog.txt', data_save,
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

        np.savetxt(working_dir+'merged-deep-catalog.txt', data_save,
            fmt='%8i %8.2f %8.2f %6.3f %6.4f %6.3f %6.4f %6.3f %6.4f %6.3f %6.4f %6.3f %6.4f %6.3f %6.4f %6.3f %6.4f %6.3f %6.4f %6.3f %6.4f %6.3f %6.4f %4i %4i %4i %4i %4i %4i %4i %4i %4i %4i %6.3f %6.3f %6.3f %6.3f %6.3f %6.1f %6.3f %3i %02i %05.2f %+03i %02i %04.1f %0.2f',
            header=head)


def read_merged_catalog(data_dir):#, center_ra, center_dec):

    dtype_comb = np.dtype([('id', int), ('x', float), ('y', float), ('V', float),
        ('Ver', float), ('B', float), ('Ber', float), ('I', float), ('Ier', float),
        ('R', float), ('Rer', float), ('U', float), ('Uer', float), ('3.6', float),
        ('3.6er', float), ('4.5', float), ('4.5er', float), ('nV', int),
        ('nB', int), ('nI', int), ('nR', int), ('nU', int), ('n3.6', int),
        ('n4.5', int), ('chi', float), ('sharp', float), ('var1', float),
        ('var2', float), ('var3', float), ('var4', float), ('var5', float),
        ('ra_h', int), ('ra_m', int), ('ra_s', float),
        ('dec_d', int), ('dec_m', int), ('dec_s', float), ('rad_dist', float)])

    data = np.loadtxt(data_dir+'/merged-deep-catalog.txt', dtype=dtype_comb)

    return data

def read_merged_catalog2(data_dir):#, center_ra, center_dec):

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

    data = np.loadtxt(data_dir+'merged-deep-catalog.txt', dtype=dtype_comb)

    return data
