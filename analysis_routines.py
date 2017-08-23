import numpy as np
import matplotlib.pyplot as mp
import variables
import glob
import read_dao
import re
from astropy.io import fits
import progressbar
from time import sleep


def make_mir_catalog(channels, target, dao_ids, data_dir=''):

    img_folder = data_dir+'data/'
    all_files = np.array([], dtype='S30')
    for channel in channels:
        temp = glob.glob(img_folder+'*.cal')
        all_files = np.append(all_files, temp)
        if len(temp) == 0:
            all_files = np.append(all_files, glob.glob(img_folder+'*.alf'))
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
        bar = progressbar.ProgressBar(maxval=len(file_list), \
            widgets=[progressbar.Bar('=', '[',']'), ' ', progressbar.Percentage()])
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
        first_index = num_channel+prev_num
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
        errs_no_nan = all_errs[~np.isnan(all_errs)]
        times_no_nan = all_times[~np.isnan(all_mags)]
        if len(mags_no_nan) > 0:
            kvar[ind], ivar[ind] = variables.welch_stetson_indices(mags_no_nan, errs_no_nan, times_no_nan)
        else:
            kvar[ind] = np.nan
            ivar[ind] = np.nan


    data_save = np.array(zip(dao_ids, mean_mag[:,0], mean_err[:,0],
        mean_mag[:,1], mean_err[:,1], n_obs[:,0], n_obs[:,1], kvar, ivar),
        dtype=[('c1', 'S8'), ('c2', float), ('c3', float), ('c4', float),
        ('c5', float), ('c6', int), ('c7', int), ('c8', float), ('c9', float)])
    np.savetxt(data_dir+'mir-catalog.txt', data_save,
        fmt='%8s %6.3f %5.3f %6.3f %5.3f %4i %4i %5.3f %7.2f')

def make_cmd_step2(dao_ids, phot_data, data_dir=''):

    x_new = np.zeros(len(dao_ids))
    y_new = np.zeros(len(dao_ids))
    mean_mag = np.zeros(len(dao_ids))
    mean_std = np.zeros(len(dao_ids))

    for ind, star in enumerate(dao_ids):

        x_new[ind] = np.mean(phot_data['x'][ind])
        y_new[ind] = np.mean(phot_data['y'][ind])

        all_mags = phot_data['psf_mag'][ind]
        all_errs = phot_data['psf_err'][ind]

        good_mags = all_mags[~np.isnan(all_mags)]
        good_errs = all_errs[~np.isnan(all_mags)]
        N = len(good_mags)
        if N > 0:
            mean_mag[ind] = variables.robust_weighted_mean(good_mags, good_errs)
            mean_std[ind] = np.std(all_mags[~np.isnan(all_mags)])
        else:
            mean_mag[ind] = np.nan
            mean_std[ind] = np.nan

    data_save = np.array(zip(dao_ids, x_new, y_new, mean_mag, mean_std),
        dtype=[('c1', 'S8'), ('c2', float), ('c3', float), ('c4', float),
        ('c5', float)])
    np.savetxt(data_dir+'test.txt', data_save, fmt='%8s %10.4f %10.4f %6.3f %5.3f')
