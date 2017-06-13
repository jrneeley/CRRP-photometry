import numpy as np
import re
import glob
import read_dao
from astropy.io import fits
import matplotlib.pyplot as mp
from astropy.stats import LombScargle
import plotting_utilities
from scipy import stats

def make_lcv(channels, stars, dao_ids):

    for channel in channels:
        alf_list = glob.glob('all/'+channel+'*.alf')

        phot_data = np.zeros(len(dao_ids), dtype=[('filter', 'S2', len(alf_list)),
            ('id', 'S8'), ('aor', 'i8', len(alf_list)), ('mjd', float, len(alf_list)),
            ('f_num', 'i2', len(alf_list)), ('x', float, len(alf_list)),
            ('y', float, len(alf_list)), ('psf_mag', float, len(alf_list)),
            ('psf_err', float, len(alf_list))])


        phot_data['id'] = dao_ids

        for ind in range(0,len(alf_list)):

            alf_id, x, y, alf_mag, alf_err = read_dao.read_alf(alf_list[ind])
            fits_file = re.sub(".alf",".fits",alf_list[ind])
            hdulist = fits.open(fits_file, mode='update')
            prihdr = hdulist[0].header
            mjd = prihdr['mjd_obs']

            for ind2 in range(0,len(dao_ids)):
                alf_match = np.argwhere(alf_id == dao_ids[ind2])
                trash1, aor_num, frame_num, trash2 = alf_list[ind].split('_')
                trash1, aor_num, trash2 = alf_list[ind].split('_')
                phot_data['aor'][ind2,ind] = int(aor_num)
                phot_data['f_num'][ind2, ind] = int(frame_num)
                phot_data['mjd'][ind2,ind] = mjd
                phot_data['filter'][ind2, ind] = channel
                if len(alf_match):
                    phot_data['x'][ind2,ind] = x[alf_match]
                    phot_data['y'][ind2,ind] = y[alf_match]
                    phot_data['psf_mag'][ind2, ind] = alf_mag[alf_match]
                    phot_data['psf_err'][ind2, ind] = alf_err[alf_match]
                else:
                    phot_data['x'][ind2,ind] = float('NaN')
                    phot_data['y'][ind2,ind] = float('NaN')
                    phot_data['psf_mag'][ind2,ind] = float('NaN')
                    phot_data['psf_err'][ind2,ind] = float('NaN')
        print 'Writing to file...'
        for ind in range(0,len(dao_ids)):

            data_save = np.array(zip(phot_data['filter'][ind], phot_data['aor'][ind],
                phot_data['f_num'][ind], phot_data['mjd'][ind],
                phot_data['x'][ind], phot_data['y'][ind],
                phot_data['psf_mag'][ind], phot_data['psf_err'][ind]),
                dtype=[('c1', 'S2'), ('c2', int), ('c3', int), ('c4', float),
                ('c5', float), ('c6', float), ('c7', float), ('c8', float)])
            if channel == channels[0]:
                f_handle = open('lcvs/'+stars[ind]+'.lcv', 'w')
            else:
                f_handle = open('lcvs/'+stars[ind]+'.lcv', 'a')
            np.savetxt(f_handle, data_save, comments='', fmt='%s %8i %2i %10.4f %7.3f %7.3f %6.3f %6.4f')
            f_handle.close()

def make_mosaic_lcv(channels, stars, dao_ids):

    for channel in channels:
        alf_list = glob.glob('mosaics/'+channel+'*.alf')

        phot_data = np.zeros(len(dao_ids), dtype=[('filter', 'S2', len(alf_list)),
            ('id', 'S8'), ('aor', 'i8', len(alf_list)), ('mjd', float, len(alf_list)),
            ('f_num', 'i2', len(alf_list)), ('x', float, len(alf_list)),
            ('y', float, len(alf_list)), ('psf_mag', float, len(alf_list)),
            ('psf_err', float, len(alf_list))])


        phot_data['id'] = dao_ids

        for ind in range(0,len(alf_list)):

            alf_id, x, y, alf_mag, alf_err = read_dao.read_alf(alf_list[ind])
            fits_file = re.sub(".alf",".fits",alf_list[ind])
            hdulist = fits.open(fits_file, mode='update')
            prihdr = hdulist[0].header
            mjd = prihdr['mjd_obs']

            for ind2 in range(0,len(dao_ids)):
                alf_match = np.argwhere(alf_id == dao_ids[ind2])
                trash1, aor_num, trash2 = alf_list[ind].split('_')
                phot_data['aor'][ind2,ind] = int(aor_num)
                phot_data['f_num'][ind2, ind] = 1
                phot_data['mjd'][ind2,ind] = mjd
                phot_data['filter'][ind2, ind] = channel
                if len(alf_match):
                    phot_data['x'][ind2,ind] = x[alf_match]
                    phot_data['y'][ind2,ind] = y[alf_match]
                    phot_data['psf_mag'][ind2, ind] = alf_mag[alf_match]
                    phot_data['psf_err'][ind2, ind] = alf_err[alf_match]
                else:
                    phot_data['x'][ind2,ind] = float('NaN')
                    phot_data['y'][ind2,ind] = float('NaN')
                    phot_data['psf_mag'][ind2,ind] = float('NaN')
                    phot_data['psf_err'][ind2,ind] = float('NaN')
        print 'Writing to file...'
        for ind in range(0,len(dao_ids)):

            data_save = np.array(zip(phot_data['filter'][ind], phot_data['aor'][ind],
                phot_data['f_num'][ind], phot_data['mjd'][ind],
                phot_data['x'][ind], phot_data['y'][ind],
                phot_data['psf_mag'][ind], phot_data['psf_err'][ind]),
                dtype=[('c1', 'S2'), ('c2', int), ('c3', int), ('c4', float),
                ('c5', float), ('c6', float), ('c7', float), ('c8', float)])
            if channel == channels[0]:
                f_handle = open('mosaic_lcvs/'+stars[ind]+'.lcv', 'w')
            else:
                f_handle = open('mosaic_lcvs/'+stars[ind]+'.lcv', 'a')
            np.savetxt(f_handle, data_save, comments='', fmt='%s %8i %2i %10.4f %7.3f %7.3f %6.3f %6.4f')
            f_handle.close()

def compare_lcv(lcv_file):

    star = lcv_file.split('/')[1]
    dtype1 = np.dtype([('mjd', float), ('mag', float), ('err', float)])
    data = np.loadtxt(lcv_file, dtype=dtype1, usecols=(3,6,7))

    file2 = 'Peter/NGC6121'+star
    dtype2 = np.dtype([('mag', float), ('err', float), ('band', int), ('jd', float)])
    data2 = np.loadtxt(file2, dtype=dtype2, usecols=(0,1,2,5))
    p_mag = data2['mag'][data2['band'] == 1]
    p_mjd = data2['jd'][data2['band'] == 1]+55999.5

    if len(p_mag):

        avg_me = np.nanmean(data['mag'])
        avg_p = np.nanmean(p_mag)
        mp.plot(data['mjd'], data['mag']-avg_me, 'bo')
        mp.plot(p_mjd, p_mag-avg_p, 'ro')
        mp.xlabel('MJD')
        mp.ylabel('mag - avg_mag')

        plt_file = re.sub('.lcv','.pdf', lcv_file)
        mp.savefig(plt_file)
        mp.gcf().clear()
def compare_phased_lcv(lcv_file):

#    star = lcv_file.split('/')[1]
    dtype1 = np.dtype([('phase', float), ('mag', float), ('err', float)])
    data = np.loadtxt(lcv_file, dtype=dtype1, usecols=(1,2,3))

#    file2 = 'Peter/'+star+'_corr.phased'
    file2 = re.sub('lcvs/', 'Peter/', lcv_file)
    file2 = re.sub('.phased', '_corr.phased', file2)
    dtype2 = np.dtype([('band', int), ('phase', float), ('mag', float), ('err', float)])
    data2 = np.loadtxt(file2, dtype=dtype2, usecols=(0,1,2,3))
    p_mag = data2['mag'][data2['band'] == 1]
    p_ph = data2['phase'][data2['band'] == 1]
    p_er = data2['err'][data2['band'] ==1]
    if len(p_mag):

        avg_me = np.nanmean(data['mag'])
        avg_p = np.nanmean(p_mag)
#        mp.plot(data['phase'], data['mag']-avg_me, 'bo')
        mp.errorbar(p_ph, p_mag, yerr=p_er, fmt='o')
        mp.ylim((np.max(p_mag)+0.1, np.max(p_mag)-0.5))
        mp.xlabel('phase')
        mp.ylabel('mag')

        plt_file = re.sub('.phased','.pdf', file2)
        mp.savefig(plt_file)
        mp.gcf().clear()

def phase_lcv(lcv_file, period, T0, bin=1, save=1, plot=0):

    dtype1 = np.dtype([('filter', 'S2'), ('aor', 'i8'), ('mjd', float),
        ('mag', float), ('err', float)])
    data = np.loadtxt(lcv_file, dtype=dtype1, usecols=(0,1,3,6,7))

    filters = np.unique(data['filter'])
    for filt in filters:

        mag_all = data['mag'][data['filter'] == filt]
        err_all = data['err'][data['filter'] == filt]
        mjd_all = data['mjd'][data['filter'] == filt]
        aor_all = data['aor'][data['filter'] == filt]

        if (bin == 1):
            aors = np.unique(aor_all)
            num_aors = len(aors)
            mag = np.zeros(num_aors)
            err = np.zeros(num_aors)
            mjd = np.zeros(num_aors)
            band = np.zeros(num_aors, dtype='S2')

            for ind, aor in enumerate(aors):
                band[ind] = filt
                num = len(mag_all[~np.isnan(mag_all[aor_all == aor])])
                if (num >= 2):
                    mag[ind] = np.nanmean(mag_all[aor_all == aor])
                    mjd[ind] = np.nanmean(mjd_all[aor_all == aor])
                    err[ind] = np.nanstd(mag_all[aor_all== aor])
                if (num == 1):
                    epoch_mag = mag_all[aor_all == aor]
                    epoch_err = err_all[aor_all == aor]
                    epoch_mjd = mjd_all[aor_all == aor]
                    mag[ind] = epoch_mag[~np.isnan(epoch_mag)]
                    mjd[ind] = epoch_mjd[~np.isnan(epoch_mag)]
                    err[ind] = epoch_err[~np.isnan(epoch_mag)]
                if (num == 0):
                    mag[ind] = np.nan
                    mjd[ind] = np.nan
                    err[ind] = np.nan
        else:
            mag = mag_all
            err = err_all
            mjd = mjd_all
            band = np.repeat(filt, len(mag))

        phase = np.mod((mjd - T0)/period, 1)
        if save == 1:
            phased_file = re.sub('.lcv', '.phased', lcv_file)
            if filt == filters[0]:
                f_handle = open(phased_file, 'w')
            else:
                f_handle = open(phased_file, 'a')
            data_save = np.array(zip(band, phase, mag, err), dtype=[('c1', 'S2'),
                ('c2', float), ('c3', float), ('c4', float)])
            np.savetxt(f_handle, data_save, fmt='%s %8.6f %6.3f %5.3f')
            f_handle.close()

            # remove any NaN values
            phase = phase[~np.isnan(mag)]
            err = err[~np.isnan(mag)]
            mag = mag[~np.isnan(mag)]

            mp.errorbar(phase, mag, yerr=err, fmt='o')
            mp.ylim((np.max(mag)+0.2, np.min(mag)-0.2))
            mp.xlabel('Phase')
            mp.ylabel('Mag')
            plot_file = re.sub('.phased', '_'+filt+'_ph.pdf', phased_file)
            mp.savefig(plot_file)
            mp.gcf().clear()

        if plot == 1:
            phase = phase[~np.isnan(mag)]
            err = err[~np.isnan(mag)]
            mag = mag[~np.isnan(mag)]

            mp.errorbar(phase, mag, yerr=err, fmt='o')
            mp.ylim((np.max(mag)+0.2, np.min(mag)-0.2))
            mp.xlabel('Phase')
            mp.ylabel('Mag')
            mp.show()

def plot_raw_lcv(lcv_file, dao_id):

    dtype1 = np.dtype([('filter', 'S2'), ('aor', 'i8'), ('mjd', float),
        ('mag', float), ('err', float)])
    data = np.loadtxt(lcv_file, dtype=dtype1, usecols=(0,1,3,6,7))

    filters = np.unique(data['filter'])
    for filt in filters:

        mag_all = data['mag'][data['filter'] == filt]
        err_all = data['err'][data['filter'] == filt]
        mjd_all = data['mjd'][data['filter'] == filt]
        aor_all = data['aor'][data['filter'] == filt]


        aors = np.unique(aor_all)
        num_aors = len(aors)
        mag = np.zeros(num_aors)
        err = np.zeros(num_aors)
        mjd = np.zeros(num_aors)

        for ind, aor in enumerate(aors):
            num = len(mag_all[~np.isnan(mag_all[aor_all == aor])])
            if (num >= 2):
                mag[ind] = np.nanmedian(mag_all[aor_all == aor])
                mjd[ind] = np.nanmean(mjd_all[aor_all == aor])
                err[ind] = np.nanstd(mag_all[aor_all== aor])
            if (num == 1):
                epoch_mag = mag_all[aor_all == aor]
                epoch_err = err_all[aor_all == aor]
                epoch_mjd = mjd_all[aor_all == aor]
                mag[ind] = epoch_mag[~np.isnan(epoch_mag)]
                mjd[ind] = epoch_mjd[~np.isnan(epoch_mag)]
                err[ind] = epoch_err[~np.isnan(epoch_mag)]
            if (num == 0):
                mag[ind] = np.nan
                mjd[ind] = np.nan
                err[ind] = np.nan


        if ~np.isnan(mag).any():
            mp.plot(mjd_all, mag_all, 'ro')
            mp.plot(mjd, mag, 'bo')
            mp.xlabel('MJD')
            mp.ylabel('Mag')
            mp.ylim(np.max(mag)+0.2, np.min(mag)-0.2)
            mp.title(dao_id)
            plot_file = re.sub('.lcv', '_'+filt+'_raw.pdf', lcv_file)
            mp.savefig(plot_file)
            mp.gcf().clear()

def plot_raw_mosaic_lcv(lcv_file, dao_id):

    dtype1 = np.dtype([('filter', 'S2'), ('aor', 'i8'), ('mjd', float),
        ('mag', float), ('err', float)])
    data = np.loadtxt(lcv_file, dtype=dtype1, usecols=(0,1,3,6,7))

    filters = np.unique(data['filter'])
    for filt in filters:

        mag_all = data['mag'][data['filter'] == filt]
        err_all = data['err'][data['filter'] == filt]
        mjd_all = data['mjd'][data['filter'] == filt]
        aor_all = data['aor'][data['filter'] == filt]


        if ~np.isnan(mag_all).any():
            mp.errorbar(mjd_all, mag_all, yerr=err_all, fmt='o', color='b')
            mp.xlabel('MJD')
            mp.ylabel('Mag')
            mp.ylim(np.max(mag_all)+0.2, np.min(mag_all)-0.2)
            mp.title(dao_id)
            plot_file = re.sub('\.lcv', '_'+filt+'_raw.pdf', lcv_file)
            mp.savefig(plot_file)
            mp.gcf().clear()

def read_optical_lcv(lcv_file, old=0):

    dtype1 = np.dtype([('mag', float), ('err', float), ('filter', int),
        ('year', int), ('day', float), ('source', 'S30')])
    if old == 1:
        data = np.loadtxt(lcv_file, dtype=dtype1, usecols=(0,1,2,4,5,8))
    if old == 0:
        data = np.loadtxt(lcv_file, dtype=dtype1, usecols=(0,1,2,4,5,10))

    V = np.zeros((4, len(data['filter'][data['filter'] == 1])), dtype=object)
    V[0][:] = data['mag'][data['filter'] == 1]
    V[1][:] = data['err'][data['filter'] == 1]
    V[2][:] = data['year'][data['filter'] == 1]*1000 + data['day'][data['filter'] == 1]
    V[3][:] = data['source'][data['filter'] == 1]
    B = np.zeros((4, len(data['filter'][data['filter'] == 2])), dtype=object)
    B[0][:] = data['mag'][data['filter'] == 2]
    B[1][:] = data['err'][data['filter'] == 2]
    B[2][:] = data['year'][data['filter'] == 2]*1000 + data['day'][data['filter'] == 2]
    B[3][:] = data['source'][data['filter'] == 2]
    I = np.zeros((4, len(data['filter'][data['filter'] == 3])), dtype=object)
    I[0][:] = data['mag'][data['filter'] == 3]
    I[1][:] = data['err'][data['filter'] == 3]
    I[2][:] = data['year'][data['filter'] == 3]*1000 + data['day'][data['filter'] == 3]
    I[3][:] = data['source'][data['filter'] == 3]
    R = np.zeros((4, len(data['filter'][data['filter'] == 4])), dtype=object)
    R[0][:] = data['mag'][data['filter'] == 4]
    R[1][:] = data['err'][data['filter'] == 4]
    R[2][:] = data['year'][data['filter'] == 4]*1000 + data['day'][data['filter'] == 4]
    R[3][:] = data['source'][data['filter'] == 4]
    U = np.zeros((4, len(data['filter'][data['filter'] == 5])), dtype=object)
    U[0][:] = data['mag'][data['filter'] == 5]
    U[1][:] = data['err'][data['filter'] == 5]
    U[2][:] = data['year'][data['filter'] == 5]*1000 + data['day'][data['filter'] == 5]
    U[3][:] = data['source'][data['filter'] == 5]

    return U, B, V, R, I

def plot_raw_optical_lcv(U, B, V, R, I):

    mp.errorbar(U[2], U[0], yerr = U[1], fmt='o', color='r')
    mp.errorbar(B[2], B[0]-1, yerr = B[1], fmt='o', color='b')
    mp.errorbar(V[2], V[0]-2, yerr = V[1], fmt='o', color='k')
    mp.errorbar(R[2], R[0]-3, yerr = R[1], fmt='o', color='c')
    mp.errorbar(I[2], I[0]-5, yerr = I[1], fmt='o', color='g')
    mags_all = np.append(U[0], B[0]-1)
    mags_all = np.append(mags_all, V[0]-2)
    mags_all = np.append(mags_all, R[0]-3)
    mags_all = np.append(mags_all, I[0]-5)

    mp.ylim(np.max(mags_all)+0.3, np.min(mags_all)-0.3)
    mp.show()

def plot_phased_optical_lcv(U, B, V, R, I, period, name, datasets, plot_save=0):

    fig, axs = mp.subplots(5, 1, figsize=(10,13), sharex=True)
    filters = ['U', 'B', 'V', 'R', 'I']
#    datasets = np.array(['danish', 'bond5', 'manu', 'Y1007', 'pwm', 'emmi8',
#        'wfi5', 'wfi6', 'apr97', 'wfi10', 'bond7', 'not017', 'fors20602', 'fors20605'])
#    colors = ['k', 'b', 'g', 'r', 'm', 'c', 'y', 'xkcd:coral', 'xkcd:brown',
#            'xkcd:purple', 'xkcd:lavender', 'xkcd:yellowgreen', 'xkcd:olive', 'xkcd:gold']
    for num in range(5):
        if num == 0:
            mags = U[0]
            errs = U[1]
            mjd = U[2]
            source = U[3]
        if num == 1:
            mags = B[0]
            errs = B[1]
            mjd = B[2]
            source = B[3]
        if num == 2:
            mags = V[0]
            errs = V[1]
            mjd = V[2]
            source = V[3]
        if num == 3:
            mags = R[0]
            errs = R[1]
            mjd = R[2]
            source = R[3]
        if num == 4:
            mags = I[0]
            errs = I[1]
            mjd = I[2]
            source = I[3]

        if len(source) == 0:
            continue
        sources_prefix = np.zeros(len(source), dtype='S30')
        for ind, string in enumerate(source):
            sources_prefix[ind] = string.split(':')[0]
    #    datasets = np.unique(datasets_all)

        phase = np.mod(mjd/period, 1)


        for ind, dataset in enumerate(datasets):
            ph = phase[sources_prefix == dataset]
            mag = mags[sources_prefix == dataset]
            err = errs[sources_prefix == dataset]
            if len(ph) == 0:
                continue
            axs[num].errorbar(ph, mag, yerr=err, fmt='o', color=plotting_utilities.get_color(ind))
        axs[num].set_ylim(np.max(mags)+0.2, np.min(mags)-0.2)
        axs[num].set_ylabel(filters[num])

    axs[0].set_title(name+' P = {}'.format(period))
    axs[4].set_xlabel('Phase')
    if plot_save == 1:
        mp.savefig('lcvs/'+name+'-optical.pdf')
    if plot_save == 0:
        mp.show()
#    mp.gcf().clear()
    mp.close()

def find_period(mag, error, mjd, initial_guess):

    x = np.array(mjd, dtype=float)
    y = np.array(mag, dtype=float)
    er = np.array(error, dtype=float)

    guess_frequency = 1/initial_guess
    frequency = np.linspace(guess_frequency-0.05, guess_frequency+0.05, 100)
    power = LombScargle(x, y, er).power(frequency)
    mp.plot(frequency, power)
#    mp.show()

    return 1/frequency[power == np.max(power)]#, frequency, power

def period_search(V, initial_guess, name, plot_save=0):

    x = np.array(V[2], dtype=float)
    y = np.array(V[0], dtype=float)
    er = np.array(V[1], dtype=float)

    best_period = initial_guess

    for iteration in range(3):
        if iteration == 0:
            period_offset = 0.1
            grid_num = 1000
        if iteration == 1:
            period_offset = 0.01
            grid_num = 1000
        if iteration == 2:
            period_offset = 0.001
            grid_num = 10000
        periods = np.linspace(initial_guess-period_offset/2, initial_guess+period_offset/2, grid_num)
        best_std = 99
        avg_std = np.zeros(len(periods))
        for ind, period in enumerate(periods):
            phase = np.mod(x/period, 1)
            stds, edges, bin_num = stats.binned_statistic(phase, y, statistic=np.std, bins=100)
            avg_std[ind] = np.nanmean(stds)
            if avg_std[ind] < best_std:
                best_std = avg_std[ind]
                best_period = period
        mp.plot(periods, avg_std, 'ro')
        mp.axvline(best_period)
        if plot_save == 1:
            mp.savefig('lcvs/'+name+'-period.pdf')
        if plot_save == 0:
            mp.show()
    mp.close()

    return best_period
