import numpy as np
import re
import glob
import read_dao
from astropy.io import fits
import matplotlib.pyplot as mp
from astropy.stats import LombScargle
import plotting_utilities
from scipy import stats
from scipy.interpolate import interp1d
from scipy import signal
from matplotlib import gridspec
import os.path



def make_lcv(channels, stars, dao_ids, data_dir=''):

    lcv_folder = data_dir+'lcvs/mir/'
    img_folder = data_dir+'data/'
    for channel in channels:
        file_list = glob.glob(img_folder+channel+'*.cal')
        if len(file_list) == 0:
            file_list = glob.glob(img_folder+channel+'*.alf')
        phot_data = np.zeros(len(dao_ids), dtype=[('filter', 'S2', len(file_list)),
            ('id', 'S8'), ('aor', int, len(file_list)), ('mjd', float, len(file_list)),
            ('f_num', int, len(file_list)), ('x', float, len(file_list)),
            ('y', float, len(file_list)), ('psf_mag', float, len(file_list)),
            ('psf_err', float, len(file_list))])


        phot_data['id'] = dao_ids

        for ind in range(0,len(file_list)):

            alf_id, x, y, alf_mag, alf_err = read_dao.read_alf(file_list[ind])
            fits_file = re.sub(".cal",".fits", file_list[ind])
            if fits_file == file_list[ind]:
                fits_file = re.sub(".alf", ".fits", file_list[ind])

            hdulist = fits.open(fits_file, mode='update')
            prihdr = hdulist[0].header
            mjd = prihdr['mjd_obs']

            for ind2 in range(0,len(dao_ids)):
                alf_match = np.argwhere(alf_id == dao_ids[ind2])
                trash1, aor_num, frame_num, trash2 = file_list[ind].split('_')
        #        trash1, aor_num, trash2 = alf_list[ind].split('_')
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
    #    print 'Writing to file...'
        for ind in range(0,len(dao_ids)):
            if np.all(np.isnan(phot_data['psf_mag'][ind])):
                continue

            data_save = np.array(zip(phot_data['filter'][ind], phot_data['aor'][ind],
                phot_data['f_num'][ind], phot_data['mjd'][ind],
                phot_data['x'][ind], phot_data['y'][ind],
                phot_data['psf_mag'][ind], phot_data['psf_err'][ind]),
                dtype=[('c1', 'S2'), ('c2', int), ('c3', int), ('c4', float),
                ('c5', float), ('c6', float), ('c7', float), ('c8', float)])
            if channel == channels[0]:
                f_handle = open(lcv_folder+stars[ind]+'.lcv', 'w')
            else:
                f_handle = open(lcv_folder+stars[ind]+'.lcv', 'a')
            np.savetxt(f_handle, data_save, comments='', fmt='%s %8i %2i %10.4f %7.3f %7.3f %6.3f %6.4f')
            f_handle.close()

def make_mosaic_lcv(channels, stars, dao_ids):

    for channel in channels:
        alf_list = glob.glob('mosaics/'+channel+'*.alf')

        phot_data = np.zeros(len(dao_ids), dtype=[('filter', 'S2', len(alf_list)),
            ('id', 'S8'), ('aor', int, len(alf_list)), ('mjd', float, len(alf_list)),
            ('f_num', int, len(alf_list)), ('x', float, len(alf_list)),
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


def phase_lcv(lcv_file, period, T0, bin=0, save=1, plot=0, error_threshold=0.3):

    dtype1 = np.dtype([('filter', 'S2'), ('aor', int), ('mjd', float),
        ('mag', float), ('err', float)])
    data = np.loadtxt(lcv_file, dtype=dtype1, usecols=(0,1,3,6,7))

    filters = np.unique(data['filter'])
    for filt in filters:

        mag_all = data['mag'][data['filter'] == filt]
        err_all = data['err'][data['filter'] == filt]
        mjd_all = data['mjd'][data['filter'] == filt]
        aor_all = data['aor'][data['filter'] == filt]

        # remove NaNs
        err_all = err_all[~np.isnan(mag_all)]
        mjd_all = mjd_all[~np.isnan(mag_all)]
        aor_all = aor_all[~np.isnan(mag_all)]
        mag_all = mag_all[~np.isnan(mag_all)]

        # Filter on the uncertainty
        mag_all = mag_all[err_all < error_threshold]
        mjd_all = mjd_all[err_all < error_threshold]
        aor_all = aor_all[err_all < error_threshold]
        err_all = err_all[err_all < error_threshold]


        if (bin == 1):
            aors = np.unique(aor_all)
            num_aors = len(aors)
            mag = np.zeros(num_aors)
            err = np.zeros(num_aors)
            mjd = np.zeros(num_aors)
            band = np.zeros(num_aors, dtype='S2')

            for ind, aor in enumerate(aors):
                band[ind] = filt
                dither_mags = mag_all[aor_all == aor]
                dither_errs = err_all[aor_all == aor]
                dither_mjds = mjd_all[aor_all == aor]
                num = len(dither_mags)

                if (num > 2):
                    dispersion = np.std(dither_mags)
                    mean_mag = np.mean(dither_mags)
                    mag[ind] = np.nanmean(dither_mags[abs(dither_mags - mean_mag) < 2*dispersion])
                    mjd[ind] = np.nanmean(dither_mjds[abs(dither_mags - mean_mag) < 2*dispersion])
                    err[ind] = np.nanstd(dither_mags[abs(dither_mags - mean_mag) < 2*dispersion])
                if (num ==2):
                    mag[ind] = np.nanmean(dither_mags)
                    mjd[ind] = np.nanmean(dither_mjds)
                    err[ind] = np.sqrt(dither_errs[0]**2+dither_errs[1]**2)
                if (num == 1):
                    mag[ind] = mag_all[aor_all == aor]
                    err[ind] = err_all[aor_all == aor]
                    mjd[ind] = mjd_all[aor_all == aor]
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
            phased_file = re.sub('\.lcv', '.phased', lcv_file)
            if filt == filters[0]:
                f_handle = open(phased_file, 'w')
            else:
                f_handle = open(phased_file, 'a')
            data_save = np.array(zip(band, mjd, phase, mag, err), dtype=[('c1', 'S2'),
                ('c2', float), ('c3', float), ('c4', float), ('c5', float)])
            np.savetxt(f_handle, data_save, fmt='%s %10.4f %8.6f %6.3f %5.3f')
            f_handle.close()

            # remove any NaN values
        #    phase = phase[~np.isnan(mag)]
        #    err = err[~np.isnan(mag)]
        #    mag = mag[~np.isnan(mag)]

            if filt == 'I1':
                color='b'
            if filt == 'I2':
                color='r'
            mp.errorbar(phase, mag, yerr=err, fmt='o', color=color)
            mp.ylim((np.max(mag)+0.3, np.min(mag)-0.8))
            mp.xlabel('Phase')
            mp.ylabel('Mag')
            plot_file = re.sub('.phased','-ph.pdf', phased_file)
            mp.savefig(plot_file)


        if plot == 1:
            phase = phase[~np.isnan(mag)]
            err = err[~np.isnan(mag)]
            mag = mag[~np.isnan(mag)]
            if filt == filters[0]:
                mp.subplot(211)
            else:
                mp.subplot(212)
            mp.errorbar(phase, mag, yerr=err, fmt='o')
            mp.ylim((np.max(mag)+0.2, np.min(mag)-0.2))
            mp.xlabel('Phase')
            mp.ylabel(filt+' Mag')
            mp.show()
    mp.gcf().clear()

def phase_lcv_all_bands(target, lcv, period, T0, optical_lcv=0, nir_lcv=0, mir_lcv=0, bin_mir=0, data_dir=''):

    fig = mp.figure(figsize=(8,10))
    path_to_optical = data_dir+'lcvs/optical/'
    path_to_nir = data_dir+'lcvs/nir/'
    path_to_mir = data_dir+'lcvs/mir/'

    if os.path.isfile(path_to_optical+target+lcv):
        optical_lcv = 1
    if os.path.isfile(path_to_nir+target+lcv):
        nir_lcv = 1
    if os.path.isfile(path_to_mir+lcv):
        mir_lcv = 1

    if optical_lcv == 1:
        phased_file = re.sub('\.lcv', '.phased', lcv)

        U, B, V, R, I = read_optical_lcv(path_to_optical+target+lcv)


        Uph = np.mod((U[2]-T0)/period, 1)
        Bph = np.mod((B[2]-T0)/period, 1)
        Vph = np.mod((V[2]-T0)/period, 1)
        Rph = np.mod((R[2]-T0)/period, 1)
        Iph = np.mod((I[2]-T0)/period, 1)

        ## Add to data to plot
        mp.errorbar(Uph, U[0]+1, yerr=U[1], fmt='o')
        mp.errorbar(Bph, B[0]+0.5, yerr=B[1], fmt='v')
        mp.errorbar(Vph, V[0], yerr=V[1], fmt='s')
        mp.errorbar(Rph, R[0]-0.25, yerr=R[1], fmt='p')
        mp.errorbar(Iph, I[0]-0.5, yerr=I[1], fmt='P')

        ## Add data to file
        f_handle = open(data_dir+'lcvs/'+phased_file, 'w')

        data_save = np.array(zip(np.repeat('U', len(U[2])), U[2], Uph, U[0], U[1], U[3]),
            dtype=[('c1', 'S2'), ('c2', float), ('c3', float), ('c4', float),
            ('c5', float), ('c6', 'S30')])
        np.savetxt(f_handle, data_save, fmt='%s %10.4f %8.6f %6.3f %5.3f %s')
        f_handle.close()
        f_handle = open(data_dir+'lcvs/'+phased_file, 'a')
        data_save = np.array(zip(np.repeat('B', len(B[2])), B[2], Bph, B[0], B[1], B[3]),
            dtype=[('c1', 'S2'), ('c2', float), ('c3', float), ('c4', float),
            ('c5', float), ('c6', 'S30')])
        np.savetxt(f_handle, data_save, fmt='%s %10.4f %8.6f %6.3f %5.3f %s')
        data_save = np.array(zip(np.repeat('V', len(V[2])), V[2], Vph, V[0], V[1], V[3]),
            dtype=[('c1', 'S2'), ('c2', float), ('c3', float), ('c4', float),
            ('c5', float), ('c6', 'S30')])
        np.savetxt(f_handle, data_save, fmt='%s %10.4f %8.6f %6.3f %5.3f %s')
        data_save = np.array(zip(np.repeat('R', len(R[2])), R[2], Rph, R[0], R[1], R[3]),
            dtype=[('c1', 'S2'), ('c2', float), ('c3', float), ('c4', float),
            ('c5', float), ('c6', 'S30')])
        np.savetxt(f_handle, data_save, fmt='%s %10.4f %8.6f %6.3f %5.3f %s')
        data_save = np.array(zip(np.repeat('I', len(I[2])), I[2], Iph, I[0], I[1], I[3]),
            dtype=[('c1', 'S2'), ('c2', float), ('c3', float), ('c4', float),
            ('c5', float), ('c6', 'S30')])
        np.savetxt(f_handle, data_save, fmt='%s %10.4f %8.6f %6.3f %5.3f %s')

    if mir_lcv == 1:

        dtype1 = np.dtype([('filter', 'S2'), ('aor', int), ('frame', int), ('mjd', float),
            ('mag', float), ('err', float)])
        data = np.loadtxt(path_to_mir+lcv, dtype=dtype1, usecols=(0,1,2,3,6,7))

        filters = np.unique(data['filter'])
        for filt in filters:

            mag_all = data['mag'][data['filter'] == filt]
            err_all = data['err'][data['filter'] == filt]
            mjd_all = data['mjd'][data['filter'] == filt]
            aor_all = data['aor'][data['filter'] == filt]
            frame_all = data['frame'][data['filter'] == filt]

            # remove NaNs
            err_all = err_all[~np.isnan(mag_all)]
            mjd_all = mjd_all[~np.isnan(mag_all)]
            aor_all = aor_all[~np.isnan(mag_all)]
            frame_all = frame_all[~np.isnan(mag_all)]
            mag_all = mag_all[~np.isnan(mag_all)]
            # Filter on the uncertainty
            error_threshold = 0.5
            mag_all = mag_all[err_all < error_threshold]
            mjd_all = mjd_all[err_all < error_threshold]
            aor_all = aor_all[err_all < error_threshold]
            frame_all = frame_all[err_all < error_threshold]
            err_all = err_all[err_all < error_threshold]

            if (bin_mir == 1):
                aors = np.unique(aor_all)
                num_aors = len(aors)
                mag = np.zeros(num_aors)
                err = np.zeros(num_aors)
                mjd = np.zeros(num_aors)
                band = np.zeros(num_aors, dtype='S2')
                source = np.zeros(num_aors, dtype='S15')


                for ind, aor in enumerate(aors):

                    source[ind] = 'r'+str(aor)
                    band[ind] = filt
                    dither_mags = mag_all[aor_all == aor]
                    dither_errs = err_all[aor_all == aor]
                    dither_mjds = mjd_all[aor_all == aor]
                    num = len(dither_mags)

                    if (num > 2):
                        dispersion = np.std(dither_mags)
                        mean_mag = np.mean(dither_mags)
                        mag[ind] = np.nanmean(dither_mags[abs(dither_mags - mean_mag) < 2*dispersion])
                        mjd[ind] = np.nanmean(dither_mjds[abs(dither_mags - mean_mag) < 2*dispersion])
                        err[ind] = np.nanstd(dither_mags[abs(dither_mags - mean_mag) < 2*dispersion])
                    if (num == 2):
                        mag[ind] = np.nanmean(dither_mags)
                        mjd[ind] = np.nanmean(dither_mjds)
                        err[ind] = np.sqrt(dither_errs[0]**2+dither_errs[1]**2)
                    if (num == 1):
                        mag[ind] = mag_all[aor_all == aor]
                        err[ind] = err_all[aor_all == aor]
                        mjd[ind] = mjd_all[aor_all == aor]
                    if (num == 0):
                        mag[ind] = np.nan
                        mjd[ind] = np.nan
                        err[ind] = np.nan
            else:
                mag = mag_all
                err = err_all
                mjd = mjd_all
                source = ['r'+str(aor_all[i])+':'+str(frame_all[i]) for i in range(len(aor_all))]

                band = np.repeat(filt, len(mag))

            phase = np.mod((mjd - T0)/period, 1)

            ## Add to data to plot
            if filt == 'I1':
                offset = -1.0
            if filt == 'I2':
                offset = -1.5
            mp.errorbar(phase, mag+offset, yerr=err, fmt='x')
            ## Add data to file
            data_save = np.array(zip(band, mjd, phase, mag, err, source),
                dtype=[('c1', 'S2'), ('c2', float), ('c3', float), ('c4', float),
                ('c5', float), ('c6', 'S30')])
            np.savetxt(f_handle, data_save, fmt='%2s %10.4f %8.6f %6.3f %5.3f %s')

    mp.ylim((18,10))
    mp.xlabel('Phase')
    mp.ylabel('Mag + offset')
    mp.title(lcv)
    plot_file = re.sub('\.phased', '-ph.pdf', phased_file)
    mp.savefig(data_dir+'lcvs/'+plot_file)
    mp.show()
    f_handle.close()

def plot_raw_lcv(lcv_file, dao_id):

    dtype1 = np.dtype([('filter', 'S2'), ('aor', int), ('mjd', float),
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

    dtype1 = np.dtype([('filter', 'S2'), ('aor', int), ('mjd', float),
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
        ('year', int), ('day', float), ('source', 'S35')])
    if old == 1:
        data = np.loadtxt(lcv_file, dtype=dtype1, usecols=(0,1,2,4,5,8))
    if old == 0:
        data = np.loadtxt(lcv_file, dtype=dtype1, usecols=(0,1,2,4,5,10))

    V = np.zeros((4, len(data['filter'][data['filter'] == 1])), dtype=object)
    V[0][:] = data['mag'][data['filter'] == 1]
    V[1][:] = data['err'][data['filter'] == 1]
    V[2][:] = data['year'][data['filter'] == 1]*1000 + data['day'][data['filter'] == 1]
    V[2] = V[2] - 2400000.5
    V[3][:] = data['source'][data['filter'] == 1]
    B = np.zeros((4, len(data['filter'][data['filter'] == 2])), dtype=object)
    B[0][:] = data['mag'][data['filter'] == 2]
    B[1][:] = data['err'][data['filter'] == 2]
    B[2][:] = data['year'][data['filter'] == 2]*1000 + data['day'][data['filter'] == 2]
    B[2] = B[2] - 2400000.5
    B[3][:] = data['source'][data['filter'] == 2]
    I = np.zeros((4, len(data['filter'][data['filter'] == 3])), dtype=object)
    I[0][:] = data['mag'][data['filter'] == 3]
    I[1][:] = data['err'][data['filter'] == 3]
    I[2][:] = data['year'][data['filter'] == 3]*1000 + data['day'][data['filter'] == 3]
    I[2] = I[2] - 2400000.5
    I[3][:] = data['source'][data['filter'] == 3]
    R = np.zeros((4, len(data['filter'][data['filter'] == 4])), dtype=object)
    R[0][:] = data['mag'][data['filter'] == 4]
    R[1][:] = data['err'][data['filter'] == 4]
    R[2][:] = data['year'][data['filter'] == 4]*1000 + data['day'][data['filter'] == 4]
    R[2] = R[2] - 2400000.5
    R[3][:] = data['source'][data['filter'] == 4]
    U = np.zeros((4, len(data['filter'][data['filter'] == 5])), dtype=object)
    U[0][:] = data['mag'][data['filter'] == 5]
    U[1][:] = data['err'][data['filter'] == 5]
    U[2][:] = data['year'][data['filter'] == 5]*1000 + data['day'][data['filter'] == 5]
    U[2] = U[2] - 2400000.5
    U[3][:] = data['source'][data['filter'] == 5]

    return U, B, V, R, I
def select_datasets(V, datasets):

    datasets_prefix = np.zeros(len(V[3]), dtype='S30')
    for ind, string in enumerate(V[3]):
        datasets_prefix[ind] = string.split(':')[0]
    num_match = len(V[3][datasets_prefix == dataset])

    V_new = np.zeros((4, num_match), dtype=object)

    V_new[0] = V[0][datasets_prefix == dataset]
    V_new[1] = V[1][datasets_prefix == dataset]
    V_new[2] = V[2][datasets_prefix == dataset]
    V_new[3] = V[3][datasets_prefix == dataset]
    return V_new

def plot_raw_optical_lcv(U):#, B, V, R, I):

    mags = U[0].astype(float)
    phase = U[2].astype(float)
    errs = U[1].astype(float)

    mp.errorbar(phase, mags, yerr = errs, fmt='o', color='r')
#    mp.errorbar(B[2], B[0]-1, yerr = B[1], fmt='o', color='b')
#    mp.errorbar(V[2], V[0]-2, yerr = V[1], fmt='o', color='k')
#    mp.errorbar(R[2], R[0]-3, yerr = R[1], fmt='o', color='c')
#    mp.errorbar(I[2], I[0]-5, yerr = I[1], fmt='o', color='g')
#    mags_all = np.append(U[0], B[0]-1)
#    mags_all = np.append(mags_all, V[0]-2)
#    mags_all = np.append(mags_all, R[0]-3)
#    mags_all = np.append(mags_all, I[0]-5)

    mp.ylim(np.nanmax(mags)+0.1, np.nanmin(mags)-0.1)
    mp.show()

def plot_phased_optical_lcv(U, B, V, R, I, period, name, datasets, plot_save=0,
    data_dir='', error_threshold=0, colors=None):


    fig, axs = mp.subplots(5, 1, figsize=(10,13), sharex=True)
    filters = ['U', 'B', 'V', 'R', 'I']

    for num in range(5):
        if num == 0:
            mags = U[0].astype(float)
            errs = U[1].astype(float)
            mjd = U[2].astype(float)
            source = U[3]
        if num == 1:
            mags = B[0].astype(float)
            errs = B[1].astype(float)
            mjd = B[2].astype(float)
            source = B[3]
        if num == 2:
            mags = V[0].astype(float)
            errs = V[1].astype(float)
            mjd = V[2].astype(float)
            source = V[3]
        if num == 3:
            mags = R[0].astype(float)
            errs = R[1].astype(float)
            mjd = R[2].astype(float)
            source = R[3]
        if num == 4:
            mags = I[0].astype(float)
            errs = I[1].astype(float)
            mjd = I[2].astype(float)
            source = I[3]


        if error_threshold > 0:
            mags = mags[errs < error_threshold]
            mjd = mjd[errs < error_threshold]
            source = source[errs < error_threshold]
            errs = errs[errs < error_threshold]
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
            if colors is None:
                color = plotting_utilities.get_color(ind)
            else:
                color = colors[ind]
            axs[num].errorbar(ph, mag, yerr=err, fmt='o', color=color)
        axs[num].set_ylim(np.nanmax(mags)+0.2, np.nanmin(mags)-0.2)
        axs[num].set_ylabel(filters[num])

    axs[0].set_title(name+' P = {}'.format(period))
    axs[4].set_xlabel('Phase')
    if plot_save == 1:
        mp.savefig(data_dir+'lcvs/optical/'+name+'-optical.pdf')
    if plot_save == 0:
        mp.show()
#    mp.gcf().clear()
    mp.close()

def plot_phased_optical_one_band(data, period, name, datasets, plot_save=0):

    mags = data[0]
    errs = data[1]
    mjd = data[2]
    source = data[3]

#    sources_prefix = np.zeros(len(source), dtype='S30')
#    for ind, string in enumerate(source):
#        sources_prefix[ind] = string.split(':')[0]

    phase = np.mod(mjd/period, 1)

#    for ind, dataset in enumerate(datasets):
#        ph = phase[sources_prefix == dataset]
#        mag = mags[sources_prefix == dataset]
#        err = errs[sources_prefix == dataset]
#        if len(ph) == 0:
#            continue
#        mp.errorbar(ph, mag, yerr=err, fmt='o', color=plotting_utilities.get_color(ind))
    mp.errorbar(phase, mags, yerr=errs, fmt='o', color='k')
    mp.ylim(np.max(mags)+0.05, np.min(mags)-0.05)

    mp.title(name+' P = {}'.format(period))
    mp.xlabel('Phase')
    mp.ylabel('mag')
    mp.show()
#    mp.gcf().clear()

def find_period(mag, error, mjd, initial_guess=0.5, num_investigate=5):

    x = np.array(mjd, dtype=float)
    y = np.array(mag, dtype=float)
    er = np.array(error, dtype=float)

    guess_frequency = 1/initial_guess
    min_freq = 1/(initial_guess-0.4)
    max_freq = 1/(initial_guess+0.4)
    frequency = np.linspace(min_freq, max_freq, 1000)
    power = LombScargle(x, y, er).power(frequency)
    mp.plot(frequency, power)
    order = np.argsort(power)
    candidate_periods = 1/frequency[order[-num_investigate:]]
    mp.show()

    return candidate_periods

def period_search(V, initial_guess, name, plot_save=0, error_threshold=0.05):

    x = np.array(V[2][V[1] < error_threshold], dtype=float)
    y = np.array(V[0][V[1] < error_threshold], dtype=float)
    er = np.array(V[1][V[1] < error_threshold], dtype=float)

    best_period = initial_guess
    fig, axs = mp.subplots(4, 1, figsize=(8,10))

    search_range = 0.01
    grid_num = 1000

    for iteration in range(4):
        best_std = 99
        if iteration != 0:
            search_range = search_range/(2*iteration)
            grid_num = search_range*10**(iteration+5)
        min_period = best_period - search_range/2
        max_period = best_period + search_range/2
        periods = np.linspace(min_period, max_period, num=grid_num+1)
        avg_std = np.zeros(len(periods))
        for ind, period in enumerate(periods):
            phase = np.mod(x/period, 1)
            stds, edges, bin_num = stats.binned_statistic(phase, y, statistic=np.std, bins=100)
            avg_std[ind] = np.nanmean(stds)
        order = np.argsort(avg_std)
        best_period = periods[order[0]]
        # Apply a median filter
        yy_smoothed = signal.medfilt(avg_std, 101)
        print iteration, search_range, grid_num, best_period
        axs[iteration].plot(periods, avg_std, 'ro')
        axs[iteration].plot(periods, yy_smoothed, 'b-')
        axs[iteration].axvline(best_period)
    axs[iteration].set_xlabel('Period (days)')
    if plot_save == 1:
        mp.savefig('lcvs/'+name+'-period.pdf')
    if plot_save == 0:
        mp.show()
    mp.close()

    return best_period

# Use to do one round of LombScargle and identify possible periods for RRL star
def period_search_LS(V, name, min_period = 0.2, max_period=1.0, plot_save=0, error_threshold=0.05):

    x1 = np.array(V[2][V[1] < error_threshold], dtype=float)
    y1 = np.array(V[0][V[1] < error_threshold], dtype=float)
    er1 = np.array(V[1][V[1] < error_threshold], dtype=float)

    freq_max = 1/(min_period)
    freq_min = 1/(max_period)
#    frequency = np.linspace(freq_min, freq_max, 10000)
    frequency, power = LombScargle(x1, y1, er1).autopower(minimum_frequency=freq_min,
                        maximum_frequency=freq_max )
    order = np.argsort(power)
    candidate_periods = 1/frequency[order[-10:]]
    mp.plot(1/frequency, power)
#    for period in candidate_periods:
#        mp.axvline(period, color='r', alpha=0.5)

    if plot_save == 0:
        mp.show()
    mp.close()

    best_frequency = frequency[np.argmax(power)]
    t_fit = np.linspace(0,1)
    y_fit = LombScargle(x1, y1, er1, nterms=2).model(t_fit, best_frequency)

    phase_data = np.mod(x1*best_frequency, 1)
    mp.plot(phase_data, y1, 'o')
    mp.plot(t_fit, y_fit)
    mp.show()
#    return candidate_periods
    return 1/best_frequency

def period_search_hybrid(first_band, initial_guess, name, second_band=None,
    plot_save=0, error_threshold=0.1, search_window=0.002,
    num_investigate=5, step_size=10, data_dir=''):


    x = np.array(first_band[2][first_band[1] < error_threshold], dtype=float)
    y = np.array(first_band[0][first_band[1] < error_threshold], dtype=float)
    er = np.array(first_band[1][first_band[1] < error_threshold], dtype=float)
    if second_band is not None:
        x2 = np.array(second_band[2][second_band[1] < error_threshold], dtype=float)
        y2 = np.array(second_band[0][second_band[1] < error_threshold], dtype=float)
        er2 = np.array(second_band[1][second_band[1] < error_threshold], dtype=float)

    delta_time = np.max(x) - np.min(x)
    approx_p = np.round(initial_guess, 1)
    n_cycles = np.floor(delta_time/approx_p)
    max_precision = n_cycles * approx_p / (n_cycles - 0.01) - approx_p
    order = np.ceil(np.abs(np.log10(max_precision)))
    precision = 10**order
#    print 'Max precision = 10^'+str(order)



    best_period = initial_guess
    min_period = best_period - search_window/2
    max_period = best_period + search_window/2
    if initial_guess == np.nan:
        best_period = 0.5
        period_range = 0.2
# Do initial Lomb Scargle
    freq_max = 1/min_period
    freq_min = 1/max_period
#    frequency = np.linspace(freq_min, freq_max, 1000)
    frequency, power = LombScargle(x, y, er).autopower(minimum_frequency=freq_min,
                        maximum_frequency=freq_max )
    order = np.argsort(power)

    possible_periods = 1/frequency[order[-num_investigate:]]
    bins = np.linspace(min_period, max_period, 11)
    bin_index = np.digitize(possible_periods, bins)
    bins_with_peaks = np.unique(bin_index)
    num_candidates = len(bins_with_peaks)
    candidate_periods = np.zeros(num_candidates)
    for ind, i in enumerate(bins_with_peaks):
        if i == len(bins):
            continue
        max_f = 1/bins[i-1]
        min_f = 1/bins[i]
        power_in_bin = power[(frequency <= max_f) & (frequency >= min_f)]
        max_power = np.max(power_in_bin)
        candidate_periods[ind] = 1/frequency[power == max_power]


    fig = mp.figure(figsize=(4*num_candidates, 6))
    ax1 = mp.subplot2grid((2,num_candidates), (0,0), colspan=num_candidates)
    ax1.plot(1/frequency, power)

    best_periods = np.copy(candidate_periods)
    best_stds = np.zeros(num_candidates)
    for ind, period in enumerate(candidate_periods):
        ax1.axvline(period, color='r')
        period_range = search_window

#        for iteration in range(2):

#            period_range = period_range/step_size

#            freq_max = 1/(best_periods[ind]-period_range/2)
#            freq_min = 1/(best_periods[ind]+period_range/2)
#            frequency = np.linspace(freq_min, freq_max, 1000)
#            power = LombScargle(x, y, er).power(frequency)
#            order = np.argsort(power)
#            best_periods[ind] = 1/frequency[order[-1]]

#            ax = mp.subplot2grid((4,num_candidates), (iteration+1, ind))
#            ax.plot(1/frequency, power)

    #    print best_period

        search_range = period_range/100
        grid_num = search_range*precision

        if grid_num > 100000:
            grid_num = 100000
        min_period = best_periods[ind] - search_range/2
        max_period = best_periods[ind] + search_range/2
        periods = np.linspace(min_period, max_period, num=grid_num+1)
        avg_std = np.zeros(len(periods))
        for ind2, trial_period in enumerate(periods):
            phase = np.mod(x/trial_period, 1)
            stds, edges, bin_num = stats.binned_statistic(phase, y, statistic=np.std, bins=100)
            counts, edges, bin_num = stats.binned_statistic(phase, y, statistic='count', bins=100)
            if second_band is not None:
                phase2 = np.mod(x2/trial_period, 1)
                stds2, edges2, bin_num2 = stats.binned_statistic(phase2, y2, statistic=np.std, bins=100)
                counts2, edges2, bin_num2 = stats.binned_statistic(phase2, y2, statistic='count', bins=100)
                avg_std[ind2] = np.mean(stds2[counts2 > 3]) + np.mean(stds[counts > 3])
            else:
                avg_std[ind2] = np.mean(stds[counts > 3])
        #print stds[counts == 0]
        order = np.argsort(avg_std)
        best_periods[ind] = periods[order[0]]
        best_stds[ind] = avg_std[order[0]]

        ax = mp.subplot2grid((2,num_candidates), (1, ind))
        ax.plot(periods, avg_std, 'ro')
        ax.axvline(best_periods[ind])

    if plot_save == 1:
        mp.savefig(data_dir+'lcvs/'+name+'-hybrid.pdf')
    if plot_save == 0:
        mp.show()
    mp.close()
#    print best_periods
#    print best_stds
    best_period = best_periods[best_stds == np.min(best_stds)]
    if len(best_period) > 1:
        best_period = best_period[0]
    return best_period

def gloess(phased_lcv_file, clean=0, smoothing_params=None, plot_save=0, data_dir=''):

    figtosave = mp.figure(figsize=(8,10))
    ax = figtosave.add_subplot(111)

    master_filters = np.array(['U', 'B', 'V', 'R', 'I', 'J', 'H', 'K',
        'I1', 'I2'], dtype='S2')

    dtype1 = np.dtype([('filter', 'S2'), ('mjd', float), ('phase', float), ('mag', float), ('err', float)])
    data = np.loadtxt(phased_lcv_file, dtype=dtype1, usecols=(0,1,2,3,4))
    print phased_lcv_file
    filters = np.unique(data['filter'])
    num_filters = len(filters)

    if smoothing_params is None:
        smoothing_params = np.repeat(1.0, 10)

    master_avg_mag = np.zeros(10)
    master_amp = np.zeros(10)
    master_sigma = np.zeros(10)

    for filt in filters:

        if filt == 'R' or filt == 'U':
            master_avg_mag[master_filters == filt] = np.nan
            master_amp[master_filters == filt] = np.nan
            master_sigma[master_filters == filt] = np.nan
            continue

        sigma = float(smoothing_params[master_filters == filt])

        if sigma == 1.0:
            sigma = 0.1

        phase = data['phase'][data['filter'] == filt]
        mag = data['mag'][data['filter'] == filt]
        err = data['err'][data['filter'] == filt]
        mjd = data['mjd'][data['filter'] == filt]

        phase = phase[~np.isnan(mag)]
        err = err[~np.isnan(mag)]
        mag = mag[~np.isnan(mag)]

        if clean == 1:
            # remove obvious bad data
            dispersion = np.std(mag)
            phase = phase[abs(mag - np.mean(mag)) < 2*dispersion]
            err = err[abs(mag - np.mean(mag)) < 2*dispersion]
            mag = mag[abs(mag - np.mean(mag)) < 2*dispersion]


        phase_copy = np.concatenate((phase-2, phase-1, phase, phase+1, phase+2))
        mag_copy = np.tile(mag, 5)
        err_copy = np.tile(err, 5)

        x = np.arange(0, 1, 0.001)

        n_data = len(mag_copy)
        n_fit = len(x)

        happy = 'n'
        while happy == 'n':
            smoothed_mag = np.zeros(n_fit)
            weight = np.zeros(n_data)

            for ind, step in enumerate(x):

                dist = phase_copy - step
                weight = err_copy * np.exp(dist**2/sigma**2)

                xphase = phase_copy - step

                fit = np.polyfit(phase_copy, mag_copy, 2, w=1/weight)

                smoothed_mag[ind] = fit[2] + fit[1]*step + fit[0]*step**2

            smoothed_mag_copy = np.tile(smoothed_mag, 5)
            x_copy = np.concatenate((x-2, x-1, x, x+1, x+2))

            figshow = mp.figure(figsize=(8,6))
            ax2 = figshow.add_subplot(111)
            ax2.plot(phase_copy, mag_copy, 'bo')
            ax2.plot(x_copy, smoothed_mag_copy, 'r-')
            ax2.set_ylim((np.max(mag)+0.2, np.min(mag)-0.2))
            ax2.set_xlim((-0.2, 1.2))
            ax2.set_xlabel('Phase')
            ax2.set_ylabel(filt+' mag')
            mp.show()


            if smoothing_params[master_filters == filt] != 1.0:
                happy = 'y'
                continue
            check = raw_input('Are you happy with this fit? [y/n]: ')
            if check == 'n':
                sigma = input('Enter new smoothing parameter: ')
            if check == 'y':
                happy = 'y'


        # Derive light curve parameters
        flux = 99*np.power(10,-smoothed_mag/2.5)
        average_flux = np.mean(flux)
        average_mag = -2.5*np.log10(average_flux/99)

        amplitude = np.max(smoothed_mag) - np.min(smoothed_mag)
        ph_max = x[smoothed_mag == np.min(smoothed_mag)]
        ph_min = x[smoothed_mag == np.max(smoothed_mag)]

        if plot_save == 1:

            if filt == 'U':
                offset = 0.8
                marker = 'D'
                color = 'xkcd:violet'
            if filt == 'B':
                offset = 0.4
                marker = 'v'
                color = 'xkcd:indigo'
            if filt == 'V':
                offset = 0
                marker = 'x'
                color = 'b'
            if filt == 'R':
                offset = -0.2
                marker = 'd'
                color = 'c'
            if filt == 'I':
                offset = -0.4
                marker = 's'
                color = 'g'
            if filt == 'J':
                offset = -0.6
                marker = 'h'
                color = 'xkcd:gold'
            if filt == 'H':
                offset = -0.8
                marker = '>'
                color = 'xkcd:orange'
            if filt == 'K':
                offset = -1
                marker = '^'
                color = 'r'
            if filt == 'I1':
                offset = -1
                marker = 'o'
                color = 'xkcd:maroon'
            if filt == 'I2':
                offset = -1.6
                marker = 'p'
                color = 'm'
            ax.errorbar(phase, mag+offset, yerr=err, fmt=marker, color=color)
            ax.plot(x_copy, smoothed_mag_copy+offset, 'r-')

        print filt, average_mag, amplitude
        master_avg_mag[master_filters == filt] = average_mag
        master_amp[master_filters == filt] = amplitude
        master_sigma[master_filters == filt] = sigma

    if plot_save == 1:

        ax.set_xlim((-0.2,1.2))
        ax.set_ylim((18,9))
        ax.set_xlabel('Phase')
        ax.set_ylabel('Mag + offset')
    #    ax.set_title(phased_lcv_file)
        plot_file = re.sub('\.phased', '-fit.pdf', phased_lcv_file)
        figtosave.savefig(plot_file)




    return master_filters, master_avg_mag, master_amp, master_sigma
