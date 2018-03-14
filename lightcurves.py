import numpy as np
import re
import glob
import dao
from astropy.io import fits
import matplotlib.pyplot as mp
from astropy.stats import LombScargle
import plotting_utilities
from scipy import stats
from scipy.interpolate import interp1d
from scipy import signal
from matplotlib import gridspec
import os.path
import peakutils
from astropy.stats import sigma_clip
from IPython import display
import config

def make_lcv(target, channels, stars, dao_ids):

    data_dir = config.top_dir+target

    lcv_folder = data_dir+'/lcvs/mir/'
    img_folder = data_dir+'/mosaics/'

    for channel in channels:
        file_list = glob.glob(img_folder+channel+'*.cal')
        if len(file_list) == 0:
            print 'Using uncalibrated files!'
            file_list = glob.glob(img_folder+channel+'*.alf')
        phot_data = np.zeros(len(dao_ids), dtype=[('filter', 'S2', len(file_list)),
            ('id', 'S8'), ('aor', int, len(file_list)), ('mjd', float, len(file_list)),
            ('f_num', int, len(file_list)), ('x', float, len(file_list)),
            ('y', float, len(file_list)), ('psf_mag', float, len(file_list)),
            ('psf_err', float, len(file_list))])


        phot_data['id'] = dao_ids

        for ind in range(0,len(file_list)):

            alf_data = dao.read_alf(file_list[ind])
            fits_file = re.sub(".cal",".fits", file_list[ind])
            if fits_file == file_list[ind]:
                fits_file = re.sub(".alf", ".fits", file_list[ind])

            hdulist = fits.open(fits_file, mode='update')
            prihdr = hdulist[0].header
            mjd = prihdr['mjd_obs']

            for ind2 in range(0,len(dao_ids)):
                alf_match = np.argwhere(alf_data['id'] == dao_ids[ind2])

                trash, aor_num, trash = file_list[ind].split('_')
                frame_num = 1

                phot_data['aor'][ind2,ind] = int(aor_num)
                phot_data['f_num'][ind2, ind] = int(frame_num)
                phot_data['mjd'][ind2,ind] = mjd
                phot_data['filter'][ind2, ind] = channel
                if len(alf_match):
                    phot_data['x'][ind2,ind] = alf_data['x'][alf_match]
                    phot_data['y'][ind2,ind] = alf_data['y'][alf_match]
                    phot_data['psf_mag'][ind2, ind] = alf_data['mag'][alf_match]
                    phot_data['psf_err'][ind2, ind] = alf_data['err'][alf_match]
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
                    #sigma clip bin
                    filtered_mags = sigma_clip(dither_mags, sigma=3, iters=1)
                    dispersion = np.nanstd(filtered_mags)
                    mean_mag = np.nanmean(filtered_mags)
                    mag[ind] = mean_mag
                    mjd[ind] = np.nanmean(dither_mjds)
                    err[ind] = dispersion
                #    mag[ind] = np.nanmean(dither_mags[abs(dither_mags - mean_mag) < 2*dispersion])
                #    mjd[ind] = np.nanmean(dither_mjds[abs(dither_mags - mean_mag) < 2*dispersion])
                #    err[ind] = np.nanstd(dither_mags[abs(dither_mags - mean_mag) < 2*dispersion])
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

def phase_lcv_all_bands(target, lcv, period, T0, optical_lcv=0, nir_lcv=0, mir_lcv=0, bin_mir=0, data_dir='', old=0, mosaics=0):

    fig = mp.figure(figsize=(8,10))
    path_to_optical = data_dir+'lcvs/optical/'
    path_to_nir = data_dir+'lcvs/nir/'
    if mosaics == 0: path_to_mir = data_dir+'lcvs/mir/'
    if mosaics == 1: path_to_mir = data_dir+'mosaic_lcvs/'

    if os.path.isfile(path_to_optical+target+lcv):
        optical_lcv = 1
    if os.path.isfile(path_to_nir+target+lcv):
        nir_lcv = 1
    if os.path.isfile(path_to_mir+lcv):
        mir_lcv = 1
    master_filters = np.array(['U', 'B', 'V', 'R', 'I', 'J', 'H', 'K',
        'I1', 'I2'], dtype='S2')
    master_markers = np.array(['P', 'v', 'D', '>', 'x', 'p', 'd', '^', 'o', 's'])
    master_offset = np.array([1.0, 0.5, 0.0, -0.25, -0.5, 0.2, 0.5, 0.8, -1.0, -1.5 ])
    master_colors = np.array(['xkcd:violet', 'xkcd:periwinkle', 'xkcd:sapphire',
        'xkcd:sky blue', 'xkcd:emerald', 'xkcd:avocado', 'xkcd:goldenrod',
        'xkcd:orange', 'xkcd:pink', 'xkcd:scarlet'])

    if optical_lcv == 1:
        if mosaics == 0: phased_file = re.sub('\.lcv', '.phased', lcv)
        if mosaics == 1: phased_file = re.sub('\.lcv', '.mphased', lcv)
        U, B, V, R, I = read_optical_lcv(path_to_optical+target+lcv, old=old)


        Uph = np.mod((U[2]-T0)/period, 1)
        Bph = np.mod((B[2]-T0)/period, 1)
        Vph = np.mod((V[2]-T0)/period, 1)
        Rph = np.mod((R[2]-T0)/period, 1)
        Iph = np.mod((I[2]-T0)/period, 1)

        ## Add to data to plot
        mp.errorbar(Uph, U[0]+master_offset[0], yerr=U[1], fmt=master_markers[0], color=master_colors[0])
        mp.errorbar(Bph, B[0]+master_offset[1], yerr=B[1], fmt=master_markers[1], color=master_colors[1])
        mp.errorbar(Vph, V[0]+master_offset[2], yerr=V[1], fmt=master_markers[2], color=master_colors[2])
        mp.errorbar(Rph, R[0]+master_offset[3], yerr=R[1], fmt=master_markers[3], color=master_colors[3])
        mp.errorbar(Iph, I[0]+master_offset[4], yerr=I[1], fmt=master_markers[4], color=master_colors[4])

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

    if nir_lcv == 1:
        if mosaics == 0: phased_file = re.sub('\.lcv', '.phased', lcv)
        if mosaics == 1: phased_file = re.sub('\.lcv', '.mphased', lcv)

        J, H, K = read_nir_lcv(path_to_nir+target+lcv, old=old)

        Jph = np.mod((J[2]-T0)/period, 1)
        Hph = np.mod((H[2]-T0)/period, 1)
        Kph = np.mod((K[2]-T0)/period, 1)


        ## Add to data to plot
        mp.errorbar(Jph, J[0]-master_offset[5], yerr=J[1], fmt=master_markers[5], color=master_colors[5])
        mp.errorbar(Hph, H[0]-master_offset[6], yerr=H[1], fmt=master_markers[6], color=master_colors[6])
        mp.errorbar(Kph, K[0]-master_offset[7], yerr=K[1], fmt=master_markers[7], color=master_colors[7])

        data_save = np.array(zip(np.repeat('J', len(J[2])), J[2], Jph, J[0], J[1], J[3]),
            dtype=[('c1', 'S2'), ('c2', float), ('c3', float), ('c4', float),
            ('c5', float), ('c6', 'S30')])
        np.savetxt(f_handle, data_save, fmt='%s %10.4f %8.6f %6.3f %5.3f %s')
        data_save = np.array(zip(np.repeat('H', len(H[2])), H[2], Hph, H[0], H[1], H[3]),
            dtype=[('c1', 'S2'), ('c2', float), ('c3', float), ('c4', float),
            ('c5', float), ('c6', 'S30')])
        np.savetxt(f_handle, data_save, fmt='%s %10.4f %8.6f %6.3f %5.3f %s')
        data_save = np.array(zip(np.repeat('K', len(K[2])), K[2], Kph, K[0], K[1], K[3]),
            dtype=[('c1', 'S2'), ('c2', float), ('c3', float), ('c4', float),
            ('c5', float), ('c6', 'S30')])


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
                offset = master_offset[8]
                marker = master_markers[8]
                color = master_colors[8]
            if filt == 'I2':
                offset = master_offset[9]
                marker = master_markers[9]
                color = master_colors[9]
            mp.errorbar(phase, mag+offset, yerr=err, fmt=marker, color=color)
            ## Add data to file
            data_save = np.array(zip(band, mjd, phase, mag, err, source),
                dtype=[('c1', 'S2'), ('c2', float), ('c3', float), ('c4', float),
                ('c5', float), ('c6', 'S30')])
            np.savetxt(f_handle, data_save, fmt='%2s %10.4f %8.6f %6.3f %5.3f %s')

    max_mag = np.mean(V[0])+3
    min_mag = np.mean(V[0])-5
    mp.ylim((max_mag, min_mag))
    mp.xlabel('Phase')
    mp.ylabel('Mag + offset')
    mp.title(lcv)
    if mosaics == 0: plot_file = re.sub('\.phased', '-ph.pdf', phased_file)
    if mosaics == 1: plot_file = re.sub('\.mphased', '-mph.pdf', phased_file)
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

def read_optical_lcv_old(lcv_file, old=0):

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

def read_optical_lcv(lcv_file, old=0):

    dtype1 = np.dtype([('mag', float), ('err', float), ('filter', 'a2'),
        ('year', int), ('day', float), ('source', 'S35')])

    if old == 1:
        data = np.loadtxt(lcv_file, dtype=dtype1, usecols=(0,1,2,4,5,8))
    if old == 0:
        data = np.loadtxt(lcv_file, dtype=dtype1, usecols=(0,1,2,4,5,10))

    data['filter'][data['filter'] == '1'] = 'V'
    data['filter'][data['filter'] == '2'] = 'B'
    data['filter'][data['filter'] == '3'] = 'I'
    data['filter'][data['filter'] == '4'] = 'R'
    data['filter'][data['filter'] == '5'] = 'U'
    data['day'] = data['year']*1000. + data['day'] - 2400000.5

    return data

def read_nir_lcv(lcv_file, old=0):

    dtype1 = np.dtype([('mag', float), ('err', float), ('filter', 'a1'),
        ('year', int), ('day', float), ('source', 'S35')])
    if old == 1:
        data = np.loadtxt(lcv_file, dtype=dtype1, usecols=(0,1,2,4,5,8))
    if old == 0:
        data = np.loadtxt(lcv_file, dtype=dtype1, usecols=(0,1,2,4,5,10))

    data['filter'][data['filter'] == '1'] = 'J'
    data['filter'][data['filter'] == '2'] = 'H'
    data['filter'][data['filter'] == '3'] = 'K'
    data['day'] = data['year']*1000. + data['day'] - 2400000.5

    return data

def read_mir_lcv(lcv_file, old=0):

    dtype1 = np.dtype([('filter', 'a2'), ('aor', 'a8'), ('frame', int),
        ('mjd', float), ('x', float), ('y', float), ('mag', float),
        ('err', float)])

    data = np.loadtxt(lcv_file, dtype=dtype1)

    return data

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
    data_dir='', error_threshold=0.05, colors=None):


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

def period_search(first_band, best_period, second_band=None,
    plot_save=0, error_threshold=0.1, search_window=0.0002, plot=0):

    x = np.array(first_band[2][first_band[1] < error_threshold], dtype=float)
    y = np.array(first_band[0][first_band[1] < error_threshold], dtype=float)
    er = np.array(first_band[1][first_band[1] < error_threshold], dtype=float)
    if second_band is not None:
        x2 = np.array(second_band[2][second_band[1] < error_threshold], dtype=float)
        y2 = np.array(second_band[0][second_band[1] < error_threshold], dtype=float)
        er2 = np.array(second_band[1][second_band[1] < error_threshold], dtype=float)

    # Calculate required precision
    delta_time = np.max(x) - np.min(x)
    approx_p = np.round(best_period, 1)
    n_cycles = np.floor(delta_time/approx_p)
    max_precision = n_cycles * approx_p / (n_cycles - 0.01) - approx_p
    order = np.ceil(np.abs(np.log10(max_precision)))
    precision = 10**order
    best_period = np.round(best_period, decimals=int(order))

    grid_num = search_window*precision

    if grid_num > 100000:
        grid_num = 100000
    min_period = best_period - search_window/2
    max_period = best_period + search_window/2
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
    order = np.argsort(avg_std)
    new_period = periods[order[0]]
    best_std = avg_std[order[0]]
    mp.plot(periods, avg_std, 'ro')
    mp.axvline(new_period)

    if plot != 0:
        mp.show()
    mp.close()

    return new_period

# Use to do one round of LombScargle and identify possible periods for RRL star
def period_search_LS(data, name, min_period = 0.2, max_period=1.0,
                    error_threshold=0.05, verbose=0, plot_save=0, data_dir=''):

    x1 = np.array(data[2][data[1] < error_threshold], dtype=float)
    y1 = np.array(data[0][data[1] < error_threshold], dtype=float)
    er1 = np.array(data[1][data[1] < error_threshold], dtype=float)

    freq_max = 1/(min_period)
    freq_min = 1/(max_period)
    frequency, power = LombScargle(x1, y1, er1).autopower(minimum_frequency=freq_min,
                        maximum_frequency=freq_max )

    # Calculate noise level
    median_power = np.median(power)

    # Find peaks
    indices = peakutils.indexes(power, thres=0.5, min_dist=5000 )
    candidate_periods = 1/frequency[indices]
    best_frequency = frequency[np.argmax(power)]
    best_period = 1/best_frequency

    fig = mp.figure(figsize=(12, 6))
    ax1 = mp.subplot2grid((1,2), (0,0))
    ax2 = mp.subplot2grid((1,2), (0,1))
    ax1.plot(1/frequency, power)
    ax1.plot(1/frequency[indices], power[indices], 'rx')
    alias_freq = np.array([best_frequency+1, best_frequency-1, best_frequency+2, best_frequency-2])
    alias_power = np.repeat(np.max(power), 4)
    ax1.plot(1/alias_freq, alias_power, 'kx')
    ax1.set_xlim(min_period, max_period)
    ax1.set_xlabel('Period (days)')
    ax1.set_ylabel('Power')


    # Calculate SNR of peaks
    snr = power[indices]/median_power
    snr_best = power[np.argmax(power)]/median_power
    if verbose == 1:
        print candidate_periods
        print snr
    t_fit = np.linspace(0,1)
    y_fit = LombScargle(x1, y1, er1).model(t_fit, best_frequency)

    phase_data = np.mod(x1*best_frequency, 1)
    ax2.plot(phase_data, y1, 'o')
    ax2.set_ylim(np.max(y1)+0.05, np.min(y1)-0.05)
    ax2.set_xlabel('Phase')
    ax2.set_ylabel('mag')
    if plot_save == 0 :
        mp.show()
    if plot_save == 1:
        mp.savefig(data_dir+'lcvs/optical/'+name+'-periodogram.pdf')
    mp.close()
#    return candidate_periods
    return best_period, snr_best

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

def gloess(phased_lcv_file, clean=1, smoothing_params=None, ask=0, filters='all', master_plot=0, mosaics=0):

    # set to 1 if you want to save a single figure for each star with all data
    if master_plot == 1:
        figtosave = mp.figure(figsize=(8,10))
        ax = figtosave.add_subplot(111)

    master_filters = np.array(['U', 'B', 'V', 'R', 'I', 'J', 'H', 'K',
        'I1', 'I2'], dtype='S2')
    master_markers = np.array(['P', 'v', 'D', '>', 'x', 'p', 'd', '^', 'o', 's'])
    master_offset = np.array([1.0, 0.5, 0.0, -0.25, -0.5, -0.7, -0.9, -1.1, -1.0, -1.5 ])
    master_colors = np.array(['xkcd:violet', 'xkcd:periwinkle', 'xkcd:sapphire',
        'xkcd:sky blue', 'xkcd:emerald', 'xkcd:avocado', 'xkcd:goldenrod',
        'xkcd:orange', 'xkcd:pink', 'xkcd:scarlet'])
    # read in the phased light curve file
    dtype1 = np.dtype([('filter', 'S2'), ('mjd', float), ('phase', float), ('mag', float), ('err', float)])
    data = np.loadtxt(phased_lcv_file, dtype=dtype1, usecols=(0,1,2,3,4))

    # which filters are available
    if filters == 'all': filters = np.unique(data['filter'])
    num_filters = len(filters)


    if smoothing_params is None: smoothing_params = np.repeat(1.0, 10)
        #smoothing_params = np.array([0.1, 0.1, 0.1, 0.1, 0.1, 0.2, 0.2, 0.2, 0.2, 0.2])

    master_avg_mag = np.zeros(10)
    master_amp = np.zeros(10)
    master_sigma = np.zeros(10)
    master_avg_mag_er = np.zeros(10)

    for filt in filters:

        phase = data['phase'][data['filter'] == filt]
        mag = data['mag'][data['filter'] == filt]
        err = data['err'][data['filter'] == filt]
        mjd = data['mjd'][data['filter'] == filt]

        mjd = mjd[~np.isnan(mag)]
        phase = phase[~np.isnan(mag)]
        err = err[~np.isnan(mag)]
        mag = mag[~np.isnan(mag)]

        if clean == 1:
            # remove data with large error bars
            filtered_err = sigma_clip(err, sigma=3, iters=2)
            mag = mag[~filtered_err.mask]
            phase = phase[~filtered_err.mask]
            err = err[~filtered_err.mask]

        # skip this band if we don't have enough observations or phase coverage
        n_obs = len(mag)
        if (filt == 'I1') or (filt == 'I2'): n_obs = 99
        delta_phase = np.max(phase) - np.min(phase)
        if (n_obs < 30) or (delta_phase < 0.7):
            continue
        sigma = float(smoothing_params[master_filters == filt])

        if sigma == 1.0:
            hist, bins = np.histogram(phase, bins='auto')
            #print 1./len(bins), 1./len(mag)*2
            sigma = 1./len(bins)

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

            # interpolate
            for ind, step in enumerate(x):

                dist = phase_copy - step
                closest_phase = np.min(np.abs(dist))
                if closest_phase > sigma: sigma = closest_phase*5
                weight = err_copy * np.exp(dist**2/sigma**2)
                fit = np.polyfit(phase_copy, mag_copy, 2, w=1/weight)
                smoothed_mag[ind] = fit[2] + fit[1]*step + fit[0]*step**2
            #    print step, np.min(np.abs(dist))

            smoothed_mag_copy = np.tile(smoothed_mag, 5)
            x_copy = np.concatenate((x-2, x-1, x, x+1, x+2))

            figshow = mp.figure(figsize=(8,6))
            ax2 = figshow.add_subplot(111)
            ax2.errorbar(phase_copy, mag_copy, yerr=err_copy, fmt='o', zorder=1)
            ax2.plot(x_copy, smoothed_mag_copy, 'r-')
            ax2.set_ylim((np.max(mag)+0.2, np.min(mag)-0.2))
            ax2.set_xlim((-0.2, 1.2))
            ax2.set_xlabel('Phase')
            ax2.set_ylabel(filt+' mag')
            display.display(mp.gcf())
            display.clear_output(wait=True)

        #    if smoothing_params[master_filters == filt] != 1.0:
            if ask == 0:
                happy = 'y'
                if mosaics == 0: plot_file = re.sub('.phased', '-'+filt+'fit.pdf', phased_lcv_file)
                if mosaics == 1: plot_file = re.sub('.mphased', '-'+filt+'mfit.pdf', phased_lcv_file)
                mp.savefig(plot_file, format='pdf')
                continue
            check = raw_input('Are you happy with this fit? [y/n]: ')
            if check == 'n':
                sigma = input('Enter new smoothing parameter: ')
            if check == 'y':
                happy = 'y'

        mp.close()
        # Derive light curve parameters
        flux = 99*np.power(10,-smoothed_mag/2.5)
        average_flux = np.mean(flux)
        average_mag = -2.5*np.log10(average_flux/99)

        amplitude = np.max(smoothed_mag) - np.min(smoothed_mag)
        ph_max = x[smoothed_mag == np.min(smoothed_mag)]
        ph_min = x[smoothed_mag == np.max(smoothed_mag)]

        err_fit = amplitude/np.sqrt(12*len(err))
        average_mag_err = np.sqrt(np.sum(err**2)/len(err)**2 + err_fit**2)

        # determine the epoch of maximum using the V band data
        if filt == 'V':
            T0 = ph_max

        #print filt, average_mag, amplitude

        master_avg_mag[master_filters == filt] = average_mag
        master_amp[master_filters == filt] = amplitude
        master_sigma[master_filters == filt] = sigma
        master_avg_mag_er[master_filters == filt] = average_mag_err

        if mosaics == 0: fit_file = re.sub('.phased', '.fit', phased_lcv_file)
        if mosaics == 1: fit_file = re.sub('.mphased', '.mfit', phased_lcv_file)
        if filt == filters[0]:
            f = open(fit_file, 'w')
        else:
            f = open(fit_file, 'a')
        dtype= np.dtype([('filt', 'S2'), ('x', float), ('mag', float)])
        data_save = np.array(zip(np.repeat(filt, len(x)), x, smoothed_mag), dtype=dtype)
        np.savetxt(f, data_save, fmt='%2s %.4f %2.3f')
        f.close()

        if master_plot == 1:

            offset = master_offset[master_filters == filt]
            marker = master_markers[master_filters == filt][0]
            color = master_colors[master_filters == filt][0]
            ax.errorbar(phase, mag+offset, yerr=err, fmt=marker, color=color, zorder=1)
            ax.plot(x_copy, smoothed_mag_copy+offset, 'k-')

    if master_plot == 1:

        max_mag = np.nanmean(data['mag'][data['filter'] == 'V'])+3
        min_mag = np.nanmean(data['mag'][data['filter'] == 'V'])-5
        ax.set_ylim((max_mag, min_mag))
        ax.set_xlim((-0.2, 2.0))
        ax.set_xlabel('Phase')
        ax.set_ylabel('Mag + offset')
        if mosaics == 0: plot_file = re.sub('\.phased', '-fit.pdf', phased_lcv_file)
        if mosaics == 1: plot_file = re.sub('\.mphased', '-mfit.pdf', phased_lcv_file)
        mp.savefig(plot_file)

    return master_filters, master_avg_mag, master_avg_mag_er, master_amp, master_sigma
