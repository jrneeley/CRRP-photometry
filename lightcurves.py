import numpy as np
import re
import glob
import read_dao
from astropy.io import fits
import matplotlib.pyplot as mp

def make_lcv(stars, dao_ids):

    alf_list = glob.glob('all/I1*.alf')

    phot_data = np.zeros(len(dao_ids), dtype=[('id', 'S8'), ('aor', 'i8', len(alf_list)),
        ('mjd', float, len(alf_list)), ('f_num', 'i2', len(alf_list)), ('x', float, len(alf_list)),
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
            phot_data['aor'][ind2,ind] = int(aor_num)
            phot_data['f_num'][ind2, ind] = int(frame_num)
            phot_data['mjd'][ind2,ind] = mjd
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
    print 'Writing files...'
    for ind in range(0,len(dao_ids)):

        np.savetxt(stars[ind]+'.lcv',
            np.c_[phot_data['aor'][ind], phot_data['f_num'][ind],phot_data['mjd'][ind],
            phot_data['x'][ind],phot_data['y'][ind], phot_data['psf_mag'][ind],
            phot_data['psf_err'][ind]],
            comments='', fmt='%8i %2i %10.4f %7.3f %7.3f %6.3f %6.4f')

def compare_lcv(lcv_file):

    dtype1 = np.dtype([('mjd', float), ('mag', float), ('err', float)])
    data = np.loadtxt(lcv_file, dtype=dtype1, usecols=(2,5,6))

    file2 = 'Peter/NGC6121'+lcv_file
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
