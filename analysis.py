import read_dao
import numpy as np
#import statsmodels.api as sm
import lightcurves
import glob
import sys

#dtype1 = np.dtype([('stars', 'S3'), ('dao_ids', 'i8')])
#data = np.loadtxt('PeterIDs.txt', dtype=dtype1)

#lightcurves.make_lcv(data['stars'], data['dao_ids'])

lcv_list = glob.glob('*.lcv')
peter_list = glob.glob('Peter/*.lcv')

for lcv in lcv_list:
    lightcurves.compare_lcv(lcv)

sys.exit()

ids, raw_phot = read_dao.read_raw('optical_alf.raw')

id_num = np.zeros(len(ids), dtype=int)
x = np.zeros(len(ids), dtype=float)
y = np.zeros(len(ids), dtype=float)
mags = np.zeros((len(ids), (len(raw_phot[0])-7)/2), dtype=float)
errs = np.zeros((len(ids), (len(raw_phot[0])-7)/2), dtype=float)
avgs = np.zeros(len(ids), dtype=float)
avg_err = np.zeros(len(ids), dtype=float)

for ind, star in enumerate(ids):
    id_num[ind] = raw_phot[ind][0]
    x[ind] = raw_phot[ind][1]
    y[ind] = raw_phot[ind][2]
    mags[ind] = raw_phot[ind][5:-2:2]
    errs[ind] = raw_phot[ind][6:-2:2]

    good = mags[ind] < 90
    data = mags[ind][good]
    weights = 1/errs[ind][good]**2
    weights = weights[abs(data-np.mean(data)) < 2.5*np.std(data)]
    data = data[abs(data-np.mean(data)) < 2.5*np.std(data)]
    avgs[ind] = np.average(data, weights=weights)
    avg_err[ind] = np.sqrt(1/np.sum(weights))

#    sm.robust.scale.huber(mags[ind][good])

data = np.c_[id_num, x, y, avgs, avg_err]
np.savetxt('test3.txt',data, fmt='%8i %9.3f %9.3f %6.3f %5.3f')
