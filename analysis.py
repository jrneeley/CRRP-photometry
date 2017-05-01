import read_dao
import numpy as np

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
    weights = 1/errs[ind][good]**2
    avgs[ind] = np.average(mags[ind][good], weights=weights)
    avg_err[ind] = np.sqrt(1/np.sum(weights))


data = np.c_[id_num, x, y, avgs, avg_err]
np.savetxt('test.txt',data, fmt='%8i %9.3f %9.3f %6.3f %5.3f')
