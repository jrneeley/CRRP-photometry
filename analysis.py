import read_dao
import numpy as np
#import statsmodels.api as sm
import lightcurves
import glob
import sys
import variables

target = sys.argv[1]

#dtype1 = np.dtype([('stars', 'S3'), ('dao_ids', 'i8')])
#data = np.loadtxt('PeterIDs.txt', dtype=dtype1)

#lightcurves.make_lcv(['I1'], data['stars'], data['dao_ids'])

#lcv_list = glob.glob('lcvs/*.phased')

#for lcv in lcv_list:
#    print lcv
#    lightcurves.compare_phased_lcv(lcv)

#dtype1= np.dtype([('id', 'S3'), ('period', float), ('t0', float)])
#data = np.loadtxt('M4_RRL.info', dtype=dtype1, usecols=(0,1,2))

#for ind, star in enumerate(data['id']):
#    try:
#        lightcurves.phase_lcv('lcvs/'+star+'.lcv', data['period'][ind], data['t0'][ind])

#    except:
#        print 'Star '+ star + ' not found.'
folder = '/Users/jrneeley/CRRP/OpticalCatalogs/'
variables.find_variables_by_coord(folder, target)
sys.exit()

ids, raw_phot = read_dao.read_raw('optical2_alf.raw')

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
np.savetxt('catalog.txt',data, fmt='%8i %9.3f %9.3f %6.3f %5.3f')
