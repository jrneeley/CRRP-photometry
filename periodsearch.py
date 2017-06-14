import matplotlib.pyplot as mp
import numpy as np
import optical
import lightcurves
import sys

target = sys.argv[1]

optical_folder = '/Users/jrneeley/CRRP/OpticalCatalogs/'

# Read in variable names and periods from Clement catalog
dtype1 = np.dtype([('id', int), ('period', float)])
data = np.loadtxt(target+'-clement.txt', dtype=dtype1, usecols=(0,3))

# Identify different datsets for this cluster
datasets = optical.compile_datasets(optical_folder, target, old=0)

print '\n\nStar  Period_old  Period_new'


for ind, lcv in enumerate(data['id']):

    lcv_file = optical_folder+target+'lcvs/'+target+'V'+str(lcv)+'.lcv'
    try:
        U, B, V, R, I = lightcurves.read_optical_lcv(lcv_file, old=0)
        period = lightcurves.period_search(V, data['period'][ind], 'V'+str(lcv), plot_save=1)
        print '%4s %0.4f %0.8f' % ('V'+str(lcv), data['period'][ind], period)
        lightcurves.plot_phased_optical_lcv(U, B, V, R, I, period, 'V'+str(lcv), datasets, plot_save=1)
    except:
        print '%4s %0.4f %s' % ('V'+str(lcv), data['period'][ind], 'Not found')
