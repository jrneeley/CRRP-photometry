import matplotlib.pyplot as mp
import numpy as np
import optical
import lightcurves
import sys
import os

target = sys.argv[1]

parent_dir = os.getcwd()

#optical_folder = '/Users/jrneeley/CRRP/OpticalCatalogs/'
optical_dir = parent_dir+'/OpticalCatalogs/'
target_dir = parent_dir+'/'+target+'/'

# Read in variable names and periods from Clement catalog
dtype1 = np.dtype([('id', 'S10'), ('period', float)])
data = np.loadtxt(target_dir+target+'-clement.txt', dtype=dtype1, usecols=(0,3))

# Identify different datsets for this cluster
datasets, colors = optical.compile_datasets(target_dir, old=0)

print '\n\nStar  Period_old  Period_new'

for ind, lcv in enumerate(data['id']):

# Open file to save periods
    if ind == 0:
        f_handle = open(target_dir+'periods.txt', 'w')
    else:
        f_handle = open(target_dir+'periods.txt', 'a')

    lcv_file = target_dir+'lcvs/optical/'+target+lcv+'.lcv'
    clement_period = data['period'][data['id'] == lcv]
    try:
        U, B, V, R, I = lightcurves.read_optical_lcv(lcv_file)
        new_guess = lightcurves.period_search_LS(V, lcv, plot_save=1, data_dir=target_dir)
        new_period = lightcurves.period_search(V, new_guess, lcv, second_band=B, search_window=0.00005)

        print '%10s %0.4f %0.8f' % (lcv, clement_period, new_period)
        f_handle.write('%10s %0.4f %0.8f\n' % (lcv, clement_period, new_period))
        lightcurves.plot_phased_optical_lcv(U, B, V, R, I, new_period, lcv, datasets, plot_save=1,data_dir=target_dir, colors=colors )
    except:
        new_period = np.nan
        print '%10s %0.4f %0.8f' % (lcv, clement_period, new_period)
        f_handle.write('%10s %0.4f %0.8f\n' % (lcv, clement_period, new_period))

    # Close the periods file
    f_handle.close()
