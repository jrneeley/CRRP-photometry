import matplotlib.pyplot as mp
import numpy as np
import optical
import lightcurves
import sys

target = sys.argv[1]

optical_folder = '/Users/jrneeley/CRRP/OpticalCatalogs/'

refine = 1
# Read in variable names and periods from Clement catalog
dtype1 = np.dtype([('id', 'S10'), ('period', float)])
data = np.loadtxt(target+'-clement.txt', dtype=dtype1, usecols=(0,3))

# Identify different datsets for this cluster
datasets = optical.compile_datasets(optical_folder, target, old=0)

#refine = ['V9']#, 'V12', 'V18', 'V42', 'V51', 'V52', 'V58', 'V66', 'V67', 'V80']
if refine == 0:

    print '\n\nStar  Period_old  Period_new'

    for ind, lcv in enumerate(data['id']):

    # Open file to save periods
        if ind == 0:
            f_handle = open('periods.txt', 'w')
        else:
            f_handle = open('periods.txt', 'a')

        lcv_file = optical_folder+target+'lcvs/'+target+lcv+'.lcv'
        period = data['period'][data['id'] == lcv]
        try:
            U, B, V, R, I = lightcurves.read_optical_lcv(lcv_file, old=0)
            new_period = lightcurves.period_search_hybrid(V, period, lcv, plot_save=1)
            print '%10s %0.4f %0.8f' % (lcv, period, new_period)
            f_handle.write('%10s %0.4f %0.8f\n' % (lcv, data['period'][ind], new_period))
            lightcurves.plot_phased_optical_lcv(U, B, V, R, I, new_period, lcv, datasets, plot_save=1)
        except:
            print '%10s %0.4f %s' % (lcv, period, 'Not found')
    # Close the periods file
        f_handle.close()

if refine == 1:

    print '\n\nStar  Period_old  Period_new'

#    refine = ['V109=NV14', 'V116=NV5', 'V12', 'V120=NV7', 'V129',
#        'V130', 'V141', 'V166=ZK22', 'V167', 'V172=ZK13', 'V174=ZK10', 'V21', 'V39', 'V41',
#        'V53', 'V59', 'V60', 'V67', 'V69', 'V81', 'ZK34', 'ZK74']
    refine = ['V67']
    for lcv in refine:
        #f = open('periods.txt', 'w')
        lcv_file = optical_folder+target+'lcvs/'+target+lcv+'.lcv'
        period = data['period'][data['id'] == lcv]

        U, B, V, R, I = lightcurves.read_optical_lcv(lcv_file, old=0)
        new_period = lightcurves.period_search_hybrid(V, period, lcv, plot_save=1, search_window=0.2)
        print '%10s %0.4f %0.8f' % (lcv, period, new_period)
        lightcurves.plot_phased_optical_lcv(U, B, V, R, I, new_period, lcv, datasets, plot_save=1)
