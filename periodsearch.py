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

f = open('periods.txt', 'w')
print '\n\nStar  Period_old  Period_new'


for ind, lcv in enumerate(data['id']):

    # Open file to save periods
    if ind == 0:
        f_handle = open('periods.txt', 'w')
    else:
        f_handle = open('periods.txt', 'a')

    lcv_file = optical_folder+target+'lcvs/'+target+'V'+str(lcv)+'.lcv'
#    try:
    U, B, V, R, I = lightcurves.read_optical_lcv(lcv_file, old=0)
    period = lightcurves.period_search(V, data['period'][ind], 'V'+str(lcv), plot_save=1)
    print '%4s %0.4f %0.8f' % ('V'+str(lcv), data['period'][ind], period)
#    data_save = np.array(zip([['V'+str(lcv)], data['period'][ind], [period]]),
#        dtype=[('id', 'S4'), ('clement', float), ('new_period', float)])
#    np.savetxt(f_handle, ['V'+str(lcv), data['period'][ind], period], comments='', fmt='%s %0.4f %0.8f')
    f_handle.write('%4s %0.4f %0.8f\n' % ('V'+str(lcv), data['period'][ind], period))
    lightcurves.plot_phased_optical_lcv(U, B, V, R, I, period, 'V'+str(lcv), datasets, plot_save=1)
#    except:
#        print '%4s %0.4f %s' % ('V'+str(lcv), data['period'][ind], 'Not found')
    # Close the periods file
    f_handle.close()
