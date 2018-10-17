import numpy as np
import re
import sys

def make_lst(channel, image):

    # You must choose your PSF stars from the .coo by hand
    ready = raw_input('Pick PSF stars by hand. Load the image into ds9, and select PSF stars. Save pixel coordinates as a regions file in xy format, named [channel]-psf.reg. \n Type continue when ready, or stop to exit: ')
    if ready == 'stop': sys.exit(0)
    if ready == 'continue':
        print '\nMaking the .lst file...'

    # Make .lst file for chosen PSF stars
    psf_stars = np.loadtxt(channel+'-psf.reg', dtype=np.dtype([('x', float), ('y', float)]))
    coo_file = re.sub('.fits', '.coo', image)
    lst_file = re.sub('.fits', '.lst', image)
    f = open(coo_file, 'r')
    f2 = open(lst_file, 'w')
    for ind in range(3):
        line = f.readline()
        f2.write(line)
    f.close()
    f2.close()
    dtype = np.dtype([('id', int), ('x', float), ('y', float), ('c1', float),
        ('c2', float), ('c3', float), ('c4', float)])
    all_stars = np.loadtxt(coo_file, dtype=dtype, skiprows=3)
    f_handle = open(lst_file, 'a')
    for x, y in zip(psf_stars['x'], psf_stars['y']):
        line = all_stars[all_stars['y'] == y]
        np.savetxt(f_handle, line, fmt='%8i %8.2f %8.2f %8.3f %8.3f %8.3f %8.3f')
    f_handle.close()
    print '.lst file complete.'
