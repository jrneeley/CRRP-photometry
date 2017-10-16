import numpy as np
import daophot_setup
import daophot
import re
import allstar
import os
import sys

def mosaic_phot(dao_dir, opt_dir, mosaic_fits, channel, exptime):

    curr_dir = os.getcwd()
    os.chdir(curr_dir+'/DeepMosaic')
    print 'Changed directory to DeepMosaic\n'
# headers of deep mosaics are wrong, so must input right numbers here.
    # convert deep mosaic to counts
    if channel == 'I1' : fluxconv = 0.1257
    if channel == 'I2' : fluxconv = 0.1447
    daophot_setup.spitzer_flux2dn(mosaic_fits, exptime=float(exptime), fluxconv=fluxconv)
    # copy OPT files to current dir
    daophot_setup.set_opt_files_mosaic(opt_dir, channel, exptime, warm=1)
    # Find stars and do initial aperture photometry on deep mosaic
    mosaic_dn = re.sub('.fits', '_dn.fits', mosaic_fits)
    daophot.deep_mosaic_phot(dao_dir, mosaic_dn)
    # You must choose your PSF stars from the .coo by hand
    ready = raw_input('Pick PSF stars by hand, type continue when ready, or stop to exit: ')
    if ready == 'stop': sys.exit()
    if ready == 'continue':
        print '\nMaking the lst file...'
    # Make .lst file for chosen PSF stars
    psf_stars = np.loadtxt(channel+'-psf.reg', dtype=np.dtype([('x', float), ('y', float)]))
    coo_file = re.sub('.fits', '.coo', mosaic_dn)
    lst_file = re.sub('.fits', '.lst', mosaic_dn)
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
    # measure the PSF
    daophot.find_psf(dao_dir, mosaic_dn)
    repeat = raw_input('Do you want to remove any PSF stars? [y/n]: ')
    if repeat == 'y':
        ready = raw_input('Delete from .lst file and type continue when ready: ')
        daophot.find_psf(dao_dir, mosaic_dn)
    if repeat == 'n':
        print 'PSF ok...continuing...\n'
    print 'Running allstar on deep mosaic...'
    allstar.allstar_deep(dao_dir, mosaic_dn)
    sub_img = re.sub('dn.fits', 'dns.fits', mosaic_dn)
    print '\nRunning FIND/PHOT on subtracted image...'
    daophot.deep_mosaic_phot(dao_dir, sub_img)
    print 'Appending new stars to star list...'
    daophot.append(dao_dir, mosaic_dn)
    print 'Running PHOT on new star list...'
    daophot.deep_mosaic_phot2(dao_dir, mosaic_dn)
    print 'Running allstar on deep mosaic again...'
    allstar.allstar_deep(dao_dir, mosaic_dn)
    print 'Finished with deep mosaic. \n'

    ## need to update .lst file with new ID numbers/positions 
