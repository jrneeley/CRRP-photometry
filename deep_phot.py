import numpy as np
import matplotlib.pyplot as mp
import daophot_setup
import daophot
import re
import allstar
import os
import sys
import optical
import coordinates
#sys.path.insert(0, '/Users/Jill/python/daophot-tools/')
import daomatch
import daomaster
import read_dao

def mosaic_phot(mosaic_fits, channel, exptime, dao_dir='/usr/local/phot/', data_dir='', opt_dir='/Volumes/Annie/CRRP/OPTfiles/'):

    curr_dir = os.getcwd()
    os.chdir(data_dir+'DeepMosaic')
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
    os.chdir('../')
    ## need to update .lst file with new ID numbers/positions

def match_optical(target, channel, data_dir='', optical_dir='/Volumes/Annie/CRRP/OpticalCatalogs/', dao_dir='/usr/local/phot/'):

    deep_mosaic_fits = target+'_'+channel+'_deep.fits'

    curr_dir = os.getcwd()
    os.chdir(data_dir+'DeepMosaic')
    ids, catalog_x, catalog_y, catalog_ra, catalog_dec = optical.read_optical_fnl(optical_dir, target)

    als_file = re.sub('.fits', '_dn.als', deep_mosaic_fits)
    dtype1 = np.dtype([('x', float), ('y', float)])
    data = np.loadtxt(als_file, dtype=dtype1, skiprows=3, usecols=(1,2))
    xmin = np.min(data['x'])
    xmax = np.max(data['x'])
    ymin = np.min(data['y'])
    ymax = np.max(data['y'])

    print "Calculating optical boundaries..."
    print xmin, xmax, ymin, ymax
    ra1, ra2, dec1, dec2 = coordinates.find_coord_window_mosaic(deep_mosaic_fits, xmin, xmax, ymin, ymax)
    print ra1, ra2, dec1, dec2
    min_x, min_y = coordinates.radec2catalogpix(ra1, dec1, catalog_x, catalog_y, catalog_ra, catalog_dec)
    max_x, max_y = coordinates.radec2catalogpix(ra2, dec2, catalog_x, catalog_y, catalog_ra, catalog_dec)
#		c1, c2, c3, c4 = coordinates.radec2pix(target, x1, x2, y1, y2, xcat, ycat, ra, dec)
    print "Xmin, Xmax, Ymin, Ymax for optical catalog:"
    print min_x, max_x, min_y, max_y

    xmin = [min_x]
    xmax = [max_x]
    ymin = [min_y]
    ymax = [max_y]
    f = ['Deep']

# Save boundary window for each field into a text file (e.g. I1-catalog-cuts.txt)
    data_save = np.array(zip(f, xmin, xmax, ymin, ymax), dtype=[('c1', 'S8'),
        ('c2', float), ('c3', float), ('c4', float), ('c5', float)])
    np.savetxt(channel+'-deep-cuts.txt', data_save, comments='', fmt='%s %0.3f %0.3f %0.3f %0.3f')

    limits = str(min_x)+','+str(max_x)+','+str(min_y)+','+str(max_y)
    image_list = ['optical:'+target+'-I.mag', als_file]
    mch_file = 'op-'+channel+'.mch'
    daomatch.daomatch(image_list, mch_file, dao_dir=dao_dir, xy_limits=limits)
    daomaster.daomaster(mch_file, dao_dir=dao_dir, frame_num='1,0.5,1')
    os.chdir(curr_dir)

def check_match(target, channel, optical_dir='/Volumes/Annie/CRRP/OpticalCatalogs/', data_dir=''):

    fig = mp.figure(figsize=(8,8))
    ax1 = fig.add_subplot(111)

    # read optical catalog and add to plots
    ids, xcat, ycat, ra, dec = optical.read_optical_fnl(optical_dir, target)
    ax1.plot(xcat, ycat, '.', color='0.25', markersize=0.75)

    # read boundaries of IRAC data
    dtype1 = np.dtype([('xmin', float), ('xmax', float), ('ymin', float), ('ymax', float)])
    cuts = np.loadtxt(data_dir+'DeepMosaic/'+channel+'-deep-cuts.txt', dtype=dtype1, usecols=(1,2,3,4))

    ax1.plot([cuts['xmin'], cuts['xmax']], [cuts['ymin'], cuts['ymin']],
        '-', color='r', linewidth=2)
    ax1.plot([cuts['xmin'], cuts['xmax']], [cuts['ymax'], cuts['ymax']],
        '-', color='r', linewidth=2)
    ax1.plot([cuts['xmin'], cuts['xmin']], [cuts['ymin'], cuts['ymax']],
        '-', color='r', linewidth=2)
    ax1.plot([cuts['xmax'], cuts['xmax']], [cuts['ymin'], cuts['ymax']],
        '-', color='r', linewidth=2)
    ax1.set_xlabel('X (pixels)')
    ax1.set_ylabel('Y (pixels)')


    # Add transformed catalogs
    data = read_dao.read_alf(data_dir+'DeepMosaic/'+target+'_'+channel+'_deep_dn.als')
    x = data['x']
    y = data['y']

    files, x_off, y_off, transform, dof = read_dao.read_mch(data_dir+'DeepMosaic/op-'+channel+'.mch')

    x_new = float(x_off[1])+float(transform[1][0])*x+float(transform[1][1])*y
    y_new = float(y_off[1])+float(transform[1][2])*x+float(transform[1][3])*y

    ax1.plot(x_new, y_new, '.', markersize=1.8, color='r')
    mp.show()
