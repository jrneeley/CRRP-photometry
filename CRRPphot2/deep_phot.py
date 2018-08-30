import numpy as np
import matplotlib.pyplot as mp
import daophot_setup
import dao
import re
import os
import sys
import optical
import coordinates
import config
import shutil
#sys.path.insert(0, '/Users/Jill/python/daophot-tools/')


def mosaic_phot(target, channel, exptime):

    data_dir = config.top_dir+target
    os.chdir(data_dir+'/DeepMosaic')
    print 'Changed directory to {}'.format(data_dir+'/DeepMosaic')

# headers of deep mosaics are wrong, so must input right numbers here.
    # convert deep mosaic to counts
    if channel == 'I1' : fluxconv = 0.1257*4
    if channel == 'I2' : fluxconv = 0.1447*4

    mosaic_fits = target+'_'+channel+'_deep.fits'
    daophot_setup.spitzer_flux2dn(mosaic_fits, exptime=float(exptime), fluxconv=fluxconv)
    print 'Copying OPT files...'
    # copy OPT files to current dir
    daophot_setup.set_opt_files(channel, exptime, warm=1, mosaic=1)
    print 'Doing initial aperture photometry...'
    # Find stars and do initial aperture photometry on deep mosaic
    mosaic_dn = target+'_'+channel+'_deep_dn.fits'
    dao.daophot(mosaic_dn, find=1, phot=1, verbose=0)
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
    dao.find_psf(mosaic_dn)
    repeat = raw_input('Do you want to remove any PSF stars? [y/n]: ')
    if repeat == 'y':
        ready = raw_input('Delete from .lst file and type continue when ready: ')
        dao.find_psf(mosaic_dn)
    if repeat == 'n':
        print 'PSF ok...continuing...\n'
    print 'Running allstar on deep mosaic...'
    dao.allstar(mosaic_dn)
    sub_img = re.sub('dn.fits', 'dns.fits', mosaic_dn)
    print '\nRunning FIND/PHOT on subtracted image...'
    dao.daophot(sub_img, find=1, phot=1)
    print 'Appending new stars to star list...'
    dao.append(mosaic_dn)
    print 'Running PHOT on new star list...'
    dao.daophot(mosaic_dn, find=0, phot=1, coo_file='.srt')
    print 'Running allstar on deep mosaic again...'
    dao.allstar(mosaic_dn)
    print 'Finished with deep mosaic. \n'
    os.chdir(data_dir)
    print 'Changed directory to {}'.format(data_dir)

def match_optical(target, channel, opt_name='None'):

    if opt_name == 'None': opt_name = target
    deep_mosaic_fits = target+'_'+channel+'_deep.fits'

    #data_dir = config.top_dir+target
    #os.chdir(data_dir+'/DeepMosaic')
    os.chdir('DeepMosaic')

    ids, catalog_x, catalog_y, catalog_ra, catalog_dec = optical.read_optical_fnl(opt_name)

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

    print 'Matching optical and MIR catalogs...'
    limits = str(min_x)+','+str(max_x)+','+str(min_y)+','+str(max_y)
    image_list = ['optical:'+opt_name+'-I.mag', als_file]
    mch_file = 'op-'+channel+'.mch'
    dao.daomatch(image_list, mch_file, xy_limits=limits)
    dao.daomaster(mch_file, frame_num='2,0.5,2', verbose=1)
    os.chdir('../')

def check_match(target, channel, opt_name='None', save=1):

    if opt_name == 'None': opt_name = target
    #data_dir = config.top_dir+target

    fig = mp.figure(figsize=(8,8))
    ax1 = fig.add_subplot(111)

    # read optical catalog and add to plots
    ids, xcat, ycat, ra, dec = optical.read_optical_fnl(opt_name)
    ax1.plot(xcat, ycat, '.', color='0.25', markersize=0.75)

    # read boundaries of IRAC data
    dtype1 = np.dtype([('xmin', float), ('xmax', float), ('ymin', float), ('ymax', float)])
    #cuts = np.loadtxt(data_dir+'/DeepMosaic/'+channel+'-deep-cuts.txt', dtype=dtype1, usecols=(1,2,3,4))
    cuts = np.loadtxt('DeepMosaic/'+channel+'-deep-cuts.txt', dtype=dtype1, usecols=(1,2,3,4))
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
    #data = dao.read_alf(data_dir+'/DeepMosaic/'+target+'_'+channel+'_deep_dn.als')
    data = dao.read_alf('DeepMosaic/'+target+'_'+channel+'_deep_dn.als')
    x = data['x']
    y = data['y']

    #files, x_off, y_off, transform, dof = dao.read_mch(data_dir+'/DeepMosaic/op-'+channel+'.mch')
    files, x_off, y_off, transform, dof = dao.read_mch('DeepMosaic/op-'+channel+'.mch')

    x_new = float(x_off[1])+float(transform[1][0])*x+float(transform[1][1])*y
    y_new = float(y_off[1])+float(transform[1][2])*x+float(transform[1][3])*y

    ax1.plot(x_new, y_new, '.', markersize=1.8, color='r')
    if save == 1:
        mp.savefig('{}-{}-match.pdf'.format(target, channel), format='pdf')
    else:
        mp.show()

# testing with known psf
def mosaic_phot2(target, channel, exptime):

    data_dir = config.top_dir+target
    os.chdir(data_dir+'/DeepMosaic')
    print 'Changed directory to {}'.format(data_dir+'/DeepMosaic')

# headers of deep mosaics are wrong, so must input right numbers here.
    # convert deep mosaic to counts
    if channel == 'I1' : fluxconv = 0.1257*4
    if channel == 'I2' : fluxconv = 0.1447*4

    mosaic_fits = target+'_'+channel+'_deep.fits'
    daophot_setup.spitzer_flux2dn(mosaic_fits, exptime=float(exptime), fluxconv=fluxconv)
    print 'Copying OPT files...'
    # copy OPT files to current dir
    daophot_setup.set_opt_files(channel, exptime, warm=1, mosaic=1)
    print 'Doing initial aperture photometry...'
    # Find stars and do initial aperture photometry on deep mosaic
    mosaic_dn = target+'_'+channel+'_deep_dn.fits'
    dao.daophot(mosaic_dn, find=1, phot=1, num_frames='60,1')

    # copy psf to frame
    psf_dir = config.top_dir+'PSF/'
    shutil.copy(psf_dir+channel+'-master.psf', target+'_'+channel+'_deep_dn.psf')

    print 'Running allstar on deep mosaic...'
    dao.allstar(mosaic_dn)
    sub_img = re.sub('dn.fits', 'dns.fits', mosaic_dn)
    print '\nRunning FIND/PHOT on subtracted image...'
    dao.daophot(sub_img, find=1, phot=1, num_frames='60,1')
    print 'Appending new stars to star list...'
    dao.append(mosaic_dn)
    print 'Running PHOT on new star list...'
    dao.daophot(mosaic_dn, find=0, phot=1, coo_file='.srt')
    print 'Running allstar on deep mosaic again...'
    dao.allstar(mosaic_dn)
    print 'Finished with deep mosaic. \n'
    os.chdir(data_dir)
    print 'Changed directory to {}'.format(data_dir)
