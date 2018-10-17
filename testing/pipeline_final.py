#!/usr/bin/env python

import sys
sys.path.insert(0,'/home/jill/python/CRRP-photometry/')
import glob
import daophot_setup
import shutil
import re
import dao
import numpy as np
import os
import config



targets = ['NGC6121', 'NGC3201', 'NGC5904', 'NGC7078', 'NGC6402']
cluster_ras = ['16:23:35.22', '10:17:36.82', '15:18:33.22', '21:29:58.33', '17:37:36.10']
cluster_decs = ['-26:31:32.7', '-46:24:44.9', '02:04:51.7', '12:10:01.2', '-03:14:45.3']

channels = ['I1', 'I2']

for target in targets:

    exptime = 30
    for channel in channels:
        print 'Working on {} in the {} band.'.format(target, channel)

        data_dir = config.top_dir+target


        # copy appropriate daophot options files to the current directory
        os.chdir(data_dir+'/newmosaics')
        print 'Changed directory to {}'.format(data_dir+'/newmosaics')
        daophot_setup.set_opt_files(channel, exptime, warm=1, mosaic=1)

        # get list of images
        file_list = glob.glob(channel+'*[0-9].fits')

        if channel == 'I1': fluxconv = 0.1257*4
        if channel == 'I2': fluxconv = 0.1447*4

        print 'Converting images to DN....'
        for image in file_list:
            daophot_setup.spitzer_flux2dn(image, "", exptime=exptime, fluxconv=fluxconv)
        print "Done.\n"


        dn_list = glob.glob(channel+'*dn.fits')
        print dn_list

        print "Starting DAOPHOT..."
        for image in dn_list:
            dao.daophot(image, find=1, phot=1, num_frames='5,1')

        print "Initial aperture photometry complete."

        # copy master PSF to every other frame
        master_psf = config.top_dir+'PSF/'+channel+'-master.psf'

        for file in dn_list:
            psfname = re.sub(".fits",".psf", file)
            shutil.copy(master_psf, psfname)

        print "Starting ALLSTAR..."
        for img in dn_list:
            dao.allstar(img)
        print "ALLSTAR complete."

        # copy deep mosaic into current directory
        #os.chdir(data_dir)
        deep = '../DeepMosaic/'+target+'_'+channel+'_deep_dn'
        deep_new = channel+'-deep'
        shutil.copy(deep+'.fits', deep_new+'.fits')
        shutil.copy(deep+'.psf', deep_new+'.psf')
        shutil.copy(deep+'.alf', deep_new+'.alf')
        als_list = glob.glob(channel+'*.als')
        als_list.insert(0, deep_new+'.alf')

        print '\nRunning DAOMATCH/DAOAMASTER between deep mosaic and individual epoch mosaics...'
        mch_file = channel+'-mosaics.mch'
        dao.daomatch(als_list, mch_file, force_scale_rot=1)
        print "DAOMATCH complete."
        dao.daomaster(mch_file, frame_num='12, 0.5, 12')
        print 'DAOMASTER complete.'

        print 'Finished {} channel {}.'.format(target, channel)
