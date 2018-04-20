#!/usr/bin/env python

import re
import shutil
from astropy.io import fits
import config


def spitzer_flux2dn(image, newname="", exptime=None, fluxconv=None):

	if exptime == 30: exptime = 23.6

	if (newname == ""):
		newname = re.sub(".fits", "_dn.fits", image)
	shutil.copy(image, newname)
	hdulist = fits.open(newname, mode='update')
	prihdr = hdulist[0].header
	scidata = hdulist[0].data
	if exptime == None : exptime = prihdr['exptime']
	if fluxconv == None : fluxconv = prihdr['fluxconv']
	scidata *= exptime/fluxconv


def set_opt_files(channel, exptime, warm=1, mosaic=1):

	opt_dir = config.opt_dir
	if warm == 1:
		opt_dir2 = opt_dir+'warm/'
	if warm == 0:
		opt_dir2 = opt_dir+'cryo/'

	if channel == 'I1':
		ch = 'ch1'
	if channel == 'I2':
		ch = 'ch2'
	if mosaic == 1:
		ext = 's-mosaic.opt'
	if mosaic == 0:
		ext = 's.opt'
	optfile = opt_dir2+ch+'-'+str(exptime)+ext

	shutil.copy(optfile, 'daophot.opt')
	shutil.copy(opt_dir+'photo-mosaic.opt', 'photo.opt')
	shutil.copy(opt_dir+'allstar-mosaic.opt', 'allstar.opt')
	shutil.copy(opt_dir+'allframe.opt', 'allframe.opt')
