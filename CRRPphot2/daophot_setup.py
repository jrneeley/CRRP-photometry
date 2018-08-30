#!/usr/bin/env python

import re
import shutil
from astropy.io import fits
import config


def spitzer_flux2dn(image, newname="", exptime=None, fluxconv=None, pixratio=1):

	# need to use exptime, not frametime
	if exptime == 2: exptime = 1.2
	if exptime == 6: exptime = 6
	if exptime == 12: exptime = 12
	if exptime == 30: exptime = 23.6
	if exptime == 100: exptime = 100

	# change default fluxconv if not in native pixel ratio
	if pixratio != 1:
		fluxconv *= pixratio

	if (newname == ""):
		newname = re.sub(".fits", "_dn.fits", image)
	shutil.copy(image, newname)
	hdulist = fits.open(newname, mode='update')
	prihdr = hdulist[0].header
	scidata = hdulist[0].data
	if exptime == None : exptime = prihdr['exptime']
	if fluxconv == None : fluxconv = prihdr['fluxconv']
	scidata *= exptime/fluxconv


def get_irac_opt_files(filters, exptime, warm=1, mosaic=1):

	opt_dir = config.opt_dir+'IRAC/'
	if warm == 1:
		opt_dir2 = opt_dir+'warm/'
	if warm == 0:
		opt_dir2 = opt_dir+'cryo/'

	shutil.copy(opt_dir+'allframe.opt', 'allframe.opt')
	shutil.copy(opt_dir+'dummy-daophot.opt', 'daophot.opt')
	if mosaic == 1:
		shutil.copy(opt_dir+'photo-mosaic.opt', 'photo.opt')
		shutil.copy(opt_dir+'allstar-mosaic.opt', 'allstar.opt')
		for filt in filters:
			opt_file = opt_dir2+filt+'-'+str(exptime)+'s-mosaic.opt'
			shutil.copy(opt_file, filt+'-daophot.opt')
	if mosaic == 0:
		shutil.copy(opt_dir+'photo.opt', 'photo.opt')
		shutil.copy(opt_dir+'allstar.opt', 'allstar.opt')
		for filt in filters:
			opt_file = opt_dir2+filt+'-'+str(exptime)+'s.opt'
			shutil.copy(opt_file, filt+'-daophot.opt')


def combine_mch_simple(mch_list, output_file='combine.mch'):

	o = open(output_file, 'w')
	for ii, mch in enumerate(mch_list):

		f = open(mch, 'r')
		lines = f.readlines()
		if ii == 0:
			o.write(lines[0])
		for jj in range(1, len(lines)):
			o.write(lines[jj])

		f.close()
	o.close()


def copy_mch_artificial(als_list, mch_file, output_file='copy.mch'):

	o = open(output_file, 'w')

	f = open(mch_file, 'r')
	mch_lines = f.readlines()

	o.write(mch_lines[0])

	for ii, als in enumerate(als_list):

		newname = als.ljust(40)
		oldline = mch_lines[1]
		#print newline[1:40]
		#newline[1:40] = newname
		newline = ' \''+newname+oldline[42:]
		o.write(newline)
	f.close()
	o.close()
