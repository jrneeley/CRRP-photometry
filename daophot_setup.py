#!/usr/bin/env python

import re
import shutil
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as mp
import sys
from matplotlib.ticker import ScalarFormatter

def folder():
	daophot_folder = raw_input("Enter path to Daophot executables: ")
	optical_folder = raw_input("Enter path to optical catalogs: ")
	return daophot_folder, optical_folder

def spitzer_flux2dn(image, newname=""):

	if (newname == ""):
		newname = re.sub(".fits", "_dn.fits", image)
	shutil.copy(image, newname)
	hdulist = fits.open(newname, mode='update')
	prihdr = hdulist[0].header
	scidata = hdulist[0].data
	exptime = prihdr['exptime']
	fluxconv = prihdr['fluxconv']
	scidata *= exptime/fluxconv


def set_opt_files(channel):
	if (channel == 'I1'):
		optfile = 'ch1.opt'
	if (channel == 'I2'):
		optfile = 'ch2.opt'
	shutil.copy(optfile, 'daophot.opt')


def find_fields(image_list, channel):
	off_list = []
	on_list = []
	center_ra = []
	center_dec = []
	for image in image_list:
		hdulist = fits.open(image)
		prihdr = hdulist[0].header
		fovid = prihdr['fovid']
# If not a map, we can determine fields directly from header
		if fovid != 81:
			if channel == 'I1':
				if fovid == 74:
					off_list.append(image)
					center_ra.append(prihdr['crval1'])
					center_dec.append(prihdr['crval2'])
				if fovid == 67:
					on_list.append(image)
					center_ra.append(prihdr['crval1'])
					center_dec.append(prihdr['crval2'])
			if channel == 'I2':
				if fovid == 67:
					off_list.append(image)
					center_ra.append(prihdr['crval1'])
					center_dec.append(prihdr['crval2'])
				if fovid == 74:
					on_list.append(image)
					center_ra.append(prihdr['crval1'])
					center_dec.append(prihdr['crval2'])
# If it is a map, we need to find the fields from the coordinates
		if fovid == 81:
			center_ra.append(prihdr['crval1'])
			center_dec.append(prihdr['crval2'])
# Plot to count number of fields
	mp.plot(center_ra, center_dec,'ro')
	mp.ylabel('Dec')
	mp.xlabel('RA')
	x_formatter = ScalarFormatter(useOffset=False)
	mp.gca().xaxis.set_major_formatter(x_formatter)
#	mp.savefig('mapping-positions.eps', format='eps')
	mp.show()
# Ask user how many separate fields there are
	num_fields = input('How many fields exist for this data set?: ')
	if fovid != 81:
		if channel == 'I1':
			field_lists = [off_list, on_list]
		if channel == 'I2':
			field_lists = [on_list, off_list]
	if fovid == 81:
		field_lists = []
		list_of_fields = raw_input('Enter file with list of fields: ')
		f = open(list_of_fields, 'r')
		lists = f.readlines()
		if len(lists) != num_fields:
			sys.exit("Wrong number of fields!")
		for lst in lists:
			lst = lst.strip()
			f2 = open(lst, 'r')
			field_images = f2.readlines()
			field_n = []
			for image in field_images:
				image = image.strip()
				field_n.append(image)
			field_lists.append(field_n)
			f2.close()
		f.close()

	return num_fields, field_lists
