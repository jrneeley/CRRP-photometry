#!/usr/bin/env python

import re
import os
import numpy as np
import pexpect
import shutil
from astropy.io import fits
import sys

def init_phot(dao_dir, target, fitsfile, mosaics=0):

	temp=re.sub(".fits","", fitsfile)

## Clean up previous runs

	extensions = ['.coo', '.lst', '.psf', '.nei', '.ap', '.als', 's.coo', 's.ap', '.srt', '.cmb', 's.fits', '.als']
    	for ext in extensions:
    		if (os.path.isfile(temp + ext)):
        		os.remove(temp+ext)
	if mosaics == 0: image = re.sub('data/',target+':', temp)
	if mosaics == 1: image = re.sub('mosaics/', target+'m:', temp)
#	print "Working on " + image

## Running daophot

	daophot = pexpect.spawn(dao_dir+'daophot')
	#daophot.logfile = sys.stdout

# attach the image
	daophot.expect("Command:")
	daophot.sendline("at " + image)

# find the stars
	daophot.expect("Command:")
	daophot.sendline("find")
	daophot.expect("Number of frames averaged, summed:")
	daophot.sendline("1,1")
	daophot.expect("File for positions")
	daophot.sendline("")
	daophot.expect("Are you happy with this?")
	daophot.sendline("y")

#	print "FIND complete"

## Aperture photometry
	daophot.expect("Command:")
	daophot.sendline("phot")
	daophot.expect("File with aperture radii")
	daophot.sendline("")
	daophot.expect("PHO>")
	daophot.sendline("")
	daophot.expect("Input position file")
	daophot.sendline(image + '.coo')
	daophot.expect("Output file")
	daophot.sendline(image + '.ap')

#	print "PHOT complete"

## Exit Daophot
	daophot.expect("Command:")
	daophot.sendline("exit")
	daophot.close(force=True)

#print "Initial aperture photometry complete."

def find_psf(dao_dir, fitsfile):

	file_stem = re.sub(".fits","", fitsfile)

## Clean up previous runs

	extensions = ['.psf', '.nei']
	for ext in extensions:
		if (os.path.isfile(file_stem + ext)):
        		os.remove(file_stem + ext)

	image = fitsfile
	print "Working on " + image

#Running daophot
	daophot = pexpect.spawn(dao_dir+'daophot')
	daophot.logfile = sys.stdout
# attach the image
	daophot.expect("Command:")
	daophot.sendline("at " + file_stem)
## PSF
	daophot.expect("Command:")
	daophot.sendline("psf")
	daophot.expect("File with aperture results")
	daophot.sendline("")
	daophot.expect("File with PSF")
	daophot.sendline("")
	daophot.expect("File for the PSF")
	daophot.sendline("")
## Exit Daophot
	daophot.expect("Command:")
	daophot.sendline("exit")
	daophot.close(force=True)
	print "PSF complete"

def deep_mosaic_phot(dao_dir, mosaic_fits):

	file_stem=re.sub(".fits","", mosaic_fits)

## Clean up previous runs
	extensions = ['.coo', '.psf', '.nei', '.ap', '.als']
	for ext in extensions:
		if (os.path.isfile(file_stem + ext)): os.remove(file_stem + ext)

	image = file_stem

	print "Working on " + image
## Running daophot
	daophot = pexpect.spawn(dao_dir+'daophot')
	#daophot.logfile = sys.stdout
# attach the image
	daophot.expect("Command:")
	daophot.sendline("at " + image)
# find the stars
	daophot.expect("Command:")
	daophot.sendline("find")
	daophot.expect("Number of frames averaged, summed:")
	daophot.sendline("1,1")
	daophot.expect("File for positions")
	daophot.sendline("")
	daophot.expect("Are you happy with this?")
	daophot.sendline("y")

## Aperture photometry
	daophot.expect("Command:")
	daophot.sendline("phot")
	daophot.expect("File with aperture radii")
	daophot.sendline("")
	daophot.expect("PHO>")
	daophot.sendline("")
	daophot.expect("Input position file")
	daophot.sendline(image + '.coo')
	daophot.expect("Output file")
	daophot.sendline(image + '.ap')

## Exit Daophot
	daophot.expect("Command:")
	daophot.sendline("exit")
	daophot.close(force=True)

def substar(dao_dir, fitsfile):

	daophot = pexpect.spawn(dao_dir+'daophot')

	daophot.expect('Command:')
	daophot.sendline('at '+fitsfile)
	daophot.expect('Command:')
	daophot.sendline('substar')
	daophot.expect('File with the PSF')
	daophot.sendline('')
	daophot.expect('File with photometry')
	daophot.sendline('.als')
	daophot.expect('stars to leave in?')
	daophot.sendline('y')
	daophot.expect('File with star list')
	daophot.sendline('')
	daophot.expect('Name for subtracted image')
	daophot.sendline('')
	daophot.expect('Command:')
	daophot.sendline('ex')
	daophot.close(force=True)

def append(dao_dir, fitsfile):

	orig_img = re.sub('.fits', '', fitsfile)
	sub_img = re.sub('dn', 'dns', orig_img)

	extensions = ['.cmb', '.srt']
	for ext in extensions:
		if (os.path.isfile(orig_img + ext)): os.remove(orig_img+ext)

	daophot = pexpect.spawn(dao_dir+'daophot')
#	daophot.logfile = sys.stdout
	daophot.expect('Command:')
	daophot.sendline('append')
	daophot.expect('First input file')
	daophot.sendline(orig_img+'.coo')
	daophot.expect('Second input file')
	daophot.sendline(sub_img+'.coo')
	daophot.expect('Output file')
	daophot.sendline('')
	daophot.expect('Command:')
	daophot.sendline('sort')
	daophot.expect('Which do you want')
	daophot.sendline('3')
	daophot.expect('Input file name')
	daophot.sendline(orig_img+'.cmb')
	daophot.expect('Output file name')
	daophot.sendline('')
	daophot.expect('stars renumbered?')
	daophot.sendline('y')
	daophot.expect('Command:')

def deep_mosaic_phot2(dao_dir, mosaic_dn):

	file_stem = re.sub('.fits', '', mosaic_dn)
	extensions = ['.ap']
	for ext in extensions:
		if (os.path.isfile(file_stem + ext)): os.rename(file_stem+ext, file_stem + '_old.ap')

	print "Working on " + file_stem

## Running daophot
	daophot = pexpect.spawn(dao_dir+'daophot')
	#daophot.logfile = sys.stdout
# attach the image
	daophot.expect("Command:")
	daophot.sendline("at " + file_stem)

## Aperture photometry
	daophot.expect("Command:")
	daophot.sendline("phot")
	daophot.expect("File with aperture radii")
	daophot.sendline("")
	daophot.expect("PHO>")
	daophot.sendline("")
	daophot.expect('Profile-fitting photometry')
	daophot.sendline('e')
	daophot.expect("Input position file")
	daophot.sendline(file_stem + '.srt')
	daophot.expect("Output file")
	daophot.sendline(file_stem + '.ap')

## Exit Daophot
	daophot.expect("Command:")
	daophot.sendline("exit")
	daophot.close(force=True)
