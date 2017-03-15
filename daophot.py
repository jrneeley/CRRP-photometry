#!/usr/bin/env python

import re
import os
import numpy as np
import pexpect
import shutil
from astropy.io import fits
import sys

def init_phot(dao_folder, target, fitsfile):

	temp=re.sub(".fits","", fitsfile)

## Clean up previous runs

	extensions = ['.coo', '.lst', '.psf', '.nei', '.ap', '.als', 's.coo', 's.ap', '.srt', '.cmb', 's.fits', '.als']
    	for ext in extensions:
    		if (os.path.isfile(temp + ext)):
        		os.remove(temp+ext)
	image = re.sub("all/",target+":", temp)

	print "Working on " + image

## Running daophot

	daophot = pexpect.spawn(dao_folder+'daophot')
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

	print "FIND complete"

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

	print "PHOT complete"

## Exit Daophot
	daophot.expect("Command:")
	daophot.sendline("exit")
	daophot.close(force=True)

#print "Initial aperture photometry complete."

def find_psf(target, fitsfile):

	temp=re.sub(".fits","", fitsfile)

## Clean up previous runs

	extensions = ['.lst', '.psf', '.nei', '.als', 's.coo', 's.ap', '.srt', '.cmb', 's.fits', '.als']
	for ext in extensions:
		if (os.path.isfile(temp + ext)):
        		os.remove(temp+ext)
	image = re.sub("all/",target+":", temp)

	print "Working on " + image

#Running daophot
	daophot = pexpect.spawn(dao_folder+'daophot')
	daophot.logfile = sys.stdout
# attach the image
	daophot.expect("Command:")
	daophot.sendline("at " + image)
## PSF
	daophot.expect("Command:")
	daophot.sendline("pick")
	daophot.expect("Input file name")
	daophot.sendline("")
	daophot.expect("Desired number")
	daophot.sendline("30,18")
	daophot.expect("Output file name")
	daophot.sendline("")
	daophot.expect("Command:")
	daophot.sendline("psf")
	daophot.expect("File with aperture results")
	daophot.sendline("")
	daophot.expect("File with PSF")
	daophot.sendline("")
	daophot.expect("File for the PSF")
	daophot.sendline("")

	print "PSF complete"

## Exit Daophot
	daophot.expect("Command:")
	daophot.sendline("exit")
	daophot.close(force=True)
