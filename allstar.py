#!/usr/bin/env python


import re
import os
import numpy as np
import shutil
import pexpect
from astropy.io import fits
import sys


def allstar_init(dao_dir, target, fitsfile):

	temp=re.sub(".fits","", fitsfile)

## Clean up previous runs

        extensions = ['.als', 's.coo']
        for ext in extensions:
                if (os.path.isfile(temp + ext)):
                        os.remove(temp+ext)
        image=re.sub("data/",target+":", temp)

    #    print "Working on " + image

## Running ALLSTAR
	allstar = pexpect.spawn(dao_dir+'allstar')
	#allstar.logfile = sys.stdout

	allstar.expect("OPT")
	allstar.sendline("")
	allstar.expect("Input image name")
	allstar.sendline(image)
	allstar.expect("File with the PSF")
	allstar.sendline("")
	allstar.expect("Input file")
	allstar.sendline("")
	allstar.expect("File for results")
	allstar.sendline("")
	allstar.expect("Name for subtracted image")
	allstar.sendline("")
	allstar.expect("stars")
	allstar.expect("Good bye")
	allstar.close()

def allstar_mosaic(dao_dir, target, fitsfile):

	temp=re.sub(".fits","", fitsfile)

## Clean up previous runs

	extensions = ['.als', 's.coo']
	for ext in extensions:
		if (os.path.isfile(temp + ext)):
			os.remove(temp+ext)
	image=re.sub("mosaics/",target+"m:", temp)

	#print "Working on " + image
## Running ALLSTAR
	allstar = pexpect.spawn(dao_folder+'allstar')
	#allstar.logfile = sys.stdout

	allstar.expect("OPT")
	allstar.sendline("fi=4.0")
	allstar.expect("OPT")
	allstar.sendline("")
	allstar.expect("Input image name")
	allstar.sendline(image)
	allstar.expect("File with the PSF")
	allstar.sendline("")
	allstar.expect("Input file")
	allstar.sendline("")
	allstar.expect("File for results")
	allstar.sendline("")
	allstar.expect("Name for subtracted image")
	allstar.sendline("")
	allstar.expect("stars")
	allstar.expect("Good bye")
	allstar.close()
