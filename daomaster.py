#!/usr/bin/env python

import sys
import pexpect
import re
import glob
import os

def daomaster_init(matchfile):

## Clean up previous runs
    magfile=re.sub(".mch", ".mag", matchfile)
    if (os.path.isfile(magfile)):
            os.remove(magfile)
    daomaster = pexpect.spawn('/Users/jrneeley/Daophot/daomaster')
    daomaster.logfile = sys.stdout

    daomaster.expect("File with list of input files")
    daomaster.sendline(matchfile)
    daomaster.expect("Minimum number")
    daomaster.sendline("2,0.5,5")
    daomaster.expect("Maximum sigma")
    daomaster.sendline("10")
    daomaster.expect("Your choice")
    daomaster.sendline("20")
    daomaster.expect("Critical match-up radius")
    daomaster.sendline("-10")
    daomaster.expect("New match-up radius")
    daomaster.sendline("8")
    daomaster.expect("New match-up radius")
    daomaster.sendline("6")
    daomaster.expect("New match-up radius")
    daomaster.sendline("4")
    daomaster.expect("New match-up radius")
    daomaster.sendline("3")
    daomaster.expect("New match-up radius")
    daomaster.sendline("3")
    daomaster.expect("New match-up radius")
    daomaster.sendline("2")
    daomaster.expect("New match-up radius")
    daomaster.sendline("2")
    daomaster.expect("New match-up radius")
    daomaster.sendline("2")
    daomaster.expect("New match-up radius")
    daomaster.sendline("1")
    daomaster.expect("New match-up radius")
    daomaster.sendline("1")
    daomaster.expect("New match-up radius")
    daomaster.sendline("1")
    daomaster.expect("New match-up radius")
    daomaster.sendline("0")
    daomaster.expect("Assign new star IDs")
    daomaster.sendline("n")
    daomaster.expect("A file with mean magnitudes")
    daomaster.sendline("y")
    daomaster.expect("Output file name")
    daomaster.sendline("")
    daomaster.expect("A file with corrected magnitudes")
    daomaster.sendline("n")
    daomaster.expect("A file with raw magnitudes")
    daomaster.sendline("n")
    daomaster.expect("A file with the new transformations")
    daomaster.sendline("y")
    daomaster.expect("Output file name")
    daomaster.sendline("")
    daomaster.expect("New output file name")
    daomaster.sendline("")
    daomaster.expect("A file with the transfer table")
    daomaster.sendline("n")
    daomaster.expect("Individual .COO files")
    daomaster.sendline("n")
    daomaster.expect("Simply transfer star IDs")
    daomaster.sendline("n")

    daomaster.close(force=True)
