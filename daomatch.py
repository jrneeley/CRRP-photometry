#!/usr/bin/env python

import sys
import pexpect
import re
import glob
import os

def daomatch_init(dao_folder, channel, target, fields, num_fields):

## Clean up previous runs
    extensions = ['*[0-9].mch']
    for ext in extensions:
        if (os.path.isfile(channel+'_field'+ ext)):
            os.remove(channel+'_field'+ ext)

    for ind in range(0,num_fields):
        img_list = list(fields[ind])
        for ii in range(0,len(img_list)):
            img_list[ii] = re.sub("all/", target+":", img_list[ii])
            img_list[ii] = re.sub('.fits', '_dn.als', img_list[ii])
## run DAOMATCH on on fields
        daomatch = pexpect.spawn(dao_folder+'daomatch')
        daomatch.logfile = sys.stdout

        daomatch.expect("Master input file")
        first_file = img_list[0]
        daomatch.sendline(first_file)
        daomatch.expect("Output file")
        daomatch.sendline(channel+"_field"+str(ind+1)+".mch")
        for image in img_list:
            if image == first_file:
                continue
            daomatch.expect("Next input file")
            daomatch.sendline(image+";") # / forces scale to be 1
        daomatch.expect("Next input file")
        daomatch.sendline("")
        daomatch.expect("Good bye")
        daomatch.close()

def daomatch_mosaic(dao_folder, channel, target, image_list):

## Clean up previous runs
    mch_file = channel+'_mosaic.mch'
    if (os.path.isfile(mch_file)):
        os.remove(mch_file)


    img_list = list(image_list)
    for ind, img in enumerate(image_list):
        img_list[ind] = re.sub("mosaics/", target+"m:", img)
        img_list[ind] = re.sub('.fits', '.als', img_list[ind])
## run DAOMATCH on on fields
    daomatch = pexpect.spawn(dao_folder+'daomatch')
    daomatch.logfile = sys.stdout

    daomatch.expect("Master input file")
    first_file = img_list[0]
    daomatch.sendline(first_file)
    daomatch.expect("Output file")
    daomatch.sendline(channel+"_mosaic.mch")
    for image in img_list:
        if image == first_file:
            continue
        daomatch.expect("Next input file")
        daomatch.sendline(image+";") # / forces scale to be 1
    daomatch.expect("Next input file")
    daomatch.sendline("")
    daomatch.expect("Good bye")
    daomatch.close()
