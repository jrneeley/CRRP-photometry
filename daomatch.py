#!/usr/bin/env python

import sys
import pexpect
import re
import glob
import os
import coordinates

def daomatch_init(dao_dir, channel, target, fields, num_fields):

## Clean up previous runs
    extensions = ['*[0-9].mch']
    for ext in extensions:
        if (os.path.isfile(channel+'_field'+ ext)):
            os.remove(channel+'_field'+ ext)

    for ind in range(0,num_fields):
        img_list = list(fields[ind])
        for ii in range(0,len(img_list)):
            img_list[ii] = re.sub("data/", target+":", img_list[ii])
            img_list[ii] = re.sub('.fits', '_dn.als', img_list[ii])
## run DAOMATCH on on fields
        daomatch = pexpect.spawn(dao_dir+'daomatch')
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

def daomatch_mosaic(dao_dir, channel, target, image_list):

## Clean up previous runs
    mch_file = channel+'_mosaic.mch'
    if (os.path.isfile(mch_file)):
        os.remove(mch_file)


    img_list = list(image_list)
    for ind, img in enumerate(image_list):
        img_list[ind] = re.sub("mosaics/", target+"m:", img)
        img_list[ind] = re.sub('.fits', '.als', img_list[ind])
## run DAOMATCH on on fields
    daomatch = pexpect.spawn(dao_dir+'daomatch')
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

def deep_mosaic(dao_dir, channel, target, bcd_list):

## Clean up previous runs

    mosaic = target+'_'+channel+'_deep_dn.fits'
    first_file = re.sub('.fits', '.als', mosaic)
    for bcd in bcd_list:
        extensions = ['.mch']
        for ext in extensions:
            if (os.path.isfile('temp'+ext)):
                os.remove('temp'+ext)
        next_bcd = re.sub("data/", target+":", bcd)
        next_bcd = re.sub('.fits', '.als', next_bcd)
        xmin, xmax, ymin, ymax = coordinates.find_deep_mos_coords(mosaic, bcd)
        limits = '{:.2f}, {:.2f}, {:.2f}, {:.2f}'.format(xmin, xmax, ymin, ymax)
        ## run DAOMATCH on this file
        daomatch = pexpect.spawn(dao_dir+'daomatch')
    #    daomatch.logfile = sys.stdout
        daomatch.expect("Master input file")
        daomatch.sendline(first_file+'*')
        daomatch.expect('Ymin, Ymax:')
        daomatch.sendline(limits)
        daomatch.expect("Output file")
        daomatch.sendline('temp.mch')
        check = daomatch.expect(["Next input file", "Write this transformation"])
        if check == 0: daomatch.sendline(next_bcd+'!') # / forces scale to be 1
        if check == 1: daomatch.sendline('y')
        daomatch.expect("Next input file")
        daomatch.sendline("")
        daomatch.expect("Good bye")
        daomatch.close()

        mosaic_mch = re.sub('.als', '.mch', first_file)
        if bcd == bcd_list [0]:
            f = open(mosaic_mch, 'w')
            f2 = open('temp.mch', 'r')
            line1 = f2.readline()
            line2 = f2.readline()
            f.write(line1)
            f.write(line2)
            f.close()
            f2.close()
        else:
            f = open(mosaic_mch, 'a')
            f2 = open('temp.mch', 'r')
            line1 = f2.readline()
            line2 = f2.readline()
            f.write(line2)
            f.close()
            f2.close()
