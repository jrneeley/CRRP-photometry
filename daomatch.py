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
## define on and off fields. *THIS ASSUMES 9 DITHERS
#    if (channel == 'I1'):
#    	on_list = glob.glob('all/'+channel+'*1[0-9]_dn.als')
#    	off_list = glob.glob('all/'+channel+'*0[0-9]_dn.als')
#    if (channel == 'I2'):
#    	on_list = glob.glob('all/'+channel+'*0[0-9]_dn.als')
#    	off_list = glob.glob('all/'+channel+'*1[0-9]_dn.als')

#    on_list = re.sub("all/", "img:", on_list)
#    off_list = re.sub("all/", "img:", off_list)

    for ind in range(0,num_fields):
        img_list = list(fields[ind])
        for ii in range(0,len(img_list)):
            img_list[ii] = re.sub("all/", target+":", img_list[ii])
            img_list[ii] = re.sub('.fits', '.als', img_list[ii])
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
## repeat for off fields
#    daomatch = pexpect.spawn('/Users/jrneeley/Daophot/daomatch')
#    daomatch.logfile = sys.stdout

#    daomatch.expect("Master input file")
#    first_file = re.sub("all/","img:",off_list[0])
#    daomatch.sendline(first_file)
#    daomatch.expect("Output file")
#    daomatch.sendline(channel+"_off.mch")
#    for line in off_list:
#        image = re.sub("all/","img:",line)
#        if image == first_file:
#            continue
#        daomatch.expect("Next input file")
#        daomatch.sendline(image+";")
#    daomatch.expect("Next input file")
#    daomatch.sendline("")
#    daomatch.expect("Good bye")
#    daomatch.close()
