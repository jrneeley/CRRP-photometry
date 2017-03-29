#!/usr/bin/env python

import calibration
import sys
import optical

target_name = sys.argv[1]
channel = sys.argv[2]

#optical_folder = raw_input("Enter path to optical catalog: ")
optical_folder = '/Users/jrneeley/CRRP/OpticalCatalogs/'
#calibration.find_cal_star_coords(optical_folder, target_name, channel)

#calibration.find_stars_in_cat(optical_folder, target_name, channel)

#calibration.find_stars_in_cat2(optical_folder, target_name, channel)

#psf_stars, num_nei = calibration.find_cal_stars(target_name, channel)

calibration.find_zp2(channel)
