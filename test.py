#!/usr/bin/env python

import calibration
import sys


target_name = sys.argv[1]
channel = sys.argv[2]

#optical_folder = raw_input("Enter path to optical catalog: ")


#calibration.find_cal_star_coords(optical_folder, target_name, channel)

psf_stars = calibration.find_cal_stars(target_name, channel)

calibration.find_zp(psf_stars, channel)
