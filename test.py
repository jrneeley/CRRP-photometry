#!/usr/bin/env python

import calibration
import sys


target_name = sys.argv[1]
channel = sys.argv[2]


psf_stars = calibration.find_cal_stars(target_name, channel)

calibration.find_zp(psf_stars, channel)
