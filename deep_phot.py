import numpy as np
import daophot_setup
import daophot
import re
import allstar

def mosaic_phot(dao_dir, opt_dir, target, channel, exptime):

    mosaic_fits = target+'_'+channel+'_deep.fits'
# headers of deep mosaics are wrong, so must input right numbers here.
    # convert deep mosaic to counts
    if channel == 'I1' : fluxconv = 0.1257
    if channel == 'I2' : fluxconv = 0.1447
    daophot_setup.spitzer_flux2dn(mosaic_fits, exptime=float(exptime), fluxconv=fluxconv)
    # copy OPT files to current dir
    daophot_setup.set_opt_files_mosaic(opt_dir, channel, exptime, warm=1)
    # Find stars and do initial aperture photometry on deep mosaic
    mosaic_dn = re.sub('.fits', '_dn.fits', mosaic_fits)
    daophot.mosaic_phot(dao_dir, mosaic_dn)
    # You must choose your PSF stars from the .coo by hand
    ready = raw_input('Pick PSF stars by hand, type continue when ready, or stop to exit: ')
    if ready == 'stop': exit()
    if ready == 'continue':

        # measure the PSF
        daophot.find_psf(dao_dir, mosaic_dn)
        allstar.allstar_deep(dao_dir, mosaic_dn)
        sub_img = re.sub('dn.fits', 'dns.fits', mosaic_dn)
        daophot.mosaic_phot(dao_dir, sub_img)
        daophot.append(dao_dir, mosaic_dn)
        daophot.mosaic_phot2(dao_dir, mosaic_dn)
        allstar.allstar_deep(dao_dir, mosaic_dn)
