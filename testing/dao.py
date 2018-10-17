import re
import os
import numpy as np
import pexpect
import sys
import matplotlib.pyplot as mp
import config

def daophot(fitsfile, threshold=None, find=1, phot=1, coo_file='',
    ap_file='', num_frames='1,1', opt_file='', verbose=0):

    dao_dir = config.dao_dir
    image = re.sub(".fits","", fitsfile)

## Running daophot
    daophot = pexpect.spawn(dao_dir+'daophot')
    if verbose == 1:
        daophot.logfile = sys.stdout
# Load appropriate opt file and edit threshold if necessary
    daophot.expect('Command:')
    daophot.sendline('op')
    daophot.expect('File with parameters')
    daophot.sendline(opt_file)
    if threshold != None:
        #daophot.expect('Command:')
        #daophot.sendline('op')
        #daophot.expect('File with parameters')
        #daophot.sendline('')
        daophot.expect('OPT>')
        daophot.sendline('th='+str(threshold))
    daophot.expect('OPT>')
    daophot.sendline('')


# ATTACH
    daophot.expect("Command:")
    daophot.sendline("at " + image)

# FIND
    if find == 1:
        daophot.expect("Command:")
        daophot.sendline("find")
        daophot.expect("Number of frames averaged, summed:")
        daophot.sendline(num_frames)
        daophot.expect("File for positions")
        daophot.sendline(coo_file)
        check = daophot.expect(['Are you happy with this?', 'OVERWRITE'])
        if check == 1:
            daophot.sendline('')
            daophot.expect("Are you happy with this?")
            daophot.sendline("y")
        if check == 0:
            daophot.sendline("y")

# PHOT
    if phot == 1:
        daophot.expect("Command:")
        daophot.sendline("phot")
        daophot.expect("File with aperture radii")
        daophot.sendline("")
        daophot.expect("PHO>")
        daophot.sendline("")
        check = daophot.expect(['Input position file', 'Profile-fitting photometry'])
        if check == 1:
            daophot.sendline('e')
            daophot.expect('Input position file')
            daophot.sendline(coo_file)
        if check == 0:
            daophot.sendline(coo_file)
        daophot.expect("Output file")
        daophot.sendline(ap_file)

## Exit Daophot
    check2 = daophot.expect(['Command:', 'OVERWRITE'])
    if check2 == 1:
        daophot.sendline('')
        daophot.expect('Command:')
        daophot.sendline('exit')
    if check2 == 0:
        daophot.sendline("exit")
    daophot.close(force=True)

#print "Initial aperture photometry complete."

def find(fitsfile, num_frames='1,1', coo_file='', opt_file='', verbose=0):

    dao_dir = config.dao_dir
    image = re.sub(".fits","", fitsfile)

## Running daophot
    daophot = pexpect.spawn(dao_dir+'daophot')
    if verbose == 1:
        daophot.logfile = sys.stdout

# Load appropriate opt file and edit threshold if necessary
    daophot.expect('Command:')
    daophot.sendline('op')
    daophot.expect('File with parameters')
    daophot.sendline(opt_file)
    daophot.expect('OPT>')
    daophot.sendline('')

# ATTACH
    daophot.expect("Command:")
    daophot.sendline("at " + image)

    daophot.expect("Command:")
    daophot.sendline("find")
    daophot.expect("Number of frames averaged, summed:")
    daophot.sendline(num_frames)
    daophot.expect("File for positions")
    daophot.sendline(coo_file)
    check = daophot.expect(['Are you happy with this?', 'OVERWRITE'])
    if check == 1:
        daophot.sendline('')
        daophot.expect("Are you happy with this?")
        daophot.sendline("y")
    if check == 0:
        daophot.sendline("y")
    daophot.expect('Command:')
    daophot.sendline('exit')
    daophot.close(force=True)

def phot(fitsfile, phot_file='', coo_file='', ap_file='', opt_file='', verbose=1):

    dao_dir = config.dao_dir
    image = re.sub(".fits","", fitsfile)

## Running daophot
    daophot = pexpect.spawn(dao_dir+'daophot')
    if verbose == 1:
        daophot.logfile = sys.stdout

# Load appropriate opt file and edit threshold if necessary
    daophot.expect('Command:')
    daophot.sendline('op')
    daophot.expect('File with parameters')
    daophot.sendline(opt_file)
    daophot.expect('OPT>')
    daophot.sendline('')

# ATTACH
    daophot.expect("Command:")
    daophot.sendline("at " + image)

# PHOT
    daophot.expect("Command:")
    daophot.sendline("phot")
    daophot.expect("File with aperture radii")
    daophot.sendline("")
    daophot.expect("PHO>")
    daophot.sendline("")
    check = daophot.expect(['Input position file', 'Profile-fitting photometry'])
    if check == 1:
        daophot.sendline('e')
        daophot.expect('Input position file')
        daophot.sendline(coo_file)
    if check == 0:
        daophot.sendline(coo_file)
    daophot.expect("Output file")
    daophot.sendline(ap_file)

## Exit Daophot
    check2 = daophot.expect(['Command:', 'OVERWRITE'])
    if check2 == 1:
        daophot.sendline('')
        daophot.expect('Command:')
        daophot.sendline('exit')
    if check2 == 0:
        daophot.sendline("exit")
    daophot.close(force=True)

def find_psf(fitsfile, opt_file=''):

    file_stem = re.sub(".fits","", fitsfile)

## Clean up previous runs

    extensions = ['.psf', '.nei']
    for ext in extensions:
        if (os.path.isfile(file_stem + ext)):
            os.remove(file_stem + ext)

    image = fitsfile
    print "Working on " + image

#Running daophot
    daophot = pexpect.spawn(config.dao_dir+'daophot')
    daophot.logfile = sys.stdout

# Load appropriate opt file and edit threshold if necessary
    daophot.expect('Command:')
    daophot.sendline('op')
    daophot.expect('File with parameters')
    daophot.sendline(opt_file)
    daophot.expect('OPT>')
    daophot.sendline('')

# attach the image
    daophot.expect("Command:")
    daophot.sendline("at " + file_stem)
## PSF
    daophot.expect("Command:")
    daophot.sendline("psf")
    daophot.expect("File with aperture results")
    daophot.sendline("")
    daophot.expect("File with PSF")
    daophot.sendline("")
    daophot.expect("File for the PSF")
    daophot.sendline("")
## Exit Daophot
    daophot.expect("Command:")
    daophot.sendline("exit")
    daophot.close(force=True)
    print "PSF complete"


def substar(fitsfile):

	daophot = pexpect.spawn(config.dao_dir+'daophot')

	daophot.expect('Command:')
	daophot.sendline('at '+fitsfile)
	daophot.expect('Command:')
	daophot.sendline('substar')
	daophot.expect('File with the PSF')
	daophot.sendline('')
	daophot.expect('File with photometry')
	daophot.sendline('.als')
	daophot.expect('stars to leave in?')
	daophot.sendline('y')
	daophot.expect('File with star list')
	daophot.sendline('')
	daophot.expect('Name for subtracted image')
	daophot.sendline('')
	daophot.expect('Command:')
	daophot.sendline('ex')
	daophot.close(force=True)

def offset(filename, id_offset=0, x_offset=0.0, y_offset=0.0, mag_offset=0.0):

    daophot = pexpect.spawn(config.dao_dir+'daophot')

    daophot.expect('Command')
    daophot.sendline('off')
    daophot.expect('Input file name')
    daophot.sendline(filename)
    daophot.expect('Additive offsets')
    off_string = '{} {} {} {}'.format(id_offset, x_offset, y_offset, mag_offset)
    daophot.sendline(off_string)
    daophot.expect('Output file name')
    daophot.sendline('')
    check = daophot.expect(['Command', 'OVERWRITE'])
    if check == 1:
        daophot.sendline('')
        daophot.expect('Command')
        daophot.sendline('ex')
    if check == 0:
        daophot.sendline('ex')
    daophot.close(force=True)

def append_old(fitsfile, file1=None, file2=None, renumber='n', verbose=0):

    orig_stem = re.sub('.fits', '', fitsfile)
    sub_stem = re.sub('dn', 'dns', orig_stem)
    if file1 == None:
        file1 = orig_stem+'.coo'
    if file2 == None:
        file2 = sub_stem+'.coo'

    daophot = pexpect.spawn(config.dao_dir+'daophot')
    if verbose == 1:
        daophot.logfile = sys.stdout
	daophot.expect('Command:')
	daophot.sendline('append')
	daophot.expect('First input file')
	daophot.sendline(file1)
	daophot.expect('Second input file')
	daophot.sendline(file2)
	daophot.expect('Output file')
	daophot.sendline('')
    check = daophot.expect(['Command:', 'OVERWRITE'])
    if check == 1:
        daophot.sendline('')
        daophot.expect('Command:')
        daophot.sendline('sort')
    if check == 0:
        daophot.sendline('sort')

    daophot.expect('Which do you want')
    daophot.sendline('3')
    daophot.expect('Input file name')
    daophot.sendline(orig_stem+'.cmb')
    daophot.expect('Output file name')
    daophot.sendline('')
    check = daophot.expect(['stars renumbered?', 'OVERWRITE'])
    if check == 1:
        daophot.sendline('')
        daophot.expect('stars renumbered?')
        daophot.sendline(renumber)
    if check == 0:
        daophot.sendline(renumber)
    #print 'made it'
    daophot.expect('Command:')
    daophot.sendline('exit')
    daophot.close(force=True)

def append(file1, file2, out_file='', verbose=0):

    daophot = pexpect.spawn(config.dao_dir+'daophot')
    if verbose == 1:
        daophot.logfile = sys.stdout
	daophot.expect('Command:')
	daophot.sendline('append')
	daophot.expect('First input file')
	daophot.sendline(file1)
	daophot.expect('Second input file')
	daophot.sendline(file2)
	daophot.expect('Output file')
	daophot.sendline(out_file)
    check = daophot.expect(['Command:', 'OVERWRITE'])
    if check == 1:
        daophot.sendline('')
        daophot.expect('Command:')
    daophot.sendline('exit')
    daophot.close(force=True)

def sort(in_file, out_file='', sort_option='3', renumber='y', verbose=0):

    # sort_option +- 1 -> increasing/decreasing ID number
    # sort option +- 2 -> increasing/decreasing X
    # sort option +- 3 -> increasing/decreasing Y
    # sort option +- 4 -> increasing/decreasing magnitude

    daophot = pexpect.spawn(config.dao_dir+'daophot')
    if verbose == 1:
        daophot.logfile = sys.stdout
	daophot.expect('Command:')
	daophot.sendline('sort')
    daophot.expect('Which do you want')
    daophot.sendline(sort_option)
    daophot.expect('Input file name')
    daophot.sendline(in_file)
    daophot.expect('Output file name')
    daophot.sendline(out_file)
    check = daophot.expect(['stars renumbered?', 'OVERWRITE'])
    if check == 1:
        daophot.sendline('')
        daophot.expect('stars renumbered?')
        daophot.sendline(renumber)
    if check == 0:
        daophot.sendline(renumber)
    #print 'made it'
    daophot.expect('Command:')
    daophot.sendline('exit')
    daophot.close(force=True)

def addstar(image, file_stem='fake', num_images = 1, seed=5, gain=999, star_list=None,
    min_mag=12, max_mag=18, num_stars=50, opt_file='', verbose=0):

    daophot = pexpect.spawn(config.dao_dir+'daophot')
    if verbose == 1: daophot.logfile = sys.stdout

    # Make sure we are using the appropriate options file
    daophot.expect('Command:')
    daophot.sendline('opt')
    daophot.expect('File with parameters')
    daophot.sendline(opt_file)
    daophot.expect('OPT')
    daophot.sendline('')

    # First attach original image
    daophot.expect('Command:')
    daophot.sendline('at '+image)
    # start addstar
    daophot.expect('Command:')
    daophot.sendline('ad')
    daophot.expect('File with the PSF')
    daophot.sendline('')
    daophot.expect('Seed')
    daophot.sendline(str(seed))
    daophot.expect('Photons per ADU')
    daophot.sendline(str(gain))
    daophot.expect('Input data file')
    if star_list != None:
        daophot.sendline(star_list)
        daophot.expect('Output picture name')
        daophot.sendline('')
        daophot.expect('Command')
        daophot.sendline('ex')
        daophot.close(force=True)
    else:
        daophot.sendline('')
        daophot.expect('Minimum, maximum magnitudes')
        daophot.sendline('{} {}'.format(min_mag, max_mag))
        daophot.expect('Number of stars to add')
        daophot.sendline(str(num_stars))
        daophot.expect('Number of new frames')
        daophot.sendline(str(num_images))
        daophot.expect('File-name stem')
        daophot.sendline(file_stem)
        check = daophot.expect(['Command', 'OVERWRITE'])
        if check == 0:
            daophot.sendline('ex')
            daophot.close(force=True)
        if check == 1:
            daophot.sendline('')
            for ii in range(num_images-1):
                daophot.expect('OVERWRITE')
                daophot.sendline('')
            daophot.expect('Command')
            daophot.sendline('ex')
            daophot.close(force=True)



def allstar(fitsfile, new_options=0, verbose=0):

    file_stem = re.sub(".fits","", fitsfile)

## Running ALLSTAR
    allstar = pexpect.spawn(config.dao_dir+'allstar', timeout=240)

    if verbose == 1:
        allstar.logfile = sys.stdout

    allstar.expect("OPT")
    if new_options == 0:
        allstar.sendline("")
    else:
        n_changes=len(new_options)
        for ii in range(n_changes):
            allstar.sendline(new_options[ii])
            allstar.expect("OPT")
        allstar.sendline("")
    allstar.expect("Input image name")
    allstar.sendline(file_stem)
    allstar.expect("File with the PSF")
    allstar.sendline("")
    allstar.expect("Input file")
    allstar.sendline("")
    allstar.expect("File for results")
    allstar.sendline("")
    check = allstar.expect(["Name for subtracted image", "OVERWRITE"])
    if check == 1:
        allstar.sendline("")
        allstar.expect("Name for subtracted image")
        allstar.sendline("")
        #print 'made it 1'
    if check == 0:
        allstar.sendline("")
    #print 'made it 2'
    #allstar.expect("stars")
    allstar.expect("Good bye")
    allstar.close(force=True)


# run daomatch on a list of images
def daomatch(image_list, output_file, verbose=0,
                    xy_limits=[], force_scale_rot=0, force_scale=0):

## run DAOMATCH on on fields
    daomatch = pexpect.spawn(config.dao_dir+'daomatch')
    if verbose == 1:
        daomatch.logfile = sys.stdout

    daomatch.expect("Master input file")
    first_file = image_list[0]
    if len(xy_limits) == 4:
        daomatch.sendline(first_file+'*')
        daomatch.expect('Ymin, Ymax')
        xylim = '{} {} {} {}'.format(xy_limits[0], xy_limits[1], xy_limits[2], xy_limits[3])
        daomatch.sendline(xylim)
    elif len(xy_limits) == 0:
        daomatch.sendline(first_file)
    else:
        daomatch.close(force=True)
        print 'Must provide 4 limits! (xmin, xmax, ymin, ymax)'
    daomatch.expect("Output file")
    daomatch.sendline(output_file)
    check = daomatch.expect(["Next input file", "OVERWRITE"])
    if check == 1:
        daomatch.sendline("")
        daomatch.expect("Next input file")
    # define appropriate suffix for images
    suffix = ''
    if force_scale_rot == 1:
        suffix = ';'
    if force_scale != 0:
        suffix = '/'+str(force_scale)
    if len(xy_limits) != 0:
        suffix += '!'

    for image in image_list:
        if image == first_file:
            continue
#            if check == 0:
        daomatch.sendline(image+suffix)
        check = daomatch.expect(["Next input file", "Write this transformation"])
        if check == 1:
            daomatch.sendline("y")
            daomatch.expect("Next input file")
        if check == 0:
            continue

    daomatch.sendline("")
    daomatch.expect("Good bye")
    daomatch.close()


def check_daomatch(mch_file):

    img_list, x_offsets, y_offsets, transform, dof = read_dao.read_mch(mch_file)

    n_imgs = len(img_list)
    master_frame = read_dao.read_alf(img_list[0])

    master_order = np.argsort(master_frame['mag'])
    master_brightest = master_order[0:100]

    for ind in range(1,n_imgs):
        fig = mp.figure(figsize=(10,8))
        ax1 = fig.add_subplot(111)
        ax1.plot(master_frame['x'][master_brightest], master_frame['y'][master_brightest], 'b.', alpha=0.5, markersize=15)
        data = read_dao.read_alf(img_list[ind])
        data_order = np.argsort(data['mag'])
        data_brightest = data_order[0:100]
        x_new = float(x_offsets[ind]) + float(transform[ind][0])*data['x'] + float(transform[ind][1])*data['y']
        y_new = float(y_offsets[ind]) + float(transform[ind][2])*data['x'] + float(transform[ind][3])*data['y']
        ax1.plot(x_new[data_brightest], y_new[data_brightest], 'r.', alpha=0.5, markersize=15)
        ax1.set_title(img_list[ind])
        ax1.set_xlabel('x')
        ax1.set_ylabel('y')
        mp.show()

def daomaster(matchfile, frame_num='12, 0.5, 12', sigma='5',
            transformation='20', new_id='n', mag_file='n', corr_file='n',
            raw_file='n', new_trans='y', verbose=0,
            match_radii=[-4, 3, 3, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1]):
# IN PROGRESS - FIX SO ALL OPTIONS WORK


    daomaster = pexpect.spawn(config.dao_dir+'daomaster')
    if verbose == 1:
        daomaster.logfile = sys.stdout

    daomaster.expect("File with list of input files")
    daomaster.sendline(matchfile)
    daomaster.expect("Minimum number")
    daomaster.sendline(frame_num)
    daomaster.expect("Maximum sigma")
    daomaster.sendline(sigma)
    daomaster.expect("Your choice")
    daomaster.sendline(transformation)
    daomaster.expect("Critical match-up radius")
    daomaster.sendline(str(match_radii[0]))
    for radii in match_radii:
        if radii == match_radii[0]: continue
        daomaster.expect("New match-up radius")
        daomaster.sendline(str(radii))
    daomaster.expect("New match-up radius")
    daomaster.sendline("0")
    daomaster.expect("Assign new star IDs")
    daomaster.sendline(new_id)
    daomaster.expect("A file with mean magnitudes")
    daomaster.sendline(mag_file)
    check = daomaster.expect(['Output file name', 'A file with corrected magnitudes'])
    if check == 0:
        daomaster.sendline('')
        check2 = daomaster.expect(['OVERWRITE', 'A file with corrected magnitudes'])
        if check2 == 0:
            daomaster.sendline('')
            daomaster.expect('A file with corrected magnitudes')
    daomaster.sendline(corr_file)
    daomaster.expect("A file with raw magnitudes")
    daomaster.sendline(raw_file)
    check = daomaster.expect(['Output file name', 'A file with the new transformations'])
    if check == 0:
        daomaster.sendline('')
        check2 = daomaster.expect(['OVERWRITE', 'A file with the new transformations'])
        if check2 == 0:
            daomaster.sendline('')
            daomaster.expect('A file with the new transformations')
    daomaster.sendline(new_trans)
    daomaster.expect("Output file name")
    daomaster.sendline("")
    daomaster.expect("New output file name")
    daomaster.sendline("")
    daomaster.expect("A file with the transfer table")
    daomaster.sendline("e")

    daomaster.close(force=True)

def allframe(image_list, star_list, verbose=0):

    # need very long timeout
    allframe = pexpect.spawn(config.dao_dir+'allframe')
    if verbose == 1:
        allframe.logfile = sys.stdout

    allframe.expect('OPT')
    allframe.sendline('')
    allframe.expect('File with list of images')
    allframe.sendline(image_list)
    allframe.expect('File with list of stars')
    allframe.sendline(star_list)
    allframe.expect('Good bye.', timeout=None)
    allframe.close(force=True)

def read_raw(raw_file):

    f = open(raw_file, 'r')
    lines = f.readlines()
    nstars = 0
    mags_all=[['header']]

    for line in lines:
        temp = line.split()
        if len(temp) == 15:
            nstars += 1
            mags_all.append(temp)
        if len(temp) < 15:
            mags_all[nstars]=mags_all[nstars]+temp

    star_ids = [item[0] for item in mags_all]
    star_ids.remove(star_ids[0])
    mags_all.remove(mags_all[0])

    return star_ids, mags_all

def read_ap(ap_file):

    dtype = np.dtype([('id', float), ('x', float), ('y', float), ('mag', float)])
    data = np.loadtxt(ap_file, dtype=dtype, skiprows=3)

    ids = data['id'][0::2].astype(int)
    x = data['x'][0::2]
    y = data['y'][0::2]
    mags = data['mag'][0::2]
    err = data['mag'][1::2]
    return ids, mags, err

def read_mag(mag_file):

    dtype1 = np.dtype([('id', int), ('x', float), ('y', float), ('mag', float),
        ('err', float)])
    data = np.loadtxt(mag_file, dtype=dtype1, usecols=(0,1,2,3,4), skiprows=3)

    return data

def read_mch(mch_file):

    f = open(mch_file, 'r')
    lines = f.readlines()
    n_lines = len(lines)

    dt = np.dtype([('file', 'S30'), ('x_offset', float), ('y_offset', float), ('transform_matrix', float, (4,))])
    simple_transform = np.zeros(n_lines, dtype=dt)

    file_list = np.zeros(n_lines, dtype='S30')
    x_offset = np.zeros(n_lines)
    y_offset = np.zeros(n_lines)
    transform = np.zeros((n_lines, 4))

    for ii, line in enumerate(lines):
        temp = line.split()
        file_name = temp[0].replace('\'','')
        file_name = file_name.split(':')
        if len(file_name) == 1:
            file_list[ii] = file_name[0]
        if len(file_name) != 1:
            file_list[ii] = file_name[1]
        x_offset[ii] = temp[2]
        y_offset[ii] = temp[3]
        transform[ii] = temp[4:8]

    simple_transform['file'] = file_list
    simple_transform['x_offset'] = x_offset
    simple_transform['y_offset'] = y_offset
    simple_transform['transform_matrix'] = transform

    return simple_transform

def read_nmg(nmg_file):

    data = ascii.read(nmg_file, delimiter=' ', data_start=2)

    id_num = np.array(data['col1'])
    x = np.array(data['col2'])
    y = np.array(data['col3'])
    mag = np.array(data['col4'])
    err = np.array(data['col5'])
    chi = np.array(data['col8'])
    sharp = np.array(data['col9'])

    return id_num, x, y, mag, err, chi, sharp

def read_lst(lst_file):

    dtype1 = np.dtype([('id', int), ('x', float), ('y', float)])
    data = np.loadtxt(lst_file, dtype = dtype1, usecols=(0,1,2), skiprows=3)

    ids = data['id']
    x = data['x']
    y = data['y']


    return ids, x, y

def read_alf(alf_file):

    dtype1 = np.dtype([('id', int), ('x', float), ('y', float), ('mag', float), ('err', float)])
    data = np.loadtxt(alf_file, dtype=dtype1, usecols=(0,1,2,3,4), skiprows=3)

    ids = data['id']
    x = data['x']
    y = data['y']
    mag = data['mag']
    err = data['err']


    return data

def read_coo_new(coo_file):

    dtype1 = np.dtype([('id', int), ('x', float), ('y', float), ('mag1', float),
        ('err', float), ('mag2', float)])
    data = np.loadtxt(coo_file, dtype=dtype1, usecols=(0,1,2,3,4,5), skiprows=3)

    return data

def read_add(add_file):

    dtype = np.dtype([('id', int), ('x', float), ('y', float), ('mag', float)])
    data = np.loadtxt(add_file, dtype=dtype, skiprows=3)

    return data

def write_mag(ids, x, y, mags, out_file, errs=None, err1=None, frames=None, chi=None,
    sharp=None, var=None, blunder=None):

    header_line1 = ' NL    NX    NY  LOWBAD HIGHBAD  THRESH     AP1  PH/ADU  RNOISE    FRAD\n'
    header_line2 = '  1  4013  4766  -500.0 32766.5   0.000   0.000   0.000   0.000   0.000\n\n'
    head = header_line1+'\n'+header_line2+'\n'

    f = open(out_file, 'w')
    f.write(header_line1)
    f.write(header_line2)
    f.close()

    f = open(out_file, 'a')

    n_stars = len(ids)

    if errs == None:
        errs = np.repeat(0.0001, n_stars)
    if err1 == None:
        err1 = np.repeat(0.0000, n_stars)
    if frames == None:
        frames = np.repeat(1., n_stars)
    if chi == None:
        chi = np.repeat(0.000, n_stars)
    if sharp == None:
        sharp = np.repeat(0.000, n_stars)
    if var == None:
        var = np.repeat(0.00, n_stars)
    if blunder == None:
        blunder = np.repeat(1.000, n_stars)
    #dtype = np.dtype([('1', int), ('2', float), ('3', float), ('4', float), \
    #    ('5', float), ('6', float), ('7', float), ('8')])
    data_save = np.array(zip(ids, x, y, mags, errs, err1, frames, chi, sharp, \
        var, blunder))

    np.savetxt(f, data_save,
        fmt='%8i %8.3f %8.3f %8.3f %8.4f %8.4f %8.0f %8.3f %8.3f %8.2f %8.3f')
    f.close()
