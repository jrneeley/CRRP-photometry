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

def append(fitsfile, file1=None, file2=None, renumber='n', verbose=0):

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



def allstar(fitsfile, verbose=0):

    file_stem = re.sub(".fits","", fitsfile)

## Running ALLSTAR
    allstar = pexpect.spawn(config.dao_dir+'allstar', timeout=240)


    if verbose == 1:
        allstar.logfile = sys.stdout

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
    allstar.expect("stars")
    allstar.expect("Good bye")
    allstar.close(force=True)


# run daomatch on a list of images
def daomatch(image_list, output_file, verbose=0,
                    xy_limits=None, force_scale_rot=0, force_scale=0):

## run DAOMATCH on on fields
    daomatch = pexpect.spawn(config.dao_dir+'daomatch')
    if verbose == 1:
        daomatch.logfile = sys.stdout

    daomatch.expect("Master input file")
    first_file = image_list[0]
    if xy_limits != None:
        daomatch.sendline(first_file+'*')
        daomatch.expect('Ymin, Ymax')
        daomatch.sendline(xy_limits)
    if xy_limits == None:
        daomatch.sendline(first_file)
    daomatch.expect("Output file")
    daomatch.sendline(output_file)
    check = daomatch.expect(["Next input file", "New output file"])
    if check == 1:
        daomatch.sendline("")

    # define appropriate suffix for images
    suffix = ''
    if force_scale_rot == 1:
        suffix = ';'
    if force_scale != 0:
        suffix = '/'+str(force_scale)
    if xy_limits != None:
        suffix += '!'

    for image in image_list:
        if image == first_file:
            continue
#            if check == 0:
        daomatch.sendline(image+suffix)
        check = daomatch.expect(["Next input file", "Write this transformation"])

        if check == 1:
            daomatch.sendline('y')
        if check == 0:
            continue
#            daomatch.expect('Next input file')
#            daomatch.sendline(image+suffix)

#    daomatch.expect("Next input file")
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

## Clean up previous runs
    magfile=re.sub(".mch", ".mag", matchfile)
    if (os.path.isfile(magfile)):
            os.remove(magfile)
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
    daomaster.expect("A file with the new transformations")
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

    file_list = []
    x_offset = []
    y_offset = []
    transform_matrix = []
    f = open(mch_file, 'r')
    for line in f:
        temp = line.split()
        n=len(temp[0])
        file_name = temp[0].replace('\'','')
        file_name = file_name.split(':')
        if len(file_name) == 1:
            file_list.append(file_name[0])
        if len(file_name) != 1:
            file_list.append(file_name[1])
        x_offset.append(temp[2])
        y_offset.append(temp[3])
        transform_matrix.append(temp[4:-2])
    dof = len(transform_matrix[1])

    return file_list, x_offset, y_offset, transform_matrix, dof

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
