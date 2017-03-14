#!/usr/bin/env python

import sys
import numpy as np
from astropy.io import ascii
import matplotlib.pyplot as mp

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

def read_mag(mag_file):

    f = open(mag_file, 'r')
    lines = f.readlines()

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
        transform_matrix.append(temp[4:n])
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

    data = ascii.read(lst_file, delimiter=' ', data_start=2)

    ids = np.array(data['col1'])
    x = np.array(data['col2'])
    y = np.array(data['col3'])

    return ids, x, y
