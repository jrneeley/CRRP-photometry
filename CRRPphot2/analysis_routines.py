import numpy as np
import matplotlib.pyplot as mp
import os
import optical
import config


def merge_opt_deep_catalogs(target, cluster_coord=None, opt_name='None'):

    if opt_name == 'None': opt_name = target
    optical_dir = config.optical_dir
    #data_dir = config.top_dir+target

    # read in optical and nir catalogs if they exist
    include_opt = 0
    include_nir = 0
    opt_file = optical_dir+opt_name+'.fnl'
    nir_file = optical_dir+opt_name+'ir.fnl'
    if os.path.isfile(opt_file):
        if cluster_coord == None: opt_data = optical.read_fnl(opt_name)
        else:
            opt_data, rad_dist = optical.read_fnl_w_radial_dist(opt_name,
                cluster_coord[0], cluster_coord[1])
        include_opt = 1
#    if os.path.isfile(nir_file):
#        if cluster_coord == None: opt_data = optical.read_fnl(opt_name, nir=1)
#        else:
#            nir_data, rad_dist2 = optical.read_fnl_w_radial_dist(opt_name+'ir',
#                cluster_coord[0], cluster_coord[1], nir=1)
#        include_nir = 1


    dtype1 = np.dtype([('id', int), ('3.6', float), ('3.6err', float)])
    dtype2 = np.dtype([('id', int), ('4.5', float), ('4.5err', float)])
    print 'Reading MIR catalog for '+target+'...'
    irac1_file = 'Deepmosaic/I1-deep.cal'
    irac2_file = 'Deepmosaic/I2-deep.cal'
    if os.path.isfile(irac1_file) == 0:
        print 'I1 is uncalibrated!!!'
        irac1_file = 'DeepMosaic/'+target+'_I1_deep_dn.alf'
    if os.path.isfile(irac2_file) == 0:
        print 'I2 is uncalibrated!!!'
        irac2_file = 'DeepMosaic/'+target+'_I2_deep_dn.alf'

    data3p6 = np.loadtxt(irac1_file, dtype=dtype1, usecols=(0,3,4), skiprows=3)
    data4p5 = np.loadtxt(irac2_file, dtype=dtype2, usecols=(0,3,4), skiprows=3)



    if include_opt == 1:

        mir_data = np.zeros(len(opt_data['id']), dtype=([('3.6', float), ('3.6err', float),
            ('4.5', float), ('4.5err', float), ('n3.6', int), ('n4.5', int)]))
        for ind, star in enumerate(opt_data['id']):
            match3p6 = np.argwhere(data3p6['id'] == star)
            match4p5 = np.argwhere(data4p5['id'] == star)
            if len(match3p6) > 0:
                mir_data['3.6'][ind] = data3p6['3.6'][data3p6['id'] == star]
                mir_data['3.6err'][ind] = data3p6['3.6err'][data3p6['id'] == star]
                mir_data['n3.6'][ind] = 1
            else:
                mir_data['3.6'][ind] = 99.999
                mir_data['3.6err'][ind] = 9.9999
            if len(match4p5) > 0:
                mir_data['4.5'][ind] = data4p5['4.5'][data4p5['id'] == star]
                mir_data['4.5err'][ind] = data4p5['4.5err'][data4p5['id'] == star]
                mir_data['n4.5'][ind] = 1
            else:
                mir_data['4.5'][ind] = 99.999
                mir_data['4.5err'][ind] = 9.9999



        if include_nir == 0:

            print 'Writing catalog with optical and MIR data....'
            # add header
            head = ' 7 FILTERS:                 U             B             V             R             I           [3.6]         [4.5]           n    n    n    n    n    n    n    chi  sharp |---------- variability ----------|--- RA  (2000)  Dec ----'

            dtype_comb = np.dtype([('id', int), ('x', float), ('y', float), ('U', float),
                ('Uer', float), ('B', float), ('Ber', float), ('V', float), ('Ver', float),
                ('R', float), ('Rer', float), ('I', float), ('Ier', float), ('3.6', float),
                ('3.6er', float), ('4.5', float), ('4.5er', float), ('nU', int),
                ('nB', int), ('nV', int), ('nR', int), ('nI', int), ('n3.6', int),
                ('n4.5', int), ('chi', float), ('sharp', float), ('var1', float),
                ('var2', float), ('var3', float), ('var4', float), ('var5', float),
                ('ra_h', int), ('ra_m', int), ('ra_s', float),
                ('dec_d', int), ('dec_m', int), ('dec_s', float), ('rad_dist', float)])

            data_save = np.array(zip(opt_data['id'], opt_data['x'], opt_data['y'],
                opt_data['U'], opt_data['Uer'], opt_data['B'], opt_data['Ber'],
                opt_data['V'], opt_data['Ver'], opt_data['R'], opt_data['Rer'],
                opt_data['I'], opt_data['Ier'], mir_data['3.6'], mir_data['3.6err'],
                mir_data['4.5'], mir_data['4.5err'], opt_data['nU'], opt_data['nB'],
                opt_data['nV'], opt_data['nR'], opt_data['nI'], mir_data['n3.6'],
                mir_data['n4.5'], opt_data['chi'], opt_data['sharp'], opt_data['var1'],
                opt_data['var2'], opt_data['var3'], opt_data['var4'], opt_data['var5'],
                opt_data['ra_h'], opt_data['ra_m'], opt_data['ra_s'], opt_data['dec_d'],
                opt_data['dec_m'], opt_data['dec_s'], rad_dist), dtype=dtype_comb)

            np.savetxt('merged-deep-catalog.txt', data_save,
                fmt='%8i %8.2f %8.2f ' + 7*'%6.3f %6.4f ' + 7*'%4i ' + 5*'%6.3f '+'%6.1f %6.3f %3i %02i %05.2f %+03i %02i %04.1f %0.2f',
                header=head)

##### THIS IS CURRENTLY BROKEN. NEED TO MATCH NIR DATA TO OPTICAL
        if include_nir == 1:

            print 'Writing catalog with optical, NIR, and MIR data....'
            head = '10 FILTERS:                 U             B             V             R             I             J             H             K           [3.6]         [4.5]           n    n    n    n    n    n    n    chi  sharp |---------- variability ----------|--- RA  (2000)  Dec ----'

            dtype_comb = np.dtype([('id', int), ('x', float), ('y', float), ('U', float),
                ('Uer', float), ('B', float), ('Ber', float), ('V', float), ('Ver', float),
                ('R', float), ('Rer', float), ('I', float), ('Ier', float), ('J', float),
                ('Jer', float), ('H', float), ('Her', float), ('K', float), ('Ker', float),
                ('3.6', float), ('3.6er', float), ('4.5', float), ('4.5er', float),
                ('nU', int), ('nB', int), ('nV', int), ('nR', int), ('nI', int),
                ('nJ', int), ('nH', int), ('nK', int), ('n3.6', int), ('n4.5', int),
                ('chi', float), ('sharp', float), ('var1', float),
                ('var2', float), ('var3', float), ('var4', float), ('var5', float),
                ('ra_h', int), ('ra_m', int), ('ra_s', float),
                ('dec_d', int), ('dec_m', int), ('dec_s', float), ('rad_dist', float)])


            data_save = np.array(zip(opt_data['id'], opt_data['x'], opt_data['y'],
                opt_data['U'], opt_data['Uer'], opt_data['B'], opt_data['Ber'],
                opt_data['V'], opt_data['Ver'], opt_data['R'], opt_data['Rer'],
                opt_data['I'], opt_data['Ier'], nir_data['J'], nir_data['Jer'],
                nir_data['H'], nir_data['Her'], nir_data['K'], nir_data['Ker'],
                mir_data['3.6'], mir_data['3.6err'], mir_data['4.5'], mir_data['4.5err'],
                opt_data['nU'], opt_data['nB'], opt_data['nV'], opt_data['nR'], opt_data['nI'],
                nir_data['nJ'], nir_data['nH'], nir_data['nK'], mir_data['n3.6'],
                mir_data['n4.5'], opt_data['chi'], opt_data['sharp'], opt_data['var1'],
                opt_data['var2'], opt_data['var3'], opt_data['var4'], opt_data['var5'],
                opt_data['ra_h'], opt_data['ra_m'], opt_data['ra_s'], opt_data['dec_d'],
                opt_data['dec_m'], opt_data['dec_s'], rad_dist), dtype=dtype_comb)

            np.savetxt('merged-deep-catalog.txt', data_save,
                fmt='%8i %8.2f %8.2f ' + 10*'%6.3f %6.4f ' + 10*'%4i ' + 5*'%6.3f '+'%6.1f %6.3f %3i %02i %05.2f %+03i %02i %04.1f %0.2f',
                header=head)

    if (include_nir == 1) & (include_opt == 0):


        mir_data = np.zeros(len(nir_data['id']), dtype=([('3.6', float), ('3.6err', float),
            ('4.5', float), ('4.5err', float), ('n3.6', int), ('n4.5', int)]))
        for ind, star in enumerate(nir_data['id']):
            match3p6 = np.argwhere(data3p6['id'] == star)
            match4p5 = np.argwhere(data4p5['id'] == star)
            if len(match3p6) > 0:
                mir_data['3.6'][ind] = data3p6['3.6'][data3p6['id'] == star]
                mir_data['3.6err'][ind] = data3p6['3.6err'][data3p6['id'] == star]
                mir_data['n3.6'][ind] = 1
            else:
                mir_data['3.6'][ind] = 99.999
                mir_data['3.6err'][ind] = 9.9999
            if len(match4p5) > 0:
                mir_data['4.5'][ind] = data4p5['4.5'][data4p5['id'] == star]
                mir_data['4.5err'][ind] = data4p5['4.5err'][data4p5['id'] == star]
                mir_data['n4.5'][ind] = 1
            else:
                mir_data['4.5'][ind] = 99.999
                mir_data['4.5err'][ind] = 9.9999

        print 'Writing catalog with NIR and MIR data....'
        # add header
        head = ' 5 FILTERS:                 J             H             K           [3.6]         [4.5]           n    n    n    n    n    chi  sharp |---------- variability ----------|--- RA  (2000)  Dec ----'

        dtype_comb = np.dtype([('id', int), ('x', float), ('y', float), ('J', float),
            ('Jer', float), ('H', float), ('Her', float), ('K', float), ('Ker', float),
            ('3.6', float), ('3.6er', float), ('4.5', float), ('4.5er', float),
            ('nJ', int), ('nH', int), ('nK', int), ('n3.6', int), ('n4.5', int),
            ('chi', float), ('sharp', float), ('var1', float),
            ('var2', float), ('var3', float), ('var4', float), ('var5', float),
            ('ra_h', int), ('ra_m', int), ('ra_s', float),
            ('dec_d', int), ('dec_m', int), ('dec_s', float), ('rad_dist', float)])

        data_save = np.array(zip(nir_data['id'], nir_data['x'], nir_data['y'],
            nir_data['J'], nir_data['Jer'], nir_data['H'], nir_data['Her'],
            nir_data['K'], nir_data['Ker'], mir_data['3.6'], mir_data['3.6err'],
            mir_data['4.5'], mir_data['4.5err'], nir_data['nJ'], nir_data['nH'],
            nir_data['nK'], mir_data['n3.6'], mir_data['n4.5'],
            nir_data['chi'], nir_data['sharp'], nir_data['var1'],
            nir_data['var2'], nir_data['var3'], nir_data['var4'], nir_data['var5'],
            nir_data['ra_h'], nir_data['ra_m'], nir_data['ra_s'], nir_data['dec_d'],
            nir_data['dec_m'], nir_data['dec_s'], rad_dist2), dtype=dtype_comb)

        np.savetxt('merged-deep-catalog.txt', data_save,
            fmt='%8i %8.2f %8.2f ' + 5*'%6.3f %6.4f ' + 5*'%4i ' + 5*'%6.3f '+'%6.1f %6.3f %3i %02i %05.2f %+03i %02i %04.1f %0.2f',
            header=head)




def read_merged_catalog(target):#, center_ra, center_dec):

    data_dir = config.top_dir+target

    # determine number of columns to know format
    f = open(data_dir+'/merged-deep-catalog.txt', 'r')
    head = f.readline()
    head = f.readline()
    head = f.readline()
    line1 = f.readline()
    firstline = line1.split()
    n_cols = len(firstline)
    f.close()

    # 32 columns means 5 filters (NIR + MIR)
    # 38 columns means 7 filters (OPT + MIR)
    # 47 columns means 10 filters (OPT + NIR + MIR)
    if n_cols == 32:

        dtype_comb = np.dtype([('id', int), ('x', float), ('y', float), ('J', float),
            ('Jer', float), ('H', float), ('Her', float), ('K', float), ('Ker', float),
            ('3.6', float), ('3.6er', float), ('4.5', float), ('4.5er', float),
            ('nJ', int), ('nH', int), ('nK', int), ('n3.6', int), ('n4.5', int),
            ('chi', float), ('sharp', float), ('var1', float),
            ('var2', float), ('var3', float), ('var4', float), ('var5', float),
            ('ra_h', int), ('ra_m', int), ('ra_s', float),
            ('dec_d', int), ('dec_m', int), ('dec_s', float), ('rad_dist', float)])

    if n_cols == 38:

        dtype_comb = np.dtype([('id', int), ('x', float), ('y', float), ('U', float),
            ('Uer', float), ('B', float), ('Ber', float), ('V', float), ('Ver', float),
            ('R', float), ('Rer', float), ('I', float), ('Ier', float), ('3.6', float),
            ('3.6er', float), ('4.5', float), ('4.5er', float), ('nU', int),
            ('nB', int), ('nV', int), ('nR', int), ('nI', int), ('n3.6', int),
            ('n4.5', int), ('chi', float), ('sharp', float), ('var1', float),
            ('var2', float), ('var3', float), ('var4', float), ('var5', float),
            ('ra_h', int), ('ra_m', int), ('ra_s', float),
            ('dec_d', int), ('dec_m', int), ('dec_s', float), ('rad_dist', float)])

    if n_cols == 47:

        dtype_comb = np.dtype([('id', int), ('x', float), ('y', float), ('U', float),
            ('Uer', float), ('B', float), ('Ber', float), ('V', float), ('Ver', float),
            ('R', float), ('Rer', float), ('I', float), ('Ier', float), ('J', float),
            ('Jer', float), ('H', float), ('Her', float), ('K', float), ('Ker', float),
            ('3.6', float), ('3.6er', float), ('4.5', float), ('4.5er', float),
            ('nU', int), ('nB', int), ('nV', int), ('nR', int), ('nI', int),
            ('nJ', int), ('nH', int), ('nK', int), ('n3.6', int), ('n4.5', int),
            ('chi', float), ('sharp', float), ('var1', float),
            ('var2', float), ('var3', float), ('var4', float), ('var5', float),
            ('ra_h', int), ('ra_m', int), ('ra_s', float),
            ('dec_d', int), ('dec_m', int), ('dec_s', float), ('rad_dist', float)])

    data = np.loadtxt(data_dir+'/merged-deep-catalog.txt', dtype=dtype_comb)

    return data
