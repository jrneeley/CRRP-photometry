import numpy as np
import coordinates
import optical
import lightcurves

# Find variables listed in the Clement Catalog
def find_variables_by_coord(optical_folder, target):

    # read in list of Clement variables
    dtype1 = np.dtype([('id', 'S8'), ('ra', 'S10'), ('dec', 'S11'), ('period', float),
    ('var_type', 'S8')])
    data = np.loadtxt(target+'-clement.txt', dtype=dtype1)
    RRc = data[:][data['var_type'] == 'RR1']
    RRab = data[:][data['var_type'] == 'RR0']

    ra_variables, dec_variables = coordinates.radec_string2deg(data['ra'], data['dec'])

    ids, xcat, ycat, ra, dec = optical.read_optical_catalog(optical_folder, target)
    # Limit search to horizontal branch
#    ids, xcat, ycat, ra, dec = optical.read_optical_catalog(optical_folder, target+'-HB')

    ra1 = np.radians(ra)
    dec1 = np.radians(dec)
    id_match = []
    ra_match = []
    dec_match = []
    num_neighbors = []

    for obj in range(0,len(ra_variables)):


        ra2 = np.radians(ra_variables[obj])
        dec2 = np.radians(dec_variables[obj])
        x1 = np.cos(dec1)*np.cos(ra1)
        y1 = np.cos(dec1)*np.sin(ra1)
        z1 = np.sin(dec1)
        x2 = np.cos(dec2)*np.cos(ra2)
        y2 = np.cos(dec2)*np.sin(ra2)
        z2 = np.sin(dec2)

        dist = np.sqrt( (x2 - x1)**2 + (y2 - y1)**2 + (z2 - z1)**2)
        dist = np.degrees(dist)*3600.

        matches = ids[dist < 3]
        code = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h']
        print 'V'+data['id'][obj]+' - '+str(len(matches))+' matches'
        for ind, star in enumerate(matches):
            if ind == 0:
                variable_star = 'V'+data['id'][obj]
                lightcurves.make_lcv(['I1'], [variable_star], [star])
            if ind <= 8 and ind > 0:
                variable_star = 'V'+data['id'][obj]+code[ind-1]
                lightcurves.make_lcv(['I1'], [variable_star], [star])
        #    lightcurves.phase_lcv('lcvs/'+variable_star+'.lcv', data['period'][obj], 50000, bin=0)
            lightcurves.plot_raw_lcv('lcvs/'+variable_star+'.lcv', star)

def find_variables_by_coord_mosaic(optical_folder, target):

    # read in list of Clement variables
    dtype1 = np.dtype([('id', 'S8'), ('ra', 'S10'), ('dec', 'S11'), ('period', float),
    ('var_type', 'S8')])
    data = np.loadtxt(target+'-clement.txt', dtype=dtype1)

    RRc = data[:][data['var_type'] == 'RR1']
    RRab = data[:][data['var_type'] == 'RR0']

    ra_variables, dec_variables = coordinates.radec_string2deg(data['ra'], data['dec'])

    ids, xcat, ycat, ra, dec = optical.read_optical_fnl(optical_folder, target)
    # Limit search to horizontal branch
#    ids, xcat, ycat, ra, dec = optical.read_optical_catalog(optical_folder, target+'-HB')

    ra1 = np.radians(ra)
    dec1 = np.radians(dec)
    id_match = []
    ra_match = []
    dec_match = []
    num_neighbors = []

    for obj in range(0,len(ra_variables)):


        ra2 = np.radians(ra_variables[obj])
        dec2 = np.radians(dec_variables[obj])
        x1 = np.cos(dec1)*np.cos(ra1)
        y1 = np.cos(dec1)*np.sin(ra1)
        z1 = np.sin(dec1)
        x2 = np.cos(dec2)*np.cos(ra2)
        y2 = np.cos(dec2)*np.sin(ra2)
        z2 = np.sin(dec2)

        dist = np.sqrt( (x2 - x1)**2 + (y2 - y1)**2 + (z2 - z1)**2)
        dist = np.degrees(dist)*3600.

        matches = ids[dist < 3]
        code = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h']
        print 'V'+data['id'][obj]+' - '+str(len(matches))+' matches'
        for ind, star in enumerate(matches):
            if ind == 0:
                variable_star = 'V'+data['id'][obj]
                lightcurves.make_mosaic_lcv(['I1'], [variable_star], [star])
            if ind <= 8 and ind > 0:
                variable_star = 'V'+data['id'][obj]+code[ind-1]
                lightcurves.make_mosaic_lcv(['I1'], [variable_star], [star])
        #    lightcurves.phase_lcv('lcvs/'+variable_star+'.lcv', data['period'][obj], 50000, bin=0)
            lightcurves.plot_raw_mosaic_lcv('mosaic_lcvs/'+variable_star+'.lcv', star)



# Calculate the Welch-Stetson variability index
# def calculate_variability_index():
