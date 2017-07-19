import numpy as np
import coordinates
import optical
import lightcurves
from astropy.stats import LombScargle
import matplotlib.pyplot as mp

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
        code = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l']
        print 'V'+data['id'][obj]+' - '+str(len(matches))+' matches'
        for ind, star in enumerate(matches):
            if ind == 0:
                variable_star = 'V'+data['id'][obj]
                lightcurves.make_mosaic_lcv(['I1'], [variable_star], [star])
            if ind <= 12 and ind > 0:
                variable_star = 'V'+data['id'][obj]+code[ind-1]
                lightcurves.make_mosaic_lcv(['I1'], [variable_star], [star])
        #    lightcurves.phase_lcv('lcvs/'+variable_star+'.lcv', data['period'][obj], 50000, bin=0)
            lightcurves.plot_raw_mosaic_lcv('mosaic_lcvs/'+variable_star+'.lcv', star)



def candidate_variables(V, name, plot_save=0, error_threshold=0.05, min_period=0.2,
    max_period=1.0, grid_num=1000):

    x = np.array(V[2][V[1] < error_threshold], dtype=float)
    y = np.array(V[0][V[1] < error_threshold], dtype=float)
    er = np.array(V[1][V[1] < error_threshold], dtype=float)

    mp.figure(figsize=(8,5))

    freq_max = 1/(min_period)
    freq_min = 1/(max_period)
    search_window = max_period - min_period

    frequency = np.linspace(freq_min, freq_max, grid_num)

    power = LombScargle(x, y, er).power(frequency)
    order = np.argsort(power)
#    best_periods = 1/frequency[order[-10:]]
    possible_periods = 1/frequency[order[-10:]]
    bins = np.linspace(min_period, max_period, 11)
    bin_index = np.digitize(possible_periods, bins)
    bins_with_peaks = np.unique(bin_index)
    num_candidates = len(bins_with_peaks)
    candidate_periods = np.zeros(num_candidates)
    for ind, i in enumerate(bins_with_peaks):
        if i == len(bins):
            continue
        max_f = 1/bins[i-1]
        min_f = 1/bins[i]
        power_in_bin = power[(frequency <= max_f) & (frequency >= min_f)]
        max_power = np.max(power_in_bin)
        candidate_periods[ind] = 1/frequency[power == max_power]

    print candidate_periods

    mp.plot(1/frequency, power)
    for period in candidate_periods:
        mp.axvline(period, color='r', alpha = 0.5)
    mp.show()

    return candidate_periods


# Calculate the Welch-Stetson variability index
def welch_stetson_index(band1, band1_errs, band2, band2_errs):

    mean_band1 = sum(band1/band1_errs**2)/sum(1/band1_errs**2)
    mean_band2 = sum(band2/band2_errs**2)/sum(1/band2_errs**2)

    delta_band1 = (band1 - mean_band1)/band1_errs
    delta_band2 = (band2 - mean_band2)/band2_errs

    N = len(band1)
    ws_index = np.sqrt(1/(N*(N-1)))*sum(delta_band1 * delta_band2)

    return ws_index
# def calculate_variability_index():
