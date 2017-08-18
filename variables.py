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


# Calculate the Welch-Stetson variability index for a single star
def welch_stetson_indices(mags, errs, times, mags2=None, errs2=None, times2=None):


    mean_mag = np.sum(mags/errs**2)/np.sum(1/errs**2)
    N = len(mags)
    res = (mags - mean_mag)/errs
    K_index=0
    #K_index = (1/N*np.sum(np.abs(res)))/np.sqrt(1/N*np.sum(res**2))
    # Find observation pairs
    order = np.argsort(times)
    mags = mags[order]
    errs = errs[order]
    times = times[order]
    delta_t = np.diff(times)
    res1 = np.array([])
    res2 = np.array([])
    test = np.array([])
    for ind, delt in enumerate(delta_t):
        if delt < 0.5:
            res1 = np.append(res1, (mags[ind] - mean_mag)/errs[ind])
            res2 = np.append(res2, (mags[ind+1] - mean_mag)/errs[ind+1])
            test = np.append(test, delt)
    n_pairs = len(res1)
    I_index = np.sqrt(1/(n_pairs*(n_pairs-1)))*np.sum(res1*res2)



    print np.sum(res1*res2)
    print n_pairs

    return I_index, K_index

def robust_weighted_mean(mags, errs):

    errs=errs.astype(float)
    mags=mags.astype(float)
    weights = 1/errs**2
    mean_mag = np.average(mags, weights=weights)
    residuals = (mags - mean_mag)/errs

    for x in range(0,10):
        old_mean = mean_mag
        weights = weights*(1+(np.abs(residuals)/2)**2)**-1
        mean_mag = np.average(mags, weights=weights)

    return mean_mag
# def calculate_variability_index():
