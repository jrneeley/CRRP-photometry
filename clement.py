import coordinates
import numpy as np


center_ra, center_dec = coordinates.radec_string2deg('21:29:58.33', '12:10:01.2')

dtype1 = np.dtype([('id', 'S8'), ('ra', 'S11'), ('dec', 'S11')])
data = np.loadtxt('NGC7078-clement-full.txt', dtype=dtype1, usecols=(0,1,2))

for star, ra, dec in zip(data['id'], data['ra'], data['dec']):

    ra_deg, dec_deg = coordinates.radec_string2deg(ra, dec)
    dist = coordinates.radial_dist(ra_deg, dec_deg, center_ra, center_dec)

    print star, ra, dec, dist
