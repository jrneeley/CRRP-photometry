


def compare_lcv(lcv_file):

    star = lcv_file.split('/')[1]
    dtype1 = np.dtype([('mjd', float), ('mag', float), ('err', float)])
    data = np.loadtxt(lcv_file, dtype=dtype1, usecols=(3,6,7))

    file2 = 'Peter/NGC6121'+star
    dtype2 = np.dtype([('mag', float), ('err', float), ('band', int), ('jd', float)])
    data2 = np.loadtxt(file2, dtype=dtype2, usecols=(0,1,2,5))
    p_mag = data2['mag'][data2['band'] == 1]
    p_mjd = data2['jd'][data2['band'] == 1]+55999.5

    if len(p_mag):

        avg_me = np.nanmean(data['mag'])
        avg_p = np.nanmean(p_mag)
        mp.plot(data['mjd'], data['mag']-avg_me, 'bo')
        mp.plot(p_mjd, p_mag-avg_p, 'ro')
        mp.xlabel('MJD')
        mp.ylabel('mag - avg_mag')

        plt_file = re.sub('.lcv','.pdf', lcv_file)
        mp.savefig(plt_file)
        mp.gcf().clear()
def compare_phased_lcv(lcv_file):

#    star = lcv_file.split('/')[1]
    dtype1 = np.dtype([('phase', float), ('mag', float), ('err', float)])
    data = np.loadtxt(lcv_file, dtype=dtype1, usecols=(1,2,3))

#    file2 = 'Peter/'+star+'_corr.phased'
    file2 = re.sub('lcvs/', 'Peter/', lcv_file)
    file2 = re.sub('.phased', '_corr.phased', file2)
    dtype2 = np.dtype([('band', int), ('phase', float), ('mag', float), ('err', float)])
    data2 = np.loadtxt(file2, dtype=dtype2, usecols=(0,1,2,3))
    p_mag = data2['mag'][data2['band'] == 1]
    p_ph = data2['phase'][data2['band'] == 1]
    p_er = data2['err'][data2['band'] ==1]
    if len(p_mag):

        avg_me = np.nanmean(data['mag'])
        avg_p = np.nanmean(p_mag)
#        mp.plot(data['phase'], data['mag']-avg_me, 'bo')
        mp.errorbar(p_ph, p_mag, yerr=p_er, fmt='o')
        mp.ylim((np.max(p_mag)+0.1, np.max(p_mag)-0.5))
        mp.xlabel('phase')
        mp.ylabel('mag')

        plt_file = re.sub('.phased','.pdf', file2)
        mp.savefig(plt_file)
        mp.gcf().clear()
