# This script is intended to make a list of objects and ra/dec coordinates from the pre-sorted TOI lists that can be
# downloaded from ExoFOP.

import math
import numpy as np

catalog_path = '../data/exofop_tess_tois_20210112.csv'
output_path = '../data/TOIS_staralt_20210121.cat'
SG2_table_path = '../data/SG2_parameter_table_INOFFICIAL.csv'
RMexpect_path = '../data/TOIS_RMexpect_20210121.dat'

# read all fields of the ExoFOP export file
TIC_ID, TOI_ID, CTOI, prio_master, SG1A, SG1B, SG2, SG3, SG4, SG5, ACWG, TESS_Dispo, TFOP_Dispo, TESSmag, TESSmag_err, \
pipe_ID, source, ra_deg, dec_deg, PM_ra, PM_ra_err, PM_dec, PM_dec_err, T0_BJD, T0_err, P_d, P_err, dur_h, \
dur_err, depth_mmag, depth_mmag_err, depth_ppm, depth_ppm_err, Prad_rearth, Prad_err, Pinso_fluxearth, Peqtemp_K, \
Psnr, Sdist_pc, Sdist_err, Steff_K, Steff_err, Slogg_cms2, Slogg_err, Srad_rsun, Srad_err, Smetal, Smetal_err, \
sectors, comments, date_TOIalert, date_TOIupdate, date_modified = np.loadtxt(catalog_path, unpack=True,
                                                                             delimiter='","', skiprows=9, dtype=str)

# this list contains TOIs that are already sorted out
TOIs_already_sorted_out = [
    '635',  # elevation too small in Jan 2021
    '886',  # elevation too small in Jan 2021
    '522',  # elevation too small in Jan 2021
    '664',  # elevation too small in Jan 2021
    '1012',  # elevation too small in Jan 2021
    '2076',  # period uncertain
    '1898',  # period unknown, single transit
    '1835',  # period uncertain
    '887',  # period uncertain
    '1796'  # transit very short (1h) and very faint object (10.7 mag V-band)
]

# make a file for staralt input
with open(output_path, 'w') as outfile:
    for i in range(len(TOI_ID)):
        if TOI_ID[i][:-3] not in TOIs_already_sorted_out:
            outstring = TOI_ID[i][:-3] + ' ' + ra_deg[i] + ' ' + dec_deg[i] + ' J2000.0\n'
            outfile.write(outstring)

# make a file with all parameters relevant for the RM estimation
# first, there are some parameters that are not available in ExoFOP but have to be searched manually (e.g. SIMBAD)
# most of them are found in the SG2 table I got from Jana
dict_manual_params = {
    'TOI_name': ['[0]equatorial rotation', '[1]K velocity amplitude', '[2]inclination angle',
                 '[3]inclination angle error', '[4]relative size of orbit', '[5]error of size of orbit',
                 '[6]other ID', '[7]V mag', '[8]impact parameter', '[9]comment'],
    '1153.01': [10.0, 13.0, 0, 0, 9.84, 3.0, 'HD106018', 8.773, 0.97, 'possible binary'],
    '1447.01': [85.0, 8.21, 0, 0, 1.34, 4.0, 'HD117173', 8.45, 0.44, '4 R_sun; possible variability, '
                                                                     'has some really strange values'],
    '1148.01': [19.49, 462.0, 89.17, 0.75, 9.95, 0.18, 'KELT-24', 8.389, 0.134, 'known planet, already on our schedule'],
    '1719.01': [8.0, 19.84, 0, 0, 10.0, 9.9, 'HD76854', 8.29, 0.09, 'possible SB2, disputed; 4 R_sun; '
                                                                    'Hebrard and Addison data'],
    '1271.01': [4.7, 37.42, 88.75, 1.0, 7.24, 0.3, 'HD118203', 8.05, 0.11, 'known planet, HD118203'],
    '1778.01': [4.0, 2.67, 0, 0, 15.0, 4.5, 'HD77946', 8.99, 0.28, 'Prio1, HARPS-N GTO time, Howard data'],
    '1415.01': [22.0, 4.59, 0, 0, 20.6, 3.5, 'HIP71409', 8.86, 0.01, 'variable star, Hebrard data'],
    '1776.01': [1.0, 1.14, 0, 0, 11.2, 2.5, 'HD95072', 8.26, 0.52, 'Prio2, Gandolfi and Howard data'],
    '2076.01': [0, 3.73, 0, 0, 59.29, 0, 'BD+40 2790', 9.139, 0, 'Period could be 2x or 6x, possible multiple system'],
    '1665.01': [8.0, 47.07, 0, 0, 3.26, 0, 'HD39315', 6.819, 0, 'crowded field, unresolved pair, '
                                                                'Buchhave data'],
    '1799.01': [2.0, 1.09, 0, 0, 19.0, 3.0, 'HD96735', 8.98, 0.14, 'Prio2, low SNR, possible systematic, Howard data'],
    '1718.01': [4.0, 5.84, 0, 0, 14.4, 2.0, 'HD58727', 8.96, 0.09, 'Hebrard data, measured mass'],
    '1773.01': [0.3, 2.98, 83.3, 0.9, 3.05, 1.2, '55 Cnc e', 5.95, 0.41, 'known planet, Howard data, '
                                                                         'THYME validation paper'],
    '1726.01': [7.5, 1.75, 89.38, 0.43, 17.15, 0.3, 'HD63433', 6.92, 0.18, 'known RM, known planet, multiple system, '
                                                                           'Howard data'],
    '1726.02': [7.5, 1.73, 89.147, 0.069, 34.8, 1.0, 'HD63433', 6.92, 0.512, 'known planet, multiple system, '
                                                                             'Howard data'],
    '1796.01': [0.3, 11.27, 86.44, 0.17, 13.9, 0.5, 'GJ436b', 10.67, 0.736, 'known planet'],
    '1821.01': [0.5, 2.9, 89.45, 0.4, 23.9, 1.33, 'HD97658', 7.76, 0.3, 'known planet, Howard data'],
    '1898.01': [7.0, 1.05, 0, 0, 1303.45, 0, 'HD83342', 7.87, 0, 'single transit, period probably too long, '
                                                                 'Dalba and Chontos data'],
    '1835.01': [1.5, 1.53, 0, 0, 16.52, 0, 'HD110067', 8.419, 0, 'period may be wrong, possible multiple system, '
                                                                 'Nowak and Palle and Luque data'],
    '1835.02': [1.5, 0.22, 0, 0, 2419.62, 0, 'HD110067', 8.419, 0, 'period may be wrong (single event), possible '
                                                                   'multiple system, Nowak and Palle and Luque data'],
    '1831.01': [35.5, 11.87, 0, 0, 1.43, 0, 'HD116206', 7.82, 0, '2.3 R_sun, some strange values, possible binary'],
    '887.01': [0, 20.91, 0, 0, 7.18, 0, 'HD43682', 8.418, 0, 'period is 3x (15-17 days), strong SB2 hints']
}

with open(RMexpect_path, 'w') as expectfile:
    for i in range(len(TOI_ID)):
        if TOI_ID[i][:-3] not in TOIs_already_sorted_out:
            # get all relevant values from the dictionary
            v_eq_rot = str(dict_manual_params[TOI_ID[i]][0])
            K_expect = str(dict_manual_params[TOI_ID[i]][1])
            inc_measured = dict_manual_params[TOI_ID[i]][2]
            inc_measured_err = dict_manual_params[TOI_ID[i]][3]
            rel_orbit = dict_manual_params[TOI_ID[i]][4]
            rel_orbit_err = dict_manual_params[TOI_ID[i]][5]
            alternat_ID = dict_manual_params[TOI_ID[i]][6]
            V_mag = str(dict_manual_params[TOI_ID[i]][7])
            b_impact = dict_manual_params[TOI_ID[i]][8]
            comment_str = dict_manual_params[TOI_ID[i]][9]

            # convert some values to the correct type for calculations
            Prad_rearth_fl = float(Prad_rearth[i])
            Prad_err_fl = float(Prad_err[i])
            Srad_rsun_fl = float(Srad_rsun[i])
            if Srad_err[i] != '':
                Srad_err_fl = float(Srad_err[i])
            else:
                Srad_err_fl = 0.1 * Srad_rsun_fl

            # for the metallicity, print a 0 if the value is not known
            if Smetal[i] == '':
                Smetal[i] = '0'
            if Smetal_err[i] == '':
                Smetal_err[i] = '0'

            # convert the stellar radius to Earth radii units to calculate the relative size of the planet
            quot_Rp_Rstar = Prad_rearth_fl / (Srad_rsun_fl * 109076.37070600962) * 1000.0
            # calculate an error for that in a very ugly way
            quot_Rp_Rstar_err = ((Prad_err_fl / Prad_rearth_fl) + (Srad_err_fl / Srad_rsun_fl)) * quot_Rp_Rstar

            # calculate an estimate for the orbital inclination from the impact parameter (for ecc=0), in degrees
            if inc_measured == 0:
                if b_impact != 0:
                    inc_estim = math.acos(b_impact / rel_orbit) * 180 / np.pi
                    inc_err = 0.05 * inc_estim
                else:
                    inc_estim = 88.0
                    inc_err = 0.15 * inc_estim
            else:
                inc_estim = inc_measured
                if inc_measured_err == 0:
                    inc_err = 0.05 * inc_estim
                else:
                    inc_err = inc_measured_err

            # if the error of the relative orbit size is not given, make an estimation
            if rel_orbit_err == 0:
                rel_orbit_err = 0.15 * rel_orbit

            expectlist = [TOI_ID[i], v_eq_rot, Smetal[i], Smetal_err[i], K_expect, T0_BJD[i], T0_err[i], P_d[i],
                          P_err[i], str(inc_estim), str(inc_err), str(quot_Rp_Rstar), str(quot_Rp_Rstar_err),
                          str(rel_orbit), str(rel_orbit_err), dur_h[i], dur_err[i], V_mag, Steff_K[i], alternat_ID,
                          TIC_ID[i][1:], comment_str]
            delimit = '|'
            expectstring = delimit.join(expectlist)
            expectfile.write(expectstring + '\n')
