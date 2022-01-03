import os, sys
import numpy as np
from operator import itemgetter

################
redmine_project = '3204'
# HD149027: '3207'; HD189733: '3155'; HAT-P-2: '3204'
datapath = '/mnt/e/science/RM_effect/data/HAT-P-2_nocomb_fresh/'
error_meth = 'order_scatter'  # 'ccf_error'
# the error method can be either 'order_scatter' or 'ccf_error'
excl_ords = ['87', '88', '91', '96', '97']  #'95', '86', '89', '90', '98', '99', '100', '101', '105', '106', '112', '114', '127', '132']  #'84', '85',
sigma_coeff = 2.0
iter = 2
################


# function for reading the file that contains the single order results from MARMOT
def read_marmot_singlefile(redmine_id, excluded_orders):
    dates_rv_array = []
    with open(os.path.join(datapath, 'RVs_ID{}_singleord_marmot_ccf.txt'.format(redmine_id)), 'r') as in_sing_file:
        for line in in_sing_file:
            # skip lines that are commented out
            if line[0] != '#':
                line = line.strip()
                line = line.split()

                if line[-1] not in excluded_orders:
                    dates_rv_array.append(line)

    # convert that array to a useful format
    dates_rv_array = make_rv_array(dates_rv_array)

    return dates_rv_array


# function to calculate the weighted mean and save it to a file
def calc_weighted_mean(redmine_id, rvs_array, error_method='order_scatter', clipping=False, sigma_coeff=1.0):
    #all_stds = []
    #all_err_ccfs = []

    with open(os.path.join(datapath, 'RVs_ID{}_newweighted_marmot_ccf.txt'.format(redmine_id)), 'w') as out_weight_file:
        for fileid in set(rvs_array[0]):
            rvs_onedate = []
            rv_err_onedate = []
            for j in range(len(rvs_array[0])):
                # only use the rows of the array that contain RV data from one observation date
                if rvs_array[0, j] == fileid:
                    date_jd = rvs_array[1, j]

                    # remove all orders that have nan values in rv_val or rv_err
                    if not np.isnan(rvs_array[2, j]):
                        if not np.isnan(rvs_array[3, j]):
                            rvs_onedate.append(rvs_array[2, j])
                            # print(np.isnan(rvs_array[2, j]))
                            rv_err_onedate.append(rvs_array[3, j])

            # do a sigma clipping for each single frame
            old_len = len(rvs_onedate)
            if clipping:
                for inum in range(iter):
                    medi = np.median(rvs_onedate)
                    sig = np.nanstd(rvs_onedate)
                    low_bord = medi - sigma_coeff * sig
                    up_bord = medi + sigma_coeff * sig
                    for this_rv in rvs_onedate:
                        if not low_bord <= this_rv <= up_bord:
                            indi = rvs_onedate.index(this_rv)
                            rvs_onedate.remove(this_rv)
                            rv_err_onedate.pop(indi)

                    new_len = len(rvs_onedate)
                    print(fileid, old_len, new_len, inum)

            # compute the weighted average for that date
            #rv_weightmean = np.median(rvs_onedate)
            #rv_weightmean = np.average(rvs_onedate, weights=(1. / np.array(rv_err_onedate)**2))
            rv_weightmean = np.average(rvs_onedate)
            # use the median RV as zero-point correction
            #rv_weightmean = rv_weightmean - med_rv
            # compute the RV error across the orders and put it in a list of RV errors
            rv_std = np.std(rvs_onedate) / np.sqrt(len(rvs_onedate))  #  * np.sqrt(2) is wrong here!
            #all_stds.append(rv_std)
            # compute the error from the CCF error like it is done in MARMOT
            rv_err_ccf = (np.sum(1. / np.array(rv_err_onedate) ** 2)) ** (-0.5)
            #all_err_ccfs.append(rv_err_ccf)

            # write the results to file, depending on the chosen error method
            if error_method == 'order_scatter':
                result_string = str(int(fileid)) + ' ' + str(date_jd) + ' ' + str(rv_weightmean) + ' ' + \
                                str(rv_std) + '\n'
            elif error_method == 'ccf_error':
                result_string = str(int(fileid)) + ' ' + str(date_jd) + ' ' + str(rv_weightmean) + ' ' + \
                                str(rv_err_ccf) + '\n'

            out_weight_file.write(result_string)


# sort, transpose and format RVs for further use
def make_rv_array(rv_list):
    # sort all RV results by date and order
    rvs_eachdate = sorted(rv_list, key=itemgetter(0))
    # convert that array to a useful format for numpy
    rvs_eachdate = np.transpose(rvs_eachdate)
    rvs_eachdate = np.asarray(rvs_eachdate).astype(np.float)
    return rvs_eachdate


data_array = read_marmot_singlefile(redmine_project, excl_ords)
calc_weighted_mean(redmine_project, data_array, error_meth, clipping=True, sigma_coeff=sigma_coeff)

print('Saved weighted RVs!')
