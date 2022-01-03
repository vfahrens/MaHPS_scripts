import os
import datetime as dt
import calendar
import shutil
import subprocess
import astropy.io.fits as fits
import numpy as np
from operator import itemgetter
import julian
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
from matplotlib.colors import rgb2hex
import barycorrpy
from astropy.time import Time

# import statements for other python scripts
import paths_and_files as pf


# update FITS or log/comment files with rsync
def rsync_files_update(only, after, filetype):
    # if date specification is missing, skip updating
    if only is None and after is None:
        print('\n')
        print('WARNING: You did not specify any date, so I will not update any FITS or logfiles.')

    if only == 'today':
        today = dt.datetime.strftime(dt.datetime.now(), '%Y%m%d')
        sync_fits(today, filetype)
    # if an explicit date is given, use that date
    elif only is not None:
        sync_fits(only, filetype)

    # if explicit date is given for "after" option, make list of dates for syncing
    if after is not None:
        startdate = dt.datetime.strptime(after, '%Y%m%d')
        dates_lst = get_after_dateslist(startdate)
        # actually sync the requested data
        sync_fits(dates_lst, filetype)

    return


# make the rsync command line text and execute it
def sync_fits(date, filetype):
    fits_update_cmd = 'rsync -avlu'  # {}:{} {}

    # subprocess needs a list of strings, so start with the base command
    cmd_list = fits_update_cmd.split(' ')

    # distinguish between log and comment files
    if filetype == 'logs':
        years_list, directory, file1 = get_category()
        cmd_list_log = []
        cmd_list_comment = []

    # for "only" option, the date will be given as a single string
    if type(date) is str:
        if filetype == 'fits':
            cmd_list.append(pf.address_focespc + ':' + pf.fcslinks_path_focespc.format(date))
            cmd_list.append(pf.abs_path_data)
        if filetype == 'logs':
            cmd_list_log.append(pf.address_ohiaaipc + ':' + pf.log_path_ohiaaipc.format(directory[0]) +
                                '/{}/{}.{}'.format(date[:4], file1[0], date[2:]))
            cmd_list_comment.append(pf.address_ohiaaipc + ':' + pf.log_path_ohiaaipc.format(directory[1]) +
                                    '/{}/{}.{}'.format(date[:4], file1[1], date[2:]))
            cmd_list_log.append(pf.abs_path_obslog)
            cmd_list_comment.append(pf.abs_path_obslog)
            cmd_list_log = cmd_list + cmd_list_log
            cmd_list_comment = cmd_list + cmd_list_comment

    # for "after" option, sync the whole list of dates with one command
    elif type(date) is list:
        if filetype == 'fits':
            for date_str in date:
                cmd_list.append(pf.address_focespc + ':' + pf.fcslinks_path_focespc.format(date_str))
            cmd_list.append(pf.abs_path_data)
        if filetype == 'logs':
            for date_str in date:
                cmd_list_log.append(pf.address_ohiaaipc + ':' + pf.log_path_ohiaaipc.format(directory[0]) +
                                    '/{}/{}.{}'.format(date_str[:4], file1[0], date_str[2:]))
                cmd_list_comment.append(pf.address_ohiaaipc + ':' + pf.log_path_ohiaaipc.format(directory[1]) +
                                        '/{}/{}.{}'.format(date_str[:4], file1[1], date_str[2:]))
            cmd_list_log.append(pf.abs_path_obslog)
            cmd_list_comment.append(pf.abs_path_obslog)
            cmd_list_log = cmd_list + cmd_list_log
            cmd_list_comment = cmd_list + cmd_list_comment

    print('\n')
    print('I am updating the FITS files for you...')
    print('Please enter your password for ltsp01:')
    if filetype == 'fits':
        subprocess.run(cmd_list)
    if filetype == 'logs':
        subprocess.run(cmd_list_log)
        subprocess.run(cmd_list_comment)

    return


def get_after_dateslist(startdate):
    now = dt.datetime.now()
    # list of years where data is available in fcs_links
    years_data = list(range(2019, now.year + 1))

    dateslist = []

    if startdate < dt.datetime.strptime(str(20190430), '%Y%m%d'):
        print('Warning: The date you chose is before the start of automatic data collection (20190430). Expect '
              'incompatibilities and errors at all places.')

    for yr in years_data:
        # this is how to handle the year when the request starts
        if yr == startdate.year:
            # handle the rest of the starting month
            startmonth = str(startdate.year) + '{:02d}'.format(startdate.month)
            # add all days that are still left from the starting month
            if startdate.year == now.year and startdate.month == now.month:
                days = list(range(startdate.day, now.day + 1))
            else:
                end_of_month = calendar.monthrange(yr, startdate.month)[1]
                days = list(range(startdate.day, end_of_month + 1))

            # add a string to the list for all individual days of the starting month
            for single_day in days:
                dateslist.append(startmonth + '{:02d}'.format(single_day))

            # handle the other months of the starting year
            if startdate.year < now.year:
                months = list(range(startdate.month + 1, 13))
            if startdate.year == now.year and startdate.month <= now.month:
                months = list(range(startdate.month + 1, now.month + 1))
            # add a string to the list for all individual months of the starting year
            for single_month in months:
                dateslist.append(str(startdate.year) + '{:02d}'.format(single_month) + '*')

        # for all years after that, copy all folders of each year
        if yr > startdate.year:
            dateslist.append(str(yr) + '*')

    return dateslist


# distinguish between logfile and comment file paths, filenames and years
def get_category():
    now = dt.datetime.now()
    # create a list with all years from when FOCES data exist
    years_log = list(range(2016, now.year + 1))
    years_comm = list(range(2018, now.year + 1))

    dif_files = ['log', 'comments']
    years_list = []
    directory = []
    file1 = []
    for cat in dif_files:
        # check whether to use the log or comments year list
        if cat == 'log':
            years_list.append(years_log)
            directory.append('log')
            file1.append('logfile')
        if cat == 'comments':
            years_list.append(years_comm)
            directory.append('comments')
            file1.append('comments')

    return years_list, directory, file1


# extract all nights with observation data for a specific redmine project and write the dates to a file
def get_obsnights(redmine_id):
    dates_for_red = []
    # extract all observation night dates for this project
    with open(pf.grep_redID_out.format(redmine_id), 'r') as grepfile:
        with open(pf.out_obsnights.format(redmine_id), 'w') as datefile:
            for line in grepfile:
                # remove whitespaces in the beginning and end of the string
                line = line.strip()
                # remove whitespaces inside the string
                line = line.replace(' ', '')
                # split the string into its single entries
                line = line.split('|')
                if line[0][0] != '#':
                    # extract the individual observation dates from the grep results
                    file_time = dt.datetime.strptime(line[0][4:18], '%Y%m%d%H%M%S')
                    folder_date = dt.datetime.strftime(file_time, '%Y%m%d')
                    day_before = file_time - dt.timedelta(days=1)
                    str_day_before = dt.datetime.strftime(day_before, '%Y%m%d')
                    if file_time.hour > 12 and folder_date not in dates_for_red:
                        dates_for_red.append(folder_date)
                        datefile.write(folder_date + '\n')
                    elif file_time.hour <= 12 and str_day_before not in dates_for_red:
                        dates_for_red.append(str_day_before)
                        datefile.write(str_day_before + '\n')
    return


# read the dates that should be used in GAMSE from a file
def get_reductiondates(redmine_id):
    red_dates = []
    for line in open(pf.out_gamse_copy.format(redmine_id), 'r'):
        line = line.strip()
        red_dates.append(line)
    return red_dates


# automatically copy the (sorted) wavelength calibrated data to the IRAF folder
def script_copy_reduced_data(redmine_id):
    total_files_copied = 0
    # check if a folder with this redmine ID exists already in the IRAF data folder
    if not os.path.exists(pf.iraf_data_folder.format(redmine_id)):
        print('Checked, but missing!')
        os.makedirs(pf.iraf_data_folder.format(redmine_id))
    else:
        print('Checked, exists.')

    # check if a folder with this redmine ID exists already in the IRAF output folder
    if not os.path.exists(pf.iraf_output_folder.format(redmine_id)):
        print('Checked, but missing!')
        os.makedirs(pf.iraf_output_folder.format(redmine_id))
    else:
        print('Checked, exists.')

    with open(pf.frames_list.format(redmine_id, redmine_id), 'w') as framelist:
        # read the results from the grep command
        with open(pf.grep_redID_out.format(redmine_id), 'r') as grepfile:
            for line in grepfile:
                filename_used = []
                # remove whitespaces in the beginning and end of the string
                line = line.strip()
                # remove whitespaces inside the string
                line = line.replace(' ', '')
                # split the string into its single entries
                line = line.split('|')
                if line[0][0] != '#':
                    # extract the name of each file from the grep results
                    file_name = line[0]
                    filename_used.append(file_name)
                    file_time = dt.datetime.strptime(line[0][4:18], '%Y%m%d%H%M%S')

                    # get the correct date for the night of observation
                    folder_date = dt.datetime.strftime(file_time, '%Y%m%d')
                    if file_time.hour <= 12:
                        day_before = file_time - dt.timedelta(days=1)
                        str_day_before = dt.datetime.strftime(day_before, '%Y%m%d')
                        folder_date = str_day_before
                    # get the path of the .tab file corresponding to the .fits file
                    data_folder_path = os.path.join(pf.abs_path_data, folder_date)
                    tab_file_path = os.path.join(data_folder_path, file_name + '.tab')

                    # from the .tab file extract the correct link name of the raw frame
                    with open(tab_file_path, 'r') as tabfile:
                        for linex in tabfile:
                            if 'LINKNAME' in linex:
                                linex = linex.strip()
                                linex = linex.split('=')
                                raw_name = linex[1]
                                filename_used.append(raw_name)

                    # generate the file name for the reduced frame and copy it
                    red_name = raw_name[:-5] + '_ods.fits'
                    result_file_path = os.path.join(pf.gamse_results_folder.format(folder_date), red_name)
                    copy_destination_path = os.path.join(pf.iraf_data_folder.format(redmine_id), red_name)

                    try:
                        shutil.copy(result_file_path, copy_destination_path)
                        total_files_copied += 1
                        framelist.write(str(filename_used[0]) + ' ' + str(filename_used[1]) + '\n')
                    except FileNotFoundError:
                        print('WARNING: File {} does not exist in the onedspec result.'.format(raw_name))

        print('Successfully copied {} files!'.format(total_files_copied))

    return


# make a list of all single extensions of the template file
def make_template_list(fname_template, redmine_id, template_orders):
    # define a default for the complete filename of the template
    fname_temp_long = fname_template

    # list all files in the output folder to find the complete name of the template file
    fname_lst = sorted(os.listdir(pf.iraf_output_folder.format(redmine_id)))
    for fname in fname_lst:
        # check if a filename starts with the string entered as template name
        if fname[:len(fname_template)] == fname_template:
            fname_temp_long = fname

    # write a list of all extensions of the template file
    with open(pf.template_list.format(redmine_id, redmine_id), 'w') as template_file:
        for index, ordnum in enumerate(template_orders):
            template_file.write(fname_temp_long + '[{}]\n'.format(str(index + 1)))

    return fname_temp_long


# get the number of orders that are in each IRAF-converted file (e.g. if only comb orders are used)
def get_number_of_orders(redmine_id):
    # list all files that were converted to IRAF format
    fname_lst = sorted(os.listdir(pf.iraf_output_folder.format(redmine_id)))
    order_numbers_dict = {}

    for fname in fname_lst:
        phys_ords_used = []
        ext_numbers = []
        # only use fits files, ignore the rest
        if fname[-14:] != '_ods_fred.fits':
            continue

        # open the fits file
        open_file = os.path.join(pf.iraf_output_folder.format(redmine_id), fname)
        with fits.open(open_file) as hdu_list:
            # get the number of orders from the length of the HDU list (subtract the empty Primary HDU)
            num_of_orders = len(hdu_list) - 1
            for hdu_num in range(1, len(hdu_list)):
                head = hdu_list[hdu_num].header
                phys_ords_used.append(head['PHYSORD'])
                ext_numbers.append(hdu_num)
        order_numbers_dict[fname[:13] + '_num_ords'] = num_of_orders
        order_numbers_dict[fname[:13] + '_phys_ords'] = phys_ords_used

    return order_numbers_dict


# generate lists with all spectra, sorted by orders
def make_orderlists(redmine_id, used_orders_dict):
    # get a list of all files in the folder
    fname_lst = sorted(os.listdir(pf.iraf_output_folder.format(redmine_id)))

    # delete all orderlist files that are already present
    for fname in fname_lst:
        if fname[:9] == 'fxcor_ord' and fname[-4:] == '.lis':
            path_of_list = os.path.join(pf.iraf_output_folder.format(redmine_id), fname)
            if os.path.exists(path_of_list):
                os.remove(path_of_list)

    # make a list of all filenames again, now after deleting some files
    fname_lst = sorted(os.listdir(pf.iraf_output_folder.format(redmine_id)))

    # remove all files from the list that are not the FITS files that should be used for fxcor
    other_files = []
    for fname in fname_lst:
        if fname[-14:] != '_ods_fred.fits':
            other_files.append(fname)
    for ff in range(len(other_files)):
        fname_lst.remove(other_files[ff])

    # save all the used filenames also in a file
    frames_list = os.path.join(pf.iraf_output_folder.format(redmine_id), pf.all_used_frames.format(redmine_id))
    with open(frames_list, 'w') as used_files_list:
        for filename in fname_lst:
            used_files_list.write(filename + '\n')

    # generate a file for each physical order that lists all the frames containing data of that order
    orderlists_path = os.path.join(pf.iraf_output_folder.format(redmine_id), 'fxcor_ord{}.lis')
    for fname in fname_lst:
        all_used_orders = used_orders_dict[fname[:13] + '_phys_ords']
        for indx, phys_ord in enumerate(all_used_orders):
            with open(orderlists_path.format(str(phys_ord)), 'a+') as ordlis:
                # for fxcor identification use the index in the MEF, not the physical order number!
                ordlis.write(fname + '[{}]\n'.format(str(int(indx) + 1)))

    return


# make a script to generate a input cl file for fxcor
def make_script_fxcor(redmine_id, template_name, output_name, template_orders, template_harps=False):
    fxcor_script_list = os.path.join(pf.iraf_output_folder.format(redmine_id), pf.fxcor_script)

    # define the default CCF options when using a template observed with FOCES
    contin_opt = 'both'
    rebin_opt = 'template'
    rsample_opt = 'p200-1848'  # 'p150-1998'
    window_opt = 'INDEF'
    wincenter_opt = 'INDEF'

    # define the CCF options needed to use the HARPS mask templates
    if template_harps:
        contin_opt = 'object'
        rebin_opt = 'object'
        rsample_opt = '*'
        window_opt = '200.0'
        wincenter_opt = '0.0'

    with open(fxcor_script_list, 'w') as script_fxcor:
        for indx, tempord in enumerate(template_orders):
            if tempord == 86:
                script_fxcor.write(
                    'fxcor @fxcor_ord{}.lis {}[{}] output={} continuum={} rebin={} osample=p200-1848 rsample={} '
                    'function=gaussian width=15.0 window={} wincenter={} interactive=no'
                    '\n'.format(tempord, template_name, str(indx + 1), output_name, contin_opt, rebin_opt, rsample_opt,
                                window_opt, wincenter_opt))
            else:
                script_fxcor.write(
                    'fxcor @fxcor_ord{}.lis {}[{}] output={} continuum={} rebin={} osample=p200-1848 rsample={} '
                    'function=gaussian width=15.0 window={} wincenter={} interactive=no'
                    '\n'.format(tempord, template_name, str(indx + 1), output_name, contin_opt, rebin_opt, rsample_opt,
                                window_opt, wincenter_opt))

    return


# get the RVs (VHELIO, VREL) and RVerrs from the image header and fxcor result file
def get_rvs(redmine_id, fxcor_outname, template_orders, object_orders_dict):
    fname_lst = sorted(os.listdir(pf.iraf_output_folder.format(redmine_id)))
    # get all physical orders for which a template spectrum exists and therefore a CCF was calculated

    with open(pf.out_RVs_single.format(redmine_id), 'w') as outfile:
        for fname in fname_lst:
            # only use the science fiber frames
            if fname[-14:] != '_ods_fred.fits':
                continue

            # get the RV error (converted to m/s) and VREL for one input FITS file from the fxcor result file
            fxcor_output = os.path.join(pf.iraf_output_folder.format(redmine_id), fxcor_outname + '.txt')
            with open(fxcor_output, 'r') as fxfile:
                rv_err_rel_dict = {}
                for line in fxfile:
                    line = line.split()
                    # for each order of the template, get the CCF result
                    for ordnum in template_orders:
                        # find the correct index for this physical order in the object MEF file
                        try:
                            indx = object_orders_dict[fname[:13] + '_phys_ords'].index(ordnum)
                        except ValueError:
                            print('No object order exists for template order {}. '.format(ordnum))
                            continue

                        fname_ord = fname + '[{}]'.format(str(indx + 1))
                        if fname_ord in line and line[-1] != 'INDEF':
                            rv_err = float(line[-1]) * 1000.0
                            rv_rel = float(line[-3]) * 1000.0
                            rv_err_rel_dict['rv_err_{}'.format(ordnum)] = rv_err
                            rv_err_rel_dict['rv_rel_{}'.format(ordnum)] = rv_rel
                        if fname_ord in line and line[4][-3:] == 'cen':
                            rv_err = 1.0
                            rv_rel = float(line[-3]) * 1000.0
                            rv_err_rel_dict['rv_err_{}'.format(ordnum)] = rv_err
                            rv_err_rel_dict['rv_rel_{}'.format(ordnum)] = rv_rel

            # get the RV (converted to m/s) and the physical order number from the header
            open_filepath = os.path.join(pf.iraf_output_folder.format(redmine_id), fname)
            print('Extracting RVs from {}.'.format(fname))
            frame_id = fname[:8] + fname[9:13]

            with fits.open(open_filepath) as datei:
                header = datei[0].header

                for ordnum in template_orders:
                    # find the correct index for this physical order in the object MEF file
                    try:
                        indx = object_orders_dict[fname[:13] + '_phys_ords'].index(ordnum)
                    except ValueError:
                        print('No object order exists for this template order. ')

                    head_ord = datei[indx + 1].header
                    phys_ord = head_ord['PHYSORD']
                    if phys_ord != ordnum:
                        print('WARNING: Something went wrong with the order identification during RV extraction! ')
                    if 'VHELIO' in head_ord:
                        date_str = header['UTMID']
                        if len(date_str) == 19:
                            date_dt = dt.datetime.strptime(date_str, '%Y-%m-%dT%H:%M:%S')
                        elif len(date_str) > 19:
                            date_dt = dt.datetime.strptime(date_str, '%Y-%m-%dT%H:%M:%S.%f')
                        else:
                            print('Warning: Date {} has unexpected format.'.format(date_str))
                        date = julian.to_jd(date_dt, fmt='jd')
                        # get the RV corrected for heliocentric velocity from the header
                        rv_value = head_ord['VHELIO'] * 1000.0
                        # get the observed RV without heliocentric correction of the object
                        # (but with heliocentric correction of the template)
                        v_obs = head_ord['VOBS'] * 1000.0
                        # get the heliocentric julian date (old date format, now BJD)
                        hjd_head = head_ord['HJD']

                        # save all the results for the single orders to a file
                        if 'rv_err_{}'.format(ordnum) in rv_err_rel_dict \
                                and 'rv_rel_{}'.format(ordnum) in rv_err_rel_dict:
                            output_singleorders = frame_id + ' ' + str(date) + ' ' + str(hjd_head) + ' ' + \
                                                  str(rv_value) + ' ' + \
                                                  str(rv_err_rel_dict['rv_err_{}'.format(ordnum)]) + ' ' + \
                                                  str(rv_err_rel_dict['rv_rel_{}'.format(ordnum)]) + ' ' + str(v_obs) \
                                                  + ' ' + str(phys_ord) + '\n'
                            # print(output_singleorders)
                            outfile.write(output_singleorders)

    return


# get the radial velocities as object and telluric velocities
def split_rvs_tel(redmine_id, all_ords=False):
    bad_orders = []
    rvs_fromfile = []
    tellurics_fromfile = []

    # if the distinction between telluric and "science" order should be made, do the following
    if not all_ords:
        # define the orders that should be used as telluric anchor
        real_tellurics = [97, 83, 82, 79, 78, 70, 69]

        # read the orders that should be excluded from the RV analysis from the file
        with open(pf.input_tel_orders, 'r') as telluric_file:
            for linet in telluric_file:
                linet = linet.strip()
                bad_orders.append(int(linet))

        # read all RV results for each frame from the single order file
        with open(pf.out_RVs_single.format(redmine_id), 'r') as infile:
            for line in infile:
                line = line.split()
                # only use physical orders that have not been excluded in the tellurics file
                if int(line[-1]) not in bad_orders:
                    # save the whole line in a larger array with all observation dates
                    rvs_fromfile.append(line)
                if int(line[-1]) in real_tellurics:
                    tellurics_fromfile.append(line)

    # if all orders should be analyzed in the same way, choose this option
    else:
        # read all RV results for each frame from the single order file
        with open(pf.out_RVs_single.format(redmine_id), 'r') as infile:
            for line in infile:
                line = line.split()
                # save the whole line in a larger array with all observation dates
                rvs_fromfile.append(line)

    rvs_single_array = make_rv_array(rvs_fromfile)
    tels_single_array = make_rv_array(tellurics_fromfile)
    return rvs_single_array, tels_single_array


# sort, transpose and format RVs for further use
def make_rv_array(rv_list):
    # sort all RV results by date and order
    rvs_eachdate = sorted(rv_list, key=itemgetter(0))
    # convert that array to a useful format for numpy
    rvs_eachdate = np.transpose(rvs_eachdate)
    rvs_eachdate = np.asarray(rvs_eachdate).astype(np.float)
    return rvs_eachdate


# compute the median of the measured single-order RVs and RV errors
def rv_and_err_median(rvs_array, rv_type, no_bary=False):
    # for normal "science" RVs
    if rv_type != 'tel':
        # use RV with heliocentric correction by barycorrpy
        if not no_bary:
            med_rv = np.median(rvs_array[2])
            # use the fxcor error to get a median error of all RV values
            med_err = np.median(rvs_array[3])
        # or those that were not corrected, if requested
        else:
            med_rv = np.median(rvs_array[-3])
            # use the fxcor error to get a median error of all RV values
            med_err = np.median(rvs_array[-4])
    # for tellurics, use relative velocity without heliocentric correction
    else:
        med_rv = np.median(rvs_array[-3])
        # use the fxcor error to get a median error of all RV values
        med_err = np.median(rvs_array[-4])

    return med_rv, med_err


# get the data for one specific observation date and compute the weighted mean of the RV and the RV error
def rv_weightedmean(redmine_id, rvs_array, med_rv, med_err, rv_type, no_bary=False):
    if rv_type != 'tel':
        rv_type_long = 'object'
    else:
        rv_type_long = 'telluric'
    print('Computing weighted averages for {} orders:'.format(rv_type_long))

    all_stds = []
    error_limit = input('Please give a limit for the max. allowed RV error in m/s: (e.g.: 29.0) ')

    if rv_type != 'tel':
        output_file = pf.out_RVs_weighted.format(redmine_id)
    else:
        output_file = pf.out_tels_weighted.format(redmine_id)
    with open(output_file, 'w') as out2file:
        for fileid in set(rvs_array[0]):
            vels_onedate = []
            v_err = []
            for j in range(len(rvs_array[0])):
                # only use the rows of the array that contain RV data from one observation date
                if rvs_array[0, j] == fileid:
                    date_jd = rvs_array[1, j]
                    if rv_type != 'tel':
                        # if the barycentric correction was made, use the corrected values
                        if not no_bary:
                            vels_onedate.append(rvs_array[2, j])
                            # if the RV error given by fxcor is zero, use the median RV error instead
                            if rvs_array[3, j] != 0.0:
                                v_err.append(rvs_array[3, j])
                            else:
                                v_err.append(med_err)
                        else:
                            # for not corrected RVs, use the relative RV without heliocentric correction
                            vels_onedate.append(rvs_array[-3, j])
                            # if the RV error given by fxcor is zero, use the median RV error instead
                            if rvs_array[-4, j] != 0.0:
                                v_err.append(rvs_array[-4, j])
                            else:
                                v_err.append(med_err)
                    else:
                        # for tellurics, use the relative RV without heliocentric correction
                        vels_onedate.append(rvs_array[-3, j])
                        # if the RV error given by fxcor is zero, use the median RV error instead
                        if rvs_array[-4, j] != 0.0:
                            v_err.append(rvs_array[-4, j])
                        else:
                            v_err.append(med_err)

            # compute the weighted average for that date and use the median RV as zero-point correction
            rv_weightmean = np.average(vels_onedate, weights=(1 / np.abs(v_err)))
            if rv_type != 'tel':
                rv_weightmean = rv_weightmean - med_rv
            # compute the RV error across the orders and put it in a list of RV errors
            rv_std = np.std(vels_onedate) * np.sqrt(2) / np.sqrt(len(vels_onedate))
            all_stds.append(rv_std)

            # write the results to a file, if the data point is good,
            # which means that the error is below a certain limit
            if rv_std < float(error_limit):
                results = str(int(fileid)) + ' ' + str(date_jd) + ' ' + str(rv_weightmean) + ' ' + str(rv_std) + '\n'
                out2file.write(results)
            else:
                date_norm = julian.from_jd(rvs_array[1, j], fmt='jd')
                print('WARNING: Date {} has larger errors than {} m/s.'.format(date_norm, error_limit))

    return all_stds


# read the weighted RV results from the file again to fix missing RV error values
def fix_missing_errors(redmine_id, rv_type, all_stds):
    rv_results = []
    if rv_type != 'tel':
        input_file = pf.out_RVs_weighted.format(redmine_id)
    else:
        input_file = pf.out_tels_weighted.format(redmine_id)

    with open(input_file, 'r') as in2file:
        for line2 in in2file:
            line2 = line2.split()
            # check if the cross-order RV error has a reasonable value, this is not the case e.g. for the template
            # a value of 0.1 m/s cross-order RV error is probably never possible with FOCES
            # replace bad RV errors with the median of all other RV errors
            if float(line2[3]) < 0.1:
                line2[3] = str(np.median(all_stds))
                rv_results.append(line2)
            else:
                rv_results.append(line2)

    rv_tofile = sorted(rv_results, key=itemgetter(0))
    # save all RV results with the now corrected error to the file again
    rv_tofile = np.transpose(rv_tofile)

    if rv_type != 'tel':
        out2_filepath = pf.out_RVs_weighted.format(redmine_id)
    else:
        out2_filepath = pf.out_tels_weighted.format(redmine_id)
    with open(out2_filepath, 'w') as out2file_corr:
        for m in range(len(rv_tofile[0])):
            results_corr = str(rv_tofile[0, m]) + ' ' + str(rv_tofile[1, m]) + ' ' + str(rv_tofile[2, m]) + ' ' + \
                           str(rv_tofile[3, m]) + '\n'
            out2file_corr.write(results_corr)

    return rv_tofile


# save the RV results corrected with the telluric shift to a file
def get_tel_correction(redmine_id, rvs_fixerr, tel_fixerr):
    with open(pf.out_RVs_telcorr.format(redmine_id), 'w') as out4file_corr:
        for m in range(len(tel_fixerr[0])):
            results_corr = str(tel_fixerr[0, m]) + ' ' + str(tel_fixerr[1, m]) + ' ' + str(
                np.float(rvs_fixerr[2, m]) - np.float(tel_fixerr[2, m])) + ' ' + str(tel_fixerr[3, m]) + '\n'
            out4file_corr.write(results_corr)

    return


# plot the RV results for all orders for each frame
def plot_single_orders(redmine_id, no_bary=False):
    dates_rv_array = []
    # choose the input file depending on whether a barycentric correction was made or not
    if not no_bary:
        single_file_in = pf.out_RVs_abc_single.format(redmine_id)
    else:
        single_file_in = pf.out_RVs_single.format(redmine_id)

    # read the single order RVs from the file
    with open(single_file_in, 'r') as singleorderfile:
        for line in singleorderfile:
            line = line.split()
            dates_rv_array.append(line)

    # convert that array to a useful format
    dates_rv_array = make_rv_array(dates_rv_array)

    # make a plot for each different frame in the array
    plots_list = []
    labels_list = []
    for frameid in set(dates_rv_array[0]):
        order_list = []
        rv_list = []
        err_list = []
        for g in range(len(dates_rv_array[0])):
            if dates_rv_array[0, g] == frameid:
                order_list.append(dates_rv_array[-1, g])
                rv_list.append(dates_rv_array[-3, g])

                # the position of the RV error in the array is different depending on whether a barycentric
                # correction was made
                if no_bary:
                    err_list.append(dates_rv_array[-4, g])
                else:
                    err_list.append(dates_rv_array[-2, g])

        rv_rms = np.sqrt(np.mean(np.array(rv_list) ** 2))
        plots_list.append([order_list, rv_list, err_list])
        labels_list.append(str(frameid) + ' med: {:.4} rms: {:.4}'.format(np.median(rv_list), rv_rms))

        # onedate_norm = julian.from_jd(onedate, fmt='jd')
        # date_out = dt.datetime.strftime(onedate_norm, '%m.%d_%H:%M:%S')

    # plot the whole thing
    fig = plt.figure()
    # plt.errorbar(order_list, rv_list, yerr=err_list, fmt='o', label=label_med, alpha=0.5)

    for i in range(len(plots_list)):
        plt.errorbar(plots_list[i][0], plots_list[i][1], yerr=plots_list[i][2], fmt='o', alpha=0.5)  # , label=labels_list[i],
                     # alpha=0.5)
        plt.hlines(np.median(plots_list[i][1]), min(plots_list[i][0]), max(plots_list[i][0]), lw=2, alpha=0.5,
                   color=plt.gca().lines[-1].get_color())

    # plt.hlines(np.median(rv_list), min(order_list), max(order_list), lw=2)
    plt.xlabel('# of physical order')
    plt.ylabel('RV in m/s')
    plt.legend()
    plt.show()

    return


# extract a specific value from the logfile for plotting against RVs
def extract_nonrv_data(redmine_id, want_value, pos, filetype):
    # read the results from the grep command
    with open(pf.grep_redID_out.format(redmine_id), 'r') as grepfile:
        string_with_value = []
        filename_used = []
        for line in grepfile:
            # remove whitespaces in the beginning and end of the string
            line = line.strip()
            # remove whitespaces inside the string
            line = line.replace(' ', '')
            # split the string into its single entries
            line = line.split('|')
            if line[0][0] != '#':
                # extract the name of each file from the grep results
                file_name = line[0]

                if filetype == 'log':
                    value = line[pos]
                    if want_value == 'ra' or want_value == 'dec':
                        value = value.split(':')
                        value = float(value[0]) + float(value[1]) / 60 + float(value[2]) / 3600
                    new_entry = [file_name, value]
                    string_with_value.append(new_entry)
                if filetype == 'tab':
                    filename_used.append(file_name)
                    file_time = dt.datetime.strptime(line[0][4:18], '%Y%m%d%H%M%S')

                    # get the correct date for the night of observation
                    folder_date = dt.datetime.strftime(file_time, '%Y%m%d')
                    if file_time.hour <= 12:
                        day_before = file_time - dt.timedelta(days=1)
                        str_day_before = dt.datetime.strftime(day_before, '%Y%m%d')
                        folder_date = str_day_before
                    # get the path of the .tab file corresponding to the .fits file
                    data_folder_path = os.path.join(pf.abs_path_data, folder_date)
                    tab_file_path = os.path.join(data_folder_path, file_name + '.tab')

                    # from the .tab file extract the desired meta data
                    with open(tab_file_path, 'r') as tabfile:
                        for linex in tabfile:
                            if pos in linex:
                                linex = linex.strip()
                                linex = linex.split('=')
                                tab_data = linex[1]
                                new_entry = [file_name, tab_data]
                                string_with_value.append(new_entry)

    with open(pf.out_RVs_weighted.format(redmine_id), 'r') as rvoutfile:
        julian_dates = []
        for line in rvoutfile:
            line = line.strip()
            line = line.split(' ')
            julian_dates.append(line[0])

    string_with_value = sorted(string_with_value, key=itemgetter(0))
    julian_dates = sorted(julian_dates, key=itemgetter(0))

    with open(pf.out_nonRV_data.format(redmine_id), 'w') as nonrv_out:
        for i in range(len(julian_dates)):
            string_out = julian_dates[i] + ' ' + str(string_with_value[i][1]) + '\n'
            nonrv_out.write(string_out)

    return


# ask which type of non-RV data should be used for plotting
def get_nonrv_type():
    want_value = input('What kind of data do you want to extract from the logfile? '
                       '("list" for overview of supported data types) ')

    # all available data types and their descriptions
    dict_nonrv_types = {'ra': 'right ascension',
                        'dec': 'declination',
                        'azi': 'telescope azimuth',
                        'alt': 'telescope altitude',
                        'airmass': 'airmass at beginning of observation',
                        'posangle': 'position angle of the object on the sky',
                        'exptime': 'exposure time [s]',
                        'ut': 'UT timestamp when observation was started',
                        'temp_m1': 'M1 mirror temperature [°C]',
                        'temp_rod': 'telescope tuberod temperature [°C]',
                        'temp_m2': 'M2 mirror temperature [°C]',
                        'temp_m3': 'M3 mirror temperature [°C]',
                        'temp_fork': 'telescope fork temperature [°C]',
                        'hex_x': 'hexapod x position',
                        'hex_y': 'hexapod y position',
                        'hex_z': 'hexapod z position',
                        'hex_u': 'hexapod u position',
                        'hex_v': 'hexapod v position',
                        'dero': 'derotator absolute position (encoder)',
                        'dero_dist': 'derotator target distance',
                        'out_temp': 'meteo 5 min median temperature',
                        'out_press': 'meteo 5 min median pressure',
                        'dero_cur_off': 'current derotator offset',
                        'focus_cur_off': 'current focus offset',
                        'airmass_long': 'airmass at start of observation (higher precision)',
                        'refraction': 'pointing model correction due to atmospheric refraction',
                        'dero_off': 'derotator offset (instrumental)',
                        'temp_m1side': 'M1 mirror temperature, measured under the mirror [°C]',
                        'temp_m1bend1': 'Bending sensor 1 at M1 mirror mount [V]',
                        'temp_m1bend2': 'Bending sensor 2 at M1 mirror mount [V]',
                        'temp_m1bend3': 'Bending sensor 3 at M1 mirror mount [V]',
                        'temp_m1bend4': 'Bending sensor 4 at M1 mirror mount [V]'}

    # all data types that can be extracted from the logfiles of  one night
    dict_nonrv_logfile = {'ra': 2, 'dec': 3, 'azi': 4, 'alt': 5, 'airmass': 6, 'posangle': 7, 'exptime': 10, 'ut': 9,
                          'temp_m1': 11, 'temp_rod': 12, 'temp_m2': 13, 'temp_m3': 14, 'temp_fork': 15, 'hex_x': 16,
                          'hex_y': 17, 'hex_z': 18, 'hex_u': 19, 'hex_v': 20}

    dict_nonrv_tabfile = {'dero': 'POSITION.INSTRUMENTAL.DEROTATOR[2].REALPOS',
                          'dero_dist': 'POSITION.INSTRUMENTAL.DEROTATOR[2].TARGETDISTANCE',
                          'out_temp': 'TELESCOPE.CONFIG.ENVIRONMENT.TEMPERATURE',
                          'out_press': 'TELESCOPE.CONFIG.ENVIRONMENT.PRESSURE',
                          'dero_cur_off': 'CURRENT.DEROTATOR_OFFSET',
                          'focus_cur_off': 'CURRENT.FOCUS_OFFSET',
                          'airmass_long': 'CURRENT.OBJECT.HORIZONTAL.AIR_MASS',
                          'refraction': 'CURRENT.OBJECT.HORIZONTAL.REFRACTION',
                          'dero_off': 'POSITION.INSTRUMENTAL.DEROTATOR[2].OFFSET',
                          'temp_m1side': 'AUXILIARY.SENSOR[4].VALUE',
                          'temp_m1bend1': 'AUXILIARY.SENSOR[15].VALUE',
                          'temp_m1bend2': 'AUXILIARY.SENSOR[16].VALUE',
                          'temp_m1bend3': 'AUXILIARY.SENSOR[17].VALUE',
                          'temp_m1bend4': 'AUXILIARY.SENSOR[18].VALUE'}

    if want_value == 'list':
        for d in dict_nonrv_types:
            print(d + ': ' + dict_nonrv_types[d])
        want_value = input('What kind of data do you want to extract from the logfile? ')

    if want_value in dict_nonrv_logfile:
        pos = dict_nonrv_logfile[want_value]
        filetype = 'log'
    if want_value in dict_nonrv_tabfile:
        pos = dict_nonrv_tabfile[want_value]
        filetype = 'tab'

    return want_value, pos, filetype


# convert a date in ISO string format to a datetime object
def date_iso_to_dt(date_iso):
    date_dt = dt.datetime.fromisoformat(date_iso)

    return date_dt


# convert a date in ISO string format to JD
def date_iso_to_jd(date_iso):
    # convert the string to an astropy Time object
    date_ast = Time(date_iso, format='isot', scale='utc')
    # convert that to JD format with numpy.longdouble precision
    date_jd = date_ast.to_value('jd', 'long')

    return date_jd


# # extract the values needed for the barycentric correction with barycorrpy from the header
# def extract_head_barycorr(redmine_id, template_orders):
#
#     # define the folder for reading the data
#     path = pf.iraf_output_folder.format(redmine_id)
#     fname_lst = sorted(os.listdir(path))
#
#     vrels_iraf = []
#     # read the results that were extracted from the IRAF results file
#     with open(pf.out_RVs_single.format(redmine_id), 'r') as vrelfile:
#         for line in vrelfile:
#             line = line.strip()
#             line = line.split()
#             # put the observation date (JD), the VREL value and the order number in a list
#             datavect = [line[0], line[3], line[-1]]
#             vrels_iraf.append(datavect)
#
#     print(vrels_iraf)
#
#     results_strings = []
#
#     for fname in fname_lst:
#
#         # only use fits files
#         if fname[-5:] != '.fits':
#             continue
#
#         print(fname)
#
#         file_in = os.path.join(path, fname)
#         # open one of the fits files
#         with fits.open(file_in) as datei:
#             # read the header of the file
#             heady = datei[0].header
#
#             # get the observation timestamp from the header
#             date_mid_iso = heady['UTMID']
#             # exptime = heady['EXPTIME']
#             # # convert the timestamp to a datetime object and calculate mid-time of observation
#             # date_mid_dt = date_iso_to_dt(date_mid_iso)
#             # date_mid_dt = date_start_dt + 0.5 * dt.timedelta(seconds=exptime)
#             # # convert the timestamps to JD
#             # date_start_jd = date_iso_to_jd(date_start_iso)
#             date_mid_jd = date_iso_to_jd(date_mid_iso)
#
#             for indx, ordnum in enumerate(template_orders):
#                 head_ord = datei[indx + 1].header
#                 phys_ord = head_ord['PHYSORD']
#                 if phys_ord != ordnum:
#                     print('WARNING: Something went wrong with the order identification during RV extraction! ')
#                 else:
#                     order = phys_ord
#
#                     # # get the systemic (radial) velocity of the object from the header in km/s
#                     # sys_vel = heady['HIERARCH ESO TEL TARG RADVEL']
#                     # # get the IRAF VOBS for that date and order in km/s
#                     # vobs_iraf_h = heady['VOBS']
#                     # # get the HARPS barycentric RV value in km/s
#                     # harps_berv = heady['HIERARCH ESO DRS BERV']
#
#                     # get the relative velocity of the observation in m/s (IRAF VREL)
#                     for i in range(len(vrels_iraf)):
#                         date_diff = np.abs(float(date_mid_jd) - float(vrels_iraf[i][0]))
#                         if date_diff < 0.000000002:
#                             if int(vrels_iraf[i][-1]) == order:
#                                 vrel_iraf_t = vrels_iraf[i][1]
#
#                     results = str(date_mid_jd) + ' ' + str(vrel_iraf_t) + ' ' + str(order) + '\n'
#
#                     results_strings.append(results)
#
#     with open(pf.params_barycorr.format(redmine_id, redmine_id), 'w') as savefile:
#         for res in range(len(results_strings)):
#             savefile.write(results_strings[res])
#
#     print('Finished extracting headers for barycorrpy.')
#
#     return


# do the barycentric correction of the single order RVs with barycorrpy
def do_barycorr(redmine_id, RVs_single_array):
    # define the geographic position of the observatory (Wendelstein 2m)
    wst_lat = 47.7036388889
    wst_lon = 12.0120555556
    wst_alt = 1838

    # ask which object is analyzed
    obj_name = input('Please give the name of the object for the barycentric correction: ')
    # define a catalog of objects for which the HIP ID numbers are known
    obj_catalog = {'ups And': 7513, '51 Peg': 113357, 'HD 189733': 98505, 'KELT-24': 52796, 'HD 210027': 109176}
    # use the HIP ID of the object instead of the string name if it is in this catalog
    if obj_name in obj_catalog:
        obj_name = obj_catalog[obj_name]

    all_bc_pars = RVs_single_array
    rvs_bc_out = []
    rvs_bc_out_strings = []

    print('Doing the barycentric correction.')

    for k in range(len(all_bc_pars[0])):
        file_id = all_bc_pars[0][k]
        date = float(all_bc_pars[1][k])
        vrel = float(all_bc_pars[-3][k])
        verr = float(all_bc_pars[-4][k])
        order = all_bc_pars[-1][k]

        if isinstance(obj_name, int):
            result = barycorrpy.get_BC_vel(JDUTC=date, hip_id=obj_name, lat=wst_lat, longi=wst_lon, alt=wst_alt,
                                           zmeas=(vrel / 299792458), leap_update=False)
        elif isinstance(obj_name, str):
            result = barycorrpy.get_BC_vel(JDUTC=date, starname=obj_name, lat=wst_lat, longi=wst_lon, alt=wst_alt,
                                           zmeas=(vrel / 299792458), leap_update=False)
        else:
            print('WARNING: Unexpected format of object ID. Please check the input for the object name. ')
        #  ephemeris='de430',

        rvs_bc_out.append([int(file_id), date, result[0][0], verr, int(order)])
        rv_bc_corr = str(int(file_id)) + ' ' + str(date) + ' ' + str(result[0][0]) + ' ' + str(verr) + ' ' + \
                     str(int(order)) + '\n'
        rvs_bc_out_strings.append(rv_bc_corr)

    # save the corrected radial velocities to a file
    with open(pf.out_RVs_abc_single.format(redmine_id), 'w') as outbcfile:
        for single_res in range(len(rvs_bc_out_strings)):
            outbcfile.write(rvs_bc_out_strings[single_res])

    print('Results of barycentric correction written to {}'.format(pf.out_RVs_abc_single.format(redmine_id)))

    # convert the results to a numpy array for further use
    rvs_bc_out_arr = make_rv_array(rvs_bc_out)

    return rvs_bc_out_arr


def plot_weighted_RVs(redmine_id):
    dates_rv_array = []
    # read the weighted RVs from the file
    with open(pf.out_RVs_weighted.format(redmine_id), 'r') as weightfile:
        for line in weightfile:
            # do not use lines that were commented out in the file of weighted averages
            if line[0] != '#':
                line = line.split()
                dates_rv_array.append(line)

    # convert that array to a useful format
    dates_rv_array = make_rv_array(dates_rv_array)

    # make a plot with all dates in the array
    dates_list = []
    rv_list = []
    err_list = []
    for g in range(len(dates_rv_array[0])):
        dates_list.append(dates_rv_array[1, g])
        rv_list.append(dates_rv_array[2, g])
        err_list.append(dates_rv_array[3, g])

    # onedate_norm = julian.from_jd(onedate, fmt='jd')
    # date_out = dt.datetime.strftime(onedate_norm, '%m.%d_%H:%M:%S')
    std_rvs = np.std(rv_list)
    med_rvs = np.median(rv_list)
    mean_rvs = np.mean(rv_list)
    label1 = redmine_id + ' std: {:.4} med: {:.4} mean: {:.4}'.format(std_rvs, med_rvs, mean_rvs)

    # plot the whole thing
    fig = plt.figure()
    plt.errorbar(dates_list, rv_list, yerr=err_list, fmt='o', label=label1, alpha=0.5)
    plt.hlines(np.median(rv_list), min(dates_list), max(dates_list), lw=2)
    plt.xlabel('time of observation')
    plt.ylabel('RV in m/s')
    plt.legend()
    plt.show()

    return


def plot_histograms(redmine_id):
    rvs_weight_array = []
    # read the weighted RVs from the file
    with open(pf.out_RVs_weighted.format(redmine_id), 'r') as weightfile:
        for line in weightfile:
            line = line.split()
            rvs_weight_array.append(line)

    # convert that array to a useful format
    rvs_weight_array = make_rv_array(rvs_weight_array)

    # make a plot with all dates in the array
    rv_list = []
    for g in range(len(rvs_weight_array[0])):
        rv_list.append(rvs_weight_array[2, g])

    std_rvs = np.std(rv_list)
    med_rvs = np.median(rv_list)
    mean_rvs = np.mean(rv_list)
    label1 = redmine_id + ' std: {:.4} med: {:.4} mean: {:.4}'.format(std_rvs, med_rvs, mean_rvs)

    # bins = 25
    bins = range(int(min(rv_list) - 1), int(max(rv_list) + 1) + 1, 1)

    # plot the whole thing
    fig1 = plt.figure()
    values, bins_w, rects = plt.hist(rv_list, bins, label=label1)

    with open(os.path.join(pf.abs_path_rvplots, '51Peg_simpoiss_weighted_100_tempSN200.txt'), 'w') as datwsave:
        for m in range(len(values)):
            datwsave.write('{} {}\n'.format(values[m], bins_w[m]))

    plt.vlines(np.median(rv_list), 0, max(values), lw=2)
    plt.xlabel('RV in m/s')
    plt.ylabel('number of results')
    plt.legend()
    fig1.savefig(os.path.join(pf.abs_path_rvplots, '51Peg_simpoiss_weighted_100_tempSN200.png'))
    plt.show()

    # here the plot of the dingle orders starts
    rvs_single_array = []
    # read the single order RVs from the file
    with open(pf.out_RVs_single.format(redmine_id), 'r') as singleorderfile:
        for line in singleorderfile:
            line = line.split()
            rvs_single_array.append(line)

    # convert that array to a useful format
    rvs_single_array = make_rv_array(rvs_single_array)

    # make a plot for each different order
    for physord in set(rvs_single_array[-1]):
        frames_list = []
        rv_single_list = []
        for g in range(len(rvs_single_array[-1])):
            if rvs_single_array[-1, g] == physord:
                frames_list.append(rvs_single_array[0, g])
                rv_single_list.append(rvs_single_array[-3, g])

        rv_rms = np.sqrt(np.mean(np.array(rv_single_list) ** 2))
        label_ord = 'Order: {} med: {:.4} rms: {:.4}'.format(int(physord), np.median(rv_single_list), rv_rms)

        bins_single = range(int(min(rv_single_list) - 1), int(max(rv_single_list) + 1) + 1, 1)

        # plot the whole thing
        fig2 = plt.figure()

        values_sing, bins_s, rects_s = plt.hist(rv_single_list, bins_single, label=label_ord, alpha=0.8)

        with open(os.path.join(pf.abs_path_rvplots, '51Peg_simpoiss_ord{}_100_tempSN200.txt'.format(int(physord))),
                  'w') as datsave:
            for n in range(len(values_sing)):
                datsave.write('{} {}\n'.format(values_sing[n], bins_s[n]))

        plt.vlines(np.median(rv_single_list), 0, max(values_sing), lw=2)
        plt.xlabel('RV in m/s')
        plt.legend()
        fig2.savefig(os.path.join(pf.abs_path_rvplots, '51Peg_simpoiss_ord{}_100_tempSN200.png'.format(int(physord))))
        # plt.show()

    return


# plot the RV results for each single order for each observation date chronologically
def plot_single_orders_chronologic(redmine_id, barycorr=False):
    dates_rv_array = []
    # read the single order RVs from the file
    # if the barycentric correction should be used, take the abc single orders file
    if barycorr:
        infiles_singleords = pf.out_RVs_abc_single.format(redmine_id)
    # if no barycentric correction is desired, read the plain single orders file
    else:
        infiles_singleords = pf.out_RVs_single.format(redmine_id)

    with open(infiles_singleords, 'r') as singleorderfile:
        for line in singleorderfile:
            line = line.split()
            dates_rv_array.append(line)

    dates_rv_array = np.asarray(dates_rv_array).astype(np.float)
    # sort all RV results by order, to have small (red) orders first
    rv_list_sort = sorted(dates_rv_array, key=itemgetter(-1))
    # convert that array to a useful format
    dates_rv_array = make_rv_array(dates_rv_array)


    # make a plot for each different order in the array
    plots_list = []
    # labels_list = []
    for orderid in sorted(set(dates_rv_array[-1])):
        frames_list = []
        dates_list = []
        rv_list = []
        err_list = []
        for g in range(len(dates_rv_array[-1])):
            if dates_rv_array[-1, g] == orderid:
                frames_list.append(dates_rv_array[0, g])
                dates_list.append(dates_rv_array[1, g])
                rv_list.append(dates_rv_array[-3, g])
                if barycorr:
                    err_list.append(dates_rv_array[-2, g])
                else:
                    err_list.append(dates_rv_array[-4, g])

                # onedate_norm = julian.from_jd(dates_rv_array[1, g], fmt='jd')
                # date_out = dt.datetime.strftime(onedate_norm, '%m.%d_%H:%M:%S')

        # rv_rms = np.sqrt(np.mean(np.array(rv_list) ** 2))
        plots_list.append([dates_list, rv_list, err_list, orderid])
        # labels_list.append(str(frameid) + ' med: {:.4} rms: {:.4}'.format(np.median(rv_list), rv_rms))


    # plot the whole thing
    fig = plt.figure()
    cmap = cm.get_cmap('jet')
    # plt.errorbar(order_list, rv_list, yerr=err_list, fmt='o', label=label_med, alpha=0.5)

    for i in range(len(plots_list)):
        color = rgb2hex(cmap(float(len(plots_list) - i) / float(len(plots_list))))
        plt.scatter(plots_list[i][0], plots_list[i][1], color=color, alpha=0.5)  # , fmt='o', # alpha=0.5)

    plt.xlabel('date of observation')
    plt.ylabel('RV in m/s')
    plt.legend()
    plt.show()

    return

#
# # copy the files that were automatically reduced by GAMSE on alioth to the IRAF data input folder
# def copy_reduced_from_alioth(redmine_id):
#     fits_update_cmd = 'rsync -avlu'  # {}:{} {}
#
#     # subprocess needs a list of strings, so start with the base command
#     cmd_list = fits_update_cmd.split(' ')
#
#     with open(pf.grep_redID_out.format(redmine_id), 'r') as grep_res:
#         for line in
#
#     # for "only" option, the date will be given as a single string
#     if type(date) is str:
#         if filetype == 'fits':
#             cmd_list.append(pf.address_focespc + ':' + pf.fcslinks_path_focespc.format(date))
#             cmd_list.append(pf.abs_path_data)
#         if filetype == 'logs':
#             cmd_list_log.append(pf.address_ohiaaipc + ':' + pf.log_path_ohiaaipc.format(directory[0]) +
#                                 '/{}/{}.{}'.format(date[:4], file1[0], date[2:]))
#             cmd_list_comment.append(pf.address_ohiaaipc + ':' + pf.log_path_ohiaaipc.format(directory[1]) +
#                                     '/{}/{}.{}'.format(date[:4], file1[1], date[2:]))
#             cmd_list_log.append(pf.abs_path_obslog)
#             cmd_list_comment.append(pf.abs_path_obslog)
#             cmd_list_log = cmd_list + cmd_list_log
#             cmd_list_comment = cmd_list + cmd_list_comment
#
#     # for "after" option, sync the whole list of dates with one command
#     elif type(date) is list:
#         if filetype == 'fits':
#             for date_str in date:
#                 cmd_list.append(pf.address_focespc + ':' + pf.fcslinks_path_focespc.format(date_str))
#             cmd_list.append(pf.abs_path_data)
#         if filetype == 'logs':
#             for date_str in date:
#                 cmd_list_log.append(pf.address_ohiaaipc + ':' + pf.log_path_ohiaaipc.format(directory[0]) +
#                                     '/{}/{}.{}'.format(date_str[:4], file1[0], date_str[2:]))
#                 cmd_list_comment.append(pf.address_ohiaaipc + ':' + pf.log_path_ohiaaipc.format(directory[1]) +
#                                         '/{}/{}.{}'.format(date_str[:4], file1[1], date_str[2:]))
#             cmd_list_log.append(pf.abs_path_obslog)
#             cmd_list_comment.append(pf.abs_path_obslog)
#             cmd_list_log = cmd_list + cmd_list_log
#             cmd_list_comment = cmd_list + cmd_list_comment
#
#     print('\n')
#     print('I am updating the FITS files for you...')
#     print('Please enter your password for ltsp01:')
#     if filetype == 'fits':
#         subprocess.run(cmd_list)
#     if filetype == 'logs':
#         subprocess.run(cmd_list_log)
#         subprocess.run(cmd_list_comment)
#
#     return

# plot_single_orders_chronologic('2864_thar', barycorr=True)
