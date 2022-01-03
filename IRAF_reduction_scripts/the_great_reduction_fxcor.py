#!/usr/bin/env python3
# makes file executable from shell without the "python3" command in front

import os
import subprocess
import re  # module for regular expressions
import datetime as dt
import shutil
import argparse

# import statements for other python scripts
import shell_script_maker as shesm
import paths_and_files as pf
import small_functions as sf
from gamse_convert import iraf_converter
from gamse_convert_comb import iraf_converter_comb
from radvel_make_conffile import make_radvel_conffile

##################################################################
# general variables and option defaults

# in_date = dt.datetime.strftime(dt.datetime.now(), '%Y%m%d')
# in_opt1 = '-o'
# in_opt2 = '-ld'

###################################################################
# argument parser definitions

# create command line argument parser
parser = argparse.ArgumentParser(description='RV calculation for FOCES with fxcor: This program takes raw frames for '
                                             'a certain object of the FOCES spectrograph, performs a wavelength '
                                             'calibration with the GAMSE software, converts the files to MEF files and '
                                             'computes radial velocities by using IRAF\'s fxcor package. There are '
                                             'different options to choose which parts of the program should be '
                                             'executed. Some other spectrographs are also (partly) supported.')

# add optional command line arguments
parser.add_argument('-fu', '--fitsupdate', help='update the local raw FITS files', action='store_true')
parser.add_argument('-lu', '--logsupdate', help='update the local observation logfiles', action='store_true')
parser.add_argument('-id', '--redmine_id', help='redmine ID of the object/project that should be analyzed', type=str)
parser.add_argument('-g', '--gamse', help='perform the wavelength calibration with GAMSE', action='store_true')
parser.add_argument('-al', '--alioth', help='download the GAMSE wavelength calibrated frames from alioth',
                    action='store_true')
parser.add_argument('-i', '--iraf_conv', help='convert the gamse data to multi-extension FITS for use with IRAF',
                    action='store_true')
parser.add_argument('-fx', '--fxcor', help='compute RVs from wavelength calibrated data with IRAF\'s fxcor',
                    action='store_true')
parser.add_argument('-c', '--comb', help='use the comb corrected wavelength calibration', action='store_true')
parser.add_argument('-ex', '--extract', help='extract the RVs from the fxcor results and do the barycentric '
                                             'correction', action='store_true')
parser.add_argument('-ps', '--plot_single', help='plot the RV results for the single orders', action='store_true')
parser.add_argument('-pw', '--plot_weighted', help='plot the RV results of the weighted averages', action='store_true')
parser.add_argument('-ph', '--plot_histogram', help='plot histograms of the RV results of the single orders and '
                                                    'weighted averages (BETA! not fully automated yet)', action='store_true')
parser.add_argument('--no_barycorr', help='do not make a barycentric correction of the RV results',
                    action='store_true')
parser.add_argument('--all_orders', help='analyze all orders in the same way, do not automatically select orders, '
                                         'e.g. based on telluric contamination or comb calibration', action='store_true')
parser.add_argument('--harps', help='use this flag when working with HARPS data', action='store_true')  # not in use yet
parser.add_argument('-th', '--template_harps', help='use one of the HARPS masks as template for the CCF',
                    action='store_true')
# create a group of options that exclude their simultaneous usage
processing_dates = parser.add_mutually_exclusive_group()
processing_dates.add_argument('-o', '--only', help='process data only for the date specified, date format: YYYYMMDD or '
                                                   '"today"', type=str)
processing_dates.add_argument('-a', '--after', help='process data for all dates starting with and after the one '
                                                    'specified, date format: YYYYMMDD', type=str)

# actually parse the arguments given in the command line in the current execution of this script
args = parser.parse_args()

if args.only is None and args.after is None:
    in_opt = '-e'
    in_date = '20200101'
if args.only == 'today':
    in_opt = '-o'
    now = dt.datetime.now()
    in_date = dt.datetime.strftime(now, '%Y%m%d')
    print(in_date)
elif args.only is not None:
    in_opt = '-o'
    in_date = str(args.only)
if args.after is not None:
    in_opt = '-a'
    in_date = str(args.after)

#####################################################################
# actual calls of functions

# Welcome message for the user
print('Hello! I see you want to reduce some data.')

# check whether an update of the data and logfiles is needed
if args.fitsupdate:
    sf.rsync_files_update(args.only, args.after, filetype='fits')
if args.logsupdate:
    sf.rsync_files_update(args.only, args.after, filetype='logs')
else:
    print('\n')
    print('No FITS or log/comment files were updated.')

# # check whether an update of the data file headers is needed
# print('\n')
# yn_header = input('Do you also want to update the data file headers? ')
# if re.match(r'^y', yn_header, re.I) or re.match(r'^j', yn_header, re.I):
#     print('\n')
#     print('Please enter a date and option:')
#     date_restrict_head = input('(yyyymmdd; option: -a = after, -o = only, -e = everything)\n')
#     # split the input into its parts
#     in_date_head = date_restrict_head[:8]
#     in_opt_head = date_restrict_head[-2:]
#
#     script_add = shesm.script_add_radec(in_date_head, in_opt_head)
#     subprocess.call(['dos2unix', str(pf.script_add)])
#     print('\n')
#     print('Starting to update the headers now...')
#     subprocess.call(['bash', str(pf.script_add)])
#
# else:
#     print('\n')
#     print('No headers were updated.')
#     print('\n')

# if args.redmine_id is not None:
#     # create and execute a script for searching the redmine ID
#     script_grep_ID = shesm.script_grep_redmineid(args.redmine_id, in_date, in_opt)
#     subprocess.run(['dos2unix', str(script_grep_ID)])
#     print('\n')
#     print('Successfully created bash script for grep. Starting to search now...')
#     subprocess.run(['bash', str(script_grep_ID)])
#
#     # write all observation dates to a file
#     sf.get_obsnights(args.redmine_id)
#
#     # create and execute a script for copying the required raw data to the reduction folders
#     shesm.script_sort_for_reduction(args.redmine_id, in_date, in_opt)
#     subprocess.run(['dos2unix', str(pf.sort_copy_cmd.format(args.redmine_id))])
#     print('\n')
#     print('Successfully created bash script for copy. Starting to copy now...')
#     subprocess.run(['bash', str(pf.sort_copy_cmd.format(args.redmine_id))])
#     print('\n')
#     print('Finished copying the files to the reduction folders.')
#
# else:
#     print('\n')
#     print('No files were searched.')

if args.alioth:
    # copy the files that were automatically reduced by GAMSE on alioth to the IRAF data input folder
    sf.copy_reduced_from_alioth(args.redmine_id)
else:
    print('\n')
    print('No files were copied from alioth to the IRAF folder.')

if args.gamse:
    # read the dates that should be reduced from the file
    use_only_these_reduction_dates = sf.get_reductiondates(args.redmine_id)
    print(use_only_these_reduction_dates)

    print('\n')
    print('Started wavelength calibration...')

    # start the wavelength calibration with the newest data folder
    use_only_these_reduction_dates.sort(reverse=True)
    for i in use_only_these_reduction_dates:
        red_folder_path = os.path.join(pf.abs_path_red_gamse, 'red_{}'.format(str(i)))
        # change the working directory to the data reduction folder, needed for GAMSE calls to work properly
        # os.chdir(str(red_folder_path))
        print('\n')
        print('Data reduction started for : {}'.format(i))
        input('Please use GAMSE in the folder "red_{}"'.format(i))
        # subprocess.run('bash -c "conda activate git-gamse"', shell=True)
        # subprocess.run(['conda activate git-gamse'], shell=True)
        # subprocess.run(['gamse', 'config'])
        # subprocess.run(['emacs', 'FOCES.*.cfg'])
        # subprocess.run(['gamse', 'list'])
        # subprocess.run(['gamse', 'reduce'])
        # subprocess.run(['conda', 'deactivate'])

    sf.script_copy_reduced_data(args.redmine_id)
else:
    print('\n')
    print('No files were reduced locally with GAMSE and copied to the IRAF folder.')

# convert the GAMSE data to IRAF readable form if required
if args.iraf_conv:
    if args.comb:
        iraf_converter_comb(pf.iraf_data_folder.format(args.redmine_id), args.redmine_id)
    else:
        iraf_converter(pf.iraf_data_folder.format(args.redmine_id), args.redmine_id)

else:
    print('\n')
    print('No files were converted to IRAF format.')

# prepare data for fxcor and give instructions for IRAF execution
if args.fxcor:
    print('List of frames for this object: {}'.format(pf.all_used_frames.format(args.redmine_id)))
    template = input('Please choose a template that should be used for the cross correlation: '
                     '(e.g.: 20190903_0114) ')
    outname = input('Please give a name for the fxcor output file: (e.g.: RVs_200723) ')

    used_orders = sf.get_number_of_orders(args.redmine_id)
    sf.make_orderlists(args.redmine_id, used_orders)

    template_orders = used_orders[template + '_phys_ords']
    templ_name_long = sf.make_template_list(template, args.redmine_id, template_orders)
    sf.make_script_fxcor(args.redmine_id, templ_name_long, outname, template_orders, args.template_harps)
    # input('Want to cancel? ')

    print('Please open a terminal now and type "xterm". Then go to the new window.')
    print('Type the following command in the xterm window: \n')
    input('$ cl')
    print('After IRAF has loaded, please execute the following commands: \n')
    print('vocl> reset obsdb=home$obsdb.dat')
    input('vocl> rv')
    print('Now navigate to the folder containing the data: \n')
    input('rv> cd {}'.format(pf.iraf_output_folder.format(args.redmine_id)))
    # print('do the heliocentric correction for the template spectra: enter the given command or open '
    #       '"epar rvcorrect" and put the template filename list as entry for '
    #       '"images": "images = @templates_ID{}.lis"\n'.format(args.redmine_id))
    # input('rv> rvcorrect images=@templates_ID{}.lis'.format(args.redmine_id))
    print('please check with "epar fxcor" if the name of the save file ({}) and the CCF settings are '
          'correct'.format(outname))
    print('finally, execute fxcor: (may take a while, please hit ENTER here when finished)\n')
    input('cl < fxcor_with_lists.cl')

if args.extract:
    if not args.fxcor:
        used_orders = sf.get_number_of_orders(args.redmine_id)
        outname = input('Please give the name of the fxcor output file: (e.g.: RVs_200723) ')
        template = input('Please choose a template that should be used for the cross correlation: '
                         '(e.g.: 20190903_0114) ')
        template_orders = used_orders[template + '_phys_ords']

    # extract the RVs from the fxcor results
    sf.get_rvs(args.redmine_id, outname, template_orders, used_orders)
    # classify orders as non-contaminated with tellurics or as "reliable" telluric orders
    RVs_single, tels_single = sf.split_rvs_tel(args.redmine_id, args.all_orders)

    # do the barycentric correction of the RV results if not explicitly requested otherwise
    if not args.no_barycorr:
        # overwrite the single order RV results with the barycentrically corrected results, because the uncorrected
        # RVs are not needed anymore after this step
        RVs_single = sf.do_barycorr(args.redmine_id, RVs_single)

    # calculate the median of all RV values and all RV errors
    RVs_single_med, RVerr_single_med = sf.rv_and_err_median(RVs_single, 'obj', args.no_barycorr)
    if len(tels_single) != 0:
        tels_single_med, telserr_single_med = sf.rv_and_err_median(tels_single, 'tel', args.no_barycorr)

    # calculate the weighted mean of all relevant orders
    RVs_stds = sf.rv_weightedmean(args.redmine_id, RVs_single, RVs_single_med, RVerr_single_med,
                                  'obj', args.no_barycorr)
    RVs_fixerr = sf.fix_missing_errors(args.redmine_id, 'obj', RVs_stds)
    if len(tels_single) != 0:
        tel_stds = sf.rv_weightedmean(args.redmine_id, tels_single, tels_single_med, telserr_single_med,
                                      'tel', args.no_barycorr)
        tel_fixerr = sf.fix_missing_errors(args.redmine_id, 'tel', tel_stds)
        sf.get_tel_correction(args.redmine_id, RVs_fixerr, tel_fixerr)

if args.plot_single:
    sf.plot_single_orders(args.redmine_id, args.no_barycorr)
    # yn_plot_singleorders = input('Do you want to plot the RV results for the single orders? ')
    # if re.match(r'^y', yn_plot_singleorders, re.I) or re.match(r'^j', yn_plot_singleorders, re.I):
    #     print('\n')
    #     sf.plot_single_orders(args.redmine_id)

if args.plot_weighted:
    sf.plot_weighted_RVs(args.redmine_id)

if args.plot_histogram:
    sf.plot_histograms(args.redmine_id)

# make a plot of the literature values compared to the FOCES data
print('\n')
yn_RV_complit = input('Do you want to plot the RV results and compare them to the literature? ')

if re.match(r'^y', yn_RV_complit, re.I) or re.match(r'^j', yn_RV_complit, re.I):
    print('\n')
    redmine_id = input('The redmine ID is needed once again: ')
    n_cand = input('Please enter the number of planet candidates: ')
    inst_list = input('Give a list of the instruments that were used (separate with space): ')
    inst_list = inst_list.split(' ')
    config_file = make_radvel_conffile(redmine_id, n_cand, inst_list)
    os.chdir(str(pf.location.resolve()))
    subprocess.call(['myradvel', 'fit', '-s', str(config_file), '-d', str(pf.abs_path_rvplots)])
    subprocess.call(['myradvel', 'plot', '-t', 'rv', '-s', str(config_file), '-d', str(pf.abs_path_rvplots)])

# make a plot of the RV values compared to other data (e.g. airmass)
print('\n')
yn_RV_nonRV = input('Do you want to plot the RV results and compare them to some other non-RV data? ')

if re.match(r'^y', yn_RV_nonRV, re.I) or re.match(r'^j', yn_RV_nonRV, re.I):
    print('\n')
    redmine_id = input('The redmine ID is needed once again: ')
    want_value, pos, filetype = sf.get_nonrv_type()
    sf.extract_nonrv_data(redmine_id, want_value, pos, filetype)
    config_file = pf.radvel_config.format(redmine_id)
    os.chdir(str(pf.location.resolve()))
    subprocess.call(['myradvel', 'fit', '-s', str(config_file), '-d', str(pf.abs_path_rvplots)])
    subprocess.call(['myradvel', 'plot', '-t', 'nonrv', '-s', str(config_file), '-d', str(pf.abs_path_rvplots)])

# do the fxcor reduction manually
print('\n')
input('Please hit ENTER when finished.')

# # plot the data produced by fxcor
# print('\n')
# yn_plotdata = input('Do you want to plot the RVs you got from fxcor? ')
#
# if re.match(r'^y', yn_wvcal_copy, re.I) or re.match(r'^j', yn_wvcal_copy, re.I):
#     print('\n')
#     redmine_id = input('I am sorry to ask for the redmine ID yet again: ')
#     plot_type = input('What kind of plot do you want to create? (single: only one dataset, '
#                       'double: comparison of two datasets)')
#
#     if re.search(r'single', plot_type, re.I):
#
#
#
