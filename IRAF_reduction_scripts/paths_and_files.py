import os
from pathlib import Path

# define the required paths
path_scripts = 'scripts/'
path_output = 'output/'
path_obslog_local = '../../../logfiles/observations'
path_data_local = '../../../FOCES_data'
path_reduce_gamse = '../../../GAMSE'
path_IRAF_datainput = 'data/'
path_gamse_results = 'onedspec'
path_gamse_reduce = 'red_{}'
path_rv_results = 'rv_results/'
path_rv_plots = 'RV_plots/'
path_smd = 'ID{}_same_min_diff/'

# definition of many filenames
file_script_USM = 'sync_obslogfiles_USM.sh'
data_script_USM = 'sync_datafiles_USM.sh'
file_script_local = 'sync_obslogfiles_local.sh'
data_script_local = 'sync_datafiles_local.sh'
add_header_script = 'add_header_entries.sh'
grep_redmineID_script = 'ID{}_grep.sh'
grep_redmineID_results = 'ID{}_grep.txt'
sort_copy_gamse_script = 'copy_ID{}_to_gamse.sh'
observed_nights_list = 'ID{}_observed_nights.txt'
copy_dates_gamse = 'copy_ID{}_to_gamse.txt'
copy_wvcal_script = 'copy_wvcal_ID{}_to_IRAF.sh'
recipe_orderlists = 'recipe_make_orderlists.sh'
recipe_cl_fxcor = 'recipe_make_fxcor_with_lists.sh'
all_used_frames = 'used_frames_ID{}.txt'
RVs_single_orders = 'RVs_ID{}_single_orders.txt'
RVs_weighted = 'RVs_ID{}_weighted.txt'
tels_weighted = 'tels_ID{}_weighted.txt'
RVs_telcorr = 'RVs_ID{}_telcorr.txt'
ords_tellurics = 'orders_with_telluric_contamination.txt'
RVs_compare = 'RVs_ID{}_telcorr_compare.txt'
literature_params = 'ID{}_literature_params.txt'
config_file_radvel = 'radvel_ID{}.py'
nonRV_data = 'nonRVs_ID{}.txt'
awk_logfiles = 'grep_redID.awk'
templates_fxcor = 'templates_ID{}.lis'
fxcor_script = 'fxcor_with_lists.cl'
barycorr_out = 'RVs_ID{}_abc_single.txt'

# initialize all paths and make them platform independent
location = Path(__file__).parent
abs_path_scripts = (location / path_scripts).resolve()
abs_path_output = (location / path_output).resolve()
abs_path_obslog = (location / path_obslog_local).resolve()
abs_path_data = (location / path_data_local).resolve()
abs_path_red_gamse = (location / path_reduce_gamse).resolve()
abs_path_IRAF = (location / path_IRAF_datainput).resolve()
abs_path_rvout = (location / path_rv_results).resolve()
abs_path_rvplots = (location / path_rv_plots).resolve()

# other paths
gamse_reduce_folder = os.path.join(abs_path_red_gamse, path_gamse_reduce)
gamse_results_folder = os.path.join(gamse_reduce_folder, path_gamse_results)
iraf_data_folder = os.path.join(abs_path_IRAF, 'ID{}')
iraf_output_folder = os.path.join(abs_path_output, 'ID{}')
iraf_smd_folder = os.path.join(abs_path_output, path_smd)

# make absolute paths to the files
script_USM = os.path.join(abs_path_scripts, file_script_USM)
script2_USM = os.path.join(abs_path_scripts, data_script_USM)
script_local = os.path.join(abs_path_scripts, file_script_local)
script2_local = os.path.join(abs_path_scripts, data_script_local)
script_add = os.path.join(abs_path_scripts, add_header_script)
grep_redID_cmd = os.path.join(abs_path_scripts, grep_redmineID_script)
grep_redID_out = os.path.join(iraf_output_folder, grep_redmineID_results)
sort_copy_cmd = os.path.join(abs_path_scripts, sort_copy_gamse_script)
out_obsnights = os.path.join(iraf_output_folder, observed_nights_list)
out_gamse_copy = os.path.join(abs_path_output, copy_dates_gamse)
copy_reduced_cmd = os.path.join(abs_path_scripts, copy_wvcal_script)
make_orderlists = os.path.join(abs_path_scripts, recipe_orderlists)
make_cl_fxcor = os.path.join(iraf_output_folder, recipe_cl_fxcor)
frames_list = os.path.join(iraf_output_folder, all_used_frames)
out_RVs_single = os.path.join(abs_path_rvout, RVs_single_orders)
out_RVs_weighted = os.path.join(abs_path_rvout, RVs_weighted)
out_tels_weighted = os.path.join(abs_path_rvout, tels_weighted)
out_RVs_telcorr = os.path.join(abs_path_rvout, RVs_telcorr)
input_tel_orders = os.path.join(abs_path_IRAF, ords_tellurics)
out_RVs_compare = os.path.join(abs_path_rvout, RVs_compare)
lit_planet_params = os.path.join(abs_path_rvout, literature_params)
radvel_config = os.path.join(abs_path_rvout, config_file_radvel)
out_nonRV_data = os.path.join(abs_path_rvout, nonRV_data)
awk_script = os.path.join(abs_path_scripts, awk_logfiles)
template_list = os.path.join(iraf_output_folder, templates_fxcor)
out_RVs_abc_single = os.path.join(abs_path_rvout, barycorr_out)

# define the remote hosts and paths for rsync
address_focespc = 'foces@195.37.68.140'
address_ohiaaipc = 'wstobserver@195.37.68.19'
address_aliothpc = 'fahrenschon@alioth.usm.uni-muenchen.de'
# fits_path_focespc = '/data/FOCES/{}'
fcslinks_path_focespc = '/data/fcs_links/{}'
log_path_ohiaaipc = '/data/3kk/{}'
onedspec_path_aliothpc = '/e/ldata/users/gamse/fcs/{}/onedspec'
