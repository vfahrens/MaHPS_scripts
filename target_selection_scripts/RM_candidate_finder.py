import os
from pathlib import Path
from astroquery.simbad import Simbad
from operator import itemgetter

catalogs_dir = '../'
work_dir = Path(__file__).parent
abs_path_cats = (work_dir / catalogs_dir).resolve()

# define filename of the TEP catalog
tepcat = 'TEPCat_obliquity_catalog.txt'
# define filename of the TEP observables catalog
tepcat_obs = 'TEPCat_observables_catalog.txt'
# define filename of the TEP catalog of stellar and planetary parameters
tepcat_sppar = 'TEPCat_starplanparams_catalog.txt'
# define filenames of the Transit Scheduler output files
transsched_out = 'RM_objects_TEPCat.txt'
transsched_other_out = 'other_objects_for_observations.txt'

tepcat_path = os.path.join(abs_path_cats, tepcat)
tepcat_obs_path = os.path.join(abs_path_cats, tepcat_obs)
tepcat_sppar_path = os.path.join(abs_path_cats, tepcat_sppar)
transsched_out_path = os.path.join(abs_path_cats, transsched_out)
transsched_other_out_path = os.path.join(abs_path_cats, transsched_other_out)

# define the string for the header of the Transit Scheduler output file
trans_sched_header = 'No\tID\tRa\tDec\tP\tRp\tDepth\tT0 - 2457000\tDuration\tVmag\tTeff\tM\tVrot\tObs\tMoon\t' \
                     'Other_Teams\tNotes\n'


# function for reding all object names from the TEP catalog
def get_tepcat_objects():
    # initialize list for all object names in the catalog
    obj_names_tep_unsorted = []

    # open the catalog and add the first entry of each line to the object name list (excluding commented lines)
    with open(tepcat_path, 'r') as in_tep:
        for line in in_tep:
            line = line.strip()
            line = line.split()
            if line[0] != '#':
                # # use this format corrections for SIMBAD object query
                # if line[0][:5] == 'Qatar':
                #     line[0] = line[0][:5] + ' ' + line[0][-1:]
                # elif line[0][:2] == 'HD' and line[0][-1] == 'b':
                #     line[0] = line[0][:-1]
                # elif line[0][:2] == 'HD' and line[0][-1] == 'c':
                #     line[0] = line[0][:-1]
                # elif line[0][:2] == 'K2' and line[0][-1] == 'b':
                #     line[0] = line[0][:-1]
                # elif line[0][:6] == '55_Cnc' and line[0][-1] == 'e':
                #     line[0] = line[0][:-2]
                # elif line[0][:6] == 'Kepler':
                #     if line[0][-1] == 'b' or line[0][-1] == 'c' or line[0][-1] == 'd':
                #         line[0] = line[0][:-1]
                # elif line[0][:8] == 'TRAPPIST':
                #     line[0] = line[0][:-1]
                # elif line[0][:4] == 'WASP' and line[0][-1] == 'b':
                #     line[0] = line[0][:-1]

                # fill the empty list with tuples of ID and T_eff values (ignore T_eff errors)
                obj_names_tep_unsorted.append((line[0], line[1]))

        # remove all duplicates of names from the list and sort it alphabetically
        obj_names_tep = sorted(set(obj_names_tep_unsorted), key=itemgetter(0))

        # # print each object name for comparison with other catalogs, e.g. photo of Jana
        # for obob in range(len(obj_names_tep)):
        #     print(obj_names_tep[obob])

    return obj_names_tep


# function for getting the observation parameters from the TEP catalog of the RM objects
def get_tepcat_observables(obj_list):
    # initialize list for all object observation parameters
    obj_observ_tep_unsorted = []
    # initialze list for the observation parameters of the RM effect objects from the obliquity catalog
    obsparams_rm_objects = []
    # initialize list for preselected targets
    preselected_obsparams_rm = []

    # open the catalog and add each line to the observation parameter list
    with open(tepcat_obs_path, 'r') as in_tep_obs:
        for line in in_tep_obs:
            line = line.strip()
            line = line.split()
            obj_observ_tep_unsorted.append(line)

    # add only the observation parameters of the obliquity catalog objects to this list and also append the T_eff value
    for rm_object in obj_list:
        for catalog_object in range(len(obj_observ_tep_unsorted)):
            if obj_observ_tep_unsorted[catalog_object][0] == rm_object[0]:
                obsparams_rm_objects.append(obj_observ_tep_unsorted[catalog_object] + [rm_object[1]])

    # do preselection that allows only objects with declination higher that -30 degrees, magnitude brighter than 9.5
    # and T_eff less than 7500 (F-type star or later type)
    for that_obj in range(len(obsparams_rm_objects)):
        if int(obsparams_rm_objects[that_obj][5]) > -30 and float(obsparams_rm_objects[that_obj][8]) < 9.5:  # \
            # and int(obsparams_rm_objects[that_obj][-1]) < 7500:
            preselected_obsparams_rm.append(obsparams_rm_objects[that_obj])

    return preselected_obsparams_rm


# function for getting the stellar mass and planetary radius from the TEP catalog
def get_tepcat_starplanparams(obj_list_preselect):
    # initialize list for stellar and planetary radius data
    starplan_unsorted = []
    starplan_rm_objects = []

    # read the catalog with the stellar and planetary parameters
    with open(tepcat_sppar_path, 'r') as in_sppar:
        for lline in in_sppar:
            lline = lline.strip()
            lline = lline.split()
            starplan_unsorted.append(lline)

    # add only the already preselected objects to a list
    for one_rm_obj in obj_list_preselect:
        for star_obj in range(len(starplan_unsorted)):
            if starplan_unsorted[star_obj][0] == one_rm_obj[0]:
                starplan_rm_objects.append(starplan_unsorted[star_obj])

    return starplan_rm_objects


# function for adding V_rot and RV amplitude values that were manually extracted from papers
def add_velocity_info():
    # dictionary with : object ID, V_rot, observed amplitude of RM effect
    vel_info_dict = {
        '55_Cnc_e': [3.3, '0.60 m/s'],
        'HAT-P-02': [20.8, '90.0 m/s'],
        'HAT-P-11': [2.09, '5.0 m/s'],
        'HAT-P-70': [100.0, 'Doppler tomography'],
        'HD_003167c': [2.0, '2.0 m/s'],
        'HD_017156': [4.1, '20.0 m/s'],
        'HD_063433b': [7.3, '6.0 m/s'],
        'HD_080606': [1.7, '12.0 m/s'],
        'HD_106315c': [12.9, 'Doppler tomography'],
        'HD_149026': [6.0, '10.0 m/s'],
        'HD_189733': [3.5, '50.0 m/s'],
        'HD_209458': [4.5, '80.0 m/s'],
        'KELT-07': [65.0, 'Doppler tomography'],
        'KELT-09': [115.0, 'Doppler tomography'],
        'KELT-17': [441.5, '300.0 m/s'],
        'KELT-20': [114.0, 'Doppler tomography'],
        'KELT-24': [20.3, '25.0 m/s'],
        'Kepler-408': [2.8, 'no spectroscopy'],
        'MASCARA-1': [109.0, 'Doppler tomography'],
        'WASP-033': [86.6, 'Doppler tomography'],
        'WASP-038': [8.3, '100.0 m/s'],
        'WASP-166': [4.6, '24.0 m/s'],
        'WASP-189': [93.1, '200.0 m/s']
    }

    return vel_info_dict


# function for getting the coordinates and magnitude of all objects
def get_simbad_parameters(obj_list):
    # define which values should be retrieved from SIMBAD
    Simbad.add_votable_fields('flux(V)')
    Simbad.add_votable_fields('flux(R)')
    Simbad.add_votable_fields('typed_id')
    # send the query command to SIMBAD
    radec_table = Simbad.query_objects(obj_list)

    # initialize a list for the preselected targets
    preselected_objects = []

    # extract the name, ra/dec coordinates and available fluxes for each object
    for each_obj in range(len(radec_table)):
        # the typed ID is a byte string, convert that to a char string
        objname = radec_table['TYPED_ID'][each_obj].decode('UTF-8')
        # use colon notation for the coordinates
        objcoord_ra = radec_table['RA'][each_obj].replace(' ', ':')
        objcoord_dec = radec_table['DEC'][each_obj].replace(' ', ':')
        # use the V-band flux if available
        if radec_table['FLUX_V'][each_obj] != '--':
            objmag = float(radec_table['FLUX_V'][each_obj])
        # if no V-band flux is available, use the R-band flux value
        elif radec_table['FLUX_R'][each_obj] != '--':
            objmag = float(radec_table['FLUX_R'][each_obj])
        # use dummy value 100.0 if no measured magnitude is availablle (then it is anyway too faint for FOCES)
        else:
            objmag = 100.0

        # use magnituse limit 9.5 and declination limit -30 degree for preselection of targets observable with FOCES
        if objmag < 9.5 and int(objcoord_dec[:3]) > -30:
            objparams = [objname, objcoord_ra, objcoord_dec, objmag]
            preselected_objects.append(objparams)

    # print(preselected_objects)
    print(len(preselected_objects))

    for bobjob in range(len(preselected_objects)):
        print(preselected_objects[bobjob])

    return


# function to write observation parameters to txt file for the Transit Scheduler
def write_transit_scheduler_file(tep_coord_list, tep_starplanparam_list, vel_info_dict):
    # define the default moon distance to be 100°
    moon_dist = 100.0

    # open the catalog and add the first entry of each line to the object name list (excluding commented lines)
    with open(transsched_out_path, 'w') as out_trasched:
        # write the header to the Transit Scheduler file
        out_trasched.write(trans_sched_header)

        # add data for each object to the Transit Scheduler file
        for one_object in range(len(tep_coord_list)):
            # do some reformatting and unit conversions where necessary
            ra_coord = tep_coord_list[one_object][2] + ':' + tep_coord_list[one_object][3] + \
                       ':' + tep_coord_list[one_object][4]
            dec_coord = tep_coord_list[one_object][5] + ':' + tep_coord_list[one_object][6] + \
                        ':' + tep_coord_list[one_object][7]
            t0_sub = float(tep_coord_list[one_object][12]) - 2457000.0
            trans_duration = float(tep_coord_list[one_object][10]) * 24.0

            # get the stellar mass and planetary radius for the current object
            for obj_thing in range(len(tep_starplanparam_list)):
                if tep_starplanparam_list[obj_thing][0] == tep_coord_list[one_object][0]:
                    star_mass_msun = tep_starplanparam_list[obj_thing][7]
                    planet_rad_rjup = tep_starplanparam_list[obj_thing][29]
                    # convert the planet radius from R_Jup t R_Earth
                    planet_rad_rearth = float(planet_rad_rjup) * 11.2

            # get the velocity information for the current object from the manual input above
            v_rot = vel_info_dict[tep_coord_list[one_object][0]][0]
            rv_meas_note = vel_info_dict[tep_coord_list[one_object][0]][1]

            string_list = [str(one_object + 1), tep_coord_list[one_object][0], ra_coord, dec_coord,
                           tep_coord_list[one_object][14], str(planet_rad_rearth), tep_coord_list[one_object][11],
                           str(t0_sub), str(trans_duration), tep_coord_list[one_object][8],
                           tep_coord_list[one_object][-1], star_mass_msun, str(v_rot), '', str(moon_dist),
                           tep_coord_list[one_object][-2], rv_meas_note]
            joined_string = '\t'.join(string_list)
            out_trasched.write(joined_string + '\n')

    return


# function to write the ra /dec data of other important observation programs for the Transit Scheduler
def write_other_objects_file(other_objs_dict):
    obj_num = 1

    with open(transsched_other_out_path, 'w') as out_othsched:
        # write the header to the Transit Scheduler file
        out_othsched.write(trans_sched_header)

        # extract the name, ra/dec coordinates and available fluxes for each object
        for objname in set(other_objs_dict):
            # use colon notation for the coordinates
            objcoord_ra = other_objs_dict[objname][0]
            objcoord_dec = other_objs_dict[objname][1]
            if ' ' in objcoord_ra:
                objcoord_ra = objcoord_ra.replace(' ', ':')
            if ' ' in objcoord_dec:
                objcoord_dec = objcoord_dec.replace(' ', ':')

            string_list = [str(obj_num), objname, objcoord_ra, objcoord_dec, '', '', '', '', '', '', '', '', '', '',
                           '', '', '']
            joined_string = '\t'.join(string_list)
            out_othsched.write(joined_string + '\n')

            obj_num += 1

    return


# # use this part to generate a list of RM effect objects for the Transit Scheduler
# olist = get_tepcat_objects()
# # get_simbad_parameters(olist)
# sortlist = get_tepcat_observables(olist)
# starplan_list = get_tepcat_starplanparams(sortlist)
# velinfo_list = add_velocity_info()
# write_transit_scheduler_file(sortlist, starplan_list, velinfo_list)

# use this part to generate a Transit Scheduler file of the other important observation projects
# (e.g. those of Matthias)
names_of_other_objects = {'SDSSJ1433': ['14:33:22.8', '+60:07:13.44'],
                          'SDSSJ1721': ['17:21:49', '+88:42:22'],
                          'A2666': ['23 50 56.2', '+27 08 41'],
                          'A2513': ['22 59 25.0', '+26 15 06'],
                          'M33': ['01 33 50.904', '+30 39 35.79'],
                          'DDO216': ['23:28:34', '+14:44:48'],
                          'M42': ['05:35:16', '-05:23:17'],
                          'M101': ['14 03 12.58', '+54 30 55.5']}
write_other_objects_file(names_of_other_objects)
