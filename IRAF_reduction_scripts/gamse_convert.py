#!/usr/bin/env python3

from pathlib import Path
import os
import numpy as np
import astropy.io.fits as fits
from datetime import datetime
import datetime as dt
from astroquery.simbad import Simbad
# import other python scripts
import paths_and_files as pf


# extract the single orders from the GAMSE result files, add header entries with wavelength calibration information
def iraf_converter(infolder, redmine_id, calib_fiber=False, raw_flux=False):
    print('Started IRAF conversion for redmine project {}...'.format(redmine_id))

    # give the readout time of the camera for the calculation of the mid-time of observation
    ccd_readtime = 87.5

    # define the folders for reading the data and saving the new files
    path_in = infolder
    path_out = pf.iraf_output_folder.format(redmine_id)

    # check if the input directory exists, give a warning if not
    if not os.path.exists(path_in):
        raise Exception('Input directory {} not found. Please check the path name '
                        'or previous reduction steps.'.format(path_in))
    # check if the output directory exists, create if it does not
    if not os.path.exists(path_out):
        os.makedirs(path_out)

    # make a list of all files in the input folder
    fname_lst = sorted(os.listdir(path_in))
    for fname in fname_lst:

        # if the file does not end with '.fits', skip it and continue with the next file in the list
        if fname[-5:] != '.fits':
            continue
        else:
            file_in = os.path.join(path_in, fname)

            # open the GAMSE fits file and read the general header and the data section
            head = fits.getheader(file_in, 0)
            data = fits.getdata(file_in)

            # calculate the midtime of observation
            # ATTENTION: the timestamp in the header is the time when the file is saved,
            # after exposure and readout time
            exp_endtime = datetime.strptime(head['FRAME'], '%Y-%m-%dT%H:%M:%S.%f')
            exp_midtime = exp_endtime - 0.5 * dt.timedelta(seconds=head['EXPOSURE']) - \
                          dt.timedelta(seconds=ccd_readtime)

            # add ra/dec coordinates to header if not present
            if 'PERA2000' not in head or 'PEDE2000' not in head:
                objname = input('Please give the name of the object for a coordinate search with SIMBAD: ')
                radec_table = Simbad.query_object(objname)
                objcoord_ra = radec_table['RA'][0].replace(' ', ':')
                objcoord_dec = radec_table['DEC'][0].replace(' ', ':')
                head['PERA2000'] = (objcoord_ra, 'Right ascension coordinate')
                head['PEDE2000'] = (objcoord_dec, 'Declination coordinate')
                print('Ra/Dec coordinates added to header.')

            head['EQUINOX'] = (2000.0, 'Epoch of observation')
            head['UTMID'] = (datetime.isoformat(exp_midtime), 'UT of midpoint of observation')

            # save the original GAMSE header as primary header of the MEF file (primary data section is empty)
            empty_primary = fits.PrimaryHDU(header=head)
            hdu_list_fred = fits.HDUList([empty_primary])
            if raw_flux:
                hdu_list_fraw = fits.HDUList([empty_primary])

            # check if the data are single or multi fiber: it is required that the reduction was made with
            # double-fiber configuration, even if the observation was without simultaneous calibration; this is
            # only the case, if the keyword 'fiber' is present in the data frame
            if 'fiber' in data.dtype.names:
                # treat the results of fiber A and B separately and independently
                for fiber in np.unique(data['fiber']):
                    # if the data from the calibration fiber is not needed, skip it and continue with the next dataset
                    if calib_fiber is False and fiber == 'B':
                        continue

                    else:
                        # make a new (empty) header for this extension, add the physical order number
                        ext_head = fits.Header()

                        # add all the header entries required for the CCF calculation
                        ext_head['PERA2000'] = head['PERA2000']
                        ext_head['PEDE2000'] = head['PEDE2000']
                        ext_head['FRAME'] = head['FRAME']
                        ext_head['UTMID'] = head['UTMID']
                        ext_head['EXPOSURE'] = head['EXPOSURE']
                        ext_head['EQUINOX'] = head['EQUINOX']

                        # add the required standard header entries for IRAF wavelength calibration
                        ext_head['WCSDIM'] = 2
                        ext_head['CTYPE1'] = 'MULTISPE'
                        ext_head['CTYPE2'] = 'MULTISPE'
                        ext_head['CDELT1'] = 1.
                        ext_head['CDELT2'] = 1.
                        ext_head['CD1_1'] = 1.
                        ext_head['CD2_2'] = 1.
                        ext_head['LTM1_1'] = 1.
                        ext_head['LTM2_2'] = 1.
                        ext_head['WAXMAP01'] = '1 0 0 0 '
                        ext_head['WAT0_001'] = 'system=multispec'
                        # only this way of writing works with fxcor!
                        ext_head['WAT1_001'] = 'wtype=multispec label=Wavelength units=angstroms'
                        # this way of writing is also ok for splot, but NOT for fxcor
                        # ext_head['WAT1_001'] = 'wtype=multispec label=Wavelength units=Wavelength'

                        # make a mask that contains all apertures which contain light of either fiber A or B
                        mask = data['fiber'] == fiber

                        # read the data row by row (aperture by aperture)
                        for row in data[mask]:
                            # each order has to be stored in a separate extension, so in each extension the
                            # data is in aperture 1
                            aperture = 1
                            order = row['order']
                            wave = row['wavelength']  # all wavelength values for the current order (aperture)
                            flux_reduced = row['flux_sum']  # all flux values as produced by GAMSE
                            if raw_flux:
                                flux_raw = row['flux_raw']  # raw flux values without flat/background/etc subtraction

                            # add the physical order number to the header
                            ext_head['PHYSORD'] = (order, 'Physical order of aperture')

                            # definition of other parameters that IRAF needs for correct interpretation
                            # of the wavelength calibration
                            dtype = 2  # non-linear dispersion function
                            wave1 = wave[0]  # wavelength coordinate of the first pixel
                            delta_wave = (wave[-1] - wave[0]) / len(wave)  # average dispersion interval per pixel
                            num_pix = len(wave)  # number of valid pixels
                            z_corr = 0.  # Doppler correction factor
                            aplow = 0.0  # dummy value for the lower aperture extraction limit
                            aphigh = 0.0  # dummy value for the upper aperture extraction limit
                            weight_i = 1.  # weight of this dispersion function
                            zero_off_i = 0.  # zero point offset of this dispersion function
                            ftype_i = 5  # pixel coordinate array is used as dispersion information

                            # put the whole string together
                            separ = ' '
                            wave_str = separ.join(str(wl) for wl in wave)  # converts the wave array to a string
                            wlcalib_paramlst = [aperture, order, dtype, wave1, delta_wave, num_pix, z_corr, aplow,
                                                aphigh, weight_i, zero_off_i, ftype_i, num_pix, wave_str]
                            wlcalib_str = separ.join(str(item) for item in wlcalib_paramlst)

                            longstring = 'wtype=multispec spec{} = "{}" '.format(aperture, wlcalib_str)

                            # define the keyword for each header entry and fill it with the corresponding
                            # part of the long string
                            head_key = 'WAT2_{:03d}'
                            i = 0
                            new_str_start = 0
                            for x in range(len(longstring)):
                                i += 1
                                head_key_num = head_key.format(i)
                                # calculate the maximum length of the string in the header
                                longstr_len = 81 - len(head_key_num) - 5

                                # find the right part of the string to write into the header
                                if new_str_start <= len(longstring) \
                                        and len(longstring) - new_str_start > 81 - 5 - len(head_key.format(i + 1)):
                                    string_part = longstring[new_str_start:new_str_start + longstr_len]
                                    new_str_start = new_str_start + longstr_len
                                    ext_head[head_key_num] = string_part

                                # at the end of the string, use the rest and leave the for loop
                                else:
                                    string_part = longstring[new_str_start:]
                                    ext_head[head_key_num] = string_part
                                    break

                            # convert the flux values to arrays for saving in a FITS file
                            flux_np = np.array(flux_reduced)
                            if raw_flux:
                                flux_raw_np = np.array(flux_raw)

                            # create multi-extension FITS files, one for the reduced and one for the raw flux
                            image_hdu_fred = fits.ImageHDU(data=flux_np, header=ext_head)
                            hdu_list_fred.append(image_hdu_fred)

                            if raw_flux:
                                image_hdu_fraw = fits.ImageHDU(data=flux_raw_np, header=ext_head)
                                hdu_list_fraw.append(image_hdu_fraw)

                    # save the new IRAF compatible multi-extension FITS files with the reduced and raw flux
                    # fname_fred = '{}_ods_fred.fits'.format(fname[:13])
                    fname_fred = '{}_fred.fits'.format(fname[:-5])
                    out_fred = os.path.join(path_out, fname_fred)
                    hdu_list_fred.writeto(out_fred, overwrite=True)

                    if raw_flux:
                        fname_fraw = '{}_ods_fraw.fits'.format(fname[:13])
                        out_fraw = os.path.join(path_out, fname_fraw)
                        hdu_list_fraw.writeto(out_fraw, overwrite=True)

                    print('Converted double fiber file: {}.'.format(fname))

            else:
                # make a new (empty) header for this extension, add the physical order number
                ext_head = fits.Header()

                # add all the header entries required for the CCF calculation
                ext_head['PERA2000'] = head['PERA2000']
                ext_head['PEDE2000'] = head['PEDE2000']
                ext_head['FRAME'] = head['FRAME']
                ext_head['UTMID'] = head['UTMID']
                ext_head['EXPOSURE'] = head['EXPOSURE']
                ext_head['EQUINOX'] = head['EQUINOX']

                # add the required standard header entries for IRAF wavelength calibration
                ext_head['WCSDIM'] = 2
                ext_head['CTYPE1'] = 'MULTISPE'
                ext_head['CTYPE2'] = 'MULTISPE'
                ext_head['CDELT1'] = 1.
                ext_head['CDELT2'] = 1.
                ext_head['CD1_1'] = 1.
                ext_head['CD2_2'] = 1.
                ext_head['LTM1_1'] = 1.
                ext_head['LTM2_2'] = 1.
                ext_head['WAXMAP01'] = '1 0 0 0 '
                ext_head['WAT0_001'] = 'system=multispec'
                # only this way of writing works with fxcor!
                ext_head['WAT1_001'] = 'wtype=multispec label=Wavelength units=angstroms'
                # this way of writing is also ok for splot, but NOT for fxcor
                # head['WAT1_001'] = 'wtype=multispec label=Wavelength units=Wavelength'

                # read the data row by row (aperture by aperture)
                for row in data:
                    # each order has to be stored in a separate extension, so in each extension the
                    # data is in aperture 1
                    aperture = 1
                    order = row['order']
                    wave = row['wavelength']  # all wavelength values for the current order (aperture)
                    flux_reduced = row['flux_sum']  # all flux values as produced by GAMSE
                    if raw_flux:
                        flux_raw = row['flux_raw']  # raw flux values without flat/background/etc subtraction

                    # add the physical order number to the header
                    ext_head['PHYSORD'] = (order, 'Physical order of aperture')

                    # definition of other parameters that IRAF needs for correct interpretation
                    # of the wavelength calibration
                    dtype = 2  # non-linear dispersion function
                    wave1 = wave[0]  # wavelength coordinate of the first pixel
                    delta_wave = (wave[-1] - wave[0]) / len(wave)  # average dispersion interval per pixel
                    num_pix = len(wave)  # number of valid pixels
                    z_corr = 0.  # Doppler correction factor
                    aplow = 0.0  # dummy value for the lower aperture extraction limit
                    aphigh = 0.0  # dummy value for the upper aperture extraction limit
                    weight_i = 1.  # weight of this dispersion function
                    zero_off_i = 0.  # zero point offset of this dispersion function
                    ftype_i = 5  # pixel coordinate array is used as dispersion information

                    # put the whole string together
                    separ = ' '
                    wave_str = separ.join(str(wl) for wl in wave)  # converts the wave array to a string
                    wlcalib_paramlst = [aperture, order, dtype, wave1, delta_wave, num_pix, z_corr, aplow,
                                        aphigh, weight_i, zero_off_i, ftype_i, num_pix, wave_str]
                    wlcalib_str = separ.join(str(item) for item in wlcalib_paramlst)

                    longstring = 'wtype=multispec spec{} = "{}" '.format(aperture, wlcalib_str)

                    # define the keyword for each header entry and fill it with the corresponding
                    # part of the long string
                    head_key = 'WAT2_{:03d}'
                    i = 0
                    new_str_start = 0
                    for x in range(len(longstring)):
                        i += 1
                        head_key_num = head_key.format(i)
                        # calculate the maximum length of the string in the header
                        longstr_len = 81 - len(head_key_num) - 5

                        # find the right part of the string to write into the header
                        if new_str_start <= len(longstring) \
                                and len(longstring) - new_str_start > 81 - 5 - len(head_key.format(i + 1)):
                            string_part = longstring[new_str_start:new_str_start + longstr_len]
                            new_str_start = new_str_start + longstr_len
                            ext_head[head_key_num] = string_part

                        # at the end of the string, use the rest and leave the for loop
                        else:
                            string_part = longstring[new_str_start:]
                            ext_head[head_key_num] = string_part
                            break

                    # convert the flux values to arrays for saving in a FITS file
                    flux_np = np.array(flux_reduced)
                    if raw_flux:
                        flux_raw_np = np.array(flux_raw)

                    # create multi-extension FITS files, one for the reduced and one for the raw flux
                    image_hdu_fred = fits.ImageHDU(data=flux_np, header=ext_head)
                    hdu_list_fred.append(image_hdu_fred)

                    if raw_flux:
                        image_hdu_fraw = fits.ImageHDU(data=flux_raw_np, header=ext_head)
                        hdu_list_fraw.append(image_hdu_fraw)

                # save the new IRAF compatible multi-extension FITS files with the reduced and raw flux
                # fname_fred = '{}_ods_fred.fits'.format(fname[:13])
                fname_fred = '{}_fred.fits'.format(fname[:-5])
                out_fred = os.path.join(path_out, fname_fred)
                hdu_list_fred.writeto(out_fred, overwrite=True)

                if raw_flux:
                    fname_fraw = '{}_ods_fraw.fits'.format(fname[:13])
                    out_fraw = os.path.join(path_out, fname_fraw)
                    hdu_list_fraw.writeto(out_fraw, overwrite=True)

                print('Converted single fiber file: {}.'.format(fname))

    print('IRAF conversion completed.')
    return
