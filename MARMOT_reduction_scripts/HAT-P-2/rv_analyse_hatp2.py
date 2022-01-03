#!/usr/bin/env python3
# -*- coding: iso-8859-1 -*-

import os, sys
from marmot import RVAnalysis
import pickle
import matplotlib.cm
import matplotlib.pyplot as plt

cmap = matplotlib.cm.get_cmap('Spectral')
matplotlib.use('TkAgg')


def main(inifile=None):

    if inifile is None:
        inifile = sys.argv[1]

    ana = RVAnalysis(ini_filepath=inifile, vlevel=3)
    #tmp_order_list = ana.order_list
    #tmp_multicore = ana.multicore
    #ana = pickle.load(open(os.path.join(ana.datainf['opath'], ana.datainf['analyproj']), "rb"))
    #ana.multicore = tmp_multicore
    #ana.order_list = tmp_order_list
    #print(ana.multicore, ana.order_list)

    ######################################################
    ## Start Analysis
    ######################################################

    if True:
        ### Reading in the data ###
        ana.load_files(multicore=ana.multicore)#mode='comb')

        plt.close('all')

        if True:
            # Identify single/double frames
            ana.single_double_list()
            all_frames = ana.double_frames.copy()
            all_frames = all_frames + ana.single_frames
            # Set the specified mask regions
            ana.set_masks(path=ana.datainf['mpath'], filename=ana.datainf['maskfile'], fiber=ana.scifiber, filelist=all_frames)
            select_use_frames=ana.select_template_frames(plot=True,mode='single',quantity=7)

            use_frames=[]
            for f in ana.double_frames:
                for day in ana.days:
                    if ana.frames[f].filename.find(day)>-1:
                        use_frames.append(f)

            ana.double_frames=use_frames

            pickle.dump(ana, open(os.path.join(ana.datainf['opath'],ana.datainf['analyproj']), "wb"))

            # Create template

            #ana.create_templates(ana.scifiber, orders=ana.order_list, plot=False, use_frames=select_use_frames[0],opath=ana.datainf['opath'], ofilename=ana.datainf['star_templ'], multicore=ana.multicore)
            #ana.save_pickle_template(opath=ana.datainf['opath'], ofilename=ana.datainf['star_templ'])
            ana.load_pickle_template(ipath=ana.datainf['opath'], ifilename=ana.datainf['star_templ'])


            # -- Option 3a or 3b --#

            # 3a New wavelength calibration but single line fitted
            # 3b New wavelength calibration but Comb fitted as Comb

            # Applying the comb calibration frames
            #ana.apply_comb_wv_cal(orders=ana.order_list, pixrange=[200,1850], use_frames=ana.double_frames, which='cal', CombAsComb=False, plotwv=False, multicore=1)

            pickle.dump(ana, open(os.path.join(ana.datainf['opath'], ana.datainf['analyproj']), "wb"))
            ana.rv_ccf(fiber='A', order_list=ana.order_list, fitrange=[200, 1850], use_frames=all_frames, opath=ana.datainf['opath'], show=False, plot=False, remove_tellurics=False, multicore=ana.multicore)
            pickle.dump(ana, open(os.path.join(ana.datainf['opath'], ana.datainf['analyproj']), "wb"))

            ### single order RVs saved to file ###
            with open(os.path.join(ana.datainf['opath'], 'RVs_ID3204_singleord_marmot_ccf.txt'), 'w') as single_file:
                for gg in all_frames:
                    fileid_s = int(ana.frames[gg].filename[:8]+ana.frames[gg].filename[9:13])
                    date_jd_s = ana.frames[gg].midexpt.jd
                    for ordi in ana.order_list:
                        rv_val_ord_s = ana.frames[gg].data['A'][str(ordi)].rv_val['ccf']
                        rv_err_ord_s = ana.frames[gg].data['A'][str(ordi)].rv_err['ccf']

                        singl_string = str(fileid_s) + ' ' + str(date_jd_s) + ' ' + str(rv_val_ord_s) + ' ' + str(rv_err_ord_s) + ' ' + str(ordi) + '\n'
                        single_file.write(singl_string)

            ### single order RVs saved to file ###
            with open(os.path.join(ana.datainf['opath'], 'RVs_ID3204_singleord_marmot_ccffit.txt'), 'w') as single_file:
                for gg in all_frames:
                    fileid_s = int(ana.frames[gg].filename[:8]+ana.frames[gg].filename[9:13])
                    date_jd_s = ana.frames[gg].midexpt.jd
                    for ordi in ana.order_list:
                        rv_val_ord_s = ana.frames[gg].data['A'][str(ordi)].rv_val['ccf_fit']
                        rv_err_ord_s = ana.frames[gg].data['A'][str(ordi)].rv_err['ccf_fit']

                        singl_string = str(fileid_s) + ' ' + str(date_jd_s) + ' ' + str(rv_val_ord_s) + ' ' + str(rv_err_ord_s) + ' ' + str(ordi) + '\n'
                        single_file.write(singl_string)

            ### RV average calculation and saving to file ###
            for ff in all_frames:
                ana.frames[ff].average_rv(fiber='A', which='ccf', orders=ana.order_list)

            pickle.dump(ana, open(os.path.join(ana.datainf['opath'], ana.datainf['analyproj']), "wb"))

            with open(os.path.join(ana.datainf['opath'], 'RVs_ID3204_weighted_marmot_ccf.txt'), 'w') as rv_file:
                for fif in all_frames:
                    fileid = int(ana.frames[fif].filename[:8]+ana.frames[fif].filename[9:13])
                    date_jd = ana.frames[fif].midexpt.jd
                    rv_val_marmot = ana.frames[fif].rv_val['ccf']
                    rv_err_marmot = ana.frames[fif].rv_err['ccf']

                    outstring = str(fileid) + ' ' + str(date_jd) + ' ' + str(rv_val_marmot) + ' ' + str(rv_err_marmot) + '\n'
                    rv_file.write(outstring)

            ### RV average calculation and saving to file ###
            for ff in all_frames:
                ana.frames[ff].average_rv(fiber='A', which='ccf_fit', orders=ana.order_list)

            pickle.dump(ana, open(os.path.join(ana.datainf['opath'], ana.datainf['analyproj']), "wb"))

            with open(os.path.join(ana.datainf['opath'], 'RVs_ID3204_weighted_marmot.txt'), 'w') as rv_file:
                for fif in all_frames:
                    fileid = int(ana.frames[fif].filename[:8]+ana.frames[fif].filename[9:13])
                    date_jd = ana.frames[fif].midexpt.jd
                    rv_val_marmot = ana.frames[fif].rv_val['ccf_fit']
                    rv_err_marmot = ana.frames[fif].rv_err['ccf_fit']

                    outstring = str(fileid) + ' ' + str(date_jd) + ' ' + str(rv_val_marmot) + ' ' + str(rv_err_marmot) + '\n'
                    rv_file.write(outstring)



    return ana

if __name__ == '__main__':
    ana=main()

