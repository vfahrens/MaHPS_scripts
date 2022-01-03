import os, sys
from marmot import RVAnalysis
import pickle
import matplotlib.cm
import matplotlib.pyplot as plt

cmap = matplotlib.cm.get_cmap('Spectral')
matplotlib.use('TkAgg')
import numpy as np

ana = RVAnalysis(ini_filepath='rv_analyse_hd189733.ini', vlevel=3)
ana = pickle.load(open(os.path.join(ana.datainf['opath'], ana.datainf['analyproj']), "rb"))
ana.load_pickle_template(ipath=ana.datainf['opath'], ifilename=ana.datainf['star_templ'])


fiber='A'
order_list=ana.order_list.copy()
#order_list.remove('76')
#order_list.remove('84')
#order_list.remove('86')
#order_list.remove('89')
#order_list.remove('90')
#order_list.remove('93')
#order_list.remove('94')
#order_list.remove('95')
#order_list.remove('98')
#order_list.remove('99')
#order_list.remove('100')
#order_list.remove('101')
#order_list.remove('104')
#order_list.remove('105')
#order_list.remove('106')
#order_list.remove('108')
#order_list.remove('110')
#order_list.remove('111')
#order_list.remove('112')
#order_list.remove('113')
#order_list.remove('114')
#order_list.remove('115')
#order_list.remove('118')
#order_list.remove('127')
#order_list.remove('132')
#order_list.remove('133')
show=True
saveas=None
ax=None
fold=False
dpi=200
use_frames=ana.double_frames.copy()
use_frames=use_frames + ana.single_frames
use_frames.remove(12)
T0=None
peri=None
Ampli=100
which='ccf'
which_err='ccf'
ref_frame=None
norm=False
weighted_avg=True
legend=True
mode='avg'
# mode options:
# 'avg': weighted avg
# 'mean': plain mean
# 'ref': using reference file

    
cmap = matplotlib.cm.get_cmap('jet')

def pick_marker(order):
    markers = ['+', 'x', 'o', 's', 'D', 'P', 'v', '>', '^']
    ord = int(order)
    return markers[(ord - 65) // 9]

def pick_color(order):
    colors = [cmap(i / 9) for i in range(9)]
    ord = int(order)
    return colors[(ord - 65) % 9]

print("Producing RV plot ...".format(fiber))


# calculate rvs for all frames
if not which_err:
    which_err = which

rv_ref = {}
for o in order_list:
    if mode == 'avg':
        vals = np.array([ana.frames[i].data[fiber][o].rv_val[which] for i in use_frames])
        errs = np.array([ana.frames[i].data[fiber][o].rv_err[which_err] for i in use_frames])
        rv_ref[o] = np.average(vals, weights=1/errs**2)
    if mode == 'mean':
        rv_ref[o] = np.mean([ana.frames[i].data[fiber][o].rv_val[which] for i in use_frames])
    if mode == 'ref':
        if ref_frame:
            rv_ref[o] = ana.frames[ref_frame].data[fiber][o].rv_val[which]
        else:
            self.print(FAIL+"No reference file provided!")

    #for i in use_frames:
    #    ana.frames[i].data[fiber][o].rv_rel[which] = ana.frames[i].data[fiber][o].rv_val[which] - rv_ref[o]
    #    ana.frames[i].data[fiber][o].rv_rel_err[which_err] = ana.frames[i].data[fiber][o].rv_err[which_err]
    # ana.frames[i].data[fiber][o].rv_ref[which] = rv_ref[o]

# Calculating weights for orders to be used in all frames uniformly
ord_weights = {}
for o in order_list:
    ord_weights[o] = np.sum([1/ana.frames[i].data[fiber][o].rv_err[which_err]**2 for i in use_frames])
for ij in use_frames:
    # calculate some average RVs
    print("Calculating average RV for frame at {}.".format(ana.frames[ij].midexpt))
    rv_vals = []  # rv of orders
    rv_rels = []  # rv of orders
    rv_val_errs = []  # rv error of orders
    rv_rel_errs = []  # rv error of orders
    weights = []  # weights of orders
    sn = []  # snr of orders

    for o in order_list:
        rv_vals.append(ana.frames[ij].data[fiber][o].rv_val[which])
        rv_val_errs.append(ana.frames[ij].data[fiber][o].rv_err[which_err])
        try:
            rv_rels.append(ana.frames[ij].data[fiber][o].rv_rel[which])
            rv_rel_errs.append(ana.frames[ij].data[fiber][o].rv_rel_err[which_err])
        except:
            rv_rels.append(ana.frames[ij].data[fiber][o].rv_val[which] - rv_ref[o])
            rv_rel_errs.append(ana.frames[ij].data[fiber][o].rv_err[which_err])
            #rv_rels.append(np.nan)
            print(f'No rv_rel for order {o}')
            pass
        if ana.frames[ij].data[fiber][o].rv_err[which_err] < 0:
            print(f'Found invalid error in order {o} -> mean instead of weighted average!')
            weighted_avg = False
            err = 1
        else:
            err = ana.frames[ij].data[fiber][o].rv_err[which_err]
            # if err < 0.5: err = 20.0
        weights.append(1 / err ** 2)
        sn.append(ana.frames[ij].data['A'][o].snr)

    #med = np.median(np.where(~np.isnan(rv_vals)))
    med = np.median(np.where(~np.isnan(rv_rels)))
    # print(med)
    for i,o in enumerate(order_list):
        rv_rel = rv_rels[i]
        rv_err = rv_rel_errs[i]
        #print(rv_rel-med, rv_err)
        if np.abs(rv_rel-med) > 50000*rv_err: # 5*rv_err:
            rv_vals[i]     = np.nan
            rv_rel_errs[i] = np.nan
            weights[i]     = np.nan

    if ord_weights:
        weights_2 = np.array([ord_weights[o] for o in order_list])
        weights_2 = np.ma.masked_array(weights_2, np.isnan(weights_2))

    weights = np.ma.masked_array(weights, np.isnan(weights))#np.ones_like(weights)#
    rv_vals = np.ma.masked_array(rv_vals, np.isnan(rv_vals))
    try:
        rv_rels = np.ma.masked_array(rv_rels, np.isnan(rv_rels))
    except:
        print(f'Cant do rv_rel!')
        pass

    # weighted average and error of the weighted avg
    if weighted_avg:
        ana.frames[ij].rv_val[which] = np.average(rv_vals, weights=weights)
        ana.frames[ij].rv_err[which_err] = np.sqrt(1 / np.sum(weights))
        try:
            if ord_weights:
                ana.frames[ij].rv_rel[which] = np.average(rv_rels, weights=weights_2)
            else:
                ana.frames[ij].rv_rel[which] = np.average(rv_rels, weights=weights)#np.average(np.array(rv_rels), weights=np.array(weights))
            ana.frames[ij].rv_rel_err[which_err] = np.sqrt(1 / np.sum(weights))
        except:
            pass
    # or if the measurements have no meaningful errors, just mean and std
    else:
        ana.frames[ij].rv_val[which] = np.average(rv_vals)
        ana.frames[ij].rv_err[which_err] = np.std(rv_vals) / np.sqrt(len(rv_vals))
        try:
            ana.frames[ij].rv_rel[which] = np.average(rv_rels)
            ana.frames[ij].rv_rel_err[which_err] = np.std(rv_rels) / np.sqrt(len(rv_rels))
        except:
            pass

    ana.frames[ij].snr = np.mean(sn)
    ana.frames[ij].std = np.std(np.array(rv_vals))
    

print("Calculating RVs for all frames".format(fiber))
for i in use_frames:
    print("{} +/- {}".format(ana.frames[i].rv_val[which], ana.frames[i].rv_err[which_err]))
    #print("{} +/- {}".format(ana.frames[i].rv_rel[which], ana.frames[i].rv_rel_err[which_err]))



fig, axs = plt.subplots(2, 1, dpi=dpi, figsize=(10, 6), sharex=True)

#-- Plotting every individual order for all frames --#
times = np.array([ana.frames[i].midexpt.datetime for i in use_frames])

for po in order_list:
    if not norm:
        rv_val = np.array([ana.frames[i].data[fiber][po].rv_val[which]     for i in use_frames])
        rv_err = np.array([ana.frames[i].data[fiber][po].rv_err[which] for i in use_frames])
    else:
        rv_val = np.array([ana.frames[i].data[fiber][po].rv_rel[which]     for i in use_frames])
        rv_err = np.array([ana.frames[i].data[fiber][po].rv_rel_err[which] for i in use_frames])

    comb_corr = np.array([ana.frames[i].data[fiber][po].comb_corr for i in use_frames])
    comb_err  = np.array([ana.frames[i].data[fiber][po].comb_err for i in use_frames])
    comb = comb_corr-comb_corr[0]

    med = np.median(rv_val)
    std = np.std(rv_val)

    wh = np.where( np.abs(rv_val - med) < 50000 * rv_err ) # 5 * rv_err )
    
    plt.sca(axs[0])
    plt.errorbar(times[wh], rv_val[wh], rv_err[wh],   marker=pick_marker(po), color=pick_color(po),
                                             linewidth=0.75, alpha=0.3, label=po, capsize=3)
    plt.sca(axs[1])
    plt.errorbar(times[wh], comb[wh], comb_err[wh],   marker=pick_marker(po), color=pick_color(po),
                                             linewidth=0.75, alpha=0.3, label=po)

# ---- Now plotting the frame averages ----------------------
if not norm:
    rv_rels = np.array([ana.frames[i].rv_val[which] for i in use_frames])
    rv_errs = np.array([ana.frames[i].rv_err[which] for i in use_frames])
else:
    rv_rels = np.array([ana.frames[i].rv_rel[which] for i in use_frames])
    rv_errs = np.array([ana.frames[i].rv_rel_err[which] for i in use_frames])

plt.sca(axs[0])
plt.errorbar(times, rv_rels, rv_errs, fmt='o', linewidth=1.5, color='k', capsize=3)

plt.xticks(rotation=45)
if legend: plt.legend(ncol=3)
plt.tight_layout()

if saveas is not None:
    plt.savefig(saveas)
if show:
    plt.show()


