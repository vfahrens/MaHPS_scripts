import os
import numpy as np
import matplotlib.pyplot as plt
from operator import itemgetter

import paths_and_files as pf


# id_string = 'ID2908'
# id_string2 = 'ID2864_oldscatred'
# # id_string3 = '51Peg_time'
# id_string3 = 'ID2864_oldscatred'
# results_folder = os.path.join(pf.abs_path_output, id_string)
# results_folder2 = os.path.join(pf.abs_path_output, id_string2)
# results_folder3 = os.path.join(pf.abs_path_output, id_string3)
# RV_file_path = os.path.join(results_folder, 'RVs_time_weighted_200309.txt')
# RV_file_path2 = os.path.join(results_folder2, 'RVs_time_weighted_corrtels_oldscatred.txt')
# RV_file_path3 = os.path.join(results_folder3, 'RVs_time_weighted_oldscatred.txt')


RV_file_path = os.path.join(pf.abs_path_rvout, 'RVs_time_weighted_200309.txt')
RV_file_path2 = os.path.join(pf.abs_path_rvout, 'RVs_time_weighted_corrtels_oldscatred.txt')
RV_file_path3 = os.path.join(pf.abs_path_rvout, 'RVs_time_weighted_oldscatred.txt')


dates = []
RVs = []
RVerr = []
this_list = []
dates2 = []
RVs2 = []
RVerr2 = []
this_list2 = []
dates3 = []
RVs3 = []
RVerr3 = []
this_list3 = []

# read the data of ups And
with open(pf.out_RVs_weighted.format('2864')) as RVfile:
    for line in RVfile:
        line = line.strip()
        line = line.split()
        this_list.append(line)

nicely = sorted(this_list, key=itemgetter(0))
nicely2 = np.transpose(nicely)
dates = nicely2[0]
RVs = nicely2[1]
RVerr = nicely2[2]

dates_fl = []
RVs_fl = []
RVerr_fl = []
for i in dates:
    dates_fl.append(float(i))
for j in RVs:
    RVs_fl.append(float(j))
for k in RVerr:
    RVerr_fl.append(float(k))

dates_fl = np.asarray(dates_fl).astype(np.float)
RVs_fl = np.asarray(RVs_fl).astype(np.float)
RVerr_fl = np.asarray(RVerr_fl).astype(np.float)

# read the data of HD9407
with open(pf.out_tels_weighted.format('2864')) as RVfile2:
    for line in RVfile2:
        line = line.strip()
        line = line.split()
        this_list2.append(line)

nicely_2 = sorted(this_list2, key=itemgetter(0))
nicely2_2 = np.transpose(nicely_2)
dates2 = nicely2_2[0]
RVs2 = nicely2_2[1]
RVerr2 = nicely2_2[2]

dates_fl2 = []
RVs_fl2 = []
RVerr_fl2 = []
for i in dates2:
    dates_fl2.append(float(i))
for j in RVs2:
    RVs_fl2.append(float(j))
for k in RVerr2:
    RVerr_fl2.append(float(k))

dates_fl2 = np.asarray(dates_fl2).astype(np.float)
RVs_fl2 = np.asarray(RVs_fl2).astype(np.float)
RVerr_fl2 = np.asarray(RVerr_fl2).astype(np.float)

# read the data of 51Peg
with open(pf.out_RVs_telcorr.format('2864')) as RVfile3:
    for line in RVfile3:
        line = line.strip()
        line = line.split()
        this_list3.append(line)

nicely_3 = sorted(this_list3, key=itemgetter(0))
nicely2_3 = np.transpose(nicely_3)
dates3 = nicely2_3[0]
RVs3 = nicely2_3[1]
RVerr3 = nicely2_3[2]

dates_fl3 = []
RVs_fl3 = []
RVerr_fl3 = []
for i in dates3:
    dates_fl3.append(float(i))
for j in RVs3:
    RVs_fl3.append(float(j))
for k in RVerr3:
    RVerr_fl3.append(float(k))

RVs_fl = np.array(RVs_fl)
RVs_fl2 = np.array(RVs_fl2)
RVs_fl3 = np.array(RVs_fl3)
RVerr_fl2 = np.array(RVerr_fl2)
RVerr_fl3 = np.array(RVerr_fl3)
# print(np.std(RVs_fl3 - RVs_fl2))

# plot the RVs for all frames after each other
x = np.arange(len(RVs_fl))
x2 = np.arange(len(RVs_fl2))
x3 = np.arange(len(RVs_fl3))
fig = plt.figure()
plt.errorbar(x, RVs_fl, yerr=RVerr_fl, fmt='o', label='51Peg', alpha=0.5)
plt.errorbar(x2, RVs_fl2, yerr=RVerr_fl2, fmt='o', label='tellurics', alpha=0.5)
plt.errorbar(x3, RVs_fl3, yerr=RVerr_fl3, fmt='o', label='51Peg - tel', alpha=0.5)
# plt.plot(x2, RVs_fl - RVs_fl2, 'o', label='51Peg - tel', alpha=0.5)
# plt.hlines(rv_weightmean, [0], len(x), lw=2)
plt.xlabel('# of observation')
plt.ylabel('RV in m/s')
plt.legend()
plt.show()
# plt.savefig('all3_RV_notime.pdf')

# # plot the RVs at the correct observation time
# fig = plt.figure()
# # plt.errorbar(dates_fl, RVs_fl, yerr=RVerr_fl, fmt='o', label='upsAnd', alpha=0.5)
# plt.errorbar(dates_fl2, RVs_fl2, yerr=RVerr_fl2, fmt='o', label='51Peg_tels', alpha=0.5)
# plt.errorbar(dates_fl3, RVs_fl3, yerr=RVerr_fl3, fmt='o', label='51Peg', alpha=0.5)
# # plt.hlines(rv_weightmean, [0], len(x), lw=2)
# plt.xlabel('date of observation')
# plt.ylabel('RV in m/s')
# plt.legend()
# plt.show()
# # plt.savefig('all3_RV_time.pdf')
#
#


