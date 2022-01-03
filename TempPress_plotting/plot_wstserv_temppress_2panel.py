#!/usr/bin/env python
# -*- coding: utf-8 -*-

# this script is used for plotting environmental data of FOCES that were collected in the "new" format, i.e. logfiles
# on the wstserver with names like "foces_temp_tank_190419.log"

import numpy as np

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.dates import DateFormatter
from matplotlib.ticker import FormatStrFormatter
from operator import itemgetter
import julian
from datetime import datetime
import datetime as dt

yearsFmt = DateFormatter('%Y-%m-%d %H:%M')
matplotlib.rcParams['axes.formatter.useoffset'] = False

##############################
##############################
# TO BE SPECIFIED BY THE USER

path = u'/mnt/e/logfiles/meta/'
# skip = 1  # plottet jeden Datenpunkt
skip = 10  # ueberspringt die angegebende Anzahl von Datenpunkten umd das Program schneller zu machen

# define the files from which the logged values should be read
ptc10new_logs = 'ptc10new_20201110-19.log.gz'
# ptc10old_logs = 'ptc10old_20201101-1125.log.gz'
wut_room_logs = 'wut_foces-room_20201110-19.log.gz'
# cpg2500_logs = 'cpg2500_20201101-1125.log.gz'
# ptc10new_logs = 'ptc10new_20200901-1125.log.gz'
# ptc10old_logs = 'ptc10old_20200901-1125.log.gz'
# wut_room_logs = 'wut_foces-room_20200901-1125.log.gz'
# cpg2500_logs = 'cpg2500_20200901-1125.log.gz'

# Definiert den anzuzeigenden Zeitbereich
t_end = datetime.strptime('2020-11-15 06:00:00', '%Y-%m-%d %H:%M:%S')
t_delta = dt.timedelta(days=0.5)

##############################
##############################

ab = t_end - t_delta
bis = t_end

# ##############################################
# # EINLESE DES DRUCKS CPG2500
#
# time_cpg2500 = []
# press_cpg2500 = []
#
# with open(path + cpg2500_logs) as cpg_file:
#     for cline in cpg_file:
#         cline = cline.split()
#         if len(cline) == 2:  # and len(cline[0]) == 19:
#             cpgtimestamp = datetime.strptime(cline[0], '%Y-%m-%dT%H:%M:%S:')
#             if ab < cpgtimestamp < bis and cline[1] != '0' and cline[1] != '-1':
#                 time_cpg2500.append(cpgtimestamp)
#                 press_cpg2500.append(float(cline[1]))

##############################################
# EINLESE DES DRUCKS  -- WUT Baro Room

time_baro_room = []
temp_baro_room = []
humi_baro_room = []
press_baro_room = []

with open(path + wut_room_logs) as wut_room_file:
    for wline in wut_room_file:
        wline = wline.split()
        if len(wline) == 4:  # and len(line[0]) == 19:
            wuttimestamp = datetime.strptime(wline[0], '%Y-%m-%dT%H:%M:%S:')
            if ab < wuttimestamp < bis:
                time_baro_room.append(wuttimestamp)
                temp_baro_room.append(float(wline[1][-4:]))
                humi_baro_room.append(float(wline[2][-4:]))
                press_baro_room.append(float(wline[3][-5:]))

##############################################
# EINLESE DER TEMPERATUR ---  Regel-Sensoren

time_ptc10 = []
temp_ptc10_3A = []
temp_ptc10_3B = []
temp_ptc10_3C = []
temp_ptc10_3D = []
temp_ptc10_inl1 = []
temp_ptc10_inl2 = []
temp_ptc10_inl3 = []
temp_ptc10_inl4 = []
pow_ptc10_Out1 = []
pow_ptc10_Out2 = []

with open(path + ptc10new_logs) as ptc10_file:
    i = 0
    for line in ptc10_file:
        i += 1
        if i % skip == 0:
            line = line.split()
            # Skips rows where the Temp-Controller was not answering
            if len(line) == 26 and len(line[0]) == 20:
                timestamp = datetime.strptime(line[0], '%Y-%m-%dT%H:%M:%S:')
                if ab < timestamp < bis:
                    time_ptc10.append(timestamp)
                    temp_ptc10_3A.append(float(line[9]))
                    temp_ptc10_3B.append(float(line[10]))
                    temp_ptc10_3C.append(float(line[11]))
                    temp_ptc10_3D.append(float(line[12]))
                    temp_ptc10_inl1.append(float(line[13]))
                    temp_ptc10_inl2.append(float(line[14]))
                    temp_ptc10_inl3.append(float(line[15]))
                    temp_ptc10_inl4.append(float(line[16]))
                    # temp_ptc10_Out1.append(float(line[1]))
                    # temp_ptc10_Out2.append(float(line[5]))

print(len(time_ptc10))
print(len(temp_ptc10_3D))

# ##############################################
# # EINLESE DER TEMPERATUR ---  Mess-Sensoren
#
# time_messptc10 = []
# temp_messptc10_3A = []
# temp_messptc10_3B = []
# temp_messptc10_3C = []
# temp_messptc10_3D = []
# pow_messptc10_Out1 = []
# pow_messptc10_Out2 = []
#
# with open(path + ptc10old_logs) as messptc10_file:
#     im = 0
#     for mline in messptc10_file:
#         im += 1
#         if im % skip == 0:
#             mline = mline.split()
#             # Skips rows where the Temp-Controller was not answering
#             if len(mline) == 18 and len(mline[0]) == 20:
#                 messtimestamp = datetime.strptime(mline[0], '%Y-%m-%dT%H:%M:%S:')
#                 if ab < messtimestamp < bis:
#                     time_messptc10.append(messtimestamp)
#                     temp_messptc10_3A.append(float(mline[4]))
#                     temp_messptc10_3B.append(float(mline[5]))
#                     # temp_messptc10_3C.append(float(mline[6]))
#                     temp_messptc10_3D.append(float(mline[7]))
#                     # temp_messptc10_Out1.append(float(mline[1]))
#                     # temp_messptc10_Out2.append(float(mline[3]))

##############################
# ERSTELLEN DES PLOTS
##############################

# define the "sigma" greek symbol
sigma = r'$\sigma$: '

fig = plt.figure('Temp_controlled_Room_Grat_Press', figsize=[15, 10], dpi=100, facecolor='white')

##############################
# Subplot: TANK


# fig, axs = plt.subplots(4, 1, sharex=True, gridspec_kw={'hspace': 0})
# set up 3 separate axes / plots
ax = plt.subplot(211)

ax2 = ax.twinx()
ax2.patch.set_alpha(0.0)

# ax3 = plt.subplot(2, 1, 2, sharex=True, gridspec_kw={'hspace': 0})
#
# ax4 = ax3.twinx()
# ax4.patch.set_alpha(0.0)
# # ax3.spines['right'].set_position(('axes', 1.15))
#
# # ax4 = ax.twinx()
# # ax4.patch.set_alpha(0.0)
# # ax4.spines['left'].set_position(('axes', -0.15))
# #
# # # some voodoo to make a second y axis on the left work
# # def make_patch_spines_invisible(whichax):
# #     whichax.set_frame_on(True)
# #     whichax.patch.set_visible(False)
# #     for sp in whichax.spines.values():
# #         sp.set_visible(False)
# #
# # make_patch_spines_invisible(ax4)
# #
# # ax4.spines["left"].set_visible(True)
# # ax4.yaxis.set_label_position('left')
# # ax4.yaxis.set_ticks_position('left')
#
#
# # define which quantities should be plotted in comparison to the RVs
plot1_x = time_baro_room
plot1_y = temp_baro_room
# plot2_x = time_messptc10
# plot2_y = temp_messptc10_3A
# plot3_x = time_ptc10
# plot3_y = temp_ptc10_3D
# plot4_x = time_cpg2500
# plot4_y = press_cpg2500
#
# # calculate peak-to-valley and standard deviations for the legends
# pv_plot1 = np.max(plot1_y) - np.min(plot1_y)
# std_plot1 = np.std(plot1_y)
# pv_plot2 = np.max(plot2_y) - np.min(plot2_y)
# std_plot2 = np.std(plot2_y)
# pv_plot3 = np.max(plot3_y) - np.min(plot3_y)
# std_plot3 = np.std(plot3_y)
# pv_plot4 = np.max(plot4_y) - np.min(plot4_y)
# std_plot4 = np.std(plot4_y)
#
# # define the labels for the sensors and their peak-to-valley and standard deviation values
# # rv_label = 'ThAr RVs (PV: {:.1f} m/s, {}{:.1f} m/s)'.format(pv_rvs, sigma, std_rvs)  #################################
# std_grat = 'grating (air) temperature (PV: {:.4f} K, {}{:.4f} K)'.format(pv_plot1, sigma, std_plot1)
# std_coll = 'collimator (air) temperature (PV: {:.4f} K, {}{:.4f} K)'.format(pv_plot2, sigma, std_plot2)  # '°C'
# std_room = 'room temperature (PV: {:.2f} K, {}{:.2f} K)'.format(pv_plot3, sigma, std_plot3)
# std_cpg = 'tank pressure: (PV: {:.2f} hPa, {}{:.2f} hPa)'.format(pv_plot4, sigma, std_plot4)
#
# # set the labels that should be used on the x and y axes
# ax1.set_xlabel("time [h]", fontsize=20)
# ax1.set_ylabel("temperature [°C]", fontsize=20)
# ax2.set_ylabel("temperature [°C]", fontsize=20)
# ax3.set_ylabel("temperature [°C]", fontsize=20)
# ax4.set_ylabel("pressure [hPa]", fontsize=20)
#
# # make the actual plots
p1, = ax.plot_date(plot1_x, plot1_y, fmt='-', color='tab:orange', alpha=1.0, linewidth=2.5)  # , label=std_grat)
# p2, = ax2.plot_date(plot2_x, plot2_y, fmt='-', color='darkblue', alpha=0.8, linewidth=2.5, label=std_coll)
# p3, = ax3.plot_date(plot3_x, plot3_y, fmt='-', color='royalblue', alpha=0.8, linewidth=2.5, label=std_room)
# p4, = ax4.plot_date(plot4_x, plot4_y, fmt='-', color='tab:orange', alpha=0.6, linewidth=2.5, label=std_cpg)
#
# # calculate the offset of the two temp axes from each max/min range
# mid_plot1 = np.min(plot1_y) + (0.5 * pv_plot1)
# mid_plot2 = np.min(plot2_y) + (0.5 * pv_plot2)
# axis_offset_temps = mid_plot1 - mid_plot2
#
# # # format the y-axes to span the same range of values for direct comparison
# # ax.set_ylim(16.3, 16.9)
# # ax2.set_ylim(20, 23)
# # new_ax2_ylim = ax.get_ylim() - axis_offset_temps
# # ax2.set_ylim(new_ax2_ylim)
#
# # make the axis labels in the same color as the corresponding graphs
# ax1.yaxis.label.set_color(p1.get_color())
# ax2.yaxis.label.set_color(p2.get_color())
# ax3.yaxis.label.set_color(p3.get_color())
# ax4.yaxis.label.set_color(p4.get_color())
# # make the numbers on each axis also the color of the graphs and adjust the font size for better reading
# plt.setp(ax1.get_xticklabels(), fontsize=16, rotation=30, ha='right')
# plt.setp(ax1.get_yticklabels(), color=p1.get_color(), fontsize=16, alpha=1.0)
# plt.setp(ax2.get_yticklabels(), color=p2.get_color(), fontsize=16)
# plt.setp(ax3.get_yticklabels(), color=p3.get_color(), fontsize=16)
# plt.setp(ax4.get_yticklabels(), color=p4.get_color(), fontsize=16)
#
# ax1.xaxis.set_major_formatter(yearsFmt)
#
# # add a hack to show all legends in this multi-axis plot
# lines, labels = ax1.get_legend_handles_labels()
# lines2, labels2 = ax2.get_legend_handles_labels()
# lines3, labels3 = ax3.get_legend_handles_labels()
# lines4, labels4 = ax4.get_legend_handles_labels()
# plt.legend(lines + lines2 + lines3 + lines4, labels + labels2 + labels3 + labels4, loc='center', bbox_to_anchor=(0.5, 1.1), fontsize=18)
#
# # ax3.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

ax.grid(True)
ax = plt.gca()
fmt = matplotlib.ticker.ScalarFormatter(useOffset=False)
fmt.set_scientific(False)

# plt.legend(loc="center", bbox_to_anchor=(0.5, -0.2))

plt.tight_layout()

plt.draw()
# plt.savefig('/mnt/e/plotting-tools/TempPress/wstserver_single_temp_press/plots/other_shortterm_SPIE2020.pdf')
plt.show()
plt.close(fig)
