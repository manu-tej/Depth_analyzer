#!/usr/bin/env python3
import subprocess
import numpy as np
import matplotlib
matplotlib.use('Agg')
from skimage import morphology
from datetime import datetime as dt
import sys
import matplotlib.pyplot as plt
from matplotlib import ticker, cm
import csv
import matplotlib.patches as mpatches

def draw_plot(data, edge_color, fill_color,ax):
    bp = ax.boxplot(data.values(), patch_artist=True)

    for element in ['boxes', 'whiskers', 'fliers', 'means', 'medians', 'caps']:
        plt.setp(bp[element], color=edge_color)

    for patch in bp['boxes']:
        patch.set(facecolor=fill_color)

def calc_change(trial,mask,hourlyChange):

    print("Downloading the log file")
    try:
        subprocess.call(['rclone', 'copy', 'manu:BioSci-McGrath/Apps/CichlidPiData/'+trial+'/Logfile.txt', './'+ trial+ '/'])
        log_path = './'+ trial+ '/'+"Logfile.txt"
    except:
        subprocess.call(['rclone', 'copy', 'manu:BioSci-McGrath/Apps/CichlidPiData/'+trial+'/LogFile.txt', './'+ trial+ '/'])
        log_path = './'+ trial+ '/'+"LogFile.txt"

    print("Making the timestamps")
    timeStamps = []
    index = 0
    with open(log_path,'r') as f:
        for line in f:
            line = line.rstrip()
            info_type = line.split(':')[0]
            if info_type == "FrameCaptured":
                time = line.split(',')[4]
                stamp = dt.strptime(time[6:], '%Y-%m-%d %H:%M:%S.%f')
                stamp = stamp.timestamp()
                if index == 0:
                    correction = stamp
                timeStamps.append(round((stamp-correction)/3600,2))
                index += 1



    i =0



    while( i + 12 < len(timeStamps)):
        print("(" + str(i/12+1) + "/"+str(hourlyChange.shape[0])+")")

        firstIndex = i
        if timeStamps[i+12] == timeStamps[i] + 1:
            lastIndex = i+11
        elif timeStamps[i+12] > timeStamps[i] +1:
            lastIndex = i+11
        elif timeStamps[i+12] < timeStamps[i] +1:
            lastIndex = i+11
        i = lastIndex + 1



        mask_t = smoo[lastIndex] - smoo[firstIndex]
        mask_t = np.absolute(mask_t)
        mask_t[mask_t < 0.2] = 0
        mask_t[np.isnan(mask_t)] = 0
        mask[int(firstIndex/12)] = morphology.remove_small_objects(mask_t.astype(bool), 1000)
        hourlyChange[int(firstIndex/12)]  = smoo[lastIndex] - smoo[firstIndex]

    np.save('./'+ trial+ '/'+"hourlyChange.npy", hourlyChange)
    np.save('./'+ trial+ '/'+ "hourlyMask_new.npy", mask)

    print("Made the LARGE numpy arrays")




def plotter(hourlyChange,mask,length_trial,bins,count):
    x_values = {-1:[], -1.33333:[], -1.666666:[], -2:[], -2.33333:[], -2.6666666:[], -3:[], 4:[]}
    up_values = {-1:[], -1.33333:[], -1.666666:[], -2:[], -2.33333:[], -2.6666666:[], -3:[], 4:[]}
    down_values = {-1:[], -1.33333:[], -1.666666:[], -2:[], -2.33333:[], -2.6666666:[], -3:[], 4:[]}
    abs_values = {-1:[], -1.33333:[], -1.666666:[], -2:[], -2.33333:[], -2.6666666:[], -3:[], 4:[]}


    for i in range(hourlyChange.shape[0]):
        tdata = hourlyChange[i].copy()
        #print(tdata.shape)
        #print(mask[i].shape)
        tdata[mask[i] == False] = 0
        tdata[tdata == 0] = np.nan
        tdata = tdata[~np.isnan(tdata)]
        abs_data = np.absolute(tdata)


        for power in sorted(list(x_values.keys())):
            try:
                if power == 4:
                    new_thresh = 0
                else:
                    new_thresh = np.percentile(abs_data, (1-10**power)*100)
            except IndexError:
               # print("IndexError in threshold: " + str(power))
                x_values[power].append(0)
                up_values[power].append(0)
                down_values[power].append(0)
                abs_values[power].append(0)

            else:
                chan = tdata[abs_data > new_thresh]
                x_values[power].append(-np.nansum(chan))
                up_values[power].append(-np.nansum(chan[chan > 0]))
                down_values[power].append(-np.nansum(chan[chan < 0]))
                abs_values[power].append(np.nansum(chan[chan > 0]) - np.nansum(chan[chan < 0]))

        #if len(abs_data) != 0:
        bins.append([trial,0.2,1000,i,len(abs_data),-np.nansum(tdata),np.nansum(abs_data),x_values[-1][-1],abs_values[-1][-1],x_values[-1.33333][-1],abs_values[-1.33333][-1],x_values[-1.666666][-1],abs_values[-1.666666][-1],x_values[-2][-1],abs_values[-2][-1],x_values[-2.33333][-1],abs_values[-2.33333][-1],x_values[-2.6666666][-1],abs_values[-2.6666666][-1],x_values[-3][-1], abs_values[-3][-1]])


    for power in sorted(list(x_values.keys())):



        # difference in actual depth plot
        avg_change[power].append(np.mean(np.array(x_values[power])[np.array(x_values[power]) != 0]))
        min_dif[power].append(min(x_values[power]))
        max_dif[power].append(max(x_values[power]))
        if len(np.array(x_values[power])[np.array(x_values[power]) > 0]) != 0 and len(np.array(x_values[power])[np.array(x_values[power]) != 0]) != 0:
            max_height[power].append(np.sum(np.array(x_values[power])[np.array(x_values[power]) > 0])/len(np.array(x_values[power])[np.array(x_values[power]) != 0]))
        else:
            max_height[power].append(0)
        if len(np.array(x_values[power])[np.array(x_values[power]) < 0]) != 0 and len(np.array(x_values[power])[np.array(x_values[power]) != 0]) != 0:
            min_depth[power].append(np.sum(np.array(x_values[power])[np.array(x_values[power]) < 0])/len(np.array(x_values[power])[np.array(x_values[power]) != 0]))
        else:
            min_depth[power].append(0)
        #print(str(len(min_depth[power]))+ "  " + str(len(max_height[power])))
        #print(str((min_depth[power][-1]))+ "  " + str((max_height[power][-1])))
    '''
    plt.figure(figsize=(35,20))
    plt.plot(list(range(len(x_values[4]))),x_values[4], label = 'Overall Change' )
    plt.plot(list(range(len(x_values[4]))),up_values[4], label = 'Digging')
    plt.plot(list(range(len(x_values[4]))),down_values[4], label = 'Building')
    plt.plot(list(range(len(x_values[4]))),abs_values[4], label = 'Absolute')
    plt.legend()
    plt.grid(True)
    plt.xlabel("Time (h)")
    plt.ylabel("Volume Change")
    plt.savefig('./'+ trial + "_time_series.png")
    plt.close()

    print("Uploading the plot")
    subprocess.call(['rclone', 'copy', './'+ trial+ "_time_series.png",'manu:BioSci-McGrath/Apps/CichlidPiData/___Manu/Building_trials_series/'])
    '''
    if count == 67:
        count +=1
        cv_handle = mpatches.Patch(color='blue', label='CV')
        mc_handle = mpatches.Patch(color='red', label='MC')
        ti_handle = mpatches.Patch(color='green', label='TI')
        f1_handle = mpatches.Patch(color='gray', label='TIxMC F1')
        f1_mc_cv_handle = mpatches.Patch(color='black', label='MCxCV F1')
        f2_mc_cv_handle = mpatches.Patch(color='cyan', label='MCxCV F2')
        f2_ti_mc_handle = mpatches.Patch(color='purple', label='TIxMC F2')

        print("Making the final plots")

        for power in sorted(list(x_values.keys())):
            '''
            plt.figure(figsize=(12,9))
            plt.scatter(list(range(count)), min_dif[power], c = color, s = 15,  marker='x')
            plt.scatter(list(range(count)), max_dif[power], c = color, s = 15, marker='o')

            plt.legend(handles = [cv_handle,mc_handle,ti_handle,f1_handle, f1_mc_cv_handle, f2_mc_cv_handle, f2_ti_mc_handle], loc = 'best')
            plt.grid(True)
            plt.xlabel("Individual")
            plt.ylabel("Volume Change")
            plt.savefig('./'+ str(power)+"_diff.png")
            plt.close()

            '''
            my_dict = {}
            my_dict['Pit'] = max_height[power][:16] + max_height[power][27:32]
            my_dict['Castle'] = max_height[power][16:27]
            #my_dict['TI'] = max_height[power][27:32]
            my_dict['F1'] = max_height[power][32:39] + max_height[power][39:45]
            #my_dict['MCxCV F1'] = max_height[power][39:45]
            my_dict['F2'] = max_height[power][45:53] + max_height[power][53:]
            #my_dict['TIxMC F2'] = max_height[power][53:]
            fig,ax = plt.subplots(figsize=(12,9))
            draw_plot(my_dict,'red','tan',ax)

            my_dict = {}
            my_dict['Pit'] = min_depth[power][:16] + min_depth[power][27:32]
            my_dict['Castle'] = min_depth[power][16:27]
            #my_dict['TI'] = min_depth[power][27:32]
            my_dict['F1'] = min_depth[power][32:39] + min_depth[power][39:45]
            #my_dict['MCxCV F1'] = min_depth[power][39:45]
            my_dict['F2'] = min_depth[power][45:53] + min_depth[power][53:]
            #my_dict['TIxMC F2'] = min_depth[power][53:]
            draw_plot(my_dict,'blue','cyan',ax)
            ax.set_xticklabels(my_dict.keys())
            plt.savefig('./'+ str(power)+"_seg_box.png")
            plt.close()

            my_dict['Pit'] = avg_change[power][:16] + avg_change[power][27:32]
            my_dict['Castle'] = avg_change[power][16:27]
            #my_dict['TI'] = avg_change[power][27:32]
            my_dict['F1'] = avg_change[power][32:39] + avg_change[power][39:45]
            #my_dict['MCxCV F1'] = avg_change[power][39:45]
            my_dict['F2'] = avg_change[power][45:53] + avg_change[power][53:]
            #my_dict['TIxMC F2'] = avg_change[power][53:]
            fig,ax = plt.subplots()
            ax.boxplot(my_dict.values())
            ax.set_xticklabels(my_dict.keys())
            plt.savefig('./'+ str(power)+"_up_down.png")
            plt.close()
            '''

            plt.figure(figsize=(12,9))
            plt.scatter(list(range(count)), max_height[power], c = color, s = 15,  marker='x')
            plt.scatter(list(range(count)), min_depth[power], c = color, s = 15, marker='o')

            plt.legend(handles = [cv_handle,mc_handle,ti_handle,f1_handle, f1_mc_cv_handle, f2_mc_cv_handle, f2_ti_mc_handle], loc = 'best')
            plt.grid(True)
            plt.xlabel("Individual")
            plt.ylabel("Volume Change")
            plt.savefig('./'+ str(power)+"_up_down.png")
            plt.close()
            
            plt.figure(figsize=(12,9))
            plt.scatter(list(range(count)), avg_change[power], c = color, s = 15, marker='o')

            plt.legend(handles = [cv_handle,mc_handle,ti_handle,f1_handle, f1_mc_cv_handle, f2_mc_cv_handle, f2_ti_mc_handle], loc = 'best')
            plt.grid(True)
            plt.xlabel("Individual")
            plt.ylabel("Average Volume Change")
            plt.savefig('./'+ str(power)+"_avg.png")
            plt.close()
            '''
            print("Uploading the final plots to Dropbox")
            #subprocess.call(['rclone', 'copy', './'+ str(power)+ "_diff.png",'manu:BioSci-McGrath/Apps/CichlidPiData/___Manu/Building_trials_overall/'+'Diff'])
            subprocess.call(['rclone', 'copy', './'+ str(power)+ "_up_down.png",'manu:BioSci-McGrath/Apps/CichlidPiData/___Manu/Building_trials_overall/'+'UP_down'])
            subprocess.call(['rclone', 'copy', './'+ str(power)+ "_seg_box.png",'manu:BioSci-McGrath/Apps/CichlidPiData/___Manu/Building_trials_overall/'+'box_plots'])



trials = ['CV1_1', 'CV2_1', 'CV3_1_(NightRecording)', 'CV3_2', 'CV4_1_(NightRecording)', 'CV4_2', 'CV4_3', 'CV7_1', 'CV8_0','CV8_3','CV8_4', 'CV10_3', 'CV10_4', 'CV11_1', 'CV11_2', 'CV12_2', 'MC1_2', 'MC2_1','MC3_1', 'MC6_4', 'MC6_5', 'MC9_1','MC16_1','MC16_2','MC15_3','MC20_1','MC20_2','MC21_1', 'TI2_3', 'TI2_4', 'TI3_3','TI6_2','TI10_1' ,'TIxMCF1_1_1', 'TIxMCF1_1_2', 'TIxMCF1_1_3', 'TIxMCF1_2_1', 'TIxMCF1_3_2', 'TIxMCF1_4_1', 'TIxMCF1_5_1', 'MCxCVF1_4_7', 'MCxCVF1_3_3', 'MCxCVF1_12a_1', 'MCxCVF1_12b_1', 'MCxCVF1_12a_2', 'MCxCVF1_12b_2', '_newtray_MCxCVF2_2_2', '_newtray_MCxCVF2_3_1', '_newtray_MCxCVF2_3_2', '_newtray_MCxCVF2_4_2', '_newtray_MCxCVF2_5_1', '_newtray_MCxCVF2_6_1', '_newtray_MCxCVF2_7_2', '_newtray_MCxCVF2_8_1', '_newtray_TIxMCF2_1', '_newtray_TIxMCF2_10_1', '_newtray_TIxMCF2_2_2', '_newtray_TIxMCF2_2_3', '_newtray_TIxMCF2_3_2', '_newtray_TIxMCF2_4_1', '_newtray_TIxMCF2_4_2', '_newtray_TIxMCF2_5_1', '_newtray_TIxMCF2_6_1', '_newtray_TIxMCF2_7_1', '_newtray_TIxMCF2_8_1', '_newtray_TIxMCF2_9_1', '_newtray_TIxMCF2_11_1','_newtray_TIxMCF2_12_1']

max_dif = {-1:[], -1.33333:[], -1.666666:[], -2:[], -2.33333:[], -2.6666666:[], -3:[], 4:[]}
min_dif = {-1:[], -1.33333:[], -1.666666:[], -2:[], -2.33333:[], -2.6666666:[], -3:[], 4:[]}
max_height = {-1:[], -1.33333:[], -1.666666:[], -2:[], -2.33333:[], -2.6666666:[], -3:[], 4:[]}
min_depth = {-1:[], -1.33333:[], -1.666666:[], -2:[], -2.33333:[], -2.6666666:[], -3:[], 4:[]}
avg_change = {-1:[], -1.33333:[], -1.666666:[], -2:[], -2.33333:[], -2.6666666:[], -3:[], 4:[]}
count = 0
bins = []
color = ['blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue','blue','blue', 'blue', 'red', 'red', 'red', 'red', 'red','red', 'red', 'red', 'red', 'red','red','red', 'green', 'green', 'green','green','green' ,'gray', 'gray', 'gray', 'gray', 'gray', 'gray', 'gray', 'black', 'black', 'black', 'black', 'black', 'black', 'cyan', 'cyan', 'cyan', 'cyan', 'cyan', 'cyan', 'cyan', 'cyan','purple','purple', 'purple', 'purple','purple','purple', 'purple', 'purple','purple','purple', 'purple', 'purple', 'purple', 'purple']
for trial in trials:
    print("Running the trial: " + trial )
    try:
        mask = np.load("./"+trial + "/hourlyMask_new.npy")
        hourlyChange = np.load("./"+trial + "/hourlyChange.npy")
    except:

        try:
            smoo = np.load('./'+ trial+ '/'+"smoothedDepthData.npy")
        except:
            print("Downloading the smoothed Depth data")
            try:
                subprocess.call(['rclone', 'copy', 'manu:BioSci-McGrath/Apps/CichlidPiData/'+trial+'/SubAnalysis/smoothedDepthData.npy', './'+ trial+ '/'])
                smoo = np.load('./'+ trial+ '/'+"smoothedDepthData.npy")
            except:
                subprocess.call(['rclone', 'copy', 'manu:BioSci-McGrath/Apps/CichlidPiData/'+trial+'/DepthAnalysis/smoothedDepthData.npy', './'+ trial+ '/'])
                smoo = np.load('./'+ trial+ '/'+"smoothedDepthData.npy")

        mask = np.empty(shape = (int(smoo.shape[0]/12), smoo.shape[1], smoo.shape[2]), dtype = bool)
        hourlyChange = np.empty(shape = (int(smoo.shape[0]/12), smoo.shape[1], smoo.shape[2]))

        print("Creating the mask")
        calc_change(trial,mask,hourlyChange)
        del smoo
    length_trial = int(mask.shape[0])
    plotter(hourlyChange,mask,length_trial,bins,count)
    bins.sort()
    count += 1

print("writing the frame numbers to csv file")
with open("all.csv",'w') as csvfile:
    fieldnames = ['trial','min_depth','min_region', 'frame_num', 'pixels_passing', 'actual_change', 'absolute_change',str((10**-1)*100) + ' act_change',str((10**-1)*100) + ' abs_change',str((10**-1.33333)*100)+ ' act_change',str((10**-1.33333)*100)+ ' abs_change',str((10**-1.666666)*100)+ ' act_change',str((10**-1.666666)*100)+ ' abs_change',str((10**-2)*100) + ' act_change',str((10**-2)*100) + ' abs_change',str((10**-2.33333)*100) + ' act_change',str((10**-2.33333)*100) + ' abs_change',str((10**-2.6666666)*100) + ' act_change',str((10**-2.6666666)*100) + ' abs_change',str((10**-3)*100) + ' act_change',str((10**-3)*100) + ' abs_change' ]
    writer = csv.writer(csvfile, delimiter=',',
                        quotechar='|', quoting=csv.QUOTE_MINIMAL)
    writer.writerow(fieldnames)
    for i in bins:
        writer.writerow(i)
subprocess.call(['rclone', 'copy', "./all.csv",'manu:BioSci-McGrath/Apps/CichlidPiData/___Manu/Building_trials_frames/'])
