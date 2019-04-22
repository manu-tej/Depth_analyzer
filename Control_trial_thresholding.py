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

        firstIndex = i
        if timeStamps[i+12] == timeStamps[i] + 1:
            lastIndex = i+11
        elif timeStamps[i+12] > timeStamps[i] +1:
            lastIndex = i+10
        elif timeStamps[i+12] < timeStamps[i] +1:
            lastIndex = i+12
        i = lastIndex + 1



        for min_dep in range(6):
            mask_t = smoo[lastIndex] - smoo[firstIndex]
            mask_t = np.absolute(mask_t)
            mask_t[mask_t < (0.05 + min_dep* 0.01)] = 0
            mask_t[np.isnan(mask_t)] = 0
            for min_pix in range(6):
                mask[int(firstIndex/12),min_pix, min_dep] = morphology.remove_small_objects(mask_t.astype(bool), 500 + min_pix*100)
        hourlyChange[int(firstIndex/12)]  = smoo[lastIndex] - smoo[firstIndex]

    np.save('./'+ trial+ '/'+"hourlyChange.npy", hourlyChange)
    np.save('./'+ trial+ '/'+ "hourlyMask_new.npy", mask)

    print("Uploading the files to Dropbox")
    subprocess.call(['rclone', 'copy', './'+ trial+ '/'+"hourlyChange.npy",'manu:BioSci-McGrath/Apps/CichlidPiData/___Manu/Control_trials/'+trial])
    subprocess.call(['rclone', 'copy', './'+ trial+ '/'+"hourlyMask_new.npy",'manu:BioSci-McGrath/Apps/CichlidPiData/___Manu/Control_trials/'+trial])



def plotter(hourlyChange,mask,length_trial):
    x = [[0  for i in range(6)] for i in range(6)]
    y = [[0  for i in range(6)] for i in range(6)]
    z = [[0  for i in range(6)] for i in range(6)]

    X = np.linspace(500,1000,6,dtype=int)
    Y = np.linspace(0.05,0.1,6)

    X,Y = np.meshgrid(X,Y)


    for i in range(hourlyChange.shape[0]):
        for j in range(6):
            for k in range(6):
                tdata = hourlyChange[i].copy()
                tdata[mask[i,j,k] == False] = 0
                tdata[tdata == 0] = np.nan
                tdata = tdata[~np.isnan(tdata)]
                abs_data = np.absolute(tdata)

                x[k][j] += np.sum(abs_data)
                y[k][j] += len(abs_data)
                z[k][j] += np.nansum(tdata)
    x = [[i/length_trial for i in n] for n in x]
    y = [[i/length_trial for i in n] for n in y]
    z = [[i/length_trial for i in n] for n in z]
    fig, ax = plt.subplots() 

    cs = ax.contourf(X,Y,x, locator=ticker.LinearLocator(), cmap=cm.PuBu_r)
    ax.set_xlabel("Threshold Size in Pixels (px)")
    ax.set_ylabel("Threshold depth Change")
    ax.set_title("Total Change in "+ trial )
    cbar = fig.colorbar(cs)
    plt.savefig('./'+ trial + "_total_change.png")

   
    fig,ax = plt.subplots()
    cs = ax.contourf(X,Y,y, locator=ticker.LinearLocator(), cmap=cm.PuBu_r)
    ax.set_xlabel("Threshold Size in Pixels (px)")
    ax.set_ylabel("Threshold depth Change")
    ax.set_title("Pixel Change in "+ trial )
    cbar = fig.colorbar(cs)
    plt.savefig('./'+ trial + "_pixel_change.png")

    fig,ax = plt.subplots()
    cs = ax.contourf(X,Y,z, locator=ticker.LinearLocator(), cmap=cm.PuBu_r)
    ax.set_xlabel("Threshold Size in Pixels (px)")
    ax.set_ylabel("Threshold depth Change")
    ax.set_title("Actual Change in "+ trial )
    cbar = fig.colorbar(cs)
    plt.savefig('./'+ trial + "_actual_change.png")
    plt.clf()

    print("Uploading the plots")
    subprocess.call(['rclone', 'copy', './'+ trial+ "_total_change.png",'manu:BioSci-McGrath/Apps/CichlidPiData/___Manu/Control_trials_plots/'])
    subprocess.call(['rclone', 'copy', './'+ trial+ "_actual_change.png",'manu:BioSci-McGrath/Apps/CichlidPiData/___Manu/Control_trials_plots/'])
    subprocess.call(['rclone', 'copy', './'+ trial+ "_pixel_change.png",'manu:BioSci-McGrath/Apps/CichlidPiData/___Manu/Control_trials_plots/'])



trials = ['empty_con1', 'empty_con2', 'empty_con3', 'empty_con4', 'empty_con5', 'empty_con6', 'MC_male_con1', 'MC_male_con2', 'MC_male_con3', 'MC_male_con4', 'CV_male_con1','CV_male_con2', 'CV_male_con3', 'CV_male_con4', 'MC_fem_con1', 'MC_fem_con2', 'MC_fem_con3', 'CV_fem_con1', 'CV_fem_con2', 'CV_fem_con3']


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
            subprocess.call(['rclone', 'copy', 'manu:BioSci-McGrath/Apps/CichlidPiData/'+trial+'/SubAnalysis/smoothedDepthData.npy', './'+ trial+ '/'])
            smoo = np.load('./'+ trial+ '/'+"smoothedDepthData.npy")

        mask = np.empty(shape = (int(smoo.shape[0]/12),6,6, smoo.shape[1], smoo.shape[2]), dtype = bool)
        hourlyChange = np.empty(shape = (int(smoo.shape[0]/12), smoo.shape[1], smoo.shape[2]))

        print("Creating the mask")
        calc_change(trial,mask,hourlyChange)
        del smoo
    length_trial = int(mask.shape[0])
    plotter(hourlyChange,mask,length_trial)
