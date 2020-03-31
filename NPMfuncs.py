# -*- coding: utf-8 -*-
"""
Created on Thu Dec 19 10:39:01 2019

@author: Alex
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math 
from os import listdir
from bisect import bisect_left
myCycler = ['c', 'm', 'y', 'k', 'r', 'g', 'b' ]

#%% plot and optionally save a figure made from a pandas dataframe

def plot_dataFrame(dataframe, fileLocation, plotname, subplots, filetype = '.pdf', save = True):
    if subplots == True:
        plot = dataframe.plot(subplots = True, title = plotname)
        fig = plot[1].get_figure()
    if subplots == False:
        plot = dataframe.plot(subplots = False, title = plotname)
        fig = plot.get_figure()
    if save == True:
        fig.savefig(fileLocation + plotname + filetype)

#%% find all .csv files in folder location
        
def find_csv_filenames(path_to_dir, suffix=".csv" ):
    filenamelist = listdir(path_to_dir)
    return [ filename for filename in filenamelist if filename.endswith( suffix ) ]

#%% Deinterleave Frames using Triggermode 

def deinterleave(expDF, trigMode):
    # find columns containing data 
    green = expDF.columns.str.contains('1G')
    red = expDF.columns.str.contains('0R')
    datacols = red + green
    col_names = tuple(expDF.columns[datacols])
    data_df = expDF.loc[:, col_names]
    timestamp = expDF['Timestamp']    
        
    if trigMode == '3':
    # adjust nRows if can't be divided by 3
        while len(data_df.index)%3 != 0:
            data_df = data_df.iloc[:-1,:]  
        wl_1 = data_df[::3]
        wl_2 = data_df[1::3]
        wl_3 = data_df[2::3]
        wl1_col = [str(col) + '__470' for col in data_df.columns]
        wl2_col = [str(col) + '__560' for col in data_df.columns]
        wl3_col = [str(col) + '__415' for col in data_df.columns]
        split_data = pd.DataFrame(np.concatenate((wl_1, wl_2, wl_3), axis=1), columns=(wl1_col + wl2_col + wl3_col))    
        # split timestamp to sync with behavioural data 
        split_timestamp = timestamp[::3]    
        # make tuple to multiindex by: name, color, wavelength
        colnames = split_data.columns
        split_list = [list(c.split('__')) for c in colnames]
        listwavelength = [split_list[i][1] for i in range(len(colnames))]
        listrest = [split_list[i][0] for i in range(len(colnames))]
        listcolor = list(c[-2:] for c in listrest)
        listname = list(c[:-2] for c in listrest)
        tuple_names = tuple((listname, listcolor, listwavelength))
        tuple_names = list(map(list, zip(*tuple_names)))
        split_data.columns = pd.MultiIndex.from_tuples(tuple_names, names= ('Branch', 'Filter', 'Wavelength' ))     
        split_data = split_data.reindex(sorted(split_data.columns), axis=1)
    elif trigMode == 'CNST':
        split_timestamp = timestamp
        split_data = data_df
    elif trigMode == '1':
        while len(data_df.index)%2 != 0:
            data_df = data_df.iloc[:-1,:]      
        wl_1 = data_df[::2]
        wl_2 = data_df[1::2]
        wl1_col = [str(col) + '__470/560' for col in data_df.columns]
        wl2_col = [str(col) + '__415' for col in data_df.columns]
        split_data = pd.DataFrame(np.concatenate((wl_1, wl_2), axis=1), columns=(wl1_col + wl2_col))  
        # split timestamp to sync with behavioural data 
        split_timestamp = timestamp[::2]        
        # make tuple to multiindex by: name, color, wavelength
        colnames = split_data.columns
        split_list = [list(c.split('__')) for c in colnames]
        listwavelength = [split_list[i][1] for i in range(len(colnames))]
        listrest = [split_list[i][0] for i in range(len(colnames))]
        listcolor = list(c[-2:] for c in listrest)
        listname = list(c[:-2] for c in listrest)
        tuple_names = tuple((listname, listcolor, listwavelength))
        tuple_names = list(map(list, zip(*tuple_names)))
        split_data.columns = pd.MultiIndex.from_tuples(tuple_names, names= ('Branch', 'Filter', 'Wavelength' ))           
        split_data = split_data.reindex(sorted(split_data.columns), axis=1)
    elif trigMode == '2':
        while len(data_df.index)%2 != 0:
            data_df = data_df.iloc[:-1,:]           
        wl_1 = data_df[::2]
        wl_2 = data_df[1::2]
        wl1_col = [str(col) + '__470' for col in data_df.columns]
        wl2_col = [str(col) + '__560' for col in data_df.columns]
        split_data = pd.DataFrame(np.concatenate((wl_1, wl_2), axis=1), columns=(wl1_col + wl2_col))        
        # split timestamp to sync with behavioural data 
        split_timestamp = timestamp[::2]        
        # make tuple to multiindex by: name, color, wavelength
        colnames = split_data.columns
        split_list = [list(c.split('__')) for c in colnames]
        listwavelength = [split_list[i][1] for i in range(len(colnames))]
        listrest = [split_list[i][0] for i in range(len(colnames))]
        listcolor = list(c[-2:] for c in listrest)
        listname = list(c[:-2] for c in listrest)
        tuple_names = tuple((listname, listcolor, listwavelength))
        tuple_names = list(map(list, zip(*tuple_names)))
        split_data.columns = pd.MultiIndex.from_tuples(tuple_names, names= ('Branch', 'Filter', 'Wavelength' ))     
        split_data = split_data.reindex(sorted(split_data.columns), axis=1)
    return split_data, split_timestamp

#%% Make Delta F 
# mode: how to deal with boundaries # valid = only return values not at edges
def deltaFMoving(data, window, calcType = 'median', mode = 'reflect', plotOn = 1):
    # Return delta F over F for signal trace using moving median or mean. Select method to deal with boundaries
    tail = round(window/2) 
    
#% mode that only returns values far away enough from the boundaries
    if mode == 'valid':
        deltaF = pd.DataFrame()
        for iColumn in range(data.shape[1]):
            tempData = data.iloc[:,iColumn]
            med = []
            mean = []        
            for count, i in enumerate(tempData):
                if count >= tail and count < len(tempData) - tail:
                    stretchToUse = tempData[count-int(window/2):count+int(window/2)+1]
                    sortedRange = np.sort(stretchToUse)
                    med = np.append(med,sortedRange[tail])                    
                    mean = np.append(mean, np.mean(stretchToUse))
            if calcType == 'median': 
                deltaF = pd.concat([deltaF, 100*((tempData - med)/med)], axis = 1)
            elif calcType == 'mean':               
                deltaF = pd.concat([deltaF, 100*((tempData - mean)/mean)], axis = 1)
           
#% reflects the data back at both edges and uses this to calculate baseline
  #  if mode == 'reflect':
    if mode == 'reflect':
        deltaF = pd.DataFrame()
        for iColumn in range(data.shape[1]):
            tempData = data.iloc[:,iColumn]
            flipped = np.concatenate([list(reversed(tempData)), list(tempData), list(reversed(tempData))])
            med = []
            mean = []        
            for count, i in enumerate(tempData):
                if count >= tail and count < len(tempData) - tail:
                    stretchToUse = tempData[count-int(window/2):count+int(window/2)+1]
                    sortedRange = np.sort(stretchToUse)
                    med = np.append(med,sortedRange[tail])                    
                    mean = np.append(mean, np.mean(stretchToUse))
                elif count < tail or count >= len(tempData) - tail:
                    stretchToUse = list(flipped[len(tempData)+count-tail:len(tempData)+count+tail+1])
                    sortedRange = np.sort(stretchToUse)
                    med = np.append(med,sortedRange[tail])
                    mean = np.append(mean, np.mean(stretchToUse))
            if calcType == 'median': 
                deltaF = pd.concat([deltaF, 100*((tempData - med)/med)], axis = 1)
            elif calcType == 'mean':               
                deltaF = pd.concat([deltaF, 100*((tempData - mean)/mean)], axis = 1)

    # label the columns 
    deltaF.columns = data.columns

    if plotOn == 1:
        data.plot(title = 'Uncorrected Data', subplots = True)
        deltaF.plot(title = 'Delta F', subplots = True)
              
    return deltaF;

#%% split into PSTH
    
# def psNPM(data, eventLoc, preWindow, postWindow):
    
# #%% find behavioural data corresponding to photometry
#     '''Returns data split and centered around events. NB will convert edge data to NAN'''

#     cutOut = np.empty(((len(eventLoc)), preWindow+postWindow))*np.nan 
#     for count, i in enumerate(eventLoc):
#         if (i-preWindow) >= 0: 
#             stretchToUse = data[i-preWindow:i+postWindow]  
#             cutOut[count] = stretchToUse    
#     plt.figure(0)
#     plt.plot(np.transpose(cutOut), color = 'green', label = 'Delta F')
#     plt.plot(np.nanmean(cutOut,axis = 0), color = 'red', linewidth = 2, label = 'Mean Delta F')
#     plt.axvline(preWindow, color = 'black')
#     plt.title('Trace centered around events')
          
    
#%% Add in lines for single events
# Assumes these are acquired in bonsai on same computer timestamp as photometry
def markKeyDown(behaviour, zdFF, splitTimestamp):
    
    figAll, (ax1) = plt.subplots(zdFF.shape[1]) 
    linked = zdFF
    if zdFF.shape[1]>1:
        for iAxis, axis in enumerate(ax1):
            axis.set_prop_cycle(color = ['black', 'c', 'm', 'y', 'k', 'r', 'g', 'b' ])
            axis.plot(zdFF.iloc[:,iAxis], label = 'zdFF')
    elif zdFF.shape[1] == 1:
        ax1.set_prop_cycle(color = ['black', 'c', 'm', 'y', 'k', 'r', 'g', 'b' ])
        ax1.plot(zdFF, label = 'zdFF')
    
    uniqueKeys = behaviour.iloc[:,0].unique()
    for i, key in enumerate(uniqueKeys):
        boolkey = behaviour.iloc[:,0] == key
        keyData = behaviour[boolkey]
        times = keyData.iloc[:,1]
       # find frames where timestamps match split timestamps    
        positions = []
        for j, time in enumerate(times):
            pos = bisect_left(list(splitTimestamp), time)
            positions.append(pos)
  
            if zdFF.shape[1]>1:
                for iAxis, axis in enumerate(ax1):
                    if j == 0:
                        axis.axvline(x=pos, color = myCycler[i], label = key)
                    else:
                        axis.axvline(x=pos, color = myCycler[i])
            elif zdFF.shape[1] == 1:   
                 if j == 0:
                     ax1.axvline(x=pos, color = myCycler[i], label = key)
                 else:
                     ax1.axvline(x=pos, color = myCycler[i])
                    
        # make column in zdFF where each uniquekey press moment is recorded
        addCol = np.zeros([zdFF.shape[0],1])
        addCol[positions] = 1
        linked['Key: ' + key] = addCol
  
    plt.legend(loc='upper right')
    return linked, figAll

#%% Add in values for status events 
    
def behStatus(beh, zdFF, colsToUse):
# Does not assume that this was acquired on same computer clock as photometry
# but assumes that recording start and end are synchronised
    
    factor = zdFF.shape[0]/beh.shape[0]
    intbeh = pd.DataFrame(np.zeros([zdFF.shape[0], 2]))

  #%% synch by interpolating behavioural data 
    for a, iCol in enumerate(colsToUse):
        for i,j in enumerate(zdFF.iloc[:,0]):
            intbeh.iloc[i,a] = beh.iloc[math.floor(i/factor),iCol]
    
    intbeh.columns = beh.columns.values[colsToUse]       
    linked = pd.concat([intbeh, zdFF], axis=1, join='inner')    
    return linked


         