# -*- coding: utf-8 -*-
"""
# Code for David lab, Paris 2020
Created on Tue Feb 25 10:21:33 2020

Assumes: 
1) Photometry, keydown and behavioural data can be linked as they are all 
labelled with a unique n character prefix, e.g.
'CTX02': CTXO2_photometry.csv, CTX02_keyDown.csv, CTX02_behaviour.csv
Depending on experimental design this will include which animal/trial/condition.
2) Photometry .csv files are in one folder and behavioural .csv files in another
3) Uses zdFF function by Ekaterina Martianova ekaterina.martianova.1@ulaval.ca
  so if you do use this, please reference: Martianova, E., Aronson, S., Proulx, C.D. Multi-Fiber Photometry 
  to Record Neural Activity in Freely Moving Animal. J. Vis. Exp.(152), e60278, doi:10.3791/60278 (2019)

LICENCE
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>
'''

@author: Alex Leighton
"""
#%% import the libraries you need

from NPMfuncs import deinterleave, find_csv_filenames, plot_dataFrame, markKeyDown, behStatus
import pandas as pd
import photometry_functions as ph
import matplotlib.pyplot as plt

#%% user input

# set photometryLocation to where you keep your photometry files, with double \\
photometryLocation = 'C:\\Users\\Alex\\Documents\\Work\\NPM\\Post hoc analysis python\\Test data paris\\Photometry\\'
behaviourLocation = 'C:\\Users\\Alex\\Documents\\Work\\NPM\\Post hoc analysis python\\Test data paris\\Behaviour\\'
saveLocation = 'C:\\Users\\Alex\\Documents\\Work\\NPM\\Post hoc analysis python\\Test data paris\\Plots\\'

# set number of characters that define unique prefix
nChar = 5

# set triggerMode to CSNT, 1, 2 or 3
triggerMode = '1'

# set number of patch cord branches used 
nBranches = 2

# choose whether behaviour is an event (e.g. keydown) or a status (e.g mouse in open arm)
behaviourType = 'Status' 

# are there specific columns that should be used in the behavioural data (e.g Anymaze)
colsToUse = [3,4]

#%% find all .csv files in photometry location
photometryFiles = find_csv_filenames(photometryLocation)
behaviourFiles = find_csv_filenames(behaviourLocation)
print('Will be performing analysis on: ' + str(photometryFiles))

#%% Loop over each experiment
for experiment in photometryFiles:
#%% Load in  data file, plot and save raw data
    exp = pd.read_csv(photometryLocation + experiment, error_bad_lines=False)
    # find unique experiment ID
    expID = experiment[0:5]
    # plot raw data and save plot
    plot_dataFrame(exp, saveLocation, expID + '_Raw Data', True, filetype='.pdf')

#%% Deinterleave alternating frames depending on trigger mode
    splitData, splitTimestamp = deinterleave(exp, triggerMode)
    # plot result and save plot
    plot_dataFrame(splitData, saveLocation, expID + '_Deinterleaved Data', True)
    # save deinterleaved data as .csv
    splitData.to_csv(saveLocation + 'Split.csv', encoding='utf-8', index=False)

#%% Loop over branches 
    idx = pd.IndexSlice
    # these apply for gcamp & isosbestic 
    isosbestic = splitData.loc[:, idx[:, '1G', '415']]
    signalData = splitData.loc[:, idx[:, '1G', '470/560']]
    zdFF = pd.DataFrame()

    if  nBranches == 1:
        f, (ax1, ax2) = plt.subplots(1, 2)
        ax1.plot(isosbestic)
        ax2.plot(signalData)
        ax1.set(title = 'Isosbestic, Branch 1')
        ax2.set(title = 'GCaMP Data, Branch 1')
        #calc zdFF using method by Ekaterina Martianova ekaterina.martianova.1@ulaval.ca 
        zdFFBranch = ph.get_zdFF(isosbestic.squeeze(),signalData.squeeze(), remove = 0)
        f, (ax1) = plt.subplots(1)
        ax1.plot(zdFFBranch,color = 'black')
        ax1.set(title = 'GCaMP zdFF, Branch 1')        
        zdFF['Branch 1'] = zdFFBranch         
    
#%% Choose relevant data (e.g green 470 = gcamp, green 415 = isosbestic) and double check using plot       
    # plot to make sure correct signal is used    
    elif nBranches > 1:
        for iBranch in range(nBranches):
            f, (ax1, ax2) = plt.subplots(1, 2)
            ax1.plot(isosbestic.iloc[:,iBranch], label = isosbestic.columns.levels[0][iBranch])
            ax2.plot(signalData.iloc[:,iBranch], label = signalData.columns.levels[0][iBranch])
            ax1.set(title = 'Isosbestic, Branch ' + str(iBranch))
            ax2.set(title = 'GCaMP Data, Branch ' + str(iBranch))
            #calc zdFF using method by Ekaterina Martianova ekaterina.martianova.1@ulaval.ca 
            zdFFBranch = ph.get_zdFF(isosbestic.iloc[:,iBranch],signalData.iloc[:,iBranch], remove = 0)
            f, (ax1) = plt.subplots(1)
            ax1.plot(zdFFBranch,color = 'black')
            ax1.set(title = 'GCaMP zdFF, Branch ' + str(iBranch))        
            zdFF['Branch ' + str(iBranch)] = zdFFBranch            
            
#%% Find connected behavioural data and sync to recording    
    matching = [s for s in behaviourFiles if expID in s]

    if behaviourType == 'Event': 
        behaviour = pd.read_csv(behaviourLocation+matching[0], header = None, error_bad_lines=False)
        linked, figAll = markKeyDown(behaviour, zdFF,splitTimestamp)
        
    elif behaviourType == 'Status':
        behaviour = pd.read_csv(behaviourLocation+matching[0], error_bad_lines=False)
        linked = behStatus(behaviour, zdFF, colsToUse)  
        #%% plot behavioural data and zdFF together         
        figAll, (ax1) = plt.subplots(linked.shape[1])
        for iCol, colname in enumerate(linked.columns):
            ax1[iCol].plot(linked.iloc[:,iCol], color = 'black')
            ax1[iCol].set(title = colname)
    
# save zdFF data and plots and linked behavioural and photometry data 
    figAll.savefig(saveLocation + expID + 'zdFF with behaviour.pdf' )
    linked.to_csv(saveLocation+expID + '_Linked PM and Beh.csv')