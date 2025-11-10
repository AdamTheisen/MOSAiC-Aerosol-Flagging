# -*- coding: utf-8 -*-
"""
Created on Fri Dec  8 16:22:30 2023

@author: sirna

Content:
    downloading data from hdf5 files
    applying purge and quality control
    performing rolling 3 standard deviation/ rolling mean
    done for a window of 100 and 200
"""


#%% Downloading HDF5 files from local folder

'''
HDF5 files were generated using a seperate script that contain the years worth 
of data for each variable and are stored locally
This section downloads the HDF5 file from their local file and converts them to
Pandas data frames
'''


import pandas as pd

purge_set=pd.read_hdf('YOUR .h5 file path')
cpc_set=pd.read_hdf('YOUR .h5 file path')
wind_set= pd.read_hdf('YOUR .h5 file path')

#%% Applying CPC QC

'''
this section takes the particle concentration outputed by the CPC and applies the 
AOS quality control
'''

cpc_set['aos_qc_concentration'] = cpc_set['concentration'] [cpc_set['qc_concentration'] == 0]

#%% Applying Purge

'''
The purge data is a flag that contains if the purge system was on or not 
Data recorded when the purge was on needs to be removed
This section merges CPC frame and the wind frames (MET data) each with the purge frame to remove this flag
Purge has a frequency of 5 seconds
CPC and wind have frequencies of 1 second
'''


# Drops the duplicated time column from each frame
purge_set=purge_set.drop(columns='time')
wind_set=wind_set.drop(columns='time')
cpc_set=cpc_set.drop(columns='time')

# Averages each frame by 5 seconds so they can be merged
cpc_set=cpc_set.resample('5s').mean()
purge_set=purge_set.resample('5s').mean()
wind_set=wind_set.resample('5s').mean()


# Merges wind and CPC with purge set 
cpc_set_purge=pd.merge(cpc_set, purge_set, on='time')
wind_set_purge= pd.merge(wind_set, purge_set, on='time')

'''
Applies the purge flag to particle concentration, wind direction and wind speed
Variables are labeled 'scrubbed' to mean they have been purged and quality controlled
'''

cpc_set_purge['concentration_scrubbed'] = cpc_set_purge['aos_qc_concentration'] [cpc_set_purge['stack_purge_state'] == 0]
wind_set_purge['wind_direction_purged'] = wind_set_purge['wind_direction'] [wind_set_purge['stack_purge_state'] == 0]
wind_set_purge['wind_speed_purged'] = wind_set_purge['wind_speed'] [wind_set_purge['stack_purge_state'] == 0]


#%% Master Frame

'''
A master frame is created to contain all variables used
This requires more merging so all repeated columns are dropped and wind and particle
concentration frames are merged to create the master frame

'''

wind_set_purge=wind_set_purge.drop(columns='stack_purge_state')
master_frame=pd.merge(cpc_set_purge, wind_set_purge, on='time')
master_frame=master_frame.dropna()

#%% ROLLING Statistical Filtering of outliers

'''
Pre filter statistical analysis: mean and median concentration values
'''

mean_concentration_pre= master_frame['concentration_scrubbed'].mean()
pre_purge_len=len(master_frame['concentration'].dropna())

pre_filt_len= len(master_frame['concentration_scrubbed'].dropna())
pre_filt_mean= master_frame['concentration_scrubbed'].mean()
pre_filt_median= master_frame['concentration_scrubbed'].median()



'''
Filter: Rolling 3 std from mean filter, using size 100 window
Rolling data is calculated based on scrubbed concentration for a window size of 100
The values to be used to filter are merged with the master frame to be labeled as outliers
The filtering condition: if greater than outlier removed
'''

rolling_mean_100= master_frame['concentration_scrubbed'].rolling(window=100, center=True).mean()
rolling_std_100=master_frame['concentration_scrubbed'].rolling(window=100, center=True).std()
filtering_trace_100 = rolling_mean_100 + (3*rolling_std_100)

master_frame['outliers_100'] = filtering_trace_100

master_frame['filtered_concentration_100'] = master_frame['concentration_scrubbed'] [master_frame['concentration_scrubbed'] < master_frame['outliers_100']]

'''
Post 100 window size statistics calculated
'''

post_filt_len_100=len(master_frame['filtered_concentration_100'].dropna())
post_filt_mean_100= master_frame['filtered_concentration_100'].mean()
post_filt_median_100= master_frame['filtered_concentration_100'].median()


'''
The same process as outlined above, for a rolling window size of 200
'''

rolling_mean_200= master_frame['concentration_scrubbed'].rolling(window=200, center=True).mean()
rolling_std_200=master_frame['concentration_scrubbed'].rolling(window=200, center=True).std()
filtering_trace_200 = rolling_mean_200 + (3*rolling_std_200)

master_frame['outliers_200'] = filtering_trace_200

master_frame['filtered_concentration_200'] = master_frame['concentration_scrubbed'] [master_frame['concentration_scrubbed'] < master_frame['outliers_200']]

post_filt_len_200=len(master_frame['filtered_concentration_200'].dropna())
post_filt_mean_200= master_frame['filtered_concentration_200'].mean()
post_filt_median_200= master_frame['filtered_concentration_200'].median()



  
# %% Time Series Plot 

'''
the python library Matplotlib and numpy are used to generate a time series
of filtered and unfiltered traces for the rolling window size of 100
'''

import numpy as np
from matplotlib import pyplot as plt

plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["font.size"] = 20


fig, axs= plt.subplots(nrows=3, ncols=1)

fig.supylabel('Concentration [cm$^{-3}$]', fontsize = 30)
fig.tight_layout()

axs[0].set_title('(a) Unfiltered Concentration', fontname="Times New Roman" , fontsize = 20)
axs[0].plot(master_frame.index, master_frame['concentration_scrubbed'], linestyle='-', label='unfiltered', color='black')
axs[0].grid(color='grey', linestyle='--', linewidth=0.5)
axs[0].patch.set_edgecolor('black')  
axs[0].patch.set_linewidth(1)  
axs[0].set_ylim([-1000,11000])
#axs[0].set_yticks(np.arange(0,81000, step=20000))
#axs[0].set_ylabel('Concentration [cm$^{-3}$]', fontsize = 20)
axs[0].set_facecolor('white')


axs[1].set_title('(b) Both', fontname="Times New Roman" , fontsize = 20)
axs[1].plot(master_frame.index, master_frame['concentration_scrubbed'], linestyle='-', label='unfiltered', color='black')
axs[1].plot(master_frame.index, master_frame['filtered_concentration_100'], linestyle='-', label='filtered', color='deepskyblue')
axs[1].grid(color='grey', linestyle='--', linewidth=0.5)
axs[1].patch.set_edgecolor('black')  
axs[1].patch.set_linewidth(1)  
axs[1].set_ylim([-1000,11000])
#axs[1].set_yticks(np.arange(0,81000, step=20000))
#axs[1].set_ylabel('Concentration [cm$^{-3}$]', fontsize = 20)
axs[1].set_facecolor('white')
axs[1].legend(facecolor='white', framealpha=1)

axs[2].set_title('(c) Filtered Concentration', fontname="Times New Roman" , fontsize = 20)
axs[2].plot(master_frame.index, master_frame['filtered_concentration_100'], linestyle='-', label='filtered', color='deepskyblue')
axs[2].grid(color='grey', linestyle='--', linewidth=0.5)
axs[2].patch.set_edgecolor('black')  
axs[2].patch.set_linewidth(1)  
axs[2].set_ylim([-1000,11000])
#axs[1].set_yticks(np.arange(0,81000, step=20000))
#axs[1].set_ylabel('Concentration [cm$^{-3}$]', fontsize = 20)
axs[2].set_facecolor('white')



plt.subplots_adjust(left=0.1,
                    bottom=0.1,
                    right=0.9,
                    top=0.9,
                    wspace=0.4,
                    hspace=0.4)

plt.show()




# %%POLAR PLOT PROCESS

'''
The polar projection is a degree average of all the data avaiable, this function
is used to do the degree binning
'''

def degree_binning(frame, mean_med):
    
    frame=frame.sort_values([ 'wind_direction'])
    frame=frame.dropna()


    frame['wind_direction']=frame['wind_direction'].round(decimals=0)
    frame=frame.astype({'wind_direction':'int'})

    frame=frame.reset_index('time')
    
    if mean_med == 'mean':
        frame=frame.groupby([ 'wind_direction']).mean()
    elif mean_med == 'median':
        frame=frame.groupby([ 'wind_direction']).median()
    
    frame=frame.reset_index([ 'wind_direction'])
    
    return frame

'''
Unfiltered and filtered dataframes are made sepearately to be passed through 
the degree binning function
The function returns a frame that will be used for polar plotting
'''

polar_frame_scrub=pd.DataFrame()
polar_frame_scrub['concentration_scrubbed'] = master_frame['concentration_scrubbed']
polar_frame_scrub['wind_direction'] = master_frame['wind_direction_purged']

polar_frame_filt_100 = pd.DataFrame()
polar_frame_filt_100['filtered_concentration_100'] = master_frame['filtered_concentration_100']
polar_frame_filt_100['wind_direction'] = master_frame['wind_direction_purged']

polar_frame_filt_200 = pd.DataFrame()
polar_frame_filt_200 ['filtered_concentration_200'] = master_frame['filtered_concentration_200']
polar_frame_filt_200 ['wind_direction'] = master_frame['wind_direction_purged']

polar_frame_filt_100 = degree_binning(polar_frame_filt_100, 'mean')
polar_frame_filt_200 = degree_binning(polar_frame_filt_200, 'mean')
polar_frame_scrub =degree_binning(polar_frame_scrub, 'mean')


#%% Polar Plot

'''
Creates a polar projection of the filtered and unfiltered trace
'''

fig=plt.figure()
ax= fig.add_subplot(projection= 'polar')

ax.set_theta_zero_location('N')
ax.set_theta_direction(-1)

theta_scrub=np.deg2rad(polar_frame_scrub['wind_direction'])
plt.polar(theta_scrub, polar_frame_scrub['concentration_scrubbed'] , label="unfiltered", color='black')

theta_filt_100=np.deg2rad(polar_frame_filt_100['wind_direction'])
plt.polar(theta_filt_100, polar_frame_filt_100['filtered_concentration_100'] , label= "filtered", color='deepskyblue')

# theta_filt_200=np.deg2rad(polar_frame_filt_200['wind_direction'])
# plt.polar(theta_filt_200, polar_frame_filt_200['filtered_concentration_200'], label ='filtered, window size = 200', color='green')

leg=ax.legend(bbox_to_anchor=(1.1, 1.05), frameon=True, loc='upper left')
leg.get_frame().set_edgecolor('black')
leg.get_frame().set_facecolor('none')

ax.grid(color='grey', linestyle='--', linewidth=0.5)
ax.patch.set_edgecolor('black')  
ax.patch.set_linewidth(1)  
ax.set_facecolor('white')

plt.show()



#%% Percent Filtered

'''
Printing of Statistical changes in data for both 100 window and 200 window sized
rolling statstical filter
'''

print('ROLLING WINDOW 100')

print('Length of data before filtering:', pre_filt_len)
print('Length of data after filtering: ', post_filt_len_100)
print('Average particle concentration before filtering: ', pre_filt_mean)
print('Average particle concentration after filtering: ', post_filt_mean_100)

percent_lost_100= (1- (post_filt_len_100/pre_filt_len))*100
print('Percentage of Data lost: ', percent_lost_100)

delta_avg_100=pre_filt_mean - post_filt_mean_100
print('Change in Average: ', delta_avg_100)

delta_avg_percent_100= (delta_avg_100/pre_filt_mean)*100
print('Percent Change in Average: ', delta_avg_percent_100)

delta_median_100 = pre_filt_median - post_filt_median_100
delta_median_percent_100= (delta_median_100/pre_filt_median)*100

print('Percent Change in Median: ', delta_median_percent_100)

print(' ')
print(' ')

print('ROLLING WINDOW 200')

print('Length of data before filtering:', pre_filt_len)
print('Length of data after filtering: ', post_filt_len_200)
print('Average particle concentration before filtering: ', pre_filt_mean)
print('Average particle concentration after filtering: ', post_filt_mean_200)

percent_lost_200= (1- (post_filt_len_200/pre_filt_len))*100
print('Percentage of Data lost: ', percent_lost_200)

delta_avg_200=pre_filt_mean - post_filt_mean_200
print('Change in Average: ', delta_avg_200)

delta_avg_percent_200= (delta_avg_200/pre_filt_mean)*100
print('Percent Change in Average: ', delta_avg_percent_200)

delta_median_200 = pre_filt_median - post_filt_median_200
delta_median_percent_200= (delta_median_200/pre_filt_median)*100

print('Percent Change in Median: ', delta_median_percent_200)



