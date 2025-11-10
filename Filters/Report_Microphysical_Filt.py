# -*- coding: utf-8 -*-
"""
Created on Sat Dec  9 12:47:10 2023

@author: sirna
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

smps_stats_qc= pd.read_hdf('YOUR .h5 file path')

#%% Apply CPC QC

'''
this section takes the particle concentration outputed by the CPC and applies the 
AOS quality control
'''

cpc_set['aos_qc_concentration'] = cpc_set['concentration'] [cpc_set['qc_concentration'] == 0]


#%% Applying Purge

'''
The purge data is a flag that contains if the purge system was on or not 
Data recorded when the purge was on needs to be removed
This section merges CPC frame, SMPS frame and the wind frame (MET data) each with the purge frame to remove this flag
'''


'''
Purging CPC
cpc 1s freq
purge 5s freq
'''

cpc_set=cpc_set.resample('5s').mean()
purge_set=purge_set.resample('5s').mean()

cpc_set_purge=pd.merge(cpc_set, purge_set, on='time')

cpc_set_purge['concentration_scrubbed'] = cpc_set_purge['aos_qc_concentration'] [cpc_set_purge['stack_purge_state'] == 0]
pre_merge_mean= cpc_set_purge['concentration_scrubbed'].dropna().mean()

'''
Purging MET
met 1s freq
purge 5s freq
'''

wind_set=wind_set.resample('5s').mean()
#wind_set=wind_set.drop(columns='time')

wind_set_purge= pd.merge(wind_set, purge_set, on='time')
wind_set_purge ['wind_direction_purge'] = wind_set_purge['wind_direction'] [ wind_set_purge['stack_purge_state']  == 0 ]
wind_set_purge['wind_speed_purge'] = wind_set_purge['wind_speed']  [ wind_set_purge['stack_purge_state'] ==0]

'''
Purging SMPS
purge 5s freq
smps 5 min freq
'''

smps_stats_qc=smps_stats_qc.drop(columns='time')
smps_stats_qc=smps_stats_qc.resample('5min').mean()

purge_set=purge_set.resample('5min').mean()
#Done to match the the purge set to the smps set
purge_set.drop(purge_set.loc[purge_set.index < '2019-10-11 00:45:00'].index, inplace=True)
purge_set.drop(purge_set.loc[purge_set.index > '2020-10-01 13:55:00'].index, inplace=True)

smps_stats_purge=pd.merge(smps_stats_qc, purge_set, on='time')

smps_stats_purge['geometric_mean_purged'] = smps_stats_purge['geometric_mean'] [smps_stats_purge['stack_purge_state'] == 0]
smps_stats_purge['mean_purged']= smps_stats_purge['mean'] [smps_stats_purge['stack_purge_state'] == 0]

#%% Master Frame

'''
A master frame is created to contain all variables used
This requires more merging so a repeated columns are dropped and variables 
necessary for filtering are extracted and transfered to the master frame

'''


smps_stats_purge= smps_stats_purge.drop(columns='stack_purge_state')
cpc_set_purge=cpc_set_purge.resample('5min').mean()

cpc_set_purge.drop(cpc_set_purge.loc[cpc_set_purge.index < '2019-10-11 00:45:00'].index, inplace=True)
cpc_set_purge.drop(cpc_set_purge.loc[cpc_set_purge.index > '2020-10-01 13:55:00'].index, inplace=True)


wind_set_purge = wind_set_purge.drop(columns='stack_purge_state')
wind_set_purge = wind_set_purge.resample('5min').mean()

master_frame= pd.merge(smps_stats_purge,cpc_set_purge, on='time')
master_frame= master_frame.dropna()

pre_filt_len= len(master_frame['concentration_scrubbed'].dropna())
pre_filt_mean= master_frame['concentration_scrubbed'].mean()
pre_filt_median= master_frame['concentration_scrubbed'].median()

#%% Microphysical Filtering

''' 
A filter is applied to the master frame when the number concentration is greater than 3 standard deviations away and the geometric mean diameter
is smaller than 40 nm
'''


import numpy as np

mean_concentration=master_frame['concentration_scrubbed'].mean()
std_concentration = master_frame['concentration_scrubbed'].std()
three_std_away = mean_concentration + (3*std_concentration)

master_frame['filter_flag'] = np.where( ( (master_frame['concentration_scrubbed'] > three_std_away ) &   (master_frame['geometric_mean'] < 40)  ), 1,0 )
master_frame['filtered_concentration'] = master_frame['concentration_scrubbed'] [master_frame['filter_flag'] == 0]



post_filt_len = len(master_frame['filtered_concentration'].dropna())
post_filt_mean = master_frame['filtered_concentration'].mean()
post_filt_median = master_frame['filtered_concentration'].median()

#%% Time Series plot
'''
the python library Matplotlib and numpy are used to generate a time series
of filtered and unfiltered traces
'''

from matplotlib import pyplot as plt

plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["font.size"] = 20

fig,axs = plt.subplots(nrows=3, ncols=1)

fig.supylabel('Concentration [cm$^{-3}$]', fontsize = 30)
fig.tight_layout()

axs[0].set_title('(a) Unfiltered Concentration', fontname="Times New Roman" , fontsize = 25)
axs[0].plot(master_frame.index, master_frame['concentration_scrubbed'], linestyle='-', label='unfiltered', color='black')
axs[0].grid(color='grey', linestyle='--', linewidth=0.5)
axs[0].patch.set_edgecolor('black')  
axs[0].patch.set_linewidth(1)  
axs[0].set_ylim([-1000,11000])
#axs[0].set_yticks(np.arange(0,81000, step=20000))
#axs[0].set_ylabel('Concentration [cm$^{-3}$]', fontsize = 20)

master_frame['removed'] = master_frame['concentration_scrubbed'] [master_frame['filter_flag'] == 1]


axs[1].set_title('(b) Both', fontname="Times New Roman" , fontsize = 25)
axs[1].plot(master_frame.index, master_frame['concentration_scrubbed'], linestyle='-', label='unfiltered', color='black')
axs[1].plot(master_frame.index, master_frame['filtered_concentration'], linestyle='-', label='filtered', color='deepskyblue')
#axs[1].plot(master_frame.index, master_frame['removed'], linestyle='-', label='unfiltered', color='red')
axs[1].grid(color='grey', linestyle='--', linewidth=0.5)
axs[1].patch.set_edgecolor('black')  
axs[1].patch.set_linewidth(1)  
axs[1].set_ylim([-1000,11000])
#axs[1].set_yticks(np.arange(0,81000, step=20000))
#axs[1].set_ylabel('Concentration [cm$^{-3}$]', fontsize = 20)
axs[1].set_facecolor('white')
axs[1].legend(facecolor='white', framealpha=1)


axs[2].set_title('(c) Filtered Concentration', fontname="Times New Roman" , fontsize =25)
axs[2].plot(master_frame.index, master_frame['filtered_concentration'], linestyle='-', label='filtered', color='deepskyblue')
axs[2].grid(color='grey', linestyle='--', linewidth=0.5)
axs[2].patch.set_edgecolor('black')  
axs[2].patch.set_linewidth(1)  
axs[2].set_ylim([-1000,11000])
#axs[1].set_yticks(np.arange(0,81000, step=20000))
#axs[1].set_ylabel('Concentration [cm$^{-3}$]', fontsize = 20)



plt.subplots_adjust(left=0.1,
                    bottom=0.1,
                    right=0.9,
                    top=0.9,
                    wspace=0.4,
                    hspace=0.4)

plt.show()



#%% Polar Plot Processing


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

wind_set_purge.drop(wind_set_purge.loc[wind_set_purge.index < '2019-10-11 00:45:00'].index, inplace=True)
wind_set_purge.drop(wind_set_purge.loc[wind_set_purge.index > '2020-10-01 13:55:00'].index, inplace=True)




master_frame_wind= pd.merge(master_frame, wind_set_purge, on='time')

polar_scrub_frame= pd.DataFrame()
polar_scrub_frame['wind_direction'] = master_frame_wind['wind_direction_purge']
polar_scrub_frame['concentration_scrubbed'] =master_frame_wind['concentration_scrubbed']

polar_filt_frame = pd.DataFrame()
polar_filt_frame['wind_direction'] = master_frame_wind['wind_direction_purge']
polar_filt_frame['filtered_concentration'] = master_frame_wind['filtered_concentration']

polar_scrub_frame = degree_binning(polar_scrub_frame, 'mean')
polar_filt_frame = degree_binning(polar_filt_frame, 'mean')

#%% Polar Plot

fig=plt.figure()
ax= fig.add_subplot(projection= 'polar')

ax.set_theta_zero_location('N')
ax.set_theta_direction(-1)

theta_scrub= np.deg2rad(polar_scrub_frame['wind_direction'])
plt.polar(theta_scrub, polar_scrub_frame['concentration_scrubbed'], label= 'unfiltered', color='black')

theta_filt = np.deg2rad(polar_filt_frame['wind_direction'])
plt.polar(theta_filt, polar_filt_frame['filtered_concentration'], label= 'filtered', color='deepskyblue')

leg=ax.legend(frameon=True)
leg.get_frame().set_edgecolor('black')
leg.get_frame().set_facecolor('none')

ax.grid(color='grey', linestyle='--', linewidth=0.5)
ax.patch.set_edgecolor('black')  
ax.patch.set_linewidth(1)  
ax.set_facecolor('white')

plt.show()



#%% Printing Stats

print('Microphysical Filtering')
print('  ')

print('Length of data before filtering:', pre_filt_len)
print('Lenth of data after filtering: ', post_filt_len)
print('Average particle concentration before filtering: ', pre_filt_mean)
print('Average particle concentration after filtering: ', post_filt_mean)

percent_lost= (1- (post_filt_len/pre_filt_len))*100
print('Percent of data lost: ', percent_lost)

delta_avg= pre_filt_mean - post_filt_mean
print('Change in average: ', delta_avg)

avg_per_change= (delta_avg/pre_filt_mean)*100
print('Percent change in average: ', avg_per_change)

delta_median= pre_filt_median - post_filt_median
print('Change in median: ', delta_median)

median_per_change = (delta_median/ pre_filt_median)*100
print('Percent change in median: ', median_per_change)

