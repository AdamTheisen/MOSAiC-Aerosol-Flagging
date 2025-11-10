# -*- coding: utf-8 -*-
"""
Created on Mon Aug  7 12:53:00 2023

@author: sirna

Content: Generic Statistical Filtering where user can choose number of standard deviations from mean to filter out
"""


#%% Functions


def download_data_nc(path_name, file_name, variable):
    import xarray as xr
    import pandas as pd
    from pathlib import Path
    import glob
    import os
    import netCDF4
    from datetime import datetime
    
    master_frame=pd.DataFrame()

    path = path_name
    general_name= os.path.join(path,file_name)
    data_1=[]
    start = datetime.now()

    for filename in glob.glob(general_name):
        #OPEN FILE
        data_1=xr.open_mfdataset(filename)
        #SORT BY TIME
        variable_dataset= data_1.sortby('time')
        time_2=variable_dataset.variables['time']
        #CREATE DATAFRAME
        variable_frame=variable_dataset[variable].to_dataframe()
        variable_frame['time']=time_2
        #APPEND TO A MASTER FRAME
        master_frame=pd.concat([master_frame, variable_frame])
        
    end = datetime.now()
      
    total_time = (end-start).total_seconds()
    print('Data download took {} seconds.\n'.format(total_time))
        
    # return master_frame, data_1
    return master_frame

def polar_df_prep(data_frame):
    data_frame=data_frame.sort_values('wind_direction')
    data_frame=data_frame.dropna()
    data_frame['wind_direction']= data_frame['wind_direction'].round(decimals=0)
    data_frame=data_frame.astype({'wind_direction':'int'})
    data_frame=data_frame.groupby('wind_direction').mean().reset_index('wind_direction')
        
    return data_frame

#%% Libraries
import xarray as xr
import pandas as pd
from pathlib import Path
from matplotlib import pyplot as plt
import math as math 
import numpy as np

#%% Main Script

#%%% Download data and Dataframe building

'''
Downloading Data using provided function
First input the path on your computer to where the nc files are stored
Then input the generalized name of the files
Next enter in the names of the variables interested in being extracted 
leave only a space to separate the variable names if downloading multiple variables

Sample input format:
    
    C:/Users/sirna/Aerosols 201910-202005/CPCf_ALL
    mosaoscpcfM1.b1.*.*.nc
    concentration qc_concentration
    
    C:/Users/sirna/Aerosols 201910-202005/AOSMET_ ALL
    mosaosmetM1.a1.*.*.nc
    wind_direction wind_speed
    
    C:/Users/sirna/Aerosols 201910-202005/SMPS_ALL
    mosaossmpsM1.b1.*.*.nc
    geometric_mean mean median mode geometric_std

'''

print('Welcome to Generalized Statsical Filtering of Time Series AOS Data.\n')

print('Starting Downloading Data Procedure:')
path=input('Enter path to nc files:')
file=input('Enter general name for nc files (ex. mosaosmetM1.a1.*.*.nc): ')
# variables=input('Enter list of variables to be extracted from files: ')
variables = list(map(str, input("Enter a list of variables to be extracted from files: ").split()))

# data_frame, data_set=download_data_nc(path, file, variables)
data_frame=download_data_nc(path, file, variables)

del path, file, variables

'''
Code will pass variables of interest through an AOS quality control
must input the name of the quality control variable
'''

qc_bool= input('\nDo you wish to remove AOS quality control flagged data? (Y/N)')

if qc_bool== 'Y':
    variable_concern_qc= input('Input name of variable to apply quality control to: ')
    qc_stream= input('Enter name of quality control column: ')
    qc_col_name= 'aos_qc_'
    qc_col_name += variable_concern_qc
    print('Applying AOS Quality Control...')
    data_frame[qc_col_name]= data_frame[variable_concern_qc] [data_frame[qc_stream]==0.0]

'''
will resample the data to any time interval of interest
'''

average_int= None
average_bool= input('\nDo you wish to average the data over an interval? (Y/N)')

if average_bool == 'Y':
    average_int= input('Over what interval should the data be averaged? (ex. 1min)')
    data_frame=data_frame.resample(average_int).mean()

'''
Download a wind direction data set and combine two outputted dataframes
'''

polar_bool=input('\nDo you wish to plot data in polar projections? (Y/N)')

wind_frame=pd.DataFrame()
master_frame= pd.DataFrame()


if polar_bool=='Y':
    polar_path=input('Ener path to wind nc files: ')
    polar_file=input('Enter general name for nc Wind files: ')
    polar_variables=list(map(str, input("Enter a list of variables to be extracted from wind files: ").split()))
    
    wind_frame= download_data_nc(path_name=polar_path, file_name= polar_file,variable= polar_variables)
    wind_frame= wind_frame.resample(average_int).mean()


if polar_bool == 'Y':
    master_frame=pd.merge(data_frame, wind_frame, on="time")
else:
    master_frame=data_frame.copy()


#%%% Statistical Filtering

'''
Does statistical filtering to variable of concern 
Standard deviation above or below mean to be filtered inputted by user
'''

print('\nStart of Statistical Filtering: ')
variable_concern= input('Input name of variable to base filtering on: ')
qc_question=input('Has AOS qc been done on this variable? (Y/N)')

if qc_question=='Y':
    #right now assumes focus is one variable and therefore takes variable name from above aos filtering
    raw_name=variable_concern
    variable_concern=qc_col_name    

no_std= int( input('Enter number of standard deviations away from average you want to be filtered out: ') )
above_below= input('Is this above or below mean? (above/below)')

mean_data= None
filter_point= None
col_name=None
stand_dev=None

#THIS DOES NOT ACTUALLY CALC STANDARD DEVIATION

if above_below== 'above' :
    col_name='above_mean_'
    col_name += variable_concern
    mean_data= master_frame[variable_concern].mean()
    stand_dev=master_frame[variable_concern].std()
    filter_point= mean_data+(stand_dev*no_std)
    master_frame['std_flag']= np.where((master_frame[variable_concern] > filter_point), 1, 0)
elif above_below== 'below':
    col_name='below_mean_'
    col_name += variable_concern
    mean_data= master_frame[variable_concern].mean()
    stand_dev=master_frame[variable_concern].std()
    filter_point= mean_data-(no_std*stand_dev)
    master_frame['std_flag']= np.where((master_frame[variable_concern] < filter_point), 1, 0)

#%% Time Series Plotting 

'''
Creating a Time series plot 
Title and axis names inputted by user, option for log scale
'''

cleaned_stream= 'clean_'
cleaned_stream +=variable_concern

print('\nPlotting Time Series...')
log_bool=input('Do you wish to plot on a log scale? (Y/N)')
y_label=input('Enter Y label: ')
title=input('Enter figure title: ')

master_frame[cleaned_stream]= master_frame[variable_concern] [master_frame['std_flag'] == 0]
master_frame=master_frame.reset_index('time')

time=master_frame['time']
cleaned_data_stream=master_frame[cleaned_stream]
unfiltered_data_stream=master_frame[variable_concern]


fig=plt.figure()
plt.grid()

if qc_question=='Y':
    raw_question=input('Do you wish to also plot raw data stream? (Y/N)')
    if raw_question=='Y':
        raw_stream=master_frame[raw_name]
        plt.scatter(time,raw_stream, s=4, label= 'raw data stream', color='orange')

plt.scatter(time, unfiltered_data_stream, s=4, label= 'unfiltered', color='blue')
plt.scatter(time, cleaned_data_stream, s=4, label= 'filtered', color='red')



if log_bool=='Y':
    plt.yscale('log')

plt.xlabel('Time', fontname="Times New Roman" , fontsize = 15)
plt.ylabel(y_label, fontname="Times New Roman" , fontsize = 15)
plt.title(title, fontname="Times New Roman", fontweight="bold" , fontsize = 20)

leg=fig.legend(frameon=True)
leg.get_frame().set_edgecolor('black')
leg.get_frame().set_facecolor('none')

plt.show()

#%% Polar Plotting

'''
if polar data is availble and downloaded a polar projection binned and avereaged by degree is generated 
'''


if polar_bool== 'Y':
    print('\nPlotting Polar Projection Data...')
    
    polar_raw_df=pd.DataFrame()
    polar_raw_df['wind_direction']=master_frame['wind_direction']
    polar_raw_df['wind_speed']=master_frame['wind_speed']
    polar_raw_df['raw_data']=master_frame[raw_name]
    
    unfiltered='unfiltered_stream'
    unfiltered +=variable_concern
    polar_unfiltered_df=pd.DataFrame()
    polar_unfiltered_df['wind_direction']=master_frame['wind_direction']
    polar_unfiltered_df['wind_speed']=master_frame['wind_speed']
    polar_unfiltered_df[unfiltered]=master_frame[variable_concern]

    polar_filtered_df=pd.DataFrame()
    polar_filtered_df['wind_direction']=master_frame['wind_direction']
    polar_filtered_df['wind_speed']=master_frame['wind_speed']
    polar_filtered_df[cleaned_stream]=master_frame[cleaned_stream]

    polar_raw_df=polar_df_prep(polar_raw_df)
    polar_unfiltered_df=polar_df_prep(polar_unfiltered_df)
    polar_filtered_df=polar_df_prep(polar_filtered_df)
    
    fig=plt.figure()
    ax=fig.add_subplot(projection= 'polar')
    
    polar_title=input('Enter the name of the Polar Plot: ')
    
    fig.suptitle(polar_title, fontname="Times New Roman" , fontweight="bold", fontsize = 20)
    ax.set_theta_zero_location('N')
    ax.set_theta_direction(-1)
    
    r_raw=polar_raw_df['raw_data']
    theta_raw=np.deg2rad(polar_raw_df['wind_direction'])
    plt.polar(theta_raw, r_raw, label='raw data stream', color= 'orange')
    
    r_unfiltered=polar_unfiltered_df[unfiltered]
    theta_unfilt=np.deg2rad(polar_unfiltered_df['wind_direction'])
    plt.polar(theta_unfilt,r_unfiltered, label="unfiltered", color='blue')
    
    r_filt=polar_filtered_df[cleaned_stream]
    theta_filt=np.deg2rad(polar_filtered_df['wind_direction'])
    plt.polar(theta_filt,r_filt, label= "filtered", color='red')
    
    leg=ax.legend(fontsize= '12' ,frameon=True)
    leg.get_frame().set_edgecolor('black')



#%% Filtering Statistics 

'''
Calculating and reporting difference in mean and percentage lost
'''

filtered_total_data=len(master_frame[cleaned_stream].dropna())
filtered_data_mean=master_frame[cleaned_stream].mean()
   
raw_total_data=len(master_frame[raw_name].dropna())
raw_data_mean=master_frame[raw_name].mean()
raw_delta_mean= raw_data_mean-filtered_data_mean
print('Difference between raw and filtered mean: {}.\n'.format(raw_delta_mean))
   
raw_percentage_removed= (1- (filtered_total_data/raw_total_data))*100
print('Percentage of raw data removed {}.\n'.format(raw_percentage_removed))
    
'''
if there was quality control- difference in quality control mean and percentage lost are also calculated
'''
if qc_question == 'Y':
    
    qc_total_data=len(master_frame[qc_col_name].dropna())
    qc_data_mean = master_frame[qc_col_name].mean()
    qc_delta_mean = qc_data_mean - filtered_data_mean
    print('Difference between quality controlled and filtered mean: {}.\n'.format(qc_delta_mean))
    
    qc_percentage_removed= (1- (filtered_total_data/qc_total_data))*100
    print('Percentage of quality controlled data removed {}.\n'.format(qc_percentage_removed))

     



