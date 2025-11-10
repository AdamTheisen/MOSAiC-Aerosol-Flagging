# -*- coding: utf-8 -*-
"""
Created on Fri Oct 27 11:35:17 2023

@author: sirna
"""


#%% Creating Concentration Dataframe 

def download_data_nc(path_name, file_name, variable):
    import xarray as xr
    import pandas as pd
    from pathlib import Path
    import glob
    import os
    import netCDF4
    
    master_frame=pd.DataFrame()

    path = path_name
    for filename in glob.glob(os.path.join(path,file_name)):
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
         

    return master_frame,data_1


wind_set= download_data_nc('YOUR FILE PATH', file_name= 'mosaosmetM1.a1.*.*.nc', variable=['varible1', 'variable2'])


#%% Exporting to HDF5 file

wind_set[0].to_hdf('filename.h5', key='wind_set', mode='w')

