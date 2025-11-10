Filters Read me:

In this folder are 6 python scripts which will be outlined as follows:

CPCf_Stats_General

This code is structured differently than all subsequent codes: It is built to interface with the user and everything is a user input.
This code can be run as is with no changes and a simple statistical filter will be applied to the data stream of choice with user determined standard deviation cut off. 
Sample inputs are provided in the code for each prompt asked.



The following codes are not soft coded. Some editing in terms of where data is saved will need to be done for each of the following scripts. What needs editing is clear in the file.


Report_Meteorological_Filt

1. This code needs an input of an hdf5 file for both number concentration data, meteorological data. User must edit this file location.
2. If necessary a file indicating when the purge system was on could also be inputted (this is specific to mobile AOS systems)
3. These hdf5 files can be over whatever time frame desired, the code is not time specific
4. with those input files, the code will apply the number concentration quality control, a purge if necessary and create a master data frame of both particle concentration and meteorological data
5. It will then apply the filtering scheme to the data
	Data is removed in this filter when the data point has been recorded within a wind direction of 110-225 and a wind speed greater than 5 m/s
6. The code will then provide statistical data for the filtering scheme and can plot time series and polar projections


Report_Microphysical_Filt

1. This code needs an input of an hdf5 file for both number concentration data and microphysical data (particle size). User must edit this file location.
	(Meteorological data is not necessary to apply the filter however the code does use it to create polar projections)
2. If necessary a file indicating when the purge system was on could also be inputted (this is specific to mobile AOS systems)
3. These hdf5 files can be over whatever time frame desired, the code is not time specific
4. with those input files, the code will apply the number concentration quality control, a purge if necessary and create a master data frame of both particle concentration and meteorological data
5. It will then apply the filtering scheme to the data
	Data is removed in this filter when the number concentration is greater than 3 standard deviations away from the mean of the time series and the geometric mean diameter is smaller than 40 nm
6. The code will then provide statistical data for the filtering scheme and can plot time series and polar projections

Report_Statistical_3std_filt

1. This code needs an input of an hdf5 file for number concentration data
	(Only concentration is needed in this code to apply a filter, however if you want polar projection plots a meteorological dataset is necessary). User must edit this file location.
2. If necessary a file indicating when the purge system was on could also be inputted (this is specific to mobile AOS systems)
3. These hdf5 files can be over whatever time frame desired, the code is not time specific
4. with those input files, the code will apply the number concentration quality control, a purge if necessary and create a master data frame of both particle concentration and meteorological data
5. It will then apply the filtering scheme to the data
	Data is removed in this filter when the number concentration is greater than 3 standard deviations away from the mean of the time series
6. The code will then provide statistical data for the filtering scheme and can plot time series and polar projections


NOTE: The following files are also varying types of statistical filters. They therefore follow the same structure as Report_Statistical_3std_filt, their standard deviation cut off is just slightly different for each one
 
Report_Statistical_3std_Rolling_filt: This code using a rolling window to calculate mean and standard deviation so the cut off point is always changing along the time series

Report_Statistical_4std_filt: In this code, data is removed when the number concentration is greater than 4 standard deviations away from the mean of the time series

