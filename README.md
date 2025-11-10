# MOSAiC Aerosol Flagging

To use the filters:

1. Download a data set from Data discovery as a .nc file.
2. Use the HDF5 Generator code to create a condensed version of the .nc files for whatever time scale you are working with.
	Details on using the HDF5 code is in the HDF5 Read Me
3. Once an HDF5 file is generated from step 2, located a desired filter to run from the filters file. Apply that filter to the HDF5 file
	Details on the different types of filters and how to apply them can be found in the Filters Read Me.
4. Each filter will result in a Master data frame that has the filtered data stream there as well as figures plotting time series and polar projection plots
