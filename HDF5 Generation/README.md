# HDF5 Gen Read Me:

A compressed file was generated to be the input for all filters built to avoid the carrying of excess information and to streamline the coding process

The python file named "HDF5 File Generation" is a sample code on how to generate a HDF5 file.
The code requires you to input a file path to where your .nc files are stored, the name of the files to be downloaded and the variables desired to be saved.

The result of this code is an HDF5 file with the variables desired and their time stamp.

A sample HDF5 file for purge, cpc and meteorological data has been provided so that filters in the filters folder can be ran. The data is sourced from the MOSAIC expedition.
