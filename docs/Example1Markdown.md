# Example 1 - files and data processing


## Files.
### Data Files 
The images files are located in the Example1/data directory. 
These images are a subset of the 400 images extracted from the original netcdf file and saved as tiffs. The time stamps in the orginal file were arbitrary and so are here. 
The n images here show the data processing without the need for all the volumes of data. 

### Output files.
The output files from running the example are in directory Example1/ExampleAnalysis.
Each of the files in this directory is produced by the processing pipeline and is referred to in the ‘RunExample.m‘ file and in the text below. 


## Running the Example.

The file ‘RunExample.m‘ contains all the commands to produce the files in the output direcrory. 




Either run:
```
ImageAnalysis('file name', ‘profiles’) 
```
Or make the files manytimes (see ), set the switch in manytimes to ‘profiles’ and run. One file per line in the manytimes file is created by this script: ‘file name_profiles.mat’. 

The first time the script runs it will present the first image of the first data set. Use the mouse pointer to select the top left and bottom right corners of the Region of interest. The script extracts the average intensity by row across the selected box. These profiles are plotted and the data seved in a \*.mat file.

N.B. if you have selected a region of interest previously and want to change it, then before running ```manytimes``` or ```ImageAnalysis(...)``` again you will need to clear the matlab workspace using ```clear all```. This is because the region of interest coordinates are stored as a global. 


## Example: Full list of commands to get profiles from marker foil.

1. Run ```AnalysisOptions```. Experiment name e.g. 'Zn_02', location 'NSLS:X17B2', choose 'reference image' as 'ref'. Ignore the other options. Press OK. 

2. ```makemanytimes('*target directory*')``` where *target_directory* is the directory containing the data files (assuming \*.nc files). For tiff or other individual file types it might be necessary to call ```MakeManyLstFiles``` - see file documentation. This script makes a file called *manytimes.m*.

3. Edit *manytimes.m*. Change the first line from ```process = 'boxes';``` to ```process = 'profile';```. 

4. Run *manytimes*. The first time it is run you will need to select an area to make the profile from. This will make a series of *.mat and *.jpg files. 
