# Make and extract intensity profiles from the images

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