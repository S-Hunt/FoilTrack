# Make and extract intensity profiles from the images

Either run:
```
ImageAnalysis('file name', ‘profiles’) 
```
Or make the files manytimes (see ), set the switch in manytimes to ‘profiles’ and run. One file per line in the manytimes file is created by this script: ‘file name_profiles.mat’. 

The first time the script runs it will present the first image of the first data set. Use the mouse pointer to select the top left and bottom right corners of the Region of interest. The script extracts the average intensity by row across the selected box. These profiles are plotted and the data seved in a \*.mat file.

N.B. if you have selected a region of interest previously and want to change it, then before running ```manytimes``` or ```ImageAnalysis(...)``` again you will need to clear the matlab workspace using ```clear all```. This is because the region of interest coordinates are stored as a global. 
