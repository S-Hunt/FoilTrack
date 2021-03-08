# Radiography Image Analysis

This suite of Matlab scripts has been created for analysing displacement of
features within X-radiographs from experiments at
Synchrotron Large Volume press facilities. The scripts track the position
of selected regions or boxes in the images by finding the
minimum of the Sum Squared Differences (SSD) between pairs of images.

The scripts were written for batch processing of images which come either
as NetCDF images or as a series of tiff or other format image files.

Although not necessary, the batch processing assumes a file name format along the lines of:
```
Expt_Name _ [temp C] _ [load tons] _ [period s] _ [other label] _ sequence_number . ending
```
The data files are processed according to sequence number.

The suite of scripts has been tested under: Ubuntu 14.04 / Matlab 2014a. Some parts of the
scripts require some of the Matlab toolboxes. I am not sure which tool boxes I use and where.
Earlier versions have been tested under Windows 7 / Matlab 2012a.

## Installation
Add the ImageAnalysis folder and all subfolders to the matlab path:
```
addpath(genpath('.../ImageAnalysis'))
```

## Running the Scripts
There are five parts to the running of the displacement scripts
1.	Creation of batch Files
2.	Setup – AnalysisOptions
3.	Parsing the images
4.	Selection of the boxes
5.	Calculating the Displacements.

The scripts write all the files to Matlab’s current active directory. It is recommended that a
new directory is made to contain the analysis files separately from the experimental files.

### 1. Make list files and batch files.

There are two sets of batch files that may be required to process the data. The first is if the data
is a series of single image files (e.g. a series of tiff files) and the second is if there are multiple
data sets with multiple images in each (e.g. a thermal diffusivity or anelasticity experiment).

#### Batch files for series of single images

If the data set is a single series of tiff or bitmap images (e.g. a deformation experiment) a `.lst` file
is used to pass the images into the image processing code. The `.lst` file is saved in the current matlab 
directory. To make the `.lst` file run:
```
MakeSingleLstFile('image directory', opts)
```
Where ‘image directory’ is the directory containing the image files. The optional arguments for this script
are:
* `image format`: Defines the extension on the image files. Without this option the script looks for 
the most common file type in the directory and assumes this is the image format. The recognised image formats are 
'tiff', 'tif', 'bmp' and 'edf'.
* `order`: Sorts the images in ascending (‘ascend’) or descending (‘descend’) time order. The default is ascending.

For multiple series of images in separate directories `MakeManyLstFiles` makes `.lst` files for all the subdirectories in the 
‘image directory’. These files are all saved in the current Matlab directory. The options are the same as `MakeSingleLstFile`
with the additional option: 
* `subdirs`: 0 - only makes list file for image directory; 1 - only makes list file for subdirectories of the target and
2 (default) - makes list files for directory and its subdirectories.

#### Batch files for multiple series of images

For data sets that are multiple netcdf files or multiple sets of single images in subdirectories 
`MakeManyTimes` generates the batch files to run image analysis on all the images. It will make 
the `.lst` files as well if they do not already exist. The syntax for this function is: 
```
MakeManyTimes('data directory', opts)
```
`Image format` and `order` are switches to this script as above. Recognised image formats are 
‘netcdf’, ‘nc’, 'tiff', 'tif', 'bmp' and 'edf'.

The additional argument `experiment type` changes which batch files are produced. The valid 
options are ‘rheology’, ‘anelasticity’ and ‘td’ (thermal diffusivity). The latter two options 
makes batch files for fitting sinusoids to the data and fitting the thermal diffusivity of anelastic
models to the data (see sections 6 and 7). The default is ‘rheology’.


### 2. Analysis Options
This makes Expt_Name.mat which holds all the options for the data analysis. 
Running
```
AnalysisOptions
```
opens the window (illustrated below) which presents the options required by the image analysis scripts. The green highlighted boxes need to be changed for the experiment being analysed, the orange may need to be changed and the red boxes are unlikely to be changed.

![AnalysisOptions Interface](./img/menu.png)

The options are: 

* Header rows: 

Option |  Description
------ | ------
Experiment name     | This is the name of the experiment and the name of the saved options file. It has to be the leading string in the image file names, which is usually the experiment name e.g. San_306 
Experiment location | Pull down menu to select the beamline where the experiment was performed. It determines the image pre-processing (rotation, inversions) required. The 'sides' option (where present) is a 90 degree rotation of the images so that displacements of the sides of the foils can be determined.


* ImageAnalysis options:  

These are the options for the ImageAnalysis script that processes the images:-
|Option |  Description |
|------ | ------ |
|Reference ImageAnalysis	| 'Sequence' uses each image and compares it to the following image, moving the reference image as it goes (i.e. image1-image2, image2-image3,…). 'Reference' compares the middle image in the sequence to all the other images (i.e.  image5-image1, …image5-image10).|
|Search distance			| This is the distance (in pixels) that the code calculates the SSD either side of the reference position. The default for the X17B2 / 6-BM-B beamlines is ±5.|
|Minimum SSD				| The method by which the minimum in the SSD is interpolated. The options are ‘spline’, ‘polynomial’, and ‘gaussian’. ‘Spline’ is the usual for rheology experiments and ‘polynamial’ for sinusoidal data sets. |
|Output data				| What files the ImageAnalysis script saves. The ‘displacements’ is the interpolated minimum in the SSD, whilst ‘SSD’ saves all the calculated SSD in a single file. ‘Both’ saves both file types.|
|Bright spots				| If 'yes', this removes the very bright spots from the images. This is used when the images are all black but for the spots so the data cannot be seen. The option was use for the Princeton camera at the X17B2/X17B2ss beamlines which produced 8 bit data within a 12 or 16 bit file.|
|Mask background			| Tries to recognise the regions outside the foil but within the selected box and replaces them with NaN. This removes the effect of these pixels from the analysis. |



* Phases Options: 
These are only needed for anelastic or thermal diffusivity experiments.

Option |  Description
------ | ------
Length type	|Determines how the script combines multiple foils. 'Single' treats each box individually, 'Pair' aligns boxes with the same horizontal position, with more than 2 foils it combines the two outside foils works inwards in pairs. ‘Anelastic’ assumes one box per foil and combines them in series (i.e. 1-2, 2-3, 3-4…)
Image Scaling	| Size of each pixel in microns.
Discard |	If there length type is pair and there are an odd number of foils which one to discard.
Symmetry type	|How to combine the foils for the pair length type. Is the bottom foil switched round?
Background type| 	How to detrend the data before fitting the sinusoid.
Errors	| Save these or not?


### 3.	Parse the files
Computer file systems can take multiple seconds to locate and open the image files when they are first opened, especially if the files are very large (e.g. netCDF) and stored on external hard drives. Parsing the files merely opens and closes the files so when the subsequent processes are celled the operating system knows where to look for the files. This is most useful during the creation of the boxes/areas of interest in the images when there are many '.lst' or netCDF files to be interacted with.

The parsing is done by running: 
```
ImageAnalysis('file name', ‘parse’)
```
Or if batch files are being used set ‘process’ in manytimes to ‘parse’ and run. 

### 4.	Make Boxes
Either run:
```
ImageAnalysis('file name', ‘boxes’)
```
Or set the switch in manytimes to ‘boxes’ and run. Two files are created by this routine (1) ‘file name_boxes.mat’ and ‘file name_boxes.tiff’. The first is a matlab mat file containing the x,y positions of the selected boxes and the second in an image of these boxes overlaid on the reference X-radiograph.


![Select Boxes Interfaces](./img/setup_box.png)

Running ‘boxes’ option opens and plots the reference image of the left hand side of the screen. This is followed by a series of dialogue boxes which aid the location and selection of the boxes. If there are multiple sets of data the boxes are propagated from one data set to the next via a file called data_prev.mat. If data_prev.mat exists it is assumed the box selection routine optimises the box position and presents the dialogue box for 
1.	**Rotation of the image**
The rotation is anti-clockwise about the centre of the image and is propagated to all images in the data series. The boxes are selected on the rotated image. 
2.	**Select the box positions**
Select the top left and bottom right corner of the desired boxes in the image. When all the boxes of interest are selected press ENTER. If the experiment is a thermal diffusivity experiment select the full width of the foil and this area is cut into smaller boxes in the next step. 
3.	**Box size and location selection**
Dialogue boxes to automatically size and locate the boxes follow.  The first window allows for the cutting of the selected areas into many horizontal sections (for thermal diffusivity experiments) and setting the vertical height of the boxes. It is recommended that the heights of the boxes in the image are a little greater than the shadow of the foils (but this can be adjusted later).


![Select Boxes Interfaces](./img/box_position.png)

The second dialogue box is for the selection of the method for automatically position the boxes over the features in the image. The box is centred over the feature selected and the options in the list have the following meanings:

Option |  Description
------ | ------
Minimum intensity	| Finds the pixel row in the box with the minimum intensity.
Minimum intensity interpolated (poly)	| Uses a polynomial to interpolate between the rows to find the minimum intensity.
Minimum intensity interpolated (spline)	| Uses a spline to interpolate between the rows to find the minimum intensity.
Maximum slope	| Finds the pixel row which has the maximum intensity gradient. For the foils in the top half of the image it looks for the gradient on the top of the foil and in the bottom of the image in the bottom of the foil 
Maximum slope interpolated (spline)	| Uses a polynomial to interpolate between the rows to find the maximum intensity gradient. For the foils in the top half of the image it looks for the gradient on the top of the foil and in the bottom of the image in the bottom of the foil

4.	**Validate the positions**
The script locates the boxes using the setting above, plots the box locations on the image and presents the following dialogue box: 

![Select Boxes Validation Inerface](./img/box_val.png)

If the boxes are in the right place press “yes”, then the box positions are saved and second image figure appears on the right hand side of the screen with the boxes shown. This image is saved. If the boxes are not in the correct position use one of the other options:

Option            |  Description
------            | ------
Yes	              | Finds the pixel row in the box with the minimum intensity.
Move all X,Y	    | Allows the moving of all the boxes in X and Y. Used if the press was moved between data sets.
Move Single Foils	| Allows single foils to be moved manually. 
Select New Boxes	| Returns the script to the image rotation option and allows the selection of a new set of boxes.  
Edit 	            | Allows manual editing of the boxes via the command line. 
No	              | Generates and error and quits the box selection. 

5.	**Move all (X,Y)**
Text 
6.	**An option Position single foil**
Text

### 5.	Calculate the displacements / SSD values. 
Either run:
```
ImageAnalysis('file name', ‘disp’)
```
Or set the switch in manytimes to ‘disp’ and run. The output files are called ‘experiment _position_change.txt’ and ‘experiment_SSD.mat’. 

### 6.	Calculate sinusoidal fits 
run: 
```Phases(…)``` or ```PhasesSSD(…)```

### 7.	Check the sinusoidal fits 
Run: ```PhaseCheck(… options)```
	SSD_dispalcement_compare(options). 

### 8.	Fit models to data 
1.	Thermal diffusivity 
Run: 
```
manytimes3
```
or
```
kappasolve(…)
```
2.	Maxwell model for anelasticity
Run: 
```
MaxwellModelPhases('fileglob')
```
Where fileglob is the search term for the files to be processed e.g. ‘*sine_fits*’
