# FoilTrack

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
There are five parts to the running of the displacement scripts. Documentation
can be found below or by following the links.
1.	[Creation of batch Files](./01-list-files.md)
2.	[Setup – AnalysisOptions](./02-analysis.md)
3.	[Parsing the images](./03-parse.md)
4.	[Selection of the boxes](./04-make-boxes.md)
5.	Calculating the Displacements (see below)
6.	Calculate sinusoidal fits (see below)
7.	Check sinusoidal fits (see below)

The scripts write all the files to Matlab’s current active directory. It is recommended that a
new directory is made to contain the analysis files separately from the experimental files.


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
