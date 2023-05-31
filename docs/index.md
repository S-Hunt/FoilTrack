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

The code requires the 'Image processing' and 'statistics and machine learning' toolbox.
For spped of processing parts of the script are parallelised and require the 'parallel computing' toolbox.
All of these toolboxes are the mathworks versions thereof. 

## Running the Scripts
There are five parts to the running of the displacement scripts. Documentation
can be found below or by following the links.

1.	[Creation of batch Files](./01-list-files.md)
2.	[Setup – AnalysisOptions](./02-analysis.md) and [Setup - Parsing the images](./02b-parse.md)
3.	[Selection of the boxes](./03-make-boxes.md)
4.	[Calculating the Displacements](./04-displacements.md)
5.	[Calculate sinusoidal fits or positions](./05-calculate-phases-positions.md)
6.	[Display outputs](./06-display-fits.md)

The scripts write all the files to Matlab’s current active directory. It is recommended that a
new directory is made to contain the analysis files separately from the experimental files.


## Examples
Alternatively, exectute [RunExample.m](../Example1/RunExample.m) in Example1 directory. A markdown version of the RunExample.m file is found [here](./Example1Markdown.md)

