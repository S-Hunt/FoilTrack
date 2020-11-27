# FoilTrack

This suite of Matlab scripts has been created for analysing displacement of features within X-radiographs acquired during Large Volume Press (D-DIA, DT-Cup, etc.) experiments at Synchrotron facilities. The scripts track the position of selected regions or boxes in the images by finding the minimum of the Sum Squared Differences (SSD) between pairs of images in a series.

The scripts were written for batch processing of images which come either as NetCDF images or as a series of tiff or other format image files.

Documentation, such as it is, for the scripts are in the File ImageAnalysisSuiteDocumentation.docx. 

The suite of scripts has been tested under:
	Ubuntu 14.04 / Matlab 2014a
  MacOS 10.13 (High Siera) / MatLab 2017a
Some parts of the scripts require Matlab toolboxes. However when dumping the package to GitHub I am not sure which toolboxes are used and where. 


Simon Hunt, November 2020.
