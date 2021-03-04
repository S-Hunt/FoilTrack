# FoilTrack

FoilTrack is a suite of Matlab scripts that have been created for analysing displacement of features within
X-radiographs acquired during Large Volume Press (D-DIA, DT-Cup, etc.) experiments at Synchrotron facilities.
The scripts track the position of selected regions or boxes in the images by finding the minimum of the Sum
Squared Differences (SSD) between pairs of images in a series. Tracking these regions allows changes in sample
lengths to be determined and this ability has been used to enable exererimental determination of large strain
rheology, anelasticity, and thermal conductivity. The scripts are able to process most common file formats and
can work in user-supervised or automated batch processing mode.

Documentation, such as it is, for the scripts is provided in the [documentation.md](./documentation.md)
markdown file. Note that (1) This is reseach software not written for public consumption. The many bugs are
all mine. If you find any please let me know. (2) Some parts of the scripts require Matlab toolboxes. See
below for some further information. (3) The documentation is far from complete. It was written as aide-mémoire
for when I use the scripts. (4) The scripts have been tested under Ubuntu 14.04 / Matlab 2014a and MacOS 10.13
(High Siera) / MatLab 2017a. It is possible that the scripts do not work elsewhere.

## Background and citations

Text on where this comes from, and what to cite goes here.

## Licence

GNU3 (see the LICENCE file). Other tools goes here

## Installation

Works as a toolbox...

## Use

Short version
