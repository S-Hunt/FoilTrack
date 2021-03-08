# Make list files and batch files.

There are two sets of batch files that may be required to process the data. The first is if the data
is a series of single image files (e.g. a series of tiff files) and the second is if there are multiple
data sets with multiple images in each (e.g. a thermal diffusivity or anelasticity experiment).

## Batch files for series of single images

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

## Batch files for multiple series of images

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
