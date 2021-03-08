# Parse the files

Computer file systems can take multiple seconds to locate and open the image files when they are first opened, especially if the files are very large (e.g. netCDF) and stored on external hard drives. Parsing the files merely opens and closes the files so when the subsequent processes are celled the operating system knows where to look for the files. This is most useful during the creation of the boxes/areas of interest in the images when there are many '.lst' or netCDF files to be interacted with.

The parsing is done by running: 
```
ImageAnalysis('file name', ‘parse’)
```
Or if batch files are being used set ‘process’ in manytimes to ‘parse’ and run. 
