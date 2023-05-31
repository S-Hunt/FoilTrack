# Calculate sinusoidal fits or positions



## Sinusoidal deformation fits

The two data types (for the prior and refined algorithms) are processed to find the phase and amplitude of a sinusoidal displacement by running: 

```
Phases('experiment_position_change.txt')
``` 
and 

```
PhasesSSD('experiment_SSD.mat')
``` 
respectively. 

Optional arguments for these scripts are listed in the header of the script files. The most significant of which is the optional argument ```'period', n``` where ```n``` is the expected period of the deformation. If the file names follow the structure layed out [here](./index.md) the call  ```'period', 'title'``` can be used. This gets the expected period from the title of the data file. 


For large data sets with multiple data at each period it is possible to check that all the periods are consistent:

```
PhaseCheck('type', 'disp')
```

```
PhaseCheck('type', 'ssd')
```

These scripts will produce a number of files:

* ```experiment_sine_fits.txt``` - phase and ampltuide of the sample length change as defined in the ZN08.mnat file.
* ```experiment_sine_fits_single.txt``` - phase and ampltuide for the fits of each foil separately.
* ```experiment_sine_fits_FITS.mat``` - The values stored as matlab data file.
* ```experiment_SSD_sine_fits.txt``` - The values stored as matlab data file.
* ```experiment_SSD_FITS.mat``` - The phase and amplitude for each region of interest stored as matlab data file.



## Position and rate fits
These scripds calculate the postion and position change-rates for the two data methods. 

I fund it useful to set the length of the window to work over (number of images) and the degree of the polynomials used to fit the data with as varaibles. 

```
window = 7;
degree = 2;
```


%run the script:

```
PositionTime('experiment_position_change.txt',  'window', window, 'overwrite', 'degree', degree)
```

This will produce a number of files:

* experiment_window_PositionChangeRate.txt
* experiment_window_LengthChangeRate.txt

Alternatively for the refined algorithm data: 

```
PositionTimeSSD('experiment_SSD',  'window', window, 'surface coefficients', 6,degree, 'overwrite')
```

Which produces:

* experiment_SSDwindow_positions.txt
* experiment_SSDwindow_lengths.txt
* experiment_SSDwindow_PositionChangeRate.txt
* experiment_SSDwindow_LengthChangeRate.txt




### Previous step
[4. Calculate displacements](./04-displacements.md)
### Next step
[6. Display Outputs](./06-display-fits.md)