# Display Outputs



## Sinusoidal deformation fits
This script compares the phases and amplitudes of the two data fitting methods.

```
SSD_displacement_compare('experiment', 'samples', {'corundum, top','zinc','corundum, bottom'} )
```

% this will produce the remaining image files in the outputs directory. 
% these figures are equivalent to the figures in the manuscript describing this method. 

* Fits_experiment.pdf
* PhasePhaseexperiment.pdf
* AmplitudeAmplitudeexperiment.pdf
* PeriodPeriodexperiment.pdf
* PhaseLagPhaseLagexperiment.pdf




## Position and rate fits

This script compares the positions and position change rates of the two data fitting methods.

```
SSD_strainrate_compare('experiment', 'samples', {'corundum, top','zinc','corundum, bottom'} )
```

which produces: 

* StrainRateFits_experiment.pdf



### Previous step
[5. Calculate Phases and Positions](./05-calculate-phases-positions.md)