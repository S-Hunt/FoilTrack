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
(High Siera) / MatLab 2017a. It is possible (but unlikely) that the scripts do not work elsewhere.

## Background and literature

FoilTrack achieves sub-pixes resolution tracking of foil markers in X-radiographs making use of the image processing
algorithms of Pratt (1991) and Trucco and Verri (1998). The initial implementation is described in Li et al (2003) and
that approach has been refined in several steps (e.g. Hunt et al. 2011, 2012). A description of the most recent
enhancments along with details of the approach and history can be found in Hunt et al. (2021). Users of FoilTrack are
requested to cite Hunt et al. (2021) alongside other relevent publications when describing the use of the software. Further
details can be found in the CITATION file shipped with the software.

## Licence and bundled software

FoilTrack is copyright of the authors. It is released under the terms of version 3 of the GNU General Public 
License (GPLv3). This means you are free to run, study, share, and modify the software but that any 
derivative work must be distributed under the same or equivalent license terms. For full details 
please refer to the [LICENSE](./LICENSE) file distributed with FoilTrack.

For ease of installation and avoid problems with versioning we distribute a handful of files and toolboxes 
alongside FoilTrack. These were all released by their authors under the terms of copyleft of permissive licences that
are compatible with GPLv3. Specifically we distribute the following:
* The [Splinefit toolbox](https://uk.mathworks.com/matlabcentral/fileexchange/71225-splinefit)
created by Jonas Lundgren and released under the permissive license found at FoilTrack/splinefit/license.txt
* The esrf-pmedf toolbox created by Petr Mikulik and released under version 2 of the GPL
* Jacobianest.m from John D'Errico's [DERIVEST suite](https://uk.mathworks.com/matlabcentral/fileexchange/13490-adaptive-robust-numerical-differentiation) and released under the permissive license
* GetCases.m from Rody P.S. Oldenhuis' [FEX-getCases](https://github.com/rodyo/FEX-getCases) software relased under the MIT licnese.
* Functions from [Jan Gläscher](https://uk.mathworks.com/matlabcentral/fileexchange/6837-nan-suite), 
[Peder Axensten](https://uk.mathworks.com/matlabcentral/fileexchange/16177-chi2test),
[Tom O'Haver](https://uk.mathworks.com/matlabcentral/fileexchange/23611-peakfit-m),
[Brad Ridder](https://uk.mathworks.com/matlabcentral/fileexchange/17901-propagation-of-uncertainty), and
[S. Zhang & J. Jin](https://people.sc.fsu.edu/~jburkardt/f77_src/special_functions/special_functions.html)
as documented within the code.

We are greatful for the efforts of these authors.


## Installation and use

FoilTrack should be installed and used as described in the [documentation](./docs/index.md). Please report bugs and
seek support using the [issue tracker on GitHub](https://github.com/ExperimentalMineralPhysics/FoilTrack/issues).
