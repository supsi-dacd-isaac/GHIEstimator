# GHIEstimator
This code estimates the global horizontal irradiance (GHI) from the power measurements of one or more photovoltaic (PV) systems located in the same neighborhood. The method is completely unsupervised and is based on a physical model of a PV plant. It can estimate the nominal power and orientation of multiple PV fields, using only the aggregated power signal from their PV power plant. Moreover, if more than one PV power plant is available, the different signals are reconciled using outliers detection and assessing shading patterns for each PV plant.

## Installation
### pv_lib-toolbox
- Download the toolbox from https://pvpmc.sandia.gov/applications/pv_lib-toolbox/matlab/
- Open Matlab
- Go to FILE->SET PATH
- Push “Add with Subfolders” and select PV_LIB folder and press OK (this will add the PV_LIB Toolbox to your path file)

## Getting started
The `demo.m` file runs an example using 10 minutes sampling data from four PV power plants in the same neighborhood. 

The folder data contains the PV and temperature signals, as well as GHI measurements from a local pyranometer and from the free satellite based MACC-RAD service. The last two signals are only used for evaluation of the estimated GHI.

The folder src contains the main function `get_ghi.m` and the other scripts for PV model identification and clear sky period detection. 

## Acknowledgements
The authors would like to thank CTI - Commission for Technology and Innovation (CH), and SCCER-FURIES - Swiss Competence Center for Energy Research - Future Swiss Electrical Infrastructure, for their financial and technical support to this research work.
