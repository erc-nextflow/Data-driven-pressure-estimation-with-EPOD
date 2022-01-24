===========================================================================
                Algorithm: EPOD and pressure estimation
                   Author: CHEN Junwei, Stefano DISCETTI, Marco RAIOLA
               Laboratory: EAP of UC3M
                 Platform: MATLAB
                     Date: 4th November 2021
                  Version: 1.1.0
                  Contact: junwei.chen@uc3m.es
===========================================================================

DESCRIPTION
The codes run on the simulation of channel flow from John Hopkins Turbulent Database (JHTDB). With the snapshot velocity field and virtual probe data downloaded, the time-resolved velocity fields are then reconstructed by applying EPOD (Extended Proper Orthogonal Decomposition), after that, the pressure fields are integrated using iterative method. In addition, the pressure fields are also integrated from snapshot velocity field data leveraging Taylor's hypothesis to compare the result.
The domain is h x h, and extracted in different positions in the spanwise and streamwise directions of the channel exploiting statistical homogeneity. There are 20 virtual probes for the flow field, the first 10 are on the trailing edge of the domain with equal distance, and the last 10 are on the wall.

PLATFORM
The codes are tested on [1. MATLAB R2021a (Linux 64 bit), on a 4-core tiger lake PC with 64 GB memory][2. MATLAB R2020a (Windows 64 bit), on a 12-core zen 2 PC with 128 GB memory]. But they can work on recent release of MATLAB.

RUNNING
1. Running h52mat.m to change data in Fields, Fields_Testing, Probes, Probes_Testing from *.h5 to *.mat format.
2. running EPOD_from_probes.m to let the algorithm learn from the flow field, yielding EPOD_Save.mat (which can also be converted from EPOD_save.h5 by running h52mat.m).
3. Running EPOD_from_probes_pressure.m to estimate the pressure by EPOD, and running EPOD_from_probes_TH_pressure.m to estimate the pressure by Taylor's hypothesis.

FILES
Fields                             Folder of training data
Fields_Testing                     Folder of testing data
Probes                             Folder of training data
Probes_Testing                     Folder of testing data
EPOD_from_probes.m                 MATLAB Script
EPOD_from_probes_pressure.m        MATLAB Script
EPOD_from_probes_TH_pressure.m     MATLAB Script
EPOD_save.h5                       Data
h52mat.m                           MATLAB Script
PIterSolver.m                      MATLAB Function
readme.txt                         Readme
sample.avi                         Sample video

ACKNOWLEDGEMENT
We warmly acknowledge John Hopkins Turbulent Database , and wish the user of data on Zenodo cite their articles.

EXTEND READING
S. Discetti, M. Raiola, A. Ianiro, Estimation of time-resolved turbulent fields through correlation of non-time-resolved field measurements and time-resolved point measurements, Experimental Thermal and Fluid Science 93 (2018) 119–130.
