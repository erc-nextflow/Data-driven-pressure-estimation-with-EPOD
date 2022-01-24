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
The codes run on the simulation of pinball flow. With the snapshot velocity field and virtual probe data generated, the time-resolved velocity fields are then reconstructed by applying EPOD (Extended Proper Orthogonal Decomposition).
After that, the pressure fields are integrated using an iterative method. In addition, the pressure fields are also integrated from snapshot velocity field data leveraging Taylor's hypothesis to compare the result.

PLATFORM
The codes are tested on [1. MATLAB R2021a (Linux 64 bit), on a 4-core Tiger Lake PC with 64 GB memory], [2. MATLAB R2021a (Linux 64 bit), on a 64-core Zen 2 workstation with 256 GB memory] and [3. MATLAB R2020a (Windows 64 bit), on a 12-core Zen 2 PC with 128 GB memory]. They should also work on other recent releases of MATLAB.

RUNNING
1. Converting *.h5 file in DataEstimation and DataInterP to *.mat using h52mat.m.
2. Running EPOD_estimation.m to perform EPOD on the training set saved in DataInterp. (the results are already in DataEstimation)
3. Running EPOD_estimation_pressure.m to estimate the pressure by EPOD, and running pressure_TH.m to estimate the pressure by Taylor's hypothesis. (the results are already in DataEstimation)

FILES
DataEstimation                     Folder of result
DataInterp                         Folder of prepared velocity field
EPOD_estimation.m                  MATLAB Script
EPOD_estimation_pressure.m         MATLAB Script
h52mat.m                           MATLAB Script
pressure_TH.m                      MATLAB Script
readme.txt                         Readme
sample.mp4                         The video telling the procedure of EPOD and pressure integration

ACKNOWLEDGEMENT
We warmly acknowledge M. Morzyński and B. Noack for providing access to the fluidic pinball DNS code, and wish the user of data on Zenodo cite Deng et al. (2018) in the article.

EXTEND READING
N. Deng, B. R. Noack, M. Morzyński, L. R. Pastur, Low-order model for successive bifurcations of the fluidic pinball, Journal of Fluid Mechanics 884
S. Discetti, M. Raiola, A. Ianiro, Estimation of time-resolved turbulent fields through correlation of non-time-resolved field measurements and time-resolved point measurements, Experimental Thermal and Fluid Science 93 (2018) 119–130.
