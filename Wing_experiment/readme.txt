===========================================================================
                Algorithm: EPOD and pressure estimation
                   Author: CHEN Junwei, Stefano DISCETTI, Marco RAIOLA
               Laboratory: EAP of UC3M
                 Platform: MATLAB
                     Date: 3rd November 2021
                  Version: 1.1.3
                  Contact: junwei.chen@uc3m.es
===========================================================================

DESCRIPTION
The codes run on PIV data of the wake flow behind the wing model in the water tunnel of UC3M. The flow speed is low enough to allow that the original PIV is time-resolved. 
The codes simulate the snapshot PIV (i.e. non-time-resolved) by down-sampling the velocity field, and retrieve the virtual probe data from the PIV data at original sample rate.
The time-resolved velocity fields are reconstructed by EPOD (Extended Proper Orthogonal Decomposition), after that, the pressure fields are integrated using iterative method. 
In addition, the pressure fields are also integrated from snapshot PIV data leveraging Taylor's hypothesis.

note:
When having problems in converting HDF5 file, contact us to require MAT files, even original images.

FILES
CalcPOD.m                         MATLAB Function
CreateTestingDataset.m            MATLAB Function
figures                           Folder of result for the frame Wing_006032.h5
h52mat.m                          MATLAB Script
PIV_fields                        Folder of PIV data (in Zenodo)
pressure.m                        MATLAB Script
ReadDatasetTR.m                   MATLAB Script
ReadInput.m                       MATLAB Function
readme.txt                        readme

STRUCTURE OF FILES
in PIV_fields/Grid_Wing.h5
X: 70x110 streamwise position of grid points
Y: 70x110 crosswise position of grid points
in PIV_fields/Wing_[0-9]{6}.mat
U: 70x110 streamwise velocity of flow field in each frame
V: 70x110 crosswise velocity of flow field in each frame

GITHUB SITE
The codes to process the data are available on the GitHub site:
https://github.com/erc-nextflow/Data-driven-pressure-estimation-with-EPOD/tree/main/Wing_experiment

ZENODO SITE
The flow field data will be published on Zenodo, with the link
https://doi.org/10.5281/zenodo.5913979

EXTENDED READING
Discetti, S., Raiola, M., & Ianiro, A. (2018). Estimation of time-resolved turbulent fields through correlation of non-time-resolved field measurements and time-resolved point measurements. Experimental Thermal and Fluid Science, 93, 119-130. doi: 10.1016/j.expthermflusci.2017.12.011
