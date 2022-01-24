===========================================================================
                Algorithm: EPOD and pressure estimation
                   Author: CHEN Junwei, Stefano DISCETTI, Marco Raiola
               Laboratory: EAP of UC3M
                 Platform: MATLAB
                     Date: 3rd November 2021
                  Version: 1.1.0
                  Contact: junwei.chen@uc3m.es
===========================================================================

DESCRIPTION
The codes run on PIV data of the wake flow behind the wing model in the water tunnel of UC3M. The flow speed is low enough to allow that the original PIV is time-resolved. 
The codes simulate the snapshot PIV (i.e. non-time-resolved) by down-sampling the velocity field, and retrieve the virtual probe data from the PIV data at original sample rate.
The time-resolved velocity fields are reconstructed by EPOD (Extended Proper Orthogonal Decomposition), after that, the pressure fields are integrated using iterative method. 
In addition, the pressure fields are also integrated from snapshot PIV data leveraging Taylor's hypothesis.

PLATFORM
The codes are tested on [1. MATLAB R2021a (Linux 64 bit), on a 4-core tiger lake PC with 64 GB memory][2. MATLAB 2020a (Windows 64 bit), on a 12-core zen 2 PC with 128 GB] (it takes 20GB at maximum). But they can work on recent release of MATLAB.

RUNNING
1. running h52mat.m to convert *.h5 file in PIV_fields to *.mat format and save in OUT_TRPIV, with the same filename of variables.
2. running ReadDatasetTR.m for EPOD reconstruction.
3. running pressure.m for pressure estimation.

note:
1. when having problems in converting HDF5 file, contact us to require MAT files, even original images.
2. switch parfor and for in line 78 of pressure.m to accelerate calculation/unload the computer.

FILES
CalcPOD.m                         MATLAB Function
CreateTestingDataset.m            MATLAB Function
figures                           Folder of result for the frame Wing_006032.h5
h52mat.m                          MATLAB Script
PIV_fields                        Folder of PIV data
pressure.m                        MATLAB Script
ReadDatasetTR.m                   MATLAB Script
ReadInput.m                       MATLAB Function
readme.txt                        readme

EXTENDED READING
S. Discetti, M. Raiola, A. Ianiro, Estimation of time-resolved turbulent fields through correlation of non-time-resolved field measurements and time-resolved point measurements, Experimental Thermal and Fluid Science 93 (2018) 119–130.
