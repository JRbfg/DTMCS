DT-MCS V1.0 is a free R software to analyze the geometry and the dynamics of bedforms as well as to quantify the uncertainty of the results related to the setting of input parameters. 
It contains the following applications:

1_WLT - Wavelet analysis and detrending 
2_ZC - Zerocrossing analysis 
3a_CC - Cross correlation analysis
3b_CA - Centroid analysis
4_POST - Postprocessing

Step 1 is based on the software Bedforms-ATM (Guitierrez, 2018). Steps 2 & 3a are based on RhenoBT (Frings et al., 2012).

The development of DT-MCS was part of the so-called MAhyD project (morphodynamic analyses using hydroacoustic data), 
which was carried out at the Federal Institute of Hydrology and which was funded by the Federal Ministry for Digital and Transport (BMDV).

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  

....................
:   INSTRUCTIONS   : 
....................

0_DATA:
All BEPs must be stored in this folder. A subfolder needs to be created for every new project.
If repeated measurements are available, a subfolder needs to be created within the project folder for each measurement. BEPs of different profile tracks but of the same measurement campaign are stored in the same folder.
If bedform migration is to be calculated, additionally, a csv-file containing information about all considered measurement pairs has to be provided.

The BEP-files must be in txt-format, have no row and column names, and have blank spaces as separators. 
They must have the following columns: x, z, km, ID
x is the distance in meters in measuring direction.
z is the bed elevation in meters.
km corresponds to the kilometer of the river.
ID is a unique ID for the measured profile track (integer value).

The measurements-file must be named measurements.csv and must be placed within the project folder. 
It must have no row and column names, and semicolons as separators.
It consists of 3 columns time_1, time_2, dt.
time_1: description of first measurement (e.g. date and time)
time_2: description of second measurement (e.g. date and time)
dt: time difference between both measurements in hours

The folder contains two example projects based on an MBES-dataset of the Lower Rhine in Germany, which was chosen as a testcase to demonstrate the tool.  
The study area is located close to Emmerich at Rhine-kilometer 860.0 to 860.5 and covers 500 meters in length and about 200 meters in width. 
The measurements were taken on consecutive days at medium to high flow conditions with discharges ranging from 4000 to 4200 m³/s. The total area covered with bedforms was measured four times with intervals between three and 24 hours. 
Additional measurements were taken along a single profile track in the center of the bedform field at shorter intervals, allowing a more detailed analysis of bedform migration and bedload transport. 
The project 'Emmerich_TH' contains four meausrements of the total bedform field.
The project 'Emmerich_detail' contains 10 measurements of the single BEP along the center of the bedform field, which is more suitable for investigating bedform dynamics.

1_WLT:
This folder contains step 1 (detrending and wavelet analysis according to Bedforms-ATM). 
Only the file "1_WLT.R" needs to be executed by the user. All other scripts are called within the main script.
The parameters embedded by 3 #-characters must be defined by the user and then the script must be started. Results are stored in the results folder.    

2_ZC:
This folder contains step 2 (zerocrossing procedure and bedform statistics). 
First, the file "inputs.csv" must be opened. 
Here, the ranges for the input parameters of the zerocrossing procedure are specified (results of step 1 provide orientation). 
Then, the file "2_ZC.R" has to be executed by the user. All other scripts are called within the main script.
The parameters embedded by 3 #-characters must be defined by the user and then the script must be started. Results are stored in the results folder.    

3a_CC:
This folder contains step 3 (method 1 - cross correlation analysis). Note: Step 2 must be completed first.
The file "3a_CC.R" needs to be executed by the user.
The parameters embedded by 3 #-characters must be defined by the user and then the script must be started. Results are stored in the results folder.    

3b_CA:
This folder contains step 3 (method 2 - centroid analysis). Note: Step 2 must be completed first.
The file "3b_CA.R" needs to be executed by the user.
The parameters embedded by 3 #-characters must be defined by the user and then the script must be started. Results are stored in the results folder.  

4_POST:
This folder contains postprocessing routines to visualize the results. 
The file "4_GEOM.R" can be used to visualize bedform geometries while the file "4_DYN.R" can be used to visualize bedform dynamics.
The parameters embedded by 3 #-characters must be defined by the user and then the script must be started.  


....................
:     CREDITS      : 
....................

Julius Reich (Federal Institute of Hydrology): programming
Dr. Axel Winterscheid (Federal Institute of Hydrology): technical advisor


....................
: ACKNOLEDGEMENTS  : 
....................

The measurements to collect the example dataset were performed by the Federal Waterways and Shipping Administration (Wasserstraßen- und Schifffahrtsverwaltung des Bundes, WSV). 
The hydrographic processing and the provision of the individual BEPs was carried out by Felix Lorenz and Thomas Artz from the Federal Institute of Hydrology (Department for Geodesy and Remote Sensing). 


....................
:    REFERENCES    : 
....................

Frings, R.M., Gehres, N., de Jong, K., Beckhausen, C., Schüttrumpf, H. and Vollmer, S.. (2012). Rheno BT, User Manual, Institute of Hydraulic Engineering and Water Resources Management, RWTH Aachen University.
Gutierrez, R.R., Mallma, J.A., Núñez-González, F., Link, O. and Abad, J.D.. (2018). Bedforms-ATM, an open source software to analyze the scale-based hierarchies and dimensionality of natural bed forms. SoftwareX, 7, pp.184-189.


