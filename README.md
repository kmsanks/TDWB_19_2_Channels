# TDWB_19_2_Channels
#######INSERT LINK##########

## Introduction
This README file describes kmsanks/TDWB_19_2_Channels repository of scripts used to calculate changes in channel morphology and kinematics due to marsh deposition in an experimental setting. The raw experimental data can be found in the Tulane_Sediment_Dynamics_Stratigraphy_TSDS project space at: https://sead2.ncsa.illinois.edu/spaces/5825f529e4b0f3dd19c8d93a. The data used here is TDB-18-1 and TDWB-19-2-Surface-Processes. The processed data used in the code is contained herein. 

The experimental delta data come from two experiments run at the Tulane Sediment Dynamics Laboratory. Both experiments were setup identically, except for TDWB-19-2 (treatment) had marsh proxy deposition, while TDB-18-1 (control) did not. 
For more information on experimental conditions, please see the data repositories hosted at: https://sead2.ncsa.illinois.edu/spaces/5825f529e4b0f3dd19c8d93a.

## What does this repository provide?
This repository contains all of the code relevant to produce the results and figures contained in the manuscript titled "Marsh induced backwater: the influence of non-fluvial sedimentation on a delta's channel morphology and kinematics (Sanks et al., 2023)". The code contained herein was written using MATLAB R2022a. Note: all the "Core Script Files" should run out the box, as relative paths were used in creation of this script.

## Contents
The repository contains three folders:
 
 1. code
 2. data
 3. figures
 
Please clone the repository in full in order to use the repo. Then download the .mat files from FIGSHARE REPO and put into the data folder. All figure outputs will populate within the figures folder. Note: some figures were modified in Illustrator for visual purposes.

## Data: TDWB_19_2_Channels/data
Please download the following files from https://doi.org/10.6084/m9.figshare.22320811.v1 and put into this folder. 

1. ZD_18.mat - A matrix of size 796x522x560, where 560 is time. Each timestep contains a elevation data collected from the control experiment via LiDAR and post-processed into a 5mmx5mm grid. The data has 796 rows and 522 columns corresponding to basin location.
  
  *Note that the matrix has 560 timesteps, which corresponds to one LiDAR scan per hour. The first run hour is t = 1, total run time = 560 hours.
  
2. ZD_19.mat - A matrix of size 750x747x281, where 281 is time. Each timestep contains a elevation data collected from the treatment experiment via LiDAR and post-processed into a 5mmx5mm grid. The data has 750 rows and 747 columns corresponding to basin location. 
  
  *Note that the matrix only has 281 timesteps because the LiDAR data was collected every other hour. The first run hour is t = 0, total run time = 560 hours.
  
3. CM_18.mat - A matrix of size 796x522x560. The data contained herein is the same reference frame as described above for ZD_18 (control). This data is a binary matrix, where channel pixels = 1 and non-channel pixels = 0.
 
4. CM_19.mat - A matrix of size 750x747x560. The data contained herein is the same reference frame as described above for ZD_19 (treatment), but contains a channel map for each timestep. This data is a binary matrix, where channel pixels = 1 and non-channel pixels = 0.

   *Note that some timesteps do not have a channel map due to issues with dye timing or cart artifacts in the LiDAR scan. These timesteps will have an empty matrix of NaNs.

5. flowscreen18.mat - A binary matrix of 796x522x560, where 1 = flow, 0 = no-flow. 

6. flowscreen19.mat - A binary matrix of 750x747x560, where 1 = flow, 0 = no-flow.  

   *Note that some timesteps have impaired flowmaps due to missing dye or cart artifacts in the LiDAR scan. These timesteps are removed from analyses and noted in flowfraction.m lines 147-152.  

## Core Script Files: TDWB_19_2_Channels/code
1. flowfraction.m - This script calculates total, channel, and overbank flow area and fraction for both the control and treatment experiments for the terrestrial delta. 
   * This sciript produces Figure 3 and some results from Table 1. 
   * Data needed: ZD_18.mat, CM_18.mat, flowscreen18.mat, ZD_19.mat, CM_19.mat, and flowscreen19.mat.
2. channelproperties.m - This script calculates the average channel bed elevations relative to sea level, channel depth, channel width, channel aggradation, far-field aggradation, and channel in-filling rate as a function of radial distance from the apex for both the control and treatment experiments.
   * This script produces Figures 4a, 4b, 4d, 6, and some results from Table 1.
   * Data needed: ZD_18.mat, CM_18.mat, ZD_19.mat, and CM_19.mat.
3. channellength.m - This script calculates the channel length for both the control and treatment experiments.
   * This script produces Figures 4c and some results from Table 1.
   * Data needed: ZD_18.mat, CM_18.mat, ZD_19.mat, and CM_19.mat.
4. backwaterlength.m - This script calculates the backwater length for the control and treatment experiments through time. See manuscript for more details.
   * This script produces Figure 5 and some results in Table 1.
   * Data needed: ZD_18.mat, CM_18.mat, ZD_19.mat, and CM_19.mat.
5. lateralmobility.m - This script calculates lateral channel mobility, fraction of the delta that is unmodified, and channelization statistics. 
   * This script produces Figure 7a, 7b, and SI Figures B1, B2, B6, and B7, as well as some results in Table 1.
   * Data needed: ZD_18.mat, CM_18.mat, ZD_19.mat, and CM_19.mat.
6. lateralmobility_radial.m - This script calculates radial lateral channel mobility. 
   * This script produces Figure 7c and 7d.
   * Data needed: ZD_18.mat, CM_18.mat, ZD_19.mat, and CM_19.mat.
  
## Supplementary Script Files: TDWB_19_2_Channels/code
 1. planformoverlap.m - This script calculates the channel decorrelation metric based on Wickert et al.(2013).
   * This script produces SI Figure B5.
   * Data needed: ZD_18.mat, CM_18.mat, ZD_19.mat, and CM_19.mat.
 2. efoldmobility.m - This script creates the time it takes for the channels to visit one e-fold (66% for control and 67% for treatment) of the terrestrial delta top as a function of time.
   * This script produces SI Figure B8.
   * Data needed: ZD_18.mat, CM_18.mat, ZD_19.mat, and CM_19.mat. 

## Figures: TDWB_19_2_Channels/figures
The figures will be created when the scripts are run. Currently there is only a README.txt files in this folder as a place holder. Note that some figures were modified in Illustrator for publication.

## Using this repository
Clone the repository. The scripts can be run in any order, but please note that the data (.mat files) need to be downloaded from FIGSHRE REPOSITORY and placed in the data folder first.
