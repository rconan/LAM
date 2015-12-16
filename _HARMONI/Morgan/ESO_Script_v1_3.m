
% This routine reads the .mat files containing the segment PTT as a function of time.
% Three options are supported:
% 1) an instant in time can be chosen (graphic input using mouse click on time series diagram),
% 2) a time period can be chosen (graphic input using 2 mouse clicks (start and end) on time series diagram),
% 3) a time period can be chosen (graphic input of start on time series diagram, numeric input of duratin in seconds)

clear all; clc; close all,

REP_ErrorItem = 'C:\Users\mgray\AO_HARMONI\PERT_DYNAMIQUES\';%'/home/mgray/AO_HARMONI/PERT_DYNAMIQUES/';%
REP_Segmentation = 'C:\Users\mgray\AO_HARMONI\MIROIR_M1\';%'/home/mgray/HARMONI/MIROIR_M1/';%
REP_Version = 'VERSION_1_3\';% 'VERSION_1_3/'
File.ErrorItem = [REP_ErrorItem REP_Version 'M1wind2mps_hardpact_Optical.mat'];
File.Segmentation = [REP_Segmentation REP_Version 'SegmentVertices.txt'];
% specify the resolution of the phase screen in [m/pixel]
Parameters.Resolution = 38e-3;
% specify the size of the phase screen in [m]
Parameters.Size = 38.542;
% include the spider obscuration in the pupil
Parameters.Obscuration = true;
% do you like to save the phase screen
File.Save.Flag = true;
% where ?
File.Save.File = ['C:\Users\mgray\AO_HARMONI\RESULTATS\','ESO_WL2msHP_1014'];%'/home/mgray/AO_HARMONI/RESULTATS/';%
% can be 'pdf', 'png' or 'fits'
File.Save.Format = 'fits';

%%

PS = ESO_GeneratePS_v1_3(File,Parameters);



