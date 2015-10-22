%> @file footprint_comparison_konnik.m
%> @brief Setting up a NGS system to study modal projections of Zernike and SmallFootprintProjections.
%> @author Mikhail Konnik
%> @date   12 May 2015
%> 
%> @section setupngsforfootprint Set up a NGS (can be NTAO or MCAO) system and do some simulations with it.
%> Study the modal projections of Zernike defined in the metapupils (either at h=0 and in altitude). 
%> For the time being we have two functions to compute those projections: tel.footprintProjection and 
%> AnalyticalSmallFootprintExpansion.m They both compute the same result, the 
%> latter is much faster then the former which computes the projections using numerical 
%> integration instead of an analytical method. We should therefore consider adding a better 
%> coded AnalyticalSmallFootprintExpansion to our lam_utilities sort of file or 
%> added to the zernikeStats class 
%> 
%> The important point is that the methods must return the same projection matrix for Zernike
%> functions, with or without NGS stars constellation (asterism). 

clear all
close all
clc
% format long
format short

addpath('./OOMAO-Raven/OOMAOlibUpdated/')  % Only use OOMAOlibUpdated!
addpath('./OOMAO-Raven/')

% if (exist('samplingTime') == 0)
%% Setting up Parameters
samplingTime    = 0.001;% seconds
frameTime       = 0.01; % seconds
nIteration              = 100;
exposureTime            = nIteration*frameTime;
randn('state', 25);     % sets the global random state


%% Guide Stars
ngs = source; %% create a simple NGS guide star

% dAsterism               = 1.0; %0.5, 0.75, 1.0, 1.5, 2 diameter in arcmin
% guideStarWavelength     = photometry.R;
% guideStarMagnitude      = 14;  %[14,14.5,15,15.5,16,16.5,17]; % lower later
% 
% %% Asterism based on Raven CU pinholes
% if (dAsterism ~= 0)
% ngs = source('zenith',[arcsec(45.3),arcsec(35.5),arcsec(45.3)],...
%              'azimuth',[3.0309,0.7854,-1.6815],...
%              'wavelength',guideStarWavelength,'magnitude',guideStarMagnitude);
% end

%% Define the Asterism (star constellation) for the tomography
dAsterism               = 1.0; % 0.5, 0.75, 1.0, 1.5, 2 diameter in arcmin
guideStarWavelength     = photometry.R;
guideStarMagnitude      = 14;%[14,14.5,15,15.5,16,16.5,17]; % lower later

%% Guide Stars with Regular asterism
ngs = source('asterism',{[3,arcmin(dAsterism/2),0]},...
            'wavelength',guideStarWavelength,... 
            'magnitude',guideStarMagnitude);

% ngs = source('zenith',[arcsec(25.3),arcsec(35.5),arcsec(45.3)],...
%              'azimuth',[3.0309,0.7854,-1.6815],...
%              'wavelength',guideStarWavelength,'magnitude',guideStarMagnitude);

%% Atmosphere parameters
refWavelength   = 500e-9;
wavelength      = photometry.R;

altitude = [1]*1e3;
fractionalR0 = [1];
windSpeed = [10];
windDirection = [pi/2];
% % 
% altitude        = [1, 10]*1e3;
% fractionalR0    = [0.7,0.3];
% windSpeed       = [10,5];
% windDirection   = [pi/2,pi/2];


L0 = 40;  %% outer scale
r0 = 0.19; %% Fried parameter

%% Building an atmosphere object
atm = atmosphere(refWavelength,r0,L0,...
    'altitude',altitude,...
    'fractionnalR0',fractionalR0,...
    'windSpeed',windSpeed,...
    'windDirection',windDirection);
atm.wavelength = wavelength;


%% Telescope
nPx = 128; %% pixels across the telescope's pupil size
D = 8; %% teelscope diameter [m]
tel = telescope(D, ...
    'fieldOfViewInArcMin',3.5,...
    'resolution',nPx,...
    'samplingTime',samplingTime);


%%% Here we evaluate the Zernike functions on a telescope's pupil
%%% of size tel.D and pixels nPx  and thus getting the zernProj.modes
%%% that are just bunch of Zernike functions of maximum order maxRadialDegreeProj
maxRadialDegreeProj = 3;
zernModeMaxProj = zernike.nModeFromRadialOrder(maxRadialDegreeProj);
zernProj = zernike(2:zernModeMaxProj,'resolution',nPx,'pupil',tel.pupil,'D',tel.D);
%%% just in case, the zernProj.modes contains Zernike modes evaluated on
%%% MxN pupil of the telescope.


%% Combining the atmosphere and the telescope
tel = tel+atm; %%% Acording to the Profiler, this string of code takes 7.67 sec per 1 call!

ngs=ngs.*tel; % Propagation throught the atmosphere to the telescope
turbPhase = ngs.meanRmPhase; % Saving the turbulence aberrated phase

% end %%% if (exist(tel) == 0)


%% Using the SmallFootprintExpansion algorithm from Lindstrum:
%%% The function AnalyticalSmallFootprintExpansion computes the modal projection of Zernike polynomials analytically using Noll's 
%%% indexing convention onto a smaller telescope pupil, displaced by $\Delta x$ and $\Delta y$ and (possibly) rotated 
%%% (In our case, the pupil is not rotated at all).

% [projection_Analytic, projection_Analytic_Cellarray] = ...
%     AnalyticalSmallFootprintExpansion(zernProj, tel, ngs, atm);

tic
% [projection_modal_smallfootprint_alpha, projection_modal_smallfootprint_alphaCellarray] = ...
%     tool_analytical_small_footprint_expansion(zernModeMaxProj, tel, ngs, atm);

[projection_modal_smallfootprint_alpha, projection_modal_smallfootprint_alphaCellarray] = ...
    lam_utilities.tool_analytical_small_footprint_expansion(zernModeMaxProj, tel, ngs, atm);

t_analyticalfoot = toc;


%% Using the tel.footprintProjection algorithm
%%% P = tel.footprintProjection(obj,zernModeMax,src); here 
%%% the "obj" is the "tel+atm"
%%% the zernModeMax  is the number of maximal modes of Zernike
%%% This file is from /home/mkonnik/matlab/OOMAO_LAM/oomao_raven/OOMAO-Raven/OOMAOlibUpdated/telescope.m
%%% which obviously needs a documentation (conveniently missing....)
tic
    projection_modal_footprint = footprintProjection(tel, zernModeMaxProj, ngs);
t_footprint = toc;

projection_modal_footprint_matrix = cell2mat(projection_modal_footprint); %% converting the cell back to a matrix.


% break

%% Check that the both methods, namely AnalyticalSmallFootprintExpansion and footprintProjection, give the same result.
fprintf('Computational performance comparison: \n    - AnalyticalSmallFootprintExpansion takes %2.4f seconds;\n    - tel.footprintProjection takes %2.4f seconds.\n', t_analyticalfoot, t_footprint);

delta_two_methods = norm(projection_modal_footprint_matrix - projection_modal_smallfootprint_alpha);
fprintf('The norm (difference) between the two methods is %2.3e \n', delta_two_methods);