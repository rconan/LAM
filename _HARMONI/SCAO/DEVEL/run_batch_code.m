clear all
close all

%% local defs
%load('/data/HARMONI/SCAO/TESTS/ESO_PupilM1_740.mat');
%pupil = Masq_M1; clear Masq_M1;

%load('/data/HARMONI/SCAO/TESTS/MAT_FI_M4_740.mat');
%load('/data/HARMONI/SCAO/SIMUL_MORGAN/MIROIR_M4/MAT_FI_M4_740.mat')
% static maps
%load('/data/HARMONI/SCAO/TESTS/JEFF_HSFreq_740.mat');
%% SOURCE
ngs = source('wavelength',photometry.R); % R-band 

%% ATMOSPHERE

% atm = atmosphere(photometry.V0,0.1587,50,...
%     'fractionnalR0',[0.5,0.3,0.2],'altitude',[0e3,5e3,12e3],...
%     'windSpeed',[10,5,20],'windDirection',[0,pi/2,pi]);

atm = atmosphere(photometry.V0,0.1587,50,...
    'fractionnalR0',[1],'altitude',[0e3],...
    'windSpeed',[12],'windDirection',[0]);


%% TELESCOPE
nL   = 10;%60; % E-ELT size
nPx  = 10;
nRes = nL*nPx;
D    = 37;
d    = D/nL; % lenslet pitch
samplingFreq = 500;

tel = telescope(D,'resolution',nRes,...
    'obstructionRatio',0.26,'fieldOfViewInArcsec',30,'samplingTime',1/samplingFreq);



tel = tel + atm;