%% ADAPTIVE OPTICS MODELING WITH OOMAO - adaptation to test the integration of the pyramid wave fromt sensor
% Demonstrate how to build a simple closed-loop single conjugated adaptive
% optics system

%%
%close all
clear all
clc

atm = atmosphere(photometry.V,1.5,30,...
    'altitude',0,...
    'fractionnalR0',1,...
    'windSpeed',100,...
    'windDirection',0);

%% Definition of the telescope
nPx = 60;
tel = telescope(3.6,...
    'fieldOfViewInArcMin',2.5,...
    'resolution',nPx,...
    'samplingTime',1/100);

%% Definition of a calibration source
ngs = source('wavelength',photometry.J);
ngs_sh = source('wavelength',photometry.J);
%% Definition of the wavefront sensor
%Experimental Pyramid WFS, expect some rough edges and maybe some bugs
%The pyramid takes only one argument, which is the pixel resolution of the
%telescope it is associated with.
wfs=pyramid(nPx);
nLenslet = 10;
wfs_sh = shackHartmann(nLenslet,nPx,0.75);
%%
% Propagation of the calibration source to the WFS through the telescope
ngs = ngs.*tel*wfs;
ngs_sh = ngs_sh.*tel*wfs_sh;

%%
%Calibration of the sensor on the current light
%wfs.setmodulation(3)
wfs.INIT
wfs_sh.INIT
%%
% A new frame read-out and slopes computing:
+wfs;
+wfs_sh;
%%
% The WFS camera display:
figure(1)
imagesc(wfs.camera)
figure(2)
imagesc(wfs_sh.camera)
%%
% The WFS slopes display:
figure(3)
slopesDisplay(wfs)
figure(4)
slopesDisplay(wfs_sh)


%%
%Propagating through the atmosphere
% tel = tel+atm;
% 
% for ii = 1 :10
% ngs=ngs.*+tel;
% end
% ngs_sh=ngs_sh.*tel;
% 
% ngs=ngs*wfs;
% ngs_sh=ngs_sh*wfs_sh;
% 
% % A new frame read-out and slopes computing:
% +wfs;
% +wfs_sh;

%% 
zer = zernike(3,tel.D, 'resolution', nPx);
zer.c = 5;
ngs = ngs.*tel*zer*wfs;
ngs_sh = ngs_sh.*tel*zer*wfs_sh;

% A new frame read-out and slopes computing:
+wfs;
+wfs_sh;

%%
% The WFS camera display:
figure(1)
imagesc(wfs.camera)
figure(2)
imagesc(wfs_sh.camera)
%%
% The WFS slopes display:
figure(3)
slopesDisplay(wfs)
figure(4)
slopesDisplay(wfs_sh)
% 
% wfs.hilbertSlopes(tel,ngs)
% figure(5)
% slopesDisplay(wfs)
