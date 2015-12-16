
clear all
close all


%% SOURCE
ngs = source('wavelength',photometry.R); % R-band 

%% ATMOSPHERE

% atm = atmosphere(photometry.V0,0.1587,50,...
%     'fractionnalR0',[0.5,0.3,0.2],'altitude',[0e3,5e3,12e3],...
%     'windSpeed',[10,5,20],'windDirection',[0,pi/2,pi]);
r0 = 0.1587;

atm = atmosphere(photometry.V,r0,30,...
    'altitude',[0,4,10]*1e3,...
    'fractionnalR0',[0.7,0.25,0.05],...
    'windSpeed',[5,10,20],...
    'windDirection',[0,pi/4,pi]);

%% TELESCOPE
nLenslet   = 10;%60; % E-ELT size
nPx  = 6;
nRes = nLenslet*nPx;
D    = 8;
d    = D/nLenslet; % lenslet pitch
samplingFreq = 500;

tel = telescope(D,'resolution',nRes,...
    'obstructionRatio',0,'fieldOfViewInArcMin',2.5,'samplingTime',1/samplingFreq);


%% Definition of the wavefront sensor
%Experimental Pyramid WFS, expect some rough edges and maybe some bugs
%The pyramid takes only one argument, which is the pixel resolution of the
%telescope it is associated with.
pyr = pyramid(nLenslet,nRes,'modulation',2,'binning',1);

pyr.camera.readOutNoise = 0;

% Propagation of the calibration source to the WFS through the telescope
ngs = ngs.*tel*pyr;

%Calibration of the sensor on the current light
pyr.INIT

% A new frame read-out and slopes computing:
+pyr;



%% WAVEFRONT SENSOR
wfs = shackHartmann(nLenslet,nRes,0.5);

ngs = ngs.*tel*wfs;

wfs.INIT

+wfs;
figure
imagesc(wfs.camera,'parent',subplot(3,2,[1,4]))
slopesDisplay(wfs,'parent',subplot(3,2,[5,6]))
%%
% The next 2 commands allow the displays of the frame and of the slopes to
% be updated when a new frame and new slopes are computed
wfs.camera.frameListener.Enabled = true;
wfs.slopesListener.Enabled = true;
%% SH-WFS GAIN CALIBRATION
% The WFS must be calibrated such as for 1rd of tip--tilt wavefront , it will
% measured a slopes of 1rd.
% whereas the source is progressively moved off-axis
wfs.pointingDirection = zeros(2,1);

pixelScale = ngs.wavelength/...
    (2*d*wfs.lenslets.nyquistSampling);
tipStep = pixelScale/2;
nStep   = floor(nPx/3)*2;
sx      = zeros(1,nStep+1);
u       = 0:nStep;
wfs.camera.frameListener.Enabled = false;
wfs.slopesListener.Enabled = false;
warning('off','oomao:shackHartmann:relay')
for kStep=u
    ngs.zenith = -tipStep*kStep;
    +ngs;
    drawnow
    sx(kStep+1) = median(wfs.slopes(1:end/2));
end
warning('on','oomao:shackHartmann:relay')
Ox_in  = u*tipStep*constants.radian2arcsec;
Ox_out = sx*ngs.wavelength/d/2*constants.radian2arcsec;
%figure
%plot(Ox_in,Ox_out)
%grid
slopesLinCoef = polyfit(Ox_in,Ox_out,1);
wfs.slopesUnits = 1/slopesLinCoef(1);

% The source is reset on--axis and the WFS is set to always be aligned to
% the source by setting \oop{shackHartmann}{pointingDirection} to empty.
ngs.zenith = 0;
wfs.pointingDirection = [];
%% DEFORMABLE MIRROR

bifa = influenceFunction('monotonic',0.4);
figure,show(bifa,'parent',subplot(1,2,1))
title('Monototic influence function')% The markers in the figures correspond to, from left to right, the points $P_k$ from $k=0$ to 6.


dm = deformableMirror(nLenslet+1,'modes',bifa,...
    'resolution',tel.resolution,...
    'validActuator',wfs.validActuator);

% replace theoretical IFs by M4 IFs
%bifa = influenceFunction('monotonic',0.0);
%bifa.modes = MatFI_M4;
%dm = deformableMirror(nL+1,'modes',bifa);

%% PYRAMID INTERACTION MATRIX
ngs = ngs.*tel;
calib = calibration(dm,pyr,ngs,ngs.wavelength/40);
calib.nThresholded = 30;


%% INTERACTION MATRIX
% Lets switch off the display automatic update:
wfs.camera.frameListener.Enabled = false;
wfs.slopesListener.Enabled = false;

% and setup the optical path before the DM/WFS subsystem:
ngs = ngs.*tel;
calibDm = calibration(dm,wfs,ngs,ngs.wavelength,nLenslet+1,'cond',1e2);



%% WAVEFRONT ESTIMATION

tel = tel + atm;
figure
imagesc(tel)
ngs = ngs.*tel*wfs;


%% Closed--loop NGS AO systems
science = source('wavelength',photometry.H);
cam = imager();
%%
tel = tel - atm;
science = science.*tel*cam;
figure(31416)
imagesc(cam,'parent',subplot(2,1,1))
%cam.frameListener.Enabled = true;
%%
cam.referenceFrame = cam.frame;
+science;
fprintf('Strehl ratio: %4.1f\n',cam.strehl)
%%
tel = tel + atm;
+science;
fprintf('Strehl ratio: %4.1f\n',cam.strehl)

%% SOURCE OBJECTS
tel = tel + atm;
dm.coefs = 0;

ngsCombo = source('zenith',zeros(1,1),'azimuth',zeros(1,1),'magnitude',10,'wavelength',photometry.R);
ngsCombo = ngsCombo.*tel*dm*wfs;

ngsPyr = source('zenith',zeros(1,1),'azimuth',zeros(1,1),'magnitude',10,'wavelength',photometry.R);
ngsPyr = ngsPyr.*tel*dm*pyr;


scienceCombo = source('zenith',zeros(1,1),'azimuth',zeros(1,1),'wavelength',photometry.H);
%scienceCombo = scienceCombo.*tel*StaticWaveSCI*dm*cam;
scienceCombo = scienceCombo.*tel*dm*cam;

%% %% LOOP INIT

flush(cam)
cam.frame = cam.frame*0;
cam.clockRate    = 1;
exposureTime     = 100;
cam.exposureTime = exposureTime;
startDelay       = 20;
figure(31416)
imagesc(cam,'parent',subplot(2,1,1))
% cam.frameListener.Enabled = true;
subplot(2,1,2)
h = imagesc(catMeanRmPhase(scienceCombo));
axis xy equal tight
colorbar


%%
gain_cl  = 0.4;
gain_pol = 0.7;


dm.coefs = 0;
flush(cam)
set(scienceCombo,'logging',true)
set(scienceCombo,'phaseVar',[])
slmmse.wavefrontSize = [dm.nValidActuator,1];
slmmse.warmStart = true;
cam.startDelay   = startDelay;
cam.frameListener.Enabled = false;
% set(ngsCombo,'magnitude',8)
wfs.camera.photonNoise = true;
wfs.camera.readOutNoise = 1;
wfs.framePixelThreshold = wfs.camera.readOutNoise;

%% CLOSED LOOP ITERATION WITH PYRAMID
% The loop is closed for one full exposure of the science camera.
cam.frameListener.Enabled = 1;

flush(cam)
cam.frame = cam.frame*0;
dm.coefs = zeros(dm.nValidActuator,1);
nIteration = startDelay + exposureTime;
for k=1:nIteration
    % Objects update
    +tel;
    +ngsPyr;
    +scienceCombo;
    % Closed-loop controller
    dm.coefs(:,1) = dm.coefs(:,1) - gain_cl*calib.M*pyr.slopes;    
    
    % Pseudo-open-loop controller
    %dm.coefs(:,2) = (1-gain_pol)*dm.coefs(:,2) + ...
    %    gain_pol*iF*( slmmse*( wfsSlopesStack(:,2) - calibDm.D*dm.coefs(:,2) ) );
    
    % Display
     set(h,'Cdata',catMeanRmPhase(scienceCombo))
     drawnow
     
end
imagesc(cam)
set(h,'Cdata',catMeanRmPhase(scienceCombo))



%% CLOSED LOOP ITERATION 
% The loop is closed for one full exposure of the science camera.
nIteration = startDelay + exposureTime;
wfsSlopesStack = zeros(wfs.nSlope,1);
for k=1:nIteration
    % Objects update
    +tel;
    +ngsCombo;
    +scienceCombo;
    % Closed-loop controller
    dm.coefs(:,1) = dm.coefs(:,1) - gain_cl*calibDm.M*wfs.slopes;    
    
    % Pseudo-open-loop controller
    %dm.coefs(:,2) = (1-gain_pol)*dm.coefs(:,2) + ...
    %    gain_pol*iF*( slmmse*( wfsSlopesStack(:,2) - calibDm.D*dm.coefs(:,2) ) );
    
    % Display
     set(h,'Cdata',catMeanRmPhase(scienceCombo))
     drawnow
     

end
imagesc(cam)
set(h,'Cdata',catMeanRmPhase(scienceCombo))



%% PERFORMANCE ANALYSIS

var_wfe_lsq = reshape(scienceCombo(1).phaseVar(1:nIteration*2),2,[])';
wfe_lsq = sqrt(var_wfe_lsq)/scienceCombo(1).waveNumber*1e6;
var_wfe_lmmse = reshape(scienceCombo(2).phaseVar(1:nIteration*2),2,[])';
wfe_lmmse = sqrt(var_wfe_lmmse)/scienceCombo(1).waveNumber*1e6;
atm_wfe_rms = sqrt(zernikeStats.residualVariance(1,atm,tel))/ngs.waveNumber*1e6;
marechalStrehl_lsq = 1e2*exp(-mean(var_wfe_lsq(startDelay:end,2)));
marechalStrehl_lmmse = 1e2*exp(-mean(var_wfe_lmmse(startDelay:end,2)));
psfStrel = 1e2*cam.strehl

%% THEORETICAL PERFORMANCE ANALYSIS - CASE OF THE SINGLE INTEGRATOR LOOP WITH GAIN
var_fit = 0.23*(d/atm.r0)^(5/3)*(atm.wavelength/scienceCombo(1).wavelength)^2;
var_alias = 0.07*(d/atm.r0)^(5/3)*(atm.wavelength/scienceCombo(1).wavelength)^2;
var_tempo = phaseStats.closedLoopVariance(atm, tel.samplingTime,0.001,gain_cl)*(atm.wavelength/scienceCombo(1).wavelength)^2;

marechalStrehl_lsq_theoretical = exp(-var_fit-var_alias-var_tempo)