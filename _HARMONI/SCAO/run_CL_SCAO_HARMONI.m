clear all
%close all

%% LOCAL DEFSs
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
nL   = 74;%60; % E-ELT size
nPx  = 4;
nRes = nL*nPx*2;
D    = 37;
d    = D/nL; % lenslet pitch
samplingFreq = 500;

tel = telescope(D,'resolution',nRes,...
    'obstructionRatio',0.26,'fieldOfViewInArcsec',30,'samplingTime',1/samplingFreq);

%tel.pupil = pupil;


%% WAVEFRONT SENSOR
wfs = shackHartmann(nL,nPx*nL,0.5); %0.97

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

couplingCoeff = 0.4;

% CASE 1: Fried topology (Certesian regular grid)
bifa = influenceFunction('monotonic',couplingCoeff);

dm = deformableMirror(nL+1,'modes',bifa,...
    'resolution',tel.resolution,...
    'validActuator',wfs.validActuator);

% CASE 2: M4 actuator locations
bifM4 = influenceFunction('monotonic',couplingCoeff); 
m4 = load('../Coord_RepHexa');
pitch = 31.5e-3*2;
bifM4.actuatorCoord =  (m4.Centres_Act(:,1) + 1j*m4.Centres_Act(:,2))/pitch;
dm = deformableMirror(m4.nb_act,'modes',bifM4,'resolution',tel.resolution,...
    'validActuator',true(1,m4.nb_act)); 


figure,show(bifM4,'parent',subplot(1,2,1))
axis square
title('Monototic influence function')% The markers in the figures correspond to, from left to right, the points $P_k$ from $k=0$ to 6.

subplot(1,2,2)
scatter(m4.Centres_Act(:,1)/pitch, m4.Centres_Act(:,2)/pitch)
title('M4 actuator locations')
axis tight square
box on


%% INTERACTION MATRIX
% Lets switch off the display automatic update:
wfs.camera.frameListener.Enabled = false;
wfs.slopesListener.Enabled = false;

% and setup the optical path before the DM/WFS subsystem:
ngs = ngs.*tel;

calibDm = calibration(dm,wfs,ngs,ngs.wavelength/8,nL+1,'cond',1e2);

%% WAVEFRONT ESTIMATION

tel = tel + atm;
figure
imagesc(tel)
ngs = ngs.*tel*wfs;


%% Closed--loop NGS AO systems
science = source('wavelength',photometry.K);
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

%% STATIC AND QUASI STATIC ABERRATIONS

%StaticWaveNGS = {tel.pupil;JEFF_HSFOptim*ngs.waveNumber*0}; % l'amplitude complexe {amplitude;phase} 
%StaticWaveSCI = {tel.pupil;JEFF_HSFOptim*science.waveNumber*0}; % l'amplitude complexe {amplitude;phase} 
% expected loss of performance from static aberration if dm cannot fit it
%exp(-var(StaticWaveSCI{2}(tel.pupilLogical)))

%% SOURCE OBJECTS
tel = tel + atm;
dm.coefs = zeros(dm.nValidActuator,1);

ngs = source('zenith',zeros(1,1),'azimuth',zeros(1,1),'magnitude',12,'wavelength',photometry.R);
%ngsCombo = ngsCombo.*tel*StaticWaveNGS*dm*wfs;
ngs = ngs.*tel*dm*wfs;
s = source('zenith',zeros(1,1),'azimuth',zeros(1,1),'wavelength',photometry.K);
%scienceCombo = scienceCombo.*tel*StaticWaveSCI*dm*cam;
science = science.*tel*dm*cam;

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
h = imagesc(catMeanRmPhase(science));
axis xy equal tight
colorbar


%%
gain_cl  = 0.5;

dm.coefs = zeros(dm.nValidActuator,1);
flush(cam)
set(science,'logging',true)
set(science,'phaseVar',[])
cam.startDelay   = startDelay;
cam.frameListener.Enabled = true;

wfs.camera.photonNoise = false;
wfs.camera.readOutNoise = 0;
wfs.framePixelThreshold = wfs.camera.readOutNoise;

 %% CLOSED LOOP ITERATION
% The loop is closed for one full exposure of the science camera.
nIteration = startDelay + exposureTime;
wfsSlopesStack = zeros(wfs.nSlope,1);
for k=1:nIteration
    % Objects update
    +tel;
    +ngs;
    +science;
    % Closed-loop controller
    dm.coefs(:,1) = dm.coefs(:,1) - gain_cl*calibDm.M*wfsSlopesStack(:,1);    
    
    % Display
     set(h,'Cdata',catMeanRmPhase(science))
     drawnow
     
     wfsSlopesStack = wfs.slopes;
end
imagesc(cam)
set(h,'Cdata',catMeanRmPhase(science))
%% PERFORMANCE ANALYSIS FROM LOOP TELEMETRY

var_wfe_lsq = reshape(science(1).phaseVar(1:nIteration*2),2,[])';
wfe_lsq = sqrt(var_wfe_lsq)/science(1).waveNumber*1e6;
atm_wfe_rms = sqrt(zernikeStats.residualVariance(1,atm,tel))/ngs.waveNumber*1e6;
marechalStrehl_lsq = 1e2*exp(-mean(var_wfe_lsq(startDelay:end,2)));
psfStrehl = 1e2*cam.strehl
text(50,60,['PSF Strehl:' num2str(psfStrehl) '%'])

%% THEORETICAL PERFORMANCE ANALYSIS - CASE OF THE SINGLE INTEGRATOR LOOP WITH GAIN
var_fit = 0.23*(d/atm.r0)^(5/3)*(atm.wavelength/science(1).wavelength)^2;
var_alias = 0.07*(d/atm.r0)^(5/3)*(atm.wavelength/science(1).wavelength)^2;
var_tempo = phaseStats.closedLoopVariance(atm, tel.samplingTime,0.001,gain_cl)*(atm.wavelength/science(1).wavelength)^2;
marechalStrehl_lsq_theoretical = 100*exp(-var_fit-var_alias-var_tempo)
text(50,170,['Marechal approx:' num2str(100*marechalStrehl_lsq_theoretical) '%'])
