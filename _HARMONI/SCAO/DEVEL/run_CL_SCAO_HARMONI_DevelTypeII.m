clear all
close all

%% local defs
%load('/data/HARMONI/SCAO/TESTS/ESO_PupilM1_740.mat');
%pupil = Masq_M1; clear Masq_M1;

%load('/data/HARMONI/SCAO/TESTS/MAT_FI_M4_740.mat');

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
nPx  = 10;
nRes = nL*nPx;
D    = 37;
d    = D/nL; % lenslet pitch
samplingFreq = 500;

tel = telescope(D,'resolution',nRes,...
    'obstructionRatio',0.3,'fieldOfViewInArcsec',30,'samplingTime',1/samplingFreq);

%tel.pupil = pupil;


%% WAVEFRONT SENSOR
wfs = shackHartmann(nL,nRes,0.5);

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


dm = deformableMirror(nL+1,'modes',bifa,...
    'resolution',tel.resolution,...
    'validActuator',wfs.validActuator);

% replace theoretical IFs by M4 IFs
%bifa = influenceFunction('monotonic',0.0);
%bifa.modes = MatFI_M4;
%dm = deformableMirror(nL+1,'modes',bifa);


%% INTERACTION MATRIX
% Lets switch off the display automatic update:
wfs.camera.frameListener.Enabled = false;
wfs.slopesListener.Enabled = false;

% and setup the optical path before the DM/WFS subsystem:
ngs = ngs.*tel;

calibDm = calibration(dm,wfs,ngs,ngs.wavelength,nL+1,'cond',1e2);

%% WAVEFRONT ESTIMATION

tel = tel + atm;
figure
imagesc(tel)
ngs = ngs.*tel*wfs;
%%
% <latex>
% The wavefront reconstruction is done by estimating the wavefront values
% at the corner of the lenslets.
% A new telescope is defined identical to the previous one but with a lower
% resolution and a pupil defined by the map of "valid actuator".
% </latex>
telLowRes = telescope(tel.D,'resolution',nL+1,...
    'fieldOfViewInArcsec',30,'samplingTime',1/500);
telLowRes.pupil = wfs.validActuator;
%%

telLowRes= telLowRes + atm;
ngs = ngs.*telLowRes;
phase = ngs.meanRmOpd;

bifaLowRes = influenceFunction('monotonic',0.4);
dmLowRes = deformableMirror(nL+1,'modes',bifaLowRes,'resolution',nL+1,...
    'validActuator',wfs.validActuator);

slmmse = slopesLinearMMSE(wfs,tel,atm,ngs,'mmseStar',ngs);

%%
F = 2*bifaLowRes.modes(wfs.validActuator,:);
iF = pinv(full(F),1e-1);

%% Closed--loop NGS AO systems
science = source('wavelength',photometry.H);
cam = imager(tel);
%%
tel = tel - atm;
science = science.*tel*cam;
figure(31416)
imagesc(cam,'parent',subplot(2,1,1))
cam.frameListener.Enabled = true;
%%
cam.referenceFrame = cam.frame;
+science;
fprintf('Strehl ratio: %4.1f\n',cam.strehl)
%%
tel = tel + atm;
+science;
fprintf('Strehl ratio: %4.1f\n',cam.strehl)

%%
    filename = ['/data/HARMONI/SCAO/SIMUL_MORGAN/PERTURB_M1/TT/MG_WL8msTT_' num2str(1) '.fits'];
    telDisturbance = fitsread(filename);
    StaticWaveNGS = {tel.pupil;telDisturbance*ngs.waveNumber}; % l'amplitude complexe {amplitude;phase} 
    StaticWaveSCI = {tel.pupil;telDisturbance*science.waveNumber}; % l'amplitude complexe {amplitude;phase} 

%StaticWaveNGS = {tel.pupil;JEFF_HSFOptim*ngs.waveNumber*0}; % l'amplitude complexe {amplitude;phase} 
%StaticWaveSCI = {tel.pupil;JEFF_HSFOptim*science.waveNumber*0}; % l'amplitude complexe {amplitude;phase} 

%% SOURCE OBJECTS
tel = tel - atm;
dm.coefs = zeros(dm.nValidActuator,1);

ngsCombo = source('zenith',zeros(1,1),'azimuth',zeros(1,1),'magnitude',12,'wavelength',photometry.R);
ngsCombo = ngsCombo.*tel*StaticWaveNGS*dm*wfs;
%ngsCombo = ngsCombo.*tel*dm*wfs;

scienceCombo = source('zenith',zeros(1,1),'azimuth',zeros(1,1),'wavelength',photometry.H);
scienceCombo = scienceCombo.*tel*StaticWaveSCI*dm*cam;
%scienceCombo = scienceCombo.*tel*dm*cam;

% expected loss of performance from static aberration if dm cannot fit it
%exp(-var(StaticWaveSCI{2}(tel.pupilLogical)))
%% %% LOOP INIT

flush(cam)
cam.frame = cam.frame*0;
cam.clockRate    = 1;
exposureTime     = 1000;
cam.exposureTime = exposureTime;
startDelay       = 100;
figure(31416)
imagesc(cam,'parent',subplot(2,1,1))
 cam.frameListener.Enabled = true;
subplot(2,1,2)
h = imagesc(catMeanRmPhase(scienceCombo));
axis xy equal tight
colorbar


%%
gain_cl  = 0.5;
gain_pol = 0.7;



dm.coefs = zeros(dm.nValidActuator,1);
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


zern = zernike(2:3,'resolution',nL+1, 'pupil',dmLowRes.validActuator);
TTd = iF*zern.modes(dmLowRes.validActuator,:);
TTd = zern.modes(dmLowRes.validActuator,:);
DMTTRem = eye(dmLowRes.nValidActuator) - TTd*pinv(TTd);


t = ones(wfs.nValidLenslet,1);
TTs = [t 0*t; 0*t,t];
iTTs = pinv(TTs);
SlopeTTRem = eye(2*wfs.nValidLenslet) - TTs*iTTs;
gainsTypeII(1) = 0.023;
C0 = 51.2515;
C1 = -49.5088;
D1 = 0.7426;


 %% CLOSED LOOP ITERATION
% The loop is closed for one full exposure of the science camera.
nIteration = startDelay + exposureTime;
wfsSlopesStack = zeros(wfs.nSlope,2);
uTT = zeros(dmLowRes.nValidActuator,1);
dm.coefs = zeros(dm.nValidActuator,1);
tipTilt = zeros(2,2);
o1 = zeros(2,1);
o2 = o1;
o3 = o1;
for k=1:nIteration
    k
    % Objects update
    filename = ['/data/HARMONI/SCAO/SIMUL_MORGAN/PERTURB_M1/TT/MG_WL8msTT_' num2str(k) '.fits'];
    telDisturbance = fitsread(filename);
    StaticWaveNGS = {tel.pupil;telDisturbance*ngs.waveNumber}; % l'amplitude complexe {amplitude;phase} 
    StaticWaveSCI = {tel.pupil;telDisturbance*science.waveNumber}; % l'amplitude complexe {amplitude;phase} 

    +tel;
    %+ngsCombo;
    %+scienceCombo;
    
    ngsCombo = ngsCombo.*tel*StaticWaveNGS*dm*wfs;
    scienceCombo = scienceCombo.*tel*StaticWaveSCI*dm*cam;



    % Closed-loop controller
    dm.coefs(:,1) = dm.coefs(:,1) - gain_cl*calibDm.M*wfsSlopesStack(:,1);    
    
    % Remove TT from DM
    dm.coefs = DMTTRem*dm.coefs;
    
    % Tip/Tilt closed-loop measurement
    tipTilt(:,1) =  ngsCombo.wavelength/8 *2 * iTTs*wfsSlopesStack(:,1); % factor 2 is to account for stroke units the mirror
    
    
    % SINGLE INTEGRATOR
    uTT = uTT -  gain_cl * TTd * tipTilt(:, 1);
    
    
    
    % DOUBLE INTEGRATOR
    %o1 = gainsTypeII(1)*(C0*tipTilt(:,1) + C1*tipTilt(:,2)) - D1*o1;
    %o = o1 + o2;
    %o3 = o2 + o3;
    %uTT = -TTd*o3;
    
    %tipTilt(:,2) = tipTilt(:,1);
    
        % add TT controls back to dm commands
    dm.coefs = dm.coefs - uTT;

    % Pseudo-open-loop controller
%     dm.coefs(:,2) = (1-gain_pol)*dm.coefs(:,2) + ...
%         gain_pol*iF*( slmmse*( wfsSlopesStack(:,2) - calibDm.D*dm.coefs(:,2) ) );
    
    % Display
     set(h,'Cdata',catMeanRmPhase(scienceCombo))
     drawnow
     
     wfsSlopesStack = wfs.slopes;

end
imagesc(cam)
set(h,'Cdata',catMeanRmPhase(scienceCombo))
%% PERFORMANCE ANALYSIS

var_wfe_lsq = reshape(scienceCombo(1).phaseVar(1:nIteration*2),2,[])';
wfe_lsq = sqrt(var_wfe_lsq)/scienceCombo(1).waveNumber*1e6;
%var_wfe_lmmse = reshape(scienceCombo(2).phaseVar(1:nIteration*2),2,[])';
%wfe_lmmse = sqrt(var_wfe_lmmse)/scienceCombo(1).waveNumber*1e6;
atm_wfe_rms = sqrt(zernikeStats.residualVariance(1,atm,tel))/ngs.waveNumber*1e6;
marechalStrehl_lsq = 1e2*exp(-mean(var_wfe_lsq(startDelay:end,2)));
%marechalStrehl_lmmse = 1e2*exp(-mean(var_wfe_lmmse(startDelay:end,2)));
psfStrel = 1e2*cam.strehl

%% THEORETICAL PERFORMANCE ANALYSIS - CASE OF THE SINGLE INTEGRATOR LOOP WITH GAIN
var_fit = 0.23*(d/atm.r0)^(5/3)*(atm.wavelength/scienceCombo(1).wavelength)^2;
var_alias = 0.07*(d/atm.r0)^(5/3)*(atm.wavelength/scienceCombo(1).wavelength)^2;
var_tempo = phaseStats.closedLoopVariance(atm, tel.samplingTime,0.001,gain_cl)*(atm.wavelength/scienceCombo(1).wavelength)^2;

marechalStrehl_lsq_theoretical = exp(-var_fit-var_alias-var_tempo)