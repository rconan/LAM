%% INTRODUCTION
%{
This is the template file for the Harmoni-LTAO E2E simulation
%}

%% LOCAL DEFINITIONS


%% SOURCES
% on--axis souce
ngs = source;

%% ATMOSPHERE
r0 = 16e-2;             % coherence lenght in meters at 0.5microns
L0 = 30;                % Outer scale in meters

atm = atmosphere(photometry.V0,r0,L0,...
    'fractionnalR0',[0.5,0.3,0.2],'altitude',[0e3,5e3,12e3],...
    'windSpeed',[10,5,20],'windDirection',[0,pi/2,pi]);

%% TELESCOPE
nL   = 20;              % number of lenslets
nPx  = 6;              % number of pixels per lenslet
nRes = nL*nPx;          % resolution on the pupil plane (no of pixels)
D    = 8;              % telescope primary mirror diameter
d    = D/nL;            % lenslet pitch
samplingFreq = 500;     % WFS sampling time
obstructionRatio=0.3;   % central obscuration ratio
fieldOfViewInArcsec = 120; %fieldOfViewInArcsec

tel = telescope(D,'resolution',nRes,...
    'obstructionRatio',obstructionRatio,'fieldOfViewInArcsec',fieldOfViewInArcsec,'samplingTime',1/samplingFreq);

%% OIWFS 
ngs.wavelength = photometry.H;

OIWFS = shackHartmann(1,8,1);
ngs = ngs.*tel*OIWFS;

OIWFS.INIT

+OIWFS;
figure
imagesc(OIWFS.camera,'parent',subplot(3,2,[1,4]))
slopesDisplay(OIWFS,'parent',subplot(3,2,[5,6]))

%% OIWFS gain calibration
% <latex>
% The WFS must be calibrated such as for 1rd of tip--tilt wavefront , it will
% measured a slopes of 1rd.
% To do so, the \oop{shackHartmann}{pointingDirection} property is set on-axis
% </latex>
OIWFS.pointingDirection = zeros(2,1);
ngs.wavelength = photometry.K;

pixelScale = ngs.wavelength/...
    (2*D*OIWFS.lenslets.nyquistSampling);
tipStep = pixelScale/2;
nStep   = floor(nPx/3)*2;
sx      = zeros(1,nStep+1);
u       = 0:nStep;
OIWFS.camera.frameListener.Enabled = false;
OIWFS.slopesListener.Enabled = false;
warning('off','oomao:shackHartmann:relay')
for kStep=u
    ngs.zenith = -tipStep*kStep;
    +ngs;
    drawnow
    sx(kStep+1) = median(OIWFS.slopes(1:end/2));
end
warning('on','oomao:shackHartmann:relay')
Ox_in  = u*tipStep*constants.radian2arcsec;
Ox_out = sx*ngs.wavelength/D/2*constants.radian2arcsec;
%figure
plot(Ox_in,Ox_out)
%grid
slopesLinCoef = polyfit(Ox_in,Ox_out,1);
OIWFS.slopesUnits = 1/slopesLinCoef(1);

ngs.zenith = 0;
OIWFS.pointingDirection = [];

%% LGS WAVE-FRONT SENSOR
ngs.wavelength = photometry.Na;
wfs = shackHartmann(nL,nRes,0.85);

% WFS initialisation
ngs = ngs.*tel*wfs;
wfs.INIT
+wfs;
figure
imagesc(wfs.camera,'parent',subplot(3,2,[1,4]))
slopesDisplay(wfs,'parent',subplot(3,2,[5,6]))

% The next 2 commands allow the displays of the frame and of the slopes to
% be updated when a new frame and new slopes are computed
wfs.camera.frameListener.Enabled = true;
wfs.slopesListener.Enabled = true;

%% LGS WFS GAIN CALIBRATION 
% WFS must be calibrated such as for 1rd of tip--tilt wavefront , it will measured a slopes of 1rd.
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

ngs.zenith = 0;
wfs.pointingDirection = [];

%% DEFORMABLE MIRROR
dmCrossCouplingCoeff = 0.4;
bifa = influenceFunction('monotonic',dmCrossCouplingCoeff);
figure,show(bifa,'parent',subplot(1,2,1))
title('Monototic influence function')
bifb = influenceFunction('overshoot',dmCrossCouplingCoeff);
show(bifb,'parent',subplot(1,2,2))
title('Overshooted influence function')

dm = deformableMirror(nL+1,'modes',bifa,...
    'resolution',tel.resolution,...
    'validActuator',wfs.validActuator);

bifaLowRes = influenceFunction('monotonic',dmCrossCouplingCoeff);
dmLowRes = deformableMirror(nL+1,'modes',bifaLowRes,'resolution',nL+1,...
    'validActuator',wfs.validActuator);

F = 2*bifaLowRes.modes(wfs.validActuator,:);
iF = pinv(full(F),1e-1);

%% INTERACTION MATRIX
wfs.camera.frameListener.Enabled = false;
wfs.slopesListener.Enabled = false;

ngs = ngs.*tel;

calibDm = calibration(dm,wfs,ngs,ngs.wavelength/8,nL+1,'cond',1e2);

%% GENERATE ATMOSPHERE
tel = tel + atm;
figure
imagesc(tel)
ngs = ngs.*tel*wfs;
%% LOW-RES TELESCOPE TO ESTIMATE PHASE IN POLC MODE

telLowRes = telescope(tel.D,'resolution',nL+1,...
    'obstructionRatio',obstructionRatio,'fieldOfViewInArcsec',fieldOfViewInArcsec,'samplingTime',1/samplingFreq);
telLowRes.pupil = wfs.validActuator;

telLowRes= telLowRes + atm;
ngs = ngs.*telLowRes;
phase = ngs.meanRmOpd;

%% LGS SOURCES

nLGS = 6;


% lgsAst = source('asterism',{[nLGS,arcsec(25),0]},'height',90e3,'magnitude',10);
% figure, imagesc(tel, [ngs,lgsAst])


% single altitude object
lgsAst = laserGuideStar(tel.D/nL,tel.D, 90e3, [], 2e6, [],'asterism',{[nLGS,arcsec(25),0]}, 'wavelength', photometry.Na,'height',9e4);
theta = linspace(0,2*pi-2*pi/nLGS,nLGS);
for iGS = 1:nLGS
    lgsAst(iGS).viewPoint = D/2*[cos(theta(iGS)), sin(theta(iGS))];
end


wfs.camera.photonNoise = 1;
wfs.camera.readOutNoise = 10;
wfs.framePixelThreshold = wfs.camera.readOutNoise;

lgsAst_slmmse = slopesLinearMMSE(wfs,tel,atm,lgsAst,'mmseStar',ngs,'NF',1024);


% extended source object
nLgsNaPoints = 5; u = linspace(-1,1,nLgsNaPoints);
profileNa = (exp(-((u-0.35)*5).^2)+0.65*exp(-((u+0.45)*4).^2));
lgsAst = laserGuideStar(tel.D/nL,tel.D, 90e3, [], 2e6, profileNa,'asterism',{[nLGS,arcsec(25),0]}, ...
    'wavelength', photometry.Na,'height',1e3*(linspace(-5,5,nLgsNaPoints)+90));

for kGS=1:nLgsNaPoints
    for iGS = 1:nLGS
        lgsAst(1,iGS,kGS).viewPoint = D/4*[cos(theta(iGS)), sin(theta(iGS))];
    end
end

wfs.camera.readOutNoise = 1;
wfs.framePixelThreshold = wfs.camera.readOutNoise;


% example lgs with elongation
% nLgs = 21; u = linspace(-1,1,nLgs);
% p=exp(-((u-0.35)*5).^2)+0.65*exp(-((u+0.45)*4).^2);
% lgs = laserGuideStar(25/60,25,90e3,1,1e9,p,...
% 'wavelength',photometry.Na,...
% 'height',1e3*(linspace(-5,5,nLgs)+90),...
% 'viewPoint',[-37/2,0]);

%% NGS Sources
nTT = 1;
TTAst = source('asterism',{[nTT,arcsec(10),0]},'wavelength',photometry.H);
%TTAst = source('asterism',{[1,arcsec(0),0]},'wavelength',photometry.H);

TTAst = TTAst .* tel*OIWFS;
%% SCIENCE CAMERA

science = source('wavelength',photometry.K);
cam = imager();
tel = tel - atm;
science = science.*tel*cam;
figure(31416)
imagesc(cam,'parent',subplot(2,1,1))

cam.referenceFrame = cam.frame;
+science;
fprintf('Strehl ratio: %4.1f\n',cam.strehl)

tel = tel + atm;
+science;
fprintf('Strehl ratio: %4.1f\n',cam.strehl)

%% CONTROLLER GAINS

%INTEGRATOR
gain_cl  = 0.5;
uNgs = zeros(dmLowRes.nValidActuator,1);

% PSEUDO OPEN-LOOP CONTROLLER
gain_pol = 0.6;


%% LTAO LOOP 

reset(tel)
dm.coefs = zeros(dm.nValidActuator,1);

% LOOP INIT

science = science.*tel*dm*cam;
lgsAst = lgsAst.*tel*dm*wfs;
TTAst = TTAst.*tel*dm*OIWFS;

figure(31416)
imagesc(cam,'parent',subplot(2,1,1))
subplot(2,1,2)
h = imagesc(catMeanRmPhase(science));
axis xy equal tight
colorbar
wfs.camera.frameListener.Enabled = false;
wfs.slopesListener.Enabled = false;

%
flush(cam)
cam.clockRate    = 1;
exposureTime     = 120;
cam.exposureTime = exposureTime;
startDelay       = 20
flush(cam)
cam.startDelay = startDelay;
set(science,'logging',true)
set(science,'phaseVar',[])
lgsAst_slmmse.wavefrontSize = [dm.nValidActuator,1];
lgsAst_slmmse.warmStart = true;
cam.frameListener.Enabled = true;

%% piston-tip-tilt removal 
t = ones(wfs.nValidLenslet,1);
TT = [t 0*t; 0*t,t];
SlopeTTRem = eye(2*wfs.nValidLenslet) - TT*pinv(TT);

zern = zernike(2:3,'resolution',nL+1, 'pupil',dmLowRes.validActuator);
TT = iF*zern.modes(dmLowRes.validActuator,:);
DMTTRem = eye(dmLowRes.nValidActuator) - TT*pinv(TT);


%% The loop is closed for one full exposure of the science camera.
nIteration = startDelay + exposureTime;
for k=1:cam.startDelay + cam.exposureTime
    tic
    k
    % Objects update
    +tel;
    %+lgsAst;
    lgsAst = lgsAst.*tel*dm*wfs;
    +TTAst;
    +science;
    % Pseudo-open-loop controller
    dm.coefs = (1-gain_pol)*dm.coefs + ...
        gain_pol*iF*( lgsAst_slmmse*( bsxfun( @minus, SlopeTTRem*wfs.slopes, calibDm.D*dm.coefs ) ) );
   
    % Remove TT from DM
    dm.coefs = DMTTRem*dm.coefs;
    
    % Add TT from TT stars, % Closed-loop controller
    uNgs = uNgs - TTAst(1).wavelength/8*gain_cl*TT*mean(OIWFS.slopes,2); 

    % add ngs controls back to dm commands
    dm.coefs = dm.coefs - uNgs;
    
    % Display
         set(h,'Cdata',catMeanRmPhase(science))
         drawnow
    toc
end
imagesc(cam)
set(h,'Cdata',catMeanRmPhase(science))
%% STATISTICAL ANALYSIS

var_wfe_ltao = reshape(science.phaseVar(1:nIteration*2),2,[])';
wfe_ltao = sqrt(var_wfe_ltao)/science.waveNumber*1e6;
marechalStrehl_ltao = 1e2*exp(-mean(var_wfe_ltao(startDelay:end,2)));
psfStrehl_ltao =1e2*cam.strehl



