
% Script to initiate OOMAO setup, initiate objects, calibrate the wfs and
% deformable mirror etc.

%% Initate OOMAO objects
% Atmosphere
atm = atmosphere(lambda_wfs,r0,L0,'altitude',altitude,...
    'fractionnalR0',fractionalR0,'windSpeed',windSpeed,'windDir',windDirection);
% Natural guide star
ngs = source('wavelength',lambda_wfs,'magnitude',magnitude_wfs);
% Science source
sci = source('wavelength',lambda_sci,'magnitude',magnitude_sci);
% Telescope
tel = telescope(D,'resolution',nPx,...
    'fieldOfViewInArcmin',fieldOfView,'samplingTime',samplingTime);
% Diffractive Shack-Hartmann
wfs = shackHartmann(nLenslet,nPx,illRatio);

%% Check WFS and select valid lenslets
% Propagate guide star through system
ngs = ngs.*tel*wfs;

% Initiate: only include wfs signals from valid lenslets (>illRatio)
wfs.INIT;

% Check wfs camera
+wfs;
figure
    imagesc(wfs.camera,'parent',subplot(3,2,[1,4]))
    slopesDisplay(wfs,'parent',subplot(3,2,[5,6]))
    
wfs.camera.frameListener.Enabled = true;
wfs.slopesListener.Enabled = true;

% Geometric Shack-Hartmann (with same pupil, lenslets as diffractive case,
% for accurate comparison)
%wfsG = geomGrad(wfs.validLenslet,tel.pupilLogical,1);
%nSlopes = size(wfsG.gradOp,1);

%% Calibrate WFS slope units (check linear range)
calibrateWFS

%% Define deformable mirror
% Influence function
bif = influenceFunction('monotonic',DMcoef);

figure
    show(bif)

% Define DM
dm = deformableMirror(nLenslet+1,'modes',bif,...
    'resolution',nPx,'validActuator',wfs.validActuator);

% Calibrate the dm
% Propagate ngs through the telescope
ngs = ngs.*tel;
% Calibrate the DM by poking each individual actuator to calibrate the dm
% commands to the signals of the wfs
% Cut-off point ('cond'=100).  I.e. cut-off no. of modes below a given
% threshold.

calibDm = calibration(dm,wfs,ngs,ngs.wavelength,nLenslet+1)
calibDm.nThresholded = 10;

% Reset coefs to 0
dm.coefs(:) = 0;

% Low resolution deformable mirror (for MMSE)
bifLowRes = influenceFunction('monotonic',DMcoef);
dmLowRes = deformableMirror(nLenslet+1,'modes',bifLowRes,...
    'resolution',nLenslet+1,'validActuator',wfs.validActuator);


%% Setup science camera
% Remove atmosphere from telescope path
tel = tel - atm;

% Set up 3 cameras (1 each for low order, high order and total phase)
for n=1:3 
    % Initiate camera
    cam(n) = imager(tel);
    % Propagate science object to imager (no atmosphere)
    sci = sci.*tel*cam(n);
    cam(n).referenceFrame = cam(n).frame;
    flush(cam(n))
    % Camera rate set to wfs sampling time
    cam(n).clockRate = 1;
    % Number of frames for each exposure
    cam(n).exposureTime = exposureTime;
    cam(n).startDelay = startDelay;
    flush(cam(n))
end

%% Set noise on wfs
wfs.camera.photonNoise = photonNoise;
wfs.camera.readOutNoise = readOutNoise;





