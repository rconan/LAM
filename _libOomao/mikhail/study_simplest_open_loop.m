clc
close all
clear all

tool_inside_path_adder('matlab_code/mikhail', 'matlab_code');

samplingTime     = 0.001;% seconds
exposureTime     = 100; %% [msek], equals to the number of iterations
startDelay       = 0;


refWavelength   = 500e-9;
wavelength      = photometry.R;
altitudes       = [0,    5]*1e3;
fractionalR0    = [0.7, 0.3];
windSpeed       = [2,  5];
windDirection   = [pi, pi];


L0 = 40; %% outer scale, [m]
r0 = 0.19; %% Fried parameter, [m]

atm = atmosphere(refWavelength, r0, L0,...
    'fractionnalR0',fractionalR0, ...
    'altitude',     altitudes,...
    'windSpeed',    windSpeed,...
    'windDirection',windDirection);
atm.wavelength = wavelength;
nLayers = atm.nLayer;


%% Telescope
D   = 8;  % telescope diameter [m]
nPx = 120;

tel = telescope(D,...
    'fieldOfViewInArcMin',3.5,...
    'resolution',nPx,...
    'samplingTime',samplingTime);


%% Define the Asterism (star constellation) for the tomography
%%% Don't forget that the FOV of the telescope (meta-pupil) must cover all the guide stars
dAsterism               = 1.0; % 0.5, 0.75, 1.0, 1.5, 2 diameter in arcmin
guideStarWavelength     = photometry.R;
guideStarMagnitude      = 14;%[14,14.5,15,15.5,16,16.5,17]; % lower later

% ngs = source('asterism',{[3,arcmin(dAsterism/2),0]},...
%             'wavelength',guideStarWavelength,... 
%             'magnitude',guideStarMagnitude);

ngs =source; %create a simple on-axis guide star object


%% Wavefront sensor
nLenslet = 10;
minLightRatio = 0.75; %% large values (more than 0.8  can lead to poor reconstruction - you are truncating the phase!).
wfs = shackHartmann(nLenslet,nPx,minLightRatio);  %%% causes problems in MATLAB 2015a

%%% Propagation of the calibration source to the WFS through the telescope
ngs = ngs.*tel*wfs;
wfs.INIT;
+wfs;


%% Deformable mirror
nActuator = nLenslet + 1;
bif = influenceFunction('monotonic',40/100);

dm = deformableMirror(nActuator,...
    'modes',bif,...
    'resolution',nPx,...
    'validActuator',logical(utilities.piston(nActuator)) );

nDm = length(dm);

clear ngs
%% Command Matrix
ngs = source;
ngs.magnitude = 0;
ngs = ngs.*tel;
calibDm = calibration(dm,wfs,ngs,ngs.wavelength,nLenslet+1,'cond',1e2);


sciStars = source('wavelength',photometry.J);
cam = imager(tel);

tel = tel - atm;
sciStars = sciStars.*tel*cam;
figure(31416)
imagesc(cam,'parent',subplot(2,1,1))

cam.referenceFrame = cam.frame;
+sciStars;

tel = tel + atm;
+sciStars;




ngsCombo = source('zenith',zeros(1,2),'azimuth',zeros(1,2),'magnitude',8);
ngsCombo = ngsCombo.*tel*dm*wfs;
sciStarsCombo = source('zenith',zeros(1,1),'azimuth',zeros(1,1),'wavelength',photometry.J);
sciStarsCombo = sciStarsCombo.*tel*dm*cam;


flush(cam)
cam.clockRate    = 1;
cam.exposureTime = exposureTime;
cam.startDelay   = startDelay;

figure(31416)
imagesc(cam,'parent',subplot(2,1,1))
cam.frameListener.Enabled = true;
subplot(2,1,2)
h = imagesc(catMeanRmPhase(sciStarsCombo));
axis xy equal tight
colorbar

% The closed--loop controller is a simple integrator where the DM
gain_cl  = 0.5;


dm.coefs = zeros(dm.nValidActuator,1);
flush(cam)


nIteration = startDelay + exposureTime;
for kIteration=1:nIteration

%     fprintf('\n \n   Iteration %d out of %d \n\n ',kIteration,nIteration)

    % Objects update
    +tel;
    +ngsCombo;
    +sciStarsCombo;
    % Closed-loop controller
    dm.coefs(:,1) = dm.coefs(:,1) - gain_cl*calibDm.M*wfs.slopes(:,1);    
    % Display
%     set(h,'Cdata',catMeanRmPhase(sciStarsCombo))
%     drawnow
end
imagesc(cam)
set(h,'Cdata',catMeanRmPhase(sciStarsCombo))

cam.strehl
