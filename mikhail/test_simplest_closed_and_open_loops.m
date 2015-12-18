clc
close all
clear all
format short

tool_inside_path_adder('matlab_code/mikhail', 'matlab_code');

samplingTime     = 0.001;% seconds
frameTime          = 0.001;%[0.01,0.015,0.020,0.025]; seconds
exposureTime     = 100; %% [msek], equals to the Number of Iterations
startDelay           = 0; %% skip this number of frames

fixedLagTimeInMs = 0;

params = struct; %%% creating the Structure with additional parameters for the computeARmodel
params.show_figures = 1; %% display the diagnostic plots during the simulation.

% params.LoopType= 'open loop'; 
params.LoopType= 'closed loop'; 

% params.noise_in_wfs = 0;
params.noise_in_wfs = 1;

% randn('state', 25);     % sets the global random state
randn('state', sum(100*clock)) %Initialize randn to a different state each time


%% Atmospheric parameters
refWavelength   = 500e-9;
wavelength      = photometry.R;
altitudes       = [0,    5]*1e3;
fractionalR0    = [0.7, 0.3];
windSpeed       = [2,  5];
windDirection   = [pi, pi];

L0 = 40; %% outer scale, [m]
r0 = 0.19; %% Fried parameter, [m]

%% Creating the Atmosphere object
atm = atmosphere(refWavelength, r0, L0,...
    'fractionnalR0',fractionalR0, ...
    'altitude',     altitudes,...
    'windSpeed',    windSpeed,...
    'windDirection',windDirection);
atm.wavelength = wavelength;
nLayers = atm.nLayer;


%% Wavefront sensor parameters
nLenslet   = 10; % number of lenslets
nPx  = 12; % pixels per lenslet
nRes = nLenslet*nPx; %% total number of pixels.
minLightRatio = 0.85; %% large values (more than 0.8  can lead to poor reconstruction - you are truncating the phase!).


%% Telescope
D    = 8;  % telescope diameter [m]
d    = D/nLenslet; % lenslet pitch

tel = telescope(D,...
    'fieldOfViewInArcMin',3,...
    'resolution',nRes,...
    'samplingTime',samplingTime);



%% GuideStar Asterism (star constellation) for the tomography
%%% Don't forget that the FOV of the telescope (meta-pupil) must cover all the guide stars
dAsterism                       = 0.0; % 0.5, 0.75, 1.0, 1.5, 2 diameter in arcmin
guideStarWavelength   = photometry.R;
guideStarMagnitude      = 14;%[14,14.5,15,15.5,16,16.5,17]; % lower later

%%% Guide Stars with Regular asterism
% astNGS = source('asterism',{[3,arcmin(dAsterism/2),0]},...
%             'wavelength',guideStarWavelength,... 
%             'magnitude',guideStarMagnitude);
% 
% %%%% Only one guidestar:
% astNGS = source('asterism',{[1,arcmin(dAsterism/2),0]},...
%             'wavelength',guideStarWavelength,... 
%             'magnitude',guideStarMagnitude);

%% This is LGS for the test
        astNGS = source('asterism',{[1,arcmin(dAsterism/2),0]},'height', 90e3);

nGs = length(astNGS);


%% Science Targets
scienceObjectMag = 10;  %% Apparent magnitude
scienceObjectWavelength = photometry.H; %% Infrared
sciZenithVector = arcsec(0);
sciAzimuthVector = 0;

sciStars = source('zenith',[sciZenithVector],'azimuth',[sciAzimuthVector],...
    'wavelength',scienceObjectWavelength,'magnitude',scienceObjectMag);


%% Wavefront sensor Calibration
ngs = source; %create a simple on-axis guide star object for WFS calibration

wfs = shackHartmann(nLenslet,nRes,minLightRatio);  %%% causes problems in MATLAB 2015a

%%% Propagation of the calibration source to the WFS through the telescope
ngs = ngs.*tel*wfs;
wfs.INIT;
% When the WFS detector is read again and a new processing of the frame is done, there is no more warning as only the
+wfs;
figure
imagesc(wfs.camera,'parent',subplot(3,2,[1,4]))
slopesDisplay(wfs,'parent',subplot(3,2,[5,6]))

% The next 2 commands allow the displays of the frame and of the slopes to
% be updated when a new frame and new slopes are computed
wfs.camera.frameListener.Enabled = true;
wfs.slopesListener.Enabled = true;

% The WFS must be calibrated such as for 1rad of tip--tilt wavefront , it will
% measured a slopes of 1rad.
% To do so, the \oop{shackHartmann}{pointingDirection} property is set on-axis
% </latex>
wfs.pointingDirection = zeros(2,1);


%%
% <latex>
% whereas the source is progressively moved off-axis
% </latex>
pixelScale = ngs.wavelength/(2*d*wfs.lenslets.nyquistSampling);
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

figure, plot(Ox_in,Ox_out), grid;

slopesLinCoef = polyfit(Ox_in,Ox_out,1);
wfs.slopesUnits = 1/slopesLinCoef(1);

% % break
% The source is reset on--axis and the WFS is set to always be aligned to
% the source by setting \oop{shackHartmann}{pointingDirection} to empty.
ngs.zenith = 0;
wfs.pointingDirection = [];

wfs.camera.frameListener.Enabled = false;
wfs.slopesListener.Enabled = false;



%% Deformable mirror
nActuator = nLenslet + 1;
bif = influenceFunction('monotonic',40/100);

dm = deformableMirror(nActuator,...
    'modes',bif,...
    'resolution',nRes,...
    'validActuator',wfs.validActuator );

nDm = length(dm);




clear ngs  %% there is no need for the calibration source any more.

%% Command Matrix
ngs = source;
ngs.magnitude = 0;
ngs = ngs.*tel;
calibDm = calibration(dm,wfs,ngs,ngs.wavelength,nLenslet+1,'cond',1e2);


%% Imager / Science Camera
cam = imager(tel);

tel = tel - atm;
sciStars = sciStars.*tel*cam;
figure(31416)
imagesc(cam,'parent',subplot(2,1,1))

cam.referenceFrame = cam.frame;
+sciStars;

tel = tel + atm;
+sciStars;

flush(cam)
cam.clockRate    = 1;
cam.exposureTime = exposureTime;
cam.startDelay   = startDelay;

figure(31416)
imagesc(cam,'parent',subplot(2,1,1))
cam.frameListener.Enabled = true;
subplot(2,1,2)
h = imagesc(catMeanRmPhase(sciStars.*tel*cam));
axis xy equal tight
colorbar


% The closed--loop controller is a simple integrator where the DM
gain_cl  = 0.5;

nIteration = startDelay + exposureTime;

dm.coefs = zeros(dm.nValidActuator,1);
coefsVec = zeros(dm.nValidActuator,1,nIteration);

flush(cam)

%% Noisy or Noiseless case:
                if (params.noise_in_wfs == 0) 
                        wfs.camera.photonNoise = false;
                        wfs.camera.readOutNoise = 0;  %% turning this noise to zero may fail the phaseRecosntruction!!!
                        +wfs
                else
                       wfs.camera.photonNoise = true;
                       wfs.camera.readOutNoise = 0.2;  %% turning this noise to zero may fail the phaseRecosntruction!!!
                       wfs.framePixelThreshold = 0.2; %% increase to reduce the WF error
                        +wfs
                end

% break 

switch params.LoopType

    case 'closed loop'
        
                for kIteration=1:nIteration
             %     fprintf('\n \n   Iteration %d out of %d \n\n ',kIteration,nIteration)

                    % Objects update
                    +tel;  % Always advance atm by 1ms

                %     +ngsCombo;
                    astNGS = astNGS.*tel*dm*wfs;

                %     +sciStarsCombo;
                    sciStars = sciStars.*tel*dm*cam;

                %%%    Closed-loop controller
                    dm.coefs(:,1) = dm.coefs(:,1) - gain_cl*calibDm.M*wfs.slopes(:,1);    

                    % Display the progress frame-by-frame
                    set(h,'Cdata',catMeanRmPhase(astNGS))
                    drawnow
                end

                
    case 'open loop'

                  for kIteration=1:nIteration
                    % Objects update
                    +tel;  % Always advance atm by 1ms

                %     +ngsCombo;
                    astNGS = astNGS.*tel*wfs;

                %     +sciStarsCombo;
                    sciStars = sciStars.*tel*dm*cam;

                                if kIteration == 1
                                    dm.coefs =  -coefsVec(:,1,kIteration);
                                else
                                    dm.coefs =  -coefsVec(:,1,kIteration-1);
                                end

            
                %%%    Closed-loop controller
                    coefsVec(:,1,kIteration) = calibDm.M*wfs.slopes(:,1); 
                    
                    % Display the progress frame-by-frame
                    set(h,'Cdata',catMeanRmPhase(astNGS))
                    drawnow
                end
        
end %%% for switch case.
% imagesc(cam)
% set(h,'Cdata',catMeanRmPhase(astNGS))

cam.strehl


%% 3. \phi_z^{res} = Z^\dagger * \phi_{pix}^{res}
zPhase = reshape(sciStars.meanRmPhase, nRes*nRes,1);
zphi_res_modes = pinv(full(dm.modes.modes))*zPhase;
figure, bar(zphi_res_modes);