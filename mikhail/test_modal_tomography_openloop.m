clc
close all
clear all

tool_inside_path_adder('matlab_code/mikhail', 'matlab_code');

samplingTime     = 0.001;% seconds
frameTime          = 0.001;%[0.01,0.015,0.020,0.025]; seconds
exposureTime     = 300; %% [msek], equals to the Number of Iterations
startDelay           = 0; %% skip this number of frames

fixedLagTimeInMs = 1;

params = struct; %%% creating the Structure with additional parameters for the computeARmodel
    params.show_figures = 0; %% display the diagnostic plots during the simulation.
    params.ircsSlitWidth           = 0.14; %%% % arcsec  EEwidth for the EE estimation in Imager

    %%% Noisy case:
    params.noise_in_wfs  = 1; %% display the diagnostic plots during the simulation.
    params.noise.wfs_readOutNoise = 0.2;
    params.noise.photonNoise = true;
    params.noise.framePixelThreshold = 0;

    % No Noise  (noiseless) case:
%     params.noise_in_wfs  = 0; %% display the diagnostic plots during the simulation.
%     params.noise.wfs_readOutNoise = 0;
%     params.noise.photonNoise = false;
%     params.noise.framePixelThreshold = 0;


% randn('state', 25);     % sets the global random state
% randn('state', sum(100*clock)) %Initialize randn to a different state each time


%% Atmospheric parameters
refWavelength   = 500e-9;
wavelength      = photometry.R;
altitudes       = [0,    5]*1e3;
fractionalR0    = [0.7, 0.3];
windSpeed       = [2,  5];
windDirection   = [pi, pi];

% refWavelength   = 500e-9;
% wavelength      = photometry.R;
% altitudes = [0,5.5,11.5]*1e3;
% fractionalR0 = [0.596,0.223816,0.180184];
% windSpeed = [2.84*2,6,17];
% windDirection = [pi/2, pi, pi];

L0 = 40; %% outer scale, [m]
r0 = 0.19; %% Fried parameter, [m]

% Creating the Atmosphere object
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
    'fieldOfViewInArcMin',3.5,...
    'resolution',nRes,...
    'samplingTime',samplingTime);


%% Zernike modes:
% maxRadialDegree     = 3;  %% - setting the Max Radial Degree = how many modes we will use
% maxRadialDegreeProj = 10; %% on how many modes we will reject

% maxRadialDegree     = 6;  %% - setting the Max Radial Degree = how many modes we will use
% maxRadialDegreeProj = 20; %% on how many modes we will reject

%%%%% Raven settings
% maxRadialDegree     = 9;  %% - setting the Max Radial Degree = how many modes we will use
% maxRadialDegreeProj = 23; %% on how many modes we will reject

maxRadialDegree     = 9;  %% - setting the Max Radial Degree = how many modes we will use
maxRadialDegreeProj = 30; %% on how many modes we will reject

%% GuideStar Asterism (star constellation) for the tomography
%%% Don't forget that the FOV of the telescope (meta-pupil) must cover all the guide stars
dAsterism               = 1.0; % currently 0.1, 0.5, 1.0 diameter in arcmin
guideStarWavelength     = photometry.R;
guideStarMagnitude      = 14;%[14,14.5,15,15.5,16,16.5,17]; % lower later

%% Guide Stars with Regular asterism
astNGS = source('asterism',{[3,arcmin(dAsterism/2),0]},...
            'wavelength',guideStarWavelength,... 
            'magnitude',guideStarMagnitude);

% % %%% Only one guidestar:
% astNGS = source('asterism',{[1,arcmin(dAsterism/2),0]},...
%             'wavelength',guideStarWavelength,... 
%             'magnitude',guideStarMagnitude);
% 
% % % Asterism based on Raven CU pinholes
% % astNGS = source('zenith',[arcsec(45.3),arcsec(35.5),arcsec(45.3)],...
% %     'azimuth',[3.0309,0.7854,-1.6815],...
% %     'wavelength',guideStarWavelength,'magnitude',guideStarMagnitude);
% 
nGs = length(astNGS);



%% Science Targets
scienceObjectMag = 10;  %% Apparent magnitude
scienceObjectWavelength = photometry.H; %% Infrared
sciZenithVector = arcsec(0);
sciAzimuthVector = 0;

sciStars = source('zenith',[sciZenithVector],'azimuth',[sciAzimuthVector],...
    'wavelength',scienceObjectWavelength,'magnitude',scienceObjectMag);

 nSciobj = length(sciStars);
 

%% Wavefront sensors
ngs = source; %% make a calibration source

wfs = shackHartmann(nLenslet,nRes,minLightRatio );
ngs = ngs.*tel*wfs;

wfs.INIT
+wfs;

% The next 2 commands allow the displays of the frame and of the slopes to
% be updated when a new frame and new slopes are computed
wfs.camera.frameListener.Enabled = true;
wfs.slopesListener.Enabled = true;

% The WFS must be calibrated such as for 1rd of tip--tilt wavefront , it will
% measured a slopes of 1rd.
wfs.pointingDirection = zeros(2,1);  %%% DON'T COMMENT - important for calibration

% whereas the source is progressively moved off-axis
pixelScale = ngs.wavelength/(2*d*wfs.lenslets.nyquistSampling);
tipStep = pixelScale/2;
nStep   = floor(nPx/3)*2;
sx      = zeros(1,nStep+1);
u       = 0:nStep;

wfs.camera.frameListener.Enabled = false;
wfs.slopesListener.Enabled = false;

warning('off','oomao:shackHartmann:relay')
%%% Begin calibration on the tilt angle 
            for kStep=u
                ngs.zenith = -tipStep*kStep;  %%% changing zenith - check the NGS name!
                +ngs;
                drawnow
                sx(kStep+1) = median(wfs.slopes(1:end/2));
            end
%%% Begin calibration on the tilt angle 
warning('on','oomao:shackHartmann:relay')

Ox_in  = u*tipStep*constants.radian2arcsec;
Ox_out = sx*ngs.wavelength/d/2*constants.radian2arcsec;

% figure, plot(Ox_in,Ox_out), grid

slopesLinCoef = polyfit(Ox_in,Ox_out,1);
wfs.slopesUnits = 1/slopesLinCoef(1);

% The source is reset on--axis and the WFS is set to always be aligned to
% the source by setting \oop{shackHartmann}{pointingDirection} to empty.
wfs.pointingDirection = []; %%% DO NOT UNCOMMENT IT!
clear ngs





%% Obtaining the Noise covariance matrix and Dzcc:
ngs = source; %create a simple on-axis guide star object for WFS calibration


%% Zernike coefs to WFS slopes calibration
zernModeMax = zernike.nModeFromRadialOrder(maxRadialDegree);
% zern = zernike(2:zernModeMax,'resolution',nRes,'pupil',tel.pupil);  %%%% this increases the  diff_Covar = trace(Sigma_z) - trace(Sigma_alpha)
zern = zernike(2:zernModeMax,'resolution',nRes,'pupil',tel.pupil,'D',tel.D);

zern.c = eye(zern.nMode)*ngs.wavelength/2/pi; %% NORMALISING zernike on wavenumber
ngs=ngs.*zern*wfs;  %%% propagating the ngs to WFS via Zernike
Dzcc = wfs.slopes; %% Matrix  \Gamma of gradients that the WFS can sense (or matrix G)
clear ngs


%% Getting the Noise Covariance matrix - With noise
+wfs
ngs = source; %create a simple on-axis guide star object for WFS calibration

ngs.wavelength = photometry.R;          %%% This is to match the Guidestar wavelength
ngs.magnitude = guideStarMagnitude; %%% This is to match the Guidestar wavelength

% wfs.camera.readOutNoise = params.noise.wfs_readOutNoise;
% wfs.camera.photonNoise = params.noise.photonNoise;
% wfs.framePixelThreshold = params.noise.framePixelThreshold; %% increase to reduce the WF error

wfs.camera.readOutNoise = 0.2;
wfs.camera.photonNoise = true;
wfs.framePixelThreshold = 0; %% increase to reduce the WF error

ngs=ngs.*tel*wfs;

%% noise convariance matrix
fprintf('\n\n Initialising the noisy WFS....\n')

nMeas = 1000;
slopes = zeros(wfs.nSlope,nMeas);
for kMeas=1:nMeas
    +wfs
    slopes(:,kMeas) = wfs.slopes;
end
Cn = slopes*slopes'/nMeas;


if (params.show_figures == 1)
    figure, imagesc(Cn), title('Noise covariance matrix');
    figure, plot(diag(Cn)), title('Noise covariance - diagonal elements only');
end

fprintf('... DONE! \n\n ')
clear slopes ngs

wfs.camera.readOutNoise = 0;
wfs.camera.photonNoise = false;
wfs.framePixelThreshold = 0; %% increase to reduce the WF error
+wfs


%% Constructing matrices for the Rmv reconstruction
CnAst = tool_make_blockdiag_matrix ('Cn', nGs, Cn);  %%% create the block-diag matrix CnAst from nGs matrix of Cn
DAst = tool_make_blockdiag_matrix ('Dzcc', nGs, Dzcc);  %%% create the block-diag matrix CnAst from nGs matrix of Cn


%% Deformable mirror
nActuator = nLenslet + 1;
bif = influenceFunction('monotonic',40/100);

dm = deformableMirror(nActuator,...
    'modes',bif,...
    'resolution',nRes,...
    'validActuator',wfs.validActuator );
nDm = length(dm);


%% Command Matrix
ngs = source
wfsCal = wfs; %%% Pick up a wavefront sensor for Calibation of DM %%% MAYBE ISSUE 10!!!

bif = influenceFunction('monotonic',40/100);
dmCal = deformableMirror(nActuator,...
    'modes',bif,...
    'resolution',nRes,...
    'validActuator', wfs.validActuator);

ngs = source;
ngs.magnitude = 0;
ngs = ngs.*tel;

calibDm = calibration(dmCal,wfsCal,ngs,ngs.wavelength,nLenslet+1, 'cond', 1e1); %%% from modalMethods. cond must be small. Produces better results since it effectively doubles the commands
%%% COND parameter make a lot of difference!!

F = dm.modes.modes; %% matrix of DM modes



%% Computing the Dz matrix that 
fprintf('\n\n Computing the Z2U matrix ....\n')
    Z2U = calibDm.M*Dzcc;% Using zernike derivatives and slopes to DM command matrix
fprintf('... DONE! \n\n ')
%%% calibDm.M must be approx Z2U * Dzcc


%% Imager / Science Camera
tel = tel - atm;
cam = imager(tel);

sciStars = sciStars.*tel*cam;
% figure(31416)
% imagesc(cam,'parent',subplot(2,1,1))

cam.referenceFrame = cam.frame;
+sciStars;

tel = tel + atm;
+sciStars;

flush(cam)
cam.clockRate    = 1;
cam.exposureTime = exposureTime;
cam.startDelay   = startDelay;
cam.eeWidth = params.ircsSlitWidth; % arcsec

figure(31416)
imagesc(cam,'parent',subplot(2,1,1))
cam.frameListener.Enabled = true;
subplot(2,1,2)
h = imagesc(catMeanRmPhase(sciStars.*tel*cam));
axis xy equal tight
colorbar

                if (params.noise_in_wfs == 0) 
                    wfs.camera.readOutNoise = 0;
                    wfs.camera.photonNoise = false;
                    wfs.framePixelThreshold = 0; %% increase to reduce the WF error
                    +wfs
                else 
wfs.camera.readOutNoise = params.noise.wfs_readOutNoise;
wfs.camera.photonNoise = params.noise.photonNoise;
wfs.framePixelThreshold = params.noise.framePixelThreshold; %% increase to reduce the WF error
                    +wfs
                end

%%
%% BEGIN::::: COMPUTING MATRICES for TOMOGRAPHY ::::::::::::::::
%%

%% Create extended Zernike object for Tomographic projections
zernModeMaxProj = zernike.nModeFromRadialOrder(maxRadialDegreeProj);
% zernProj = zernike(2:zernModeMaxProj,'resolution',nRes,'pupil',tel.pupil);
zernProj = zernike(2:zernModeMaxProj,'resolution',nRes,'pupil',tel.pupil,'D',tel.D);

projFileName = strcat('Projections_GS', num2str(nGs), '_Asterism',  strrep(num2str(dAsterism), '.', '_'), '_AtmLayers', num2str(nLayers), '_MaxRadDeg',num2str(maxRadialDegree), '_MaxRadDegProj',num2str(maxRadialDegreeProj), '.mat');

if exist(projFileName,'file')

    fprintf('\n\n Loading precomputed Projection_beta and Projection_alpha ....\n')
        load(projFileName)
    fprintf('... DONE! \n\n ')

else
%     [Projection_alpha, Palphacell] = AnalyticalSmallFootprintExpansion(zernProj,tel,astNGS,atm);
%     [Projection_beta, Pbetacell]   = AnalyticalSmallFootprintExpansion(zernProj,tel,sciStars,atm);
% %     save(projFileName,'Pbeta','Palpha')

    %% Zernike projection for atmosphere layers
    fprintf('\n\n Computing Projection_alpha....\n')

    altitudes = [atm.layer.altitude];
    Projection_alpha_Cell = tool_analytical_small_footprint_expansion_with_piston(zernProj, altitudes, astNGS, tel, nLayers, nGs);
    Projection_alpha = cell2mat(Projection_alpha_Cell);
    % % %     [Palpha2, Palphacell2] = AnalyticalSmallFootprintExpansion(zernProj,tel,ast,atm); %tested, coincide with tool_analytical_small_
    fprintf('... DONE! \n\n ')


    %% Modal Projection of the science object 
    fprintf('\n\n Computing Projection_beta....\n')

    altitudes = [atm.layer.altitude];
    Projection_beta_Cell = tool_analytical_small_footprint_expansion_with_piston(zernProj, altitudes, sciStars, tel, nLayers, nSciobj);
    Projection_beta = cell2mat(Projection_beta_Cell);
    % % %     [Pbeta2, Pbetacell2] = AnalyticalSmallFootprintExpansion(zernProj,tel,sciCal,atm); %tested, coincide with tool_analytical_small_
    fprintf('... DONE! \n\n ')

    save(projFileName,'Projection_alpha','Projection_beta')
end


%% Computing the cross-correlation matrix of Zernike coefficients from the modes of the atmosphere object for EACH ATM LAYER
Sigma_phi_Cell = cell(nLayers,nLayers);  %%% composite cell array for Atm Covariance

% zernWfsi = zernike(2:zernModeMaxProj,'resolution',nRes,'pupil',tel.pupil); 
zernWfsi = zernike(2:zernModeMaxProj,'resolution',nRes,'pupil',tel.pupil,'D',tel.D);

for kLayer = 1:nLayers
    for kLayer1 = 1:nLayers
        
        if kLayer == kLayer1
            zernWfsi.D = tel.diameterAt( atm.layer(kLayer).altitude ); %% re-normalize the teelscope's diameter at each altitude.
            %%% Compute cross-covariance between layers.
            Sigma_phi_Cell{kLayer,kLayer1} = ...
                phaseStats.zernikeCovariance(zernWfsi,atm)*atm.layer(kLayer).fractionnalR0; %% %% this is from zernike.m
            
        else
            Sigma_phi_Cell{kLayer,kLayer1} = zeros(zernWfsi.nMode);
        end

    end %%% for kLayer1 = 1:nLayers
end  %%% for kLayer = 1:nLayers

Sigma_phi = cell2mat(Sigma_phi_Cell); %% matrix for the Atmosphere Covariance



%% This vector of indices is 1x(zernModeMax-1), sans piston (removing the piston mode)
idy = [];
for kL=1:nLayers
    ii =  2:zernModeMaxProj;
    idy = [idy ii+zernModeMaxProj*(kL-1)];
end

idx = [];
for kGs=1:nGs
    ii = 1:zernModeMax-1;
    idx = [idx ii+(zernModeMaxProj-1)*(kGs-1)+kGs];
end


ids = [];
for kL=1:nLayers
    ii = (1:zernModeMaxProj-1);
    ids = [ids ii+(zernWfsi.nMode)*(kL-1)];
end


%% This vector of indices is 1x(zernModeMax-1), sans piston
idp = [];
for kSci=1:nSciobj
    ii = 1:zernModeMax-1;
    idp = [idp ii+(zernModeMax-1)*(kSci-1)+kSci];
end

Sigma_alpha = Projection_alpha(idx,idy)*Sigma_phi(ids,ids)*Projection_alpha(idx,idy)'; %% Computing the projection of the covariance in the direction of Guidestars
Sigma_betaalpha = Projection_beta(idp,idy)*Sigma_phi*Projection_alpha(idx,idy)'; %% Computing the Projection for the sci object:

fprintf('\n PROJECTIONS DONE\n \n')
%% End calculating the projections


%%% This is covariance of the atmosphere modes (will be Palpha * Sigma_Phi * Palpha^T in the tomography case)
Sigma_z = zernikeStats.covariance(zern,atm); %|sigma_Z

diff_Covar = trace(Sigma_z) - trace(Sigma_alpha); %% checking the covariance difference between the initial Zernike modes and the projected (Sigma_\alpha == P_\alpha * Sigma_phi_z P_alpha^T)
diff_Covar_relative = (trace(Sigma_z) - trace(Sigma_alpha))/trace(Sigma_z); %% checking the covariance difference between the initial Zernike modes and the projected (Sigma_\alpha == P_\alpha * Sigma_phi_z P_alpha^T)





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ---- PREDICTION ALGORITHMS ---
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Reconstruction methods
% switch reconstruction_method
% 
%          case 'Static MMSE'
% Minimum Variance reconstructor (static MMSE):
    Rmv = (Sigma_betaalpha*DAst')/(DAst*Sigma_alpha*DAst' + CnAst);

    Rmv_notomo = (Sigma_z*Dzcc')/(Dzcc*Sigma_z*Dzcc' + Cn);

%%%% FIXME - potential numerical problem here, consider using Rmv2 = (Sigma_betaalpha* DAst')*pinv(DAst*Sigma_alpha*DAst' + CnAst);


% end

nIteration = startDelay + exposureTime;

coefsVec = zeros(dm.nValidActuator,1,nIteration);
dm.coefs = zeros(dm.nValidActuator,1);

flush(cam); %% clean the science camera.

nSlope = wfs(1).nSlope;
slopesStack = zeros(nSlope,nGs);
StackSlopeLarge = zeros(nSlope*nGs,1);

if (params.show_figures == 1)
    figure, imagesc(tel,[astNGS,sciStars]), title('Guidestar Asterism visualisation');
end



% invDzcc = pinv(Dzcc); %% pseudoinverse of the matrix of WFS modes
invDzcc = pinv(DAst); %% pseudoinverse of the matrix of WFS modes
Fz =  0.5*pinv(full(F))*zern.modes*ngs.wavelength/(2*pi);  %% 0.5 is because of the reflective DM (doubles the phase). BUT IT DOES NOT WORK IN THE NOISY CASE!

%% Computing the Pseudoinverse via SVD
svd_threshold = 5;
[U,S,V] = svd( full(F),0); %% % F_back = U*S*V';   economy-size SVD
S( abs(S) < svd_threshold ) = 0;
S_new = S;

%%%%% Now building an inverse matrix
Fz_svd = (V*pinv(S_new)*U' ) * zern.modes*ngs.wavelength/(2*pi); %% DON'T MULTIPLY IT ON 0.5!!!! when you do noisy simulations


% break

%% #########   Begin the Loop %%%%
for kIteration=1:nIteration

%       fprintf('\n \n   Iteration %d out of %d \n\n ',kIteration,nIteration)
% fprintf('Iteration %d  \n',kIteration)
      
    for iFrame = 1:frameTime/samplingTime
    % Objects update
    +tel;  % Always advance atm by 1ms

%% OPEN LOOP   
        astNGS = astNGS.*tel*wfs;
        
        sciStars = sciStars.*tel*dm*cam;
    
        if iFrame == fixedLagTimeInMs
            if kIteration == 1
                dm.coefs =  -coefsVec(:,1,kIteration);
            else
                dm.coefs =  -coefsVec(:,1,kIteration-1);
            end
        end %% for iFrame
            
       end %%% END generating the frames

       
%% Stacking the slopes in a single vector    
if nGs ==1      %%% This is for one GuideStar    
    StackSlopeLarge(1:nSlope,1) = wfs.slopes;
end

if nGs ==3 % This is for 3 NGSes
        slopesStack = wfs.slopes;
        StackSlopeLarge(1:nSlope,1) = slopesStack(:,1);
        StackSlopeLarge(nSlope+1:nSlope*2,1) = slopesStack(:,2);
        StackSlopeLarge(2*nSlope+1:nSlope*3,1) = slopesStack(:,3);    
end

%% Reconstruction methods
%         coefsVec(:,1,kIteration) = calibDm.M*StackSlopeLarge(:);    %%% direct method, via influence functions
%         coefsVec(:,1,kIteration) = 2*calibDm.M*StackSlopeLarge(:);    %%% direct method, via influence functions

%         coefsVec(:,1,kIteration) = Z2U*invDzcc*StackSlopeLarge(:);  %% least squares


        coefsVec(:,1,kIteration) = 2*Z2U*Rmv*StackSlopeLarge(:);  %%% this works for the NOISY case
%         coefsVec(:,1,kIteration) = Z2U*Rmv*StackSlopeLarge(:);  %%% this works for the NOISE-LESS case

%         coefsVec(:,1,kIteration) = Z2U*Rmv_notomo*StackSlopeLarge(:);


%         coefsVec(:,1,kIteration) =  -Fz*invDzcc*StackSlopeLarge(:);  %%% Modal
%         coefsVec(:,1,kIteration) =  -Fz_svd*invDzcc*StackSlopeLarge(:);  %%% Modal


% % % % %     end
    
% % %     Display the progress frame-by-frame
    set(h,'Cdata',catMeanRmPhase(sciStars))
    drawnow
   
end

cam.strehl

%% 3. \phi_z^{res} = Z^\dagger * \phi_{pix}^{res}
zphi_res_modes = pinv(zern.p)*reshape(sciStars.meanRmPhase, nRes*nRes,1);
figure, bar(zphi_res_modes);