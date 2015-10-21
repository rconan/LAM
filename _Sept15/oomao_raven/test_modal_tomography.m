clc
close all
clear all
% warning off
%%% using the modalTomographyHowto.m example for setting up 3NGSs
%%% Based on demo_modal_tomography in OOMAO/mikhail.

addpath('./OOMAO-Raven/OOMAOlibUpdated/')  % Only use OOMAOlibUpdated!
addpath('./OOMAO-Raven/')

format short

params = struct; %%% creating the Structure with additional parameters for the computeARmodel
params.show_figures = 1; %% display the diagnostic plots during the simulation.

% randn('state', 25);     % sets the global random state
randn('state', sum(100*clock)) %Initialize randn to a different state each time

nIteration       = 100; %% Number of frames to simulate and take

samplingTime     = 0.001;% seconds
frameTime        = 0.01;%[0.01,0.015,0.020,0.025]; seconds
exposureTime     = nIteration*frameTime;
fixedLagTimeInMs = 0;

ircsSlitWidth    = 0.14;

reconstruction_method = 'Static MMSE';
% reconstruction_method = 'AR1';
% reconstruction_method = 'AmmseLag';


RecVec = [1]; %%% will be more methods
number_of_methods2copmare = length(RecVec); %% add different methods like MMSE and Ammse [TBD, reserved but not currently used]

%% \subsubsection{The atmosphere}
% refWavelength   = 500e-9;
% wavelength      = photometry.R;
% altitudes       = [0,    5,   10]*1e3;
% fractionalR0    = [0.5, 0.3, 0.2];
% windSpeed       = [10,  5,   5];
% windDirection   = [pi/2, pi, pi];

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
nLayers = atm.nLayer;

%% Zernike modes:
maxRadialDegree     = 3;  %% - setting the Max Radial Degree = how many modes we will use
maxRadialDegreeProj = 10; %% on how many modes we will roject


%% Define the Asterism (star constellation) for the tomography
%%% Don't forget that the FOV of the telescope (meta-pupil) must cover all the guide stars
dAsterism               = 1.0; % 0.5, 0.75, 1.0, 1.5, 2 diameter in arcmin
guideStarWavelength     = photometry.R;
guideStarMagnitude      = 14;%[14,14.5,15,15.5,16,16.5,17]; % lower later

%%% Guide Stars with Regular asterism
ast = source('asterism',{[3,arcmin(dAsterism/2),0]},...
            'wavelength',guideStarWavelength,... 
            'magnitude',guideStarMagnitude);

% ast = source('asterism',{[1,arcmin(0),0]},...
%             'wavelength',guideStarWavelength,... 
%             'magnitude',guideStarMagnitude);

nGs = length(ast);


%% Science Targets
scienceObjectMag = 10;  %% Apparent magnitude
scienceObjectWavelength = photometry.H; %% Infrared
sciZenithVector = arcsec(0);
sciAzimuthVector = 0;

sciCal = source('zenith',0*cougarConstants.arcsec2radian,'wavelength',photometry.R); %%% this is probably Projection Disk (like the telescope pupil)
nSciobj = length(sciCal);



%% Adding science objects for different reconstruction methods
sciObjectVec = cell(1,1); %% if you want to compare two or more methods - increase the number of cells, like (1,N)

%%%%% MMSE reconstructor
    sciStatic = source('zenith',[sciZenithVector(1)],'azimuth',[sciZenithVector(1)],...
        'wavelength',scienceObjectWavelength,'magnitude',scienceObjectMag);
    sciObjectVec{1} = sciStatic;
    sciObjectVec{2} = sciStatic;


%% Telescope
D   = 8;  % telescope diameter [m]
nPx = 120;

tel = telescope(D,...
    'fieldOfViewInArcMin',3,...
    'resolution',nPx,...
    'samplingTime',samplingTime);


%% Wavefront sensor
nLenslet = 10;
minLightRatio = 0.75; %% large values (more than 0.8  can lead to poor reconstruction - you are truncating the phase!).
wfs = shackHartmann(nLenslet,nPx,minLightRatio);  %%% causes problems in MATLAB 2015a

%%% Propagation of the calibration source to the WFS through the telescope
ngs = source;
ngs = ngs.*tel*wfs;

wfs.INIT;
wfs.camera.photonNoise = false;
%%%%% % out = piston(Npx) Computes a piston on Npx pixel across the diameter
wfs.validLenslet = utilities.piston(nLenslet); % used to be 10, probably for nLenslet ??????????????
+wfs;

% % % Saving the turbulence aberrated phase
% turbPhase = ngs.meanRmPhase; %% this is a turbulent phase


%% Zernike coefs to WFS slopes calibration
% From modal interaction with model of WFS
zernModeMax = zernike.nModeFromRadialOrder(maxRadialDegree);
% zern = zernike(2:zernModeMax,'resolution',nPx,'pupil',tel.pupil);  %%%% this increases the  diff_Covar = trace(Sigma_z) - trace(Sigma_alpha)
zern = zernike(2:zernModeMax,'resolution',nPx,'pupil',tel.pupil,'D',tel.D);
zern.lex = false;
zern.c = eye(zern.nMode);

    ngs=ngs.*zern*wfs;  %%% propagating the ngs to WFS via Zernike

Dzcc = wfs.slopes; %% Matrix \Gamma of gradients
zern.lex = true; %% WHY?


%% With noise
ngs.wavelength = photometry.R;
ngs.magnitude = guideStarMagnitude;
wfs.camera.readOutNoise = 0.2;
wfs.camera.photonNoise = true;
wfs.framePixelThreshold = 0;

    ngs=ngs.*tel*wfs;

%% noise convariance matrix
fprintf('\n\n Initialising the noisy WFS....\n')
nMeasurements = 2000; %% nMeasurements With noise
slopes = zeros(wfs.nSlope,nMeasurements);

for kMeas=1:nMeasurements
    +wfs;
    slopes(:,nMeasurements) = wfs.slopes;
end

Cn = slopes*slopes'/nMeasurements;

% break %%% To check the calibration slopes / covariance measurements

clear ngs wfs 
fprintf('... DONE! \n\n ')


%% Constructing matrices for the Rmv reconstruction
eval_command_tmp = 'Cn';
for ii = 2:nGs
    eval_command_tmp = strcat(eval_command_tmp, ', Cn');
end
eval_command = strcat('CnAst  = blkdiag(', eval_command_tmp, ');');
eval(eval_command);
% % % CnAst = blkdiag( Cn, Cn, Cn ); %%% block-diagonal Covariance matrix, the size of N Ngs

eval_command_tmp = 'Dzcc';
for ii = 2:nGs
    eval_command_tmp = strcat(eval_command_tmp, ', Dzcc');
end
eval_command = strcat('DAst  = blkdiag(', eval_command_tmp, ');');
eval(eval_command);
% % % % DAst  = blkdiag( Dzcc, Dzcc, Dzcc);


%% Multiple WFSes
fprintf('\n\n Initialising the WFSes....\n')
for iGs = 1:nGs
    wfs(iGs) = shackHartmann(nLenslet,nPx,minLightRatio);
    ast(iGs) = ast(iGs).*tel*wfs(iGs);
    wfs(iGs).INIT;
    wfs(iGs).validLenslet = utilities.piston(nLenslet);
    +wfs(iGs);
end
wfsCal = wfs(1); %%% Pick up a wavefront sensor for Calibation of DM
fprintf('... DONE! \n\n ')



%% Deformable mirror
nActuator = nLenslet + 1;
bif = influenceFunction('monotonic',40/100);

dm = deformableMirror(nActuator,...
    'modes',bif,...
    'resolution',nPx,...
    'validActuator',logical(utilities.piston(nActuator)) );

nDm = length(dm);


dmObjectVec = cell(1,number_of_methods2copmare);

    for kPd = 1:nSciobj
        bif = influenceFunction('monotonic',40/100);
        
        dmStatic(kPd) = deformableMirror(nActuator,...
            'modes',bif,...
            'resolution',nPx,...
            'validActuator',logical(utilities.piston(nActuator)));
    end
    dmObjectVec{1} = dmStatic;


%% Command Matrix
ngs = source;
ngs.magnitude = 0;
ngs.magnitude = guideStarMagnitude;
ngs = ngs.*tel;

    dmWfsCalib = calibration(dm,wfsCal,ngs,ngs.wavelength/2); %% why half amplitude?

dmWfsCalib.threshold = 5e5; %% for svd values


for iPd = 1:nSciobj
    commandMatrix(:,:,iPd) = dmWfsCalib.M;
end


ngs.magnitude = guideStarMagnitude;
dm.coefs = 0;

for iGs = 1:nGs
    wfs(iGs).camera.clockRate = frameTime/samplingTime;
    wfs(iGs).camera.photonNoise = true;
    wfs(iGs).camera.readOutNoise = 0.2;
    wfs(iGs).framePixelThreshold = 0.2;
end

%%%% "Modes to Dm commands" calibration
% % zernModeMax = zernike.nModeFromRadialOrder(maxRadialDegree);
% % zernDm = zernike(2:zernModeMax,'resolution',nPx,'pupil',tel.pupil,'D',tel.D);
%
% %%% Using modal interaction with model of DM
% 
% % src = source;
% % src.wavelength = guideStarWavelength;
% % zernDm.c = eye(zernDm.nMode);
% % src = src.*zernDm;
% % dmCal = dmCal\src;
% % src = src.*tel*dmCal;
% % Z2U = dmCal.coefs;
% % dmCal.coefs = 0;


%% Computing the Dz matrix that 
fprintf('\n\n Computing the Z2U matrix ....\n')

    Z2U = dmWfsCalib.M*Dzcc;% Using zernike derivatives and slopes to DM command matrix

fprintf('... DONE! \n\n ')



%% Imagers:
% imgrObjectVec = cell(1,number_of_methods2copmare);
imgrObjectVec = cell(1,2);

    imgrStatic = imager(tel);
    imgrStatic.exposureTime = exposureTime;
    imgrStatic.eeWidth = ircsSlitWidth; % arcsec %% - what is it???
    imgrStatic.imgLens.nyquistSampling = 1;
    % Get a reference frame
    sciStatic = sciStatic.*tel*dmStatic*imgrStatic;
    imgrStatic.flush
    imgrStatic.referenceFrame = imgrStatic.frame;
    imgrObjectVec{1} = imgrStatic;

% 
% 
%     imgrNoCorrection = imager(tel);
%     imgrNoCorrection.exposureTime = exposureTime;
%     imgrNoCorrection.eeWidth = ircsSlitWidth; % arcsec %% - what is it???
%     imgrNoCorrection.imgLens.nyquistSampling = 1;
%     % Get a reference frame
%     sciNoCorrection = sciStatic.*tel*dmStatic*imgrNoCorrection;
%     imgrNoCorrection.flush
%     imgrNoCorrection.referenceFrame = imgrNoCorrection.frame;
%     imgrObjectVec{2} = imgrNoCorrection;
% 



%%
%% BEGIN::::: COMPUTING MATRICES for TOMOGRAPHY ::::::::::::::::
%%

%% Create extended Zernike object for Tomographic projections
zernModeMaxProj = zernike.nModeFromRadialOrder(maxRadialDegreeProj);
zernProj = zernike(2:zernModeMaxProj,'resolution',nPx,'pupil',tel.pupil,'D',tel.D);




% break 

projFileName = strcat('PrecomputeProjections_MaxRadDeg',num2str(maxRadialDegree), '_MaxRadDeg',num2str(maxRadialDegreeProj), '.mat');

if exist(projFileName,'file')

    fprintf('\n\n Loading precomputed Projection_beta and Projection_alpha ....\n')
        load(projFileName)
    fprintf('... DONE! \n\n ')

else
% % %     [Pbeta, Pbetacell] = AnalyticalSmallFootprintExpansion(zernProj,tel,sciCal,atm);
% % %     [Palpha, Palphacell] = AnalyticalSmallFootprintExpansion(zernProj,tel,gs,atm);
% % %     save(projFileName,'Pbeta','Palpha')

    %% Zernike projection for atmosphere layers
    fprintf('\n\n Computing Projection_alpha....\n')

    altitudes = [atm.layer.altitude];
    Projection_alpha_Cell = tool_analytical_small_footprint_expansion_with_piston(zernProj, altitudes, ast, tel, nLayers, nGs);
    Projection_alpha = cell2mat(Projection_alpha_Cell);
    % % %     [Palpha2, Palphacell2] = AnalyticalSmallFootprintExpansion(zernProj,tel,ast,atm); %tested, coincide with tool_analytical_small_
    fprintf('... DONE! \n\n ')


    %% Modal Projection of the science object 
    fprintf('\n\n Computing Projection_beta....\n')

    altitudes = [atm.layer.altitude];
    Projection_beta_Cell = tool_analytical_small_footprint_expansion_with_piston(zernProj, altitudes, sciCal, tel, nLayers, nSciobj);
    Projection_beta = cell2mat(Projection_beta_Cell);
    % % %     [Pbeta2, Pbetacell2] = AnalyticalSmallFootprintExpansion(zernProj,tel,sciCal,atm); %tested, coincide with tool_analytical_small_
    fprintf('... DONE! \n\n ')

    save(projFileName,'Projection_alpha','Projection_beta')
end



%% Computing the cross-correlation matrix of Zernike coefficients from the modes of the atmosphere object for EACH ATM LAYER
Sigma_phi_Cell = cell(nLayers,nLayers);  %%% composite cell array for Atm Covariance

zernWfsi = zernike(2:zernModeMaxProj,'resolution',nPx,'pupil',tel.pupil);

for kLayer = 1:nLayers
    for kLayer1 = 1:nLayers
        
        if kLayer == kLayer1

            zernWfsi.D = tel.diameterAt( atm.layer(kLayer).altitude );

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






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ---- PREDICTION ALGORITHMS ---
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Reconstruction methods
switch reconstruction_method

         case 'Static MMSE'
%% Minimum Variance reconstructor (static MMSE):
    Rmv = (Sigma_betaalpha*DAst')/(DAst*Sigma_alpha*DAst' + CnAst);
%%%%% FIXME - potential numerical problem here, consider using Rmv2 = (Sigma_betaalpha* DAst')*pinv(DAst*Sigma_alpha*DAst' + CnAst);

         case 'AR1'
    params.flag_correlation_method = 'corr_via_spatio_angular';
%     params.flag_correlation_method = 'corr_via_PSD';

    params.order  = 1;  %% order of the AR model to fit into Zernike
    params.Npts = 20; %% Npts is Number of points to consider in the fitting
    params.fN = 1001; % number of points in the frequency domain - a smaller number of points leads to wrong results in the decorrelation functions cc, Jul2013
    params.freq_part = 0.45;  %% default is 200, or 0.2 of maximum freq.
    params.flag_messages = 1; %% show or hide test messages from the function

    [ARmodel] = lam_utilities.tool_compute_AR_model_via_covariance_bfgs_fitting(atm,tel,zernProj,params);

    AR1_model = diag(ARmodel.A(:));
    Sigma_betaalpha_AR1 = Projection_beta(idp,idy)*AR1_model*Sigma_phi*Projection_alpha(idx,idy)';
    RmvAR1 = (Sigma_betaalpha_AR1*DAst')/(DAst*Sigma_alpha*DAst' + CnAst);

        case 'AmmseLag'
    params.order  = 1;  %% order of the AR model to fit into Zernike
    params.Npts = 2; %% Npts is Number of points to consider in the fitting
    [AmmseLag] = lam_utilities.tool_compute_AR_model_via_spatio_angular_temporal_corr(atm,tel,zernProj,params);

    eval_command_tmp = 'AmmseLag{1,1}';
    for ii = 2:atm.nLayer
        eval_command_tmp = strcat(eval_command_tmp, ', AmmseLag{1,', num2str(ii),'}' );
    end
    eval_command = strcat('AmmseMatrix  = blkdiag(', eval_command_tmp, ');');
    eval(eval_command);


    Sigma_betaalphaAmmse = Projection_beta(idp,idy)*AmmseMatrix*Sigma_phi*Projection_alpha(idx,idy)';
    RmvAmmse = (Sigma_betaalphaAmmse*DAst')/(DAst*Sigma_alpha*DAst' + CnAst);

end




% break %%%%%%%%%%%%%%%%%%%%%




%% Begin Open Loop
tel = tel+atm;
coefsVec = zeros(dm.nValidActuator,number_of_methods2copmare,nIteration);
nSlope = wfs(1).nSlope;

for kIteration = 1:nIteration
    fprintf('Iteration %d out of %d \n',kIteration,nIteration)
    
    for iFrame = 1:frameTime/samplingTime

        +tel;  % Always advance atm by 1ms

        for iGs = 1:nGs;
            ast(iGs) = ast(iGs).*tel*wfs(iGs);
        end

            sciObjectVec{1} = sciObjectVec{1}.*tel*dmObjectVec{1}; %%% ???
            sciObjectVec{1} = sciObjectVec{1}*imgrObjectVec{1};


        if iFrame == fixedLagTimeInMs
            if kIteration == 1
                    dmObjectVec{1}.coefs = -coefsVec(:,1,kIteration);
            else
                    dmObjectVec{1}.coefs = -coefsVec(:,1,kIteration-1);
            end
        end

    end %%% END generating the frames


    for iGs = 1:nGs;
        slopesStack(:,iGs,kIteration) = wfs(iGs).slopes;
    end

         NGS(1:nSlope,kIteration) = slopesStack(:,1,kIteration);



    %% Reconstruction methods
    switch reconstruction_method

         case 'Static MMSE'
%%%% Static MMSE reconstruction

%% This is for 3 NGSes
        NGS(1:nSlope,kIteration) = slopesStack(:,1,kIteration);
        NGS(nSlope+1:nSlope*2,kIteration) = slopesStack(:,3,kIteration);
        NGS(2*nSlope+1:nSlope*3,kIteration) = slopesStack(:,2,kIteration);
             

        coefsVec(:,1,kIteration) = Z2U*Rmv*NGS(:,kIteration);     %%% The problem appears to be here: the slopes are a lot larger than the coefficientsVectors

         case 'AR1'
%%%% AR1 prediction
            coefsVec(:,1,kIteration) = Z2U*RmvAR1*NGS(:,kIteration);


         case 'AmmseLag'
            coefsVec(:,1,kIteration) = Z2U*RmvAmmse*NGS(:,kIteration);

    end

end %%% for kIteration

%% Why there is no instantaneous Strehl ratio available?
imgrObjectVec{1}.ee     %% encircled energy
imgrObjectVec{1}.strehl %% strehl ration
% ./OOMAO-Raven/OOMAOlibUpdated/imager.m
