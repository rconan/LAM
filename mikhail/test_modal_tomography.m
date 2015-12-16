clc
close all
clear all
% warning off
%%% using the modalTomographyHowto.m example for setting up 3NGSs
%%% Based on demo_modal_tomography in OOMAO/mikhail.

tool_inside_path_adder('matlab_code/mikhail', 'matlab_code');

format short

params = struct; %%% creating the Structure with additional parameters for the computeARmodel
params.show_figures = 1; %% display the diagnostic plots during the simulation.

% randn('state', 25);     % sets the global random state
randn('state', sum(100*clock)) %Initialize randn to a different state each time

nIteration       = 100; %% Number of frames to simulate and take.

samplingTime     = 0.001;% seconds
frameTime        = 0.001;%[0.01,0.015,0.020,0.025]; seconds

% % exposureTime     = nIteration*frameTime; %???? causes problems with OOMAO imager.

fixedLagTimeInMs = 0;


reconstruction_method = 'Static MMSE';
% reconstruction_method = 'AR1';
% reconstruction_method = 'AmmseLag';

% reconstruction_method = 'LQGviaAmmse';


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

% refWavelength   = 500e-9;
% wavelength      = photometry.R;
% altitudes = [0,5.5,11.5]*1e3;
% fractionalR0 = [0.596,0.223816,0.180184];
% windSpeed = [2.84*2,6,17];
% windDirection = [pi/2, pi, pi];


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


%% Zernike modes:
maxRadialDegree     = 3;  %% - setting the Max Radial Degree = how many modes we will use
maxRadialDegreeProj = 10; %% on how many modes we will roject


%% Define the Asterism (star constellation) for the tomography
%%% Don't forget that the FOV of the telescope (meta-pupil) must cover all the guide stars
dAsterism               = 1.0; % 0.5, 0.75, 1.0, 1.5, 2 diameter in arcmin
guideStarWavelength     = photometry.R;
guideStarMagnitude      = 14;%[14,14.5,15,15.5,16,16.5,17]; % lower later

%% Guide Stars with Regular asterism
ast = source('asterism',{[3,arcmin(dAsterism/2),0]},...
            'wavelength',guideStarWavelength,... 
            'magnitude',guideStarMagnitude);

%%%% Only one guidestar:
% ast = source('asterism',{[1,arcmin(0),0]},...
%             'wavelength',guideStarWavelength,... 
%             'magnitude',guideStarMagnitude);

% % Asterism based on Raven CU pinholes
% ast = source('zenith',[arcsec(45.3),arcsec(35.5),arcsec(45.3)],...
%     'azimuth',[3.0309,0.7854,-1.6815],...
%     'wavelength',guideStarWavelength,'magnitude',guideStarMagnitude);

nGs = length(ast);


%% Science Targets
scienceObjectMag = 10;  %% Apparent magnitude
scienceObjectWavelength = photometry.H; %% Infrared
sciZenithVector = arcsec(0);
sciAzimuthVector = 0;


% scienceObjectMag = 10;
% scienceObjectWavelength = photometry.H;
% sciZenithVector = arcsec(7);
% sciAzimuthVector = -0.7854;

sciCal = source('zenith',[sciZenithVector],'azimuth',[sciAzimuthVector],...
    'wavelength',scienceObjectWavelength,'magnitude',scienceObjectMag);


nSciobj = length(sciCal);


%% Adding science objects for different reconstruction methods
% sciObjectVec = cell(1,1); %% if you want to compare two or more methods - increase the number of cells, like (1,N)

%%%%% MMSE reconstructor
    sciStatic = source('zenith',[sciZenithVector],'azimuth',[sciAzimuthVector],...
        'wavelength',scienceObjectWavelength,'magnitude',scienceObjectMag);
%     sciObjectVec{1} = sciStatic;


% Calibration sources
ngs = source; %%% source must be initialised before the WFS
ngsCal = source('wavelength',guideStarWavelength,'magnitude',0);


%% Wavefront sensor
nLenslet = 10;
minLightRatio = 0.75; %% large values (more than 0.8  can lead to poor reconstruction - you are truncating the phase!).
wfs = shackHartmann(nLenslet,nPx,minLightRatio);  %%% causes problems in MATLAB 2015a

%%% Propagation of the calibration source to the WFS through the telescope
ngs = ngs.*tel*wfs;
wfs.INIT;
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

clear ngs wfs 

if (params.show_figures == 1)
    figure, imagesc(Cn), title('Noise covariance matrix');
    figure, plot(diag(Cn)), title('Noise covariance - diagonal elements only');
end

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


% dmObjectVec = cell(1,number_of_methods2copmare);

    for kPd = 1:nSciobj
        bif = influenceFunction('monotonic',40/100);
        
        dmStatic(kPd) = deformableMirror(nActuator,...
            'modes',bif,...
            'resolution',nPx,...
            'validActuator',logical(utilities.piston(nActuator)));
    end
%     dmObjectVec{1} = dmStatic;


%% Command Matrix
ngs = source;
ngs.magnitude = 0;
ngs = ngs.*tel;
dmWfsCalib = calibration(dm,wfsCal,ngs,ngs.wavelength/2); %% why half amplitude?
dmWfsCalib.threshold = 5e5; %% for svd values

for iPd = 1:nSciobj
    commandMatrix(:,:,iPd) = dmWfsCalib.M;
end


ngs.magnitude = guideStarMagnitude;
dm.coefs = 0;

for iGs = 1:nGs
%     wfs(iGs).camera.clockRate = frameTime/samplingTime;
    wfs(iGs).camera.clockRate = 1;
%     wfs(iGs).camera.exposureTime = frameTime;
    wfs(iGs).camera.photonNoise = false;
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
% imgrObjectVec = cell(1,1);
% 
%     imgrStatic = imager(tel);
% %     imgrStatic = imagerRaven(tel);
%     imgrStatic.exposureTime = exposureTime;
%     imgrStatic.eeWidth = ircsSlitWidth; % arcsec %% - what is it???
%     imgrStatic.imgLens.nyquistSampling = 1;
%     
%     % Get a reference frame
%     sciStatic = sciStatic.*tel*dmStatic*imgrStatic;
%     imgrStatic.flush
%     imgrStatic.referenceFrame = imgrStatic.frame;
%     imgrObjectVec{1} = imgrStatic;

cam = imager(tel);


% Get a reference frame
tel = tel - atm;  %% detach the atmosphere from the telescope
sciStatic = sciStatic.*tel*cam;
cam.referenceFrame = cam.frame;
+sciStatic;


flush(cam)
cam.clockRate    = 1;
% cam.exposureTime = exposureTime;
cam.startDelay   = 0;


%%
%% BEGIN::::: COMPUTING MATRICES for TOMOGRAPHY ::::::::::::::::
%%

%% Create extended Zernike object for Tomographic projections
zernModeMaxProj = zernike.nModeFromRadialOrder(maxRadialDegreeProj);
zernProj = zernike(2:zernModeMaxProj,'resolution',nPx,'pupil',tel.pupil,'D',tel.D);

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



% break


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

    [ARmodel] = tool_compute_AR_model_via_covariance_bfgs_fitting(atm,tel,zernProj,params);

    AR1_model = diag(ARmodel.A(:));
    Sigma_betaalpha_AR1 = Projection_beta(idp,idy)*AR1_model*Sigma_phi*Projection_alpha(idx,idy)';
    RmvAR1 = (Sigma_betaalpha_AR1*DAst')/(DAst*Sigma_alpha*DAst' + CnAst);

        case 'AmmseLag'
    params.order  = 1;  %% order of the AR model to fit into Zernike
    params.Npts = 2; %% Npts is Number of points to consider in the fitting
    [Ammse] = tool_compute_AR_model_via_spatio_angular_temporal_corr(atm,tel,zernProj,params);

    eval_command_tmp = 'AmmseLag{1,1}';
    for ii = 2:atm.nLayer
        eval_command_tmp = strcat(eval_command_tmp, ', AmmseLag{1,', num2str(ii),'}' );
    end
    eval_command = strcat('AmmseMatrix  = blkdiag(', eval_command_tmp, ');');
    eval(eval_command);


    Sigma_betaalphaAmmse = Projection_beta(idp,idy)*AmmseMatrix*Sigma_phi*Projection_alpha(idx,idy)';
    RmvAmmse = (Sigma_betaalphaAmmse*DAst')/(DAst*Sigma_alpha*DAst' + CnAst);
    
    
    case 'LQGviaAmmse'

    params.order  = 1;  %% order of the AR model to fit into Zernike
    params.Npts = 2; %% Npts is Number of points to consider in the fitting
    [AmmseLag] = tool_compute_AR_model_via_spatio_angular_temporal_corr(atm,tel,zernProj,params);

    eval_command_tmp = 'AmmseLag{1,1}';
    for ii = 2:atm.nLayer
        eval_command_tmp = strcat(eval_command_tmp, ', AmmseLag{1,', num2str(ii),'}' );
    end
    eval_command = strcat('AmmseMatrix  = blkdiag(', eval_command_tmp, ');');
    eval(eval_command);
   
    Sigma_epsilon = Sigma_phi - AmmseMatrix*Sigma_phi*AmmseMatrix';  %% excitation noise covariance
    Sigma_epsilon = Sigma_epsilon + 1e-4*eye(size(Sigma_epsilon));  %% prevent numerical errors - use crude diagonal loading

    R = CnAst;
    Q = (Sigma_epsilon+Sigma_epsilon')./2;
    GPalphaAmmse = DAst*Projection_alpha(idx,idy);

    %%% Solve DARE:
    Sigma_infty_Ammse = dare (AmmseMatrix, GPalphaAmmse', Q, R);
    M_infty_Ammse = (Sigma_infty_Ammse*GPalphaAmmse')/(GPalphaAmmse*Sigma_infty_Ammse*GPalphaAmmse' + R);

end




% break %%%%%%%%%%%%%%%%%%%%%




%% Begin Open Loop
tel = tel+atm;
coefsVec = zeros(dm.nValidActuator,number_of_methods2copmare,nIteration);
nSlope = wfs(1).nSlope;

if (params.show_figures == 1)
    figure, imagesc(tel,[ast,sciCal]), title('Guidestar Asterism visualisation');
end

StackSlopeLarge = zeros(nSlope*nGs,1);
slopesStack = zeros(nSlope,nGs);

for kIteration = 1:nIteration
    fprintf('\n \n   Iteration %d out of %d \n\n ',kIteration,nIteration)
    
    for iFrame = 1:frameTime/samplingTime

        +tel;  % Always advance atm by 1ms

        for iGs = 1:nGs;
            ast(iGs) = ast(iGs).*tel*wfs(iGs);
        end

%             sciObjectVec{1} = sciObjectVec{1}.*tel*dmObjectVec{1};
%             sciObjectVec{1} = sciObjectVec{1}*imgrObjectVec{1};
        sciCal = sciCal.*tel*dm*cam;

% ngsCombo = source('zenith',zeros(1,2),'azimuth',zeros(1,2),'magnitude',8);
% ngsCombo = ngsCombo.*tel*dm*wfs;
% 
% scienceCombo = source('zenith',zeros(1,2),'azimuth',zeros(1,2),'wavelength',photometry.J);
% scienceCombo = scienceCombo.*tel*dm*cam;

        if iFrame == fixedLagTimeInMs
            if kIteration == 1
%                     dmObjectVec{1}.coefs = -coefsVec(:,1,kIteration);
dm.coefs =  -coefsVec(:,1,kIteration);
            else
dm.coefs =  -coefsVec(:,1,kIteration-1);
%                     dmObjectVec{1}.coefs = -coefsVec(:,1,kIteration-1);
            end
        end %% for iFrame

    end %%% END generating the frames


    for iGs = 1:nGs;
        slopesStack(:,iGs) = wfs(iGs).slopes;
    end


% This is for 3 NGSes
        StackSlopeLarge(1:nSlope,1) = slopesStack(:,1);
        StackSlopeLarge(nSlope+1:nSlope*2,1) = slopesStack(:,3);
        StackSlopeLarge(2*nSlope+1:nSlope*3,1) = slopesStack(:,2);

% % %% This is for 3 NGSes
% %         StackSlopeLarge(1:nSlope) = slopesStack(:,1,kIteration);
% %         StackSlopeLarge(nSlope+1:nSlope*2) = slopesStack(:,2,kIteration);
% %         StackSlopeLarge(2*nSlope+1:nSlope*3) = slopesStack(:,3,kIteration);

    %% Reconstruction methods
    switch reconstruction_method

         case 'Static MMSE'
%%%% Static MMSE reconstruction             
        coefsVec(:,1,kIteration) = Z2U*Rmv*StackSlopeLarge(:);     %%% The problem appears to be here: the slopes are a lot larger than the coefficientsVectors

         case 'AR1'
%%%% AR1 prediction
            coefsVec(:,1,kIteration) = Z2U*RmvAR1*StackSlopeLarge(:);


         case 'AmmseLag'
%%%% Ammse SA prediction
            coefsVec(:,1,kIteration) = Z2U*RmvAmmse*StackSlopeLarge(:);

         case 'LQGviaAmmse'
%%%% LQG control via Ammse SA

            if kIteration == 1;
                phase_state = M_infty_Ammse*StackSlopeLarge(:);
            else
                phase_state = phase_state + M_infty_Ammse*(StackSlopeLarge(:) - GPalphaAmmse*phase_state);
            end

            phase_state_next = AmmseMatrix*phase_state;
            coefsVec(:,1,kIteration) = Z2U*Projection_beta(idp,idy)*phase_state_next;

                
    end

end %%% for kIteration

%% Why there is no instantaneous Strehl ratio available?
% imgrObjectVec{1}.strehl %% strehl ration
