clc
close all
clear all
% warning off
%%% using the modalTomographyHowto.m example for setting up 3NGSs
%%% Based on demo_modal_tomography in OOMAO/mikhail.

addpath('./OOMAO-Raven/OOMAOlibUpdated/')  % Only use OOMAOlibUpdated!
addpath('./OOMAO-Raven/')


params = struct; %%% creating the Structure with additional parameters for the computeARmodel
params.show_figures = 1; %% display the diagnostic plots during the simulation.

randn('state', 25);     % sets the global random state

nIteration       = 40; %% Number of frames to simulate and take

samplingTime     = 0.001;% seconds
frameTime        = 0.002;%[0.01,0.015,0.020,0.025]; seconds
exposureTime     = nIteration*frameTime;

ircsSlitWidth    = 0.14;
fixedLagTimeInMs = 1;

RecVec = [1]; %%% will be more methods
number_of_methods2copmare = length(RecVec); %% add different methods like MMSE and Ammse

%% \subsubsection{The atmosphere}
refWavelength   = 500e-9;
wavelength      = photometry.R;
altitudes       = [0,    5,   10]*1e3;
fractionalR0    = [0.5, 0.3, 0.2];
windSpeed       = [10,  5,   5];
windDirection   = [pi/2, pi, pi];

L0 = 40; %% outer scale, [m]
r0 = 0.19; %% Fried parameter, [m]


% if (exist('atm') == 0)
atm = atmosphere(refWavelength, r0, L0,...
    'fractionnalR0',fractionalR0, ...
    'altitude',     altitudes,...
    'windSpeed',    windSpeed,...
    'windDirection',windDirection);

nLayers = atm.nLayer;
% end %% if (exist('atm') 


%% Define the Asterism (star constellation) for the tomography
%%% Don't forget that the FOV of the telescope (meta-pupil) must cover all the guide stars
dAsterism               = 1.0; % 0.5, 0.75, 1.0, 1.5, 2 diameter in arcmin
guideStarWavelength     = photometry.R;
guideStarMagnitude      = 14;%[14,14.5,15,15.5,16,16.5,17]; % lower later

%%% Guide Stars with Regular asterism
ast = source('asterism',{[3,arcmin(dAsterism/2),0]},...
            'wavelength',guideStarWavelength,... 
            'magnitude',guideStarMagnitude);
nGs = length(ast);


%% Science Targets
scienceObjectMag = 10;
scienceObjectWavelength = photometry.H;
sciZenithVector = arcsec(7);
sciAzimuthVector = 0;

sci_Calibaration = source('zenith',0*cougarConstants.arcsec2radian,'wavelength',photometry.R); %%% this is probably Projection Disk (like the telescope pupil)
nSciobj = length(sci_Calibaration);

%% Adding science objects for different reconstruction methods
sciObjectVec = cell(1,1); %% if you want to compare two or more methods - increase the number of cells, like (1,N)

%%%%% MMSE reconstructor
    sciStatic = source('zenith',[sciZenithVector(1)],'azimuth',[sciZenithVector(1)],...
        'wavelength',scienceObjectWavelength,'magnitude',scienceObjectMag);
    sciObjectVec{1} = sciStatic;


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
wfs.validLenslet = utilities.piston(10);
+wfs;
break
% % % Saving the turbulence aberrated phase
% turbPhase = ngs.meanRmPhase; %% this is a turbulent phase



%% Zernike coefs to WFS slopes calibration
maxRadialDegree = 4;  %% - setting the Max Radial Degree = how many modes we will use

% From modal interaction with model of WFS
zern = zernike(2:zernike.nModeFromRadialOrder(maxRadialDegree),'resolution',nPx,'pupil',tel.pupil);
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
nMeasurements = 20; %% nMeasurements With noise
slopes = zeros(wfs.nSlope,nMeasurements);

for kMeas=1:nMeasurements
    +wfs;
    slopes(:,nMeasurements) = wfs.slopes;
end

Cn = slopes*slopes'/nMeasurements;
clear ngs wfs 


%% Constructing matrices for the Rmv reconstruction
CnAst = blkdiag( Cn, Cn, Cn ); %%% block-diag matrix, the size of N Ngs
DAst  = blkdiag( Dzcc, Dzcc, Dzcc);


%% Multiple WFSes
for iGs = 1:nGs
    wfs(iGs) = shackHartmann(nLenslet,nPx,minLightRatio);
    ast(iGs) = ast(iGs).*tel*wfs(iGs);
    wfs(iGs).INIT;
    wfs(iGs).validLenslet = utilities.piston(nLenslet);
    +wfs(iGs);
end

wfsCal = wfs(1); %%% Pick up a wavefront sensor for Calibation of DM

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
ngs = ngs.*tel;

    dmWfsCalib = calibration(dm,wfsCal,ngs,ngs.wavelength/2);

dmWfsCalib.threshold = 5e5; 

for iPd = 1:nSciobj
    commandMatrix(:,:,iPd) = dmWfsCalib.M; %%% dmWfsCalib is the command matrix
end

ngs.magnitude = guideStarMagnitude;
dm.coefs = 0;

for iGs = 1:nGs
    wfs(iGs).camera.clockRate = frameTime/samplingTime;
    wfs(iGs).camera.photonNoise = true;
    wfs(iGs).camera.readOutNoise = 0.2;
    wfs(iGs).framePixelThreshold = 0.2;
end

%% "Modes to Dm commands" calibration
zernModeMax = zernike.nModeFromRadialOrder(maxRadialDegree);
zernDm = zernike(2:zernModeMax,'resolution',nPx,'pupil',tel.pupil,'D',tel.D);

% Using zernike derivatives and slopes to DM command matrix
Z2U = dmWfsCalib.M*Dzcc;



%% Imagers:
imgrObjectVec = cell(1,number_of_methods2copmare);

    imgrStatic = imager(tel);
    imgrStatic.exposureTime = exposureTime;
    imgrStatic.eeWidth = ircsSlitWidth; % arcsec
    imgrStatic.imgLens.nyquistSampling = 1;

    % Get a reference frame
    sciStatic = sciStatic.*tel*dmStatic*imgrStatic;
    imgrStatic.flush
    imgrStatic.referenceFrame = imgrStatic.frame;
    imgrObjectVec{1} = imgrStatic;


%%
%% BEGIN::::: COMPUTING MATRICES for TOMOGRAPHY ::::::::::::::::
%%

%% Create extended Zernike object for Tomographic projections
maxRadialDegreeProj = 10;
zernModeMaxProj = zernike.nModeFromRadialOrder(maxRadialDegreeProj);
zernProj = zernike(2:zernModeMaxProj,'resolution',nPx,'pupil',tel.pupil,'D',tel.D);


%% Zernike projection for atmosphere layers
altitudes = [atm.layer.altitude];
Projection_alpha_Cell = tool_analytical_small_footprint_expansion_with_piston(zernProj, altitudes, ast, tel, nLayers, nGs);
Projection_alpha = cell2mat(Projection_alpha_Cell);


%% Modal Projection of the science object 
altitudes = [atm.layer.altitude];
Projection_beta_Cell = tool_analytical_small_footprint_expansion_with_piston(zernProj, altitudes, sci_Calibaration, tel, nLayers, nSciobj);
Projection_beta = cell2mat(Projection_beta_Cell);


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



%% This vector of indices is 1x(zernModeMax-1), sans piston
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


%% Computing the projection of the covariance in the direction of Guidestars
Sigma_alpha = Projection_alpha(idx,idy)*Sigma_phi(ids,ids)*Projection_alpha(idx,idy)'; %% 


%% This vector of indices is 1x(zernModeMax-1), sans piston
idp = [];
for kSci=1:nSciobj
    ii = 1:zernModeMax-1;
    idp = [idp ii+(zernModeMax-1)*(kSci-1)+kSci];
end


%% Computing the Projection for the sci object:
Sigma_betaalpha = Projection_beta(idp,idy)*Sigma_phi*Projection_alpha(idx,idy)';

fprintf('\n PROJECTIONS DONE\n \n')
%% End calculating the projections


%% Minimum Variance reconstructor (static MMSE):
Rmv = (Sigma_betaalpha*DAst')/(DAst*Sigma_alpha*DAst' + CnAst);
%%%%% FIXME - potential numerical problem here, consider using Rmv2 = (Sigma_betaalpha* DAst')*pinv(DAst*Sigma_alpha*DAst' + CnAst);



%% Begin Open Loop
tel = tel+atm;
coefsVec = zeros(dm.nValidActuator,number_of_methods2copmare,nIteration);
nSlope = wfs(1).nSlope;

for kIteration = 1:nIteration
    fprintf('Iteration %d out of %d \n',kIteration,nIteration)
    
    for iFrame = 1:frameTime/samplingTime

        +tel;  %Always advance atm by 1ms

        for iGs = 1:nGs;
            ast(iGs) = ast(iGs).*tel*wfs(iGs);
        end

        for iRec = 1:length(RecVec)
            RecID = RecVec(iRec);
            sciObjectVec{RecID} = sciObjectVec{RecID}.*tel*dmObjectVec{RecID};
        end


        for iRec = 1:length(RecVec)
            RecID = RecVec(iRec);
            sciObjectVec{RecID} = sciObjectVec{RecID}*imgrObjectVec{RecID};
        end

        if iFrame == fixedLagTimeInMs
            if kIteration == 1
                for iRec = 1:length(RecVec)
                    RecID = RecVec(iRec);
                    dmObjectVec{RecID}.coefs = -coefsVec(:,RecID,kIteration);
                end
            else
                for iRec = 1:length(RecVec)
                    RecID = RecVec(iRec);
                    dmObjectVec{RecID}.coefs = -coefsVec(:,RecID,kIteration-1);
                end
            end
        end

    end %%% END generating the frames



    for iGs = 1:nGs;
        slopesStack(:,iGs,kIteration) = wfs(iGs).slopes;
    end

        NGS(1:nSlope,kIteration) = slopesStack(:,1,kIteration);
        NGS(nSlope+1:nSlope*2,kIteration) = slopesStack(:,3,kIteration);
        NGS(2*nSlope+1:nSlope*3,kIteration) = slopesStack(:,2,kIteration);

    %%% The problem appears to be here: the slopes are a lot larger than the coefficientsVectors
%     if any(RecVec == 1)
        k1 = NGS(:,kIteration);
        coefsVec(:,1,kIteration) = Z2U*Rmv*k1(:);
%     end

end %%% for kIteration