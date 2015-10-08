%**************************************************************************
% Modal Reconstructor Comparison
%
% 1 = static MMSE
% 2 = AR1 prediction
% 3 = AR2 prediction
% 4 = Ammse prediction
% 5 = LQG with AR2 prediction
% 6 = LQG with Ammse prediction
%
%**************************************************************************
% rmpath(genpath('/Users/ccorreia/Documents/LAM/Simulations/_lib_matlab/OOMAO-Feb2015/'))
% rmpath(genpath('/Users/ccorreia/Documents/LAM/Simulations/OOMAO-v15/'))
% addpath(genpath('/Users/ccorreia/Documents/LAM/Simulations/OOMAO-Raven/'))

addpath('./OOMAO-Raven/OOMAOlibUpdated/')  % Only use OOMAOlibUpdated!
addpath('./OOMAO-Raven/')  % % % % % addpath('./OOMAO-Raven/ravenOOMAOlib/') DONT USE THIS!



%% Reconstructors to Compute
RecVec = [1];
% RecVec = [1 2 3 4 5 6];

%% Parameter definition
samplingTime    = 0.001;% seconds
frameTime       = 0.01;%[0.01,0.015,0.020,0.025]; seconds
% fixedLagTimeInMs = 3;
fixedLagTimeInMs = 0;

nIteration              = 100;
exposureTime            = nIteration*frameTime;
ircsSlitWidth           = 0.14;
randn('state', 25);     % sets the global random state
plotPhase = false;       % draw 2-D phase during OL execution

% System Atmosphere
refWavelength   = 500e-9;
wavelength      = photometry.R;
altitude = [0,5.5,11.5]*1e3;
fractionalR0 = [0.596,0.223816,0.180184];
windSpeed = [2.84*2,6,17];
windDirection = [pi/2, pi, pi];

L0 = 40;
r0 = 0.19;
windSpeedNeg = -windSpeed;

atm = atmosphere(refWavelength,r0,L0,...
    'altitude',altitude,...
    'fractionnalR0',fractionalR0,...
    'windSpeed',windSpeed,...
    'windDirection',windDirection);

atmNeg = atmosphere(refWavelength,r0,L0,...
    'altitude',altitude,...
    'fractionnalR0',fractionalR0,...
    'windSpeed',windSpeedNeg,...
    'windDirection',windDirection);

atm.wavelength = wavelength;
nLayer = atm.nLayer;

% Telescope
nPx = 120;
tel = telescope(8,...
    'fieldOfViewInArcMin',3.5,...
    'resolution',nPx,...
    'samplingTime',samplingTime);
telAR = telescope(8,'fieldOfViewInArcMin',3.5,...
    'resolution',nPx,...
    'samplingTime',frameTime);

%% Guide Stars

dAsterism               = 1.0; % 0.5, 0.75, 1.0, 1.5, 2 diameter in arcmin
guideStarWavelength     = photometry.R;
guideStarMagnitude      = 14;%[14,14.5,15,15.5,16,16.5,17]; % lower later

% Regular asterism
% gs = source('asterism',{[3,arcmin(dAsterism/2),0]},'wavelength',guideStarWavelength,...
%    'magnitude',guideStarMagnitude);

% Asterism based on Raven CU pinholes
gs = source('zenith',[arcsec(45.3),arcsec(35.5),arcsec(45.3)],...
    'azimuth',[3.0309,0.7854,-1.6815],...
    'wavelength',guideStarWavelength,'magnitude',guideStarMagnitude);
nGs = length(gs);

%% Science Targets

scienceObjectMag = 10;
scienceObjectWavelength = photometry.H;
sciZenithVector = arcsec(7);
sciAzimuthVector = -0.7854;
nSci = length(sciZenithVector);
sciObjectVec = cell(1,6);
sciCal = source('zenith',[arcsec(7)],'azimuth',[-0.7854],...
    'wavelength',scienceObjectWavelength,'magnitude',scienceObjectMag);

if any(RecVec == 1)
    sciStatic = source('zenith',[arcsec(7)],'azimuth',[-0.7854],...
        'wavelength',scienceObjectWavelength,'magnitude',scienceObjectMag);
    sciObjectVec{1} = sciStatic;
end

if any(RecVec == 2)
    sciPredAR1 = source('zenith',[arcsec(7)],'azimuth',[-0.7854],...
        'wavelength',scienceObjectWavelength,'magnitude',scienceObjectMag);
    sciObjectVec{2} = sciPredAR1;
end

if any(RecVec == 3)
    sciPredAR2 = source('zenith',[arcsec(7)],'azimuth',[-0.7854],...
        'wavelength',scienceObjectWavelength,'magnitude',scienceObjectMag);
    sciObjectVec{3} = sciPredAR2;
end

if any(RecVec == 4)
    sciPredAmmse = source('zenith',[arcsec(7)],'azimuth',[-0.7854],...
        'wavelength',scienceObjectWavelength,'magnitude',scienceObjectMag);
    sciObjectVec{4} = sciPredAmmse;
end

if any(RecVec == 5)
    sciLQGAR2 = source('zenith',[arcsec(7)],'azimuth',[-0.7854],...
        'wavelength',scienceObjectWavelength,'magnitude',scienceObjectMag);
    sciObjectVec{5} = sciLQGAR2;
end

if any(RecVec == 6)
    sciLQGAmmse = source('zenith',[arcsec(7)],'azimuth',[-0.7854],...
        'wavelength',scienceObjectWavelength,'magnitude',scienceObjectMag);
    sciObjectVec{6} = sciLQGAmmse;
end

% Calibration sources
ngs = source;
ngsCal = source('wavelength',guideStarWavelength,'magnitude',0);

%% Wavefront sensor
nLenslet = 10;
wfs = shackHartmann(nLenslet,nPx,0.75);
ngs = ngs.*tel*wfs;
wfs.INIT;
wfs.camera.photonNoise = false;
wfs.validLenslet = utilities.piston(10);
+wfs;

%% Zernike coefs to WFS slopes calibration
maxRadialDegree = 9;

% From modal interaction with model of WFS
zern = zernike(2:zernike.nModeFromRadialOrder(maxRadialDegree),'resolution',nPx,'pupil',tel.pupil);
zern.lex = false;
zern.c = eye(zern.nMode);
ngs=ngs.*zern*wfs;
Dzcc = wfs.slopes;
zern.lex = true;

% From analytical derivatives (not scaled properly)

% lensletMask = logical(wfsCal.validLenslet);
% nLenslet = size(lensletMask,1);
% modes       = 1:(maxRadialDegree+1)*(maxRadialDegree+2)/2;
% zern        = zernike(modes,'resolution',nLenslet);
% dZern       = [ zern.xDerivative(lensletMask,:) ; zern.yDerivative(lensletMask,:)  ];
% Dzcc = dZern(:,2:end); %to get rid of piston

%% With noise
ngs.wavelength = photometry.R;
ngs.magnitude = guideStarMagnitude;
wfs.camera.readOutNoise = 0.2;
wfs.camera.photonNoise = true;
wfs.framePixelThreshold = 0;
ngs=ngs.*tel*wfs;

%% noise convariance matrix

nMeas = 1000;
slopes = zeros(wfs.nSlope,nMeas);
parfor kMeas=1:nMeas
    +wfs
    slopes(:,kMeas) = wfs.slopes;
end
Cn = slopes*slopes'/nMeas;
clear ngs wfs

CnAst = blkdiag( Cn , Cn , Cn );
DAst = blkdiag( Dzcc , Dzcc , Dzcc  );
%% WFSs

nLenslet = 10;

for iGs = 1:nGs
    wfs(iGs) = shackHartmann(nLenslet,nPx,0.75);
    gs(iGs) = gs(iGs).*tel*wfs(iGs);
    wfs(iGs).INIT;
    wfs(iGs).validLenslet = utilities.piston(nLenslet);
    +wfs(iGs);
end
wfsCal = wfs(1);

%% Deformable Mirrors

nActuator = nLenslet+1;

bif = influenceFunction('monotonic',40/100);

dmCal = deformableMirror(nActuator,...
    'modes',bif,...
    'resolution',nPx,...
    'validActuator',logical(utilities.piston(nActuator)));
dmObjectVec = cell(1,6);
if any(RecVec ==1)
    for kPd = 1:nSci
        bif = influenceFunction('monotonic',40/100);
        
        dmStatic(kPd) = deformableMirror(nActuator,...
            'modes',bif,...
            'resolution',nPx,...
            'validActuator',logical(utilities.piston(nActuator)));
    end
    dmObjectVec{1} = dmStatic;
end

if any(RecVec ==2)
    for kPd = 1:nSci
        bif = influenceFunction('monotonic',40/100);
        
        dmPredAR1(kPd) = deformableMirror(nActuator,...
            'modes',bif,...
            'resolution',nPx,...
            'validActuator',logical(utilities.piston(nActuator)));
    end
    dmObjectVec{2} = dmPredAR1;
end

if any(RecVec ==3)
    for kPd = 1:nSci
        bif = influenceFunction('monotonic',40/100);
        
        dmPredAR2(kPd) = deformableMirror(nActuator,...
            'modes',bif,...
            'resolution',nPx,...
            'validActuator',logical(utilities.piston(nActuator)));
    end
    dmObjectVec{3} = dmPredAR2;
end

if any(RecVec ==4)
    for kPd = 1:nSci
        bif = influenceFunction('monotonic',40/100);
        
        dmPredAmmse(kPd) = deformableMirror(nActuator,...
            'modes',bif,...
            'resolution',nPx,...
            'validActuator',logical(utilities.piston(nActuator)));
    end
    dmObjectVec{4} = dmPredAmmse;
end

if any(RecVec ==5)
    for kPd = 1:nSci
        bif = influenceFunction('monotonic',40/100);
        
        dmLQGAR2(kPd) = deformableMirror(nActuator,...
            'modes',bif,...
            'resolution',nPx,...
            'validActuator',logical(utilities.piston(nActuator)));
    end
    dmObjectVec{5} = dmLQGAR2;
end

if any(RecVec ==6)
    for kPd = 1:nSci
        bif = influenceFunction('monotonic',40/100);
        
        dmLQGAmmse(kPd) = deformableMirror(nActuator,...
            'modes',bif,...
            'resolution',nPx,...
            'validActuator',logical(utilities.piston(nActuator)));
    end
    dmObjectVec{6} = dmLQGAmmse;
end

%% Command Matrix
ngs = source;
ngs.magnitude = 0;
ngs = ngs.*tel;
dmWfsCalib = calibration(dmCal,wfsCal,ngs,ngs.wavelength/2);
dmWfsCalib.threshold = 5e5;
for iPd = 1:nSci
    commandMatrix(:,:,iPd) = dmWfsCalib.M;
end

ngs.magnitude = guideStarMagnitude;

dmCal.coefs = 0;

for iGs = 1:nGs
    wfs(iGs).camera.clockRate = frameTime/samplingTime;
    wfs(iGs).camera.photonNoise = true;
    wfs(iGs).camera.readOutNoise = 0.2;
    wfs(iGs).framePixelThreshold = 0.2;
end

%% "Modes to Dm commands" calibration
zernModeMax = zernike.nModeFromRadialOrder(maxRadialDegree);
zernDm = zernike(2:zernModeMax,'resolution',nPx,'pupil',tel.pupil,'D',tel.D);

% Using modal interaction with model of DM

% src = source;
% src.wavelength = guideStarWavelength;
% zernDm.c = eye(zernDm.nMode);
% src = src.*zernDm;
% dmCal = dmCal\src;
% src = src.*tel*dmCal;
% Z2U = dmCal.coefs;
% dmCal.coefs = 0;

% Using zernike derivatives and slopes to DM command matrix

Z2U = dmWfsCalib.M*Dzcc;
%% Imagers

imgrObjectVec = cell(1,6);

if any(RecVec == 1)
    imgrStatic = imager(tel);
    imgrStatic.exposureTime = exposureTime;
    imgrStatic.eeWidth = ircsSlitWidth; % arcsec
    imgrStatic.imgLens.nyquistSampling = 1;
    % Get a reference frame
    sciStatic = sciStatic.*tel*dmStatic*imgrStatic;
    imgrStatic.flush
    imgrStatic.referenceFrame = imgrStatic.frame;
    imgrObjectVec{1} = imgrStatic;
    
end

if any(RecVec == 2)
    imgrPredAR1 = imager(tel);
    imgrPredAR1.exposureTime = exposureTime;
    imgrPredAR1.eeWidth = ircsSlitWidth; % arcsec
    imgrPredAR1.imgLens.nyquistSampling = 1;
    % Get a reference frame
    sciPredAR1 = sciPredAR1.*tel*...
        dmPredAR1*imgrPredAR1;
    imgrPredAR1.flush
    imgrPredAR1.referenceFrame = imgrPredAR1.frame;
    imgrObjectVec{2} = imgrPredAR1;
end

if any(RecVec == 3)
    imgrPredAR2 = imager(tel);
    imgrPredAR2.exposureTime = exposureTime;
    imgrPredAR2.eeWidth = ircsSlitWidth; % arcsec
    imgrPredAR2.imgLens.nyquistSampling = 1;
    % Get a reference frame
    sciPredAR2 = sciPredAR2.*tel*...
        dmPredAR2*imgrPredAR2;
    imgrPredAR2.flush
    imgrPredAR2.referenceFrame = imgrPredAR2.frame;
    imgrObjectVec{3} = imgrPredAR2;
end

if any(RecVec == 4)
    imgrPredAmmse = imager(tel);
    imgrPredAmmse.exposureTime = exposureTime;
    imgrPredAmmse.eeWidth = ircsSlitWidth; % arcsec
    imgrPredAmmse.imgLens.nyquistSampling = 1;
    % Get a reference frame
    sciPredAmmse = sciPredAmmse.*tel*...
        dmPredAmmse*imgrPredAmmse;
    imgrPredAmmse.flush
    imgrPredAmmse.referenceFrame = imgrPredAmmse.frame;
    imgrObjectVec{4} = imgrPredAmmse;
end

if any(RecVec == 5)
    imgrLQGAR2 = imager(tel);
    imgrLQGAR2.exposureTime = exposureTime;
    imgrLQGAR2.eeWidth = ircsSlitWidth; % arcsec
    imgrLQGAR2.imgLens.nyquistSampling = 1;
    % Get a reference frame
    sciLQGAR2 = sciLQGAR2.*tel*...
        dmLQGAR2*imgrLQGAR2;
    imgrLQGAR2.flush
    imgrLQGAR2.referenceFrame = imgrLQGAR2.frame;
    imgrObjectVec{5} = imgrLQGAR2;
end

if any(RecVec == 6)
    imgrLQGAmmse = imager(tel);
    imgrLQGAmmse.exposureTime = exposureTime;
    imgrLQGAmmse.eeWidth = ircsSlitWidth; % arcsec
    imgrLQGAmmse.imgLens.nyquistSampling = 1;
    % Get a reference frame
    sciLQGAmmse = sciLQGAmmse.*tel*...
        dmLQGAmmse*imgrLQGAmmse;
    imgrLQGAmmse.flush
    imgrLQGAmmse.referenceFrame = imgrLQGAmmse.frame;
    imgrObjectVec{6} = imgrLQGAmmse;
end

%% Create extended zernike object

maxRadialDegreeProj = 23;
zernModeMaxProj = zernike.nModeFromRadialOrder(maxRadialDegreeProj);
zernProj = zernike(2:zernModeMaxProj,'resolution',nPx,'pupil',tel.pupil,'D',tel.D);
projFileName = 'ProjectionsRO23.mat';

if exist(projFileName,'file')
    load(projFileName)
else
    [Pbeta, Pbetacell] = AnalyticalSmallFootprintExpansion(zernProj,tel,sciCal,atm);
    [Palpha, Palphacell] = AnalyticalSmallFootprintExpansion(zernProj,tel,gs,atm);
    save(projFileName,'Pbeta','Palpha')
end

nLayer = atm.nLayer;
Slayered = cell(nLayer,nLayer);
fr0 = [atm.layer.fractionnalR0];
altitudes = [atm.layer.altitude];
nModeAll = zernProj.nMode;
nMode = zernike.nModeFromRadialOrder(maxRadialDegree)-1;

for kLayer = 1:nLayer
    for kLayer1 = 1:nLayer
        if kLayer == kLayer1
            
            zernWfsi = zernike(2:zernModeMaxProj,'resolution',nPx,'pupil',tel.pupil);
            zernWfsi.D = tel.diameterAt(altitudes(kLayer));%/tel.D*2;
            Slayered{kLayer,kLayer1} = phaseStats.zernikeCovariance(zernWfsi,atm)*fr0(kLayer);
        else
            Slayered{kLayer,kLayer1} = zeros(zernWfsi.nMode);
        end
    end
end
Sphi = cell2mat(Slayered);

idy = [];
for kGs=1:nLayer
    ii = 2:zernModeMaxProj;
    idy = [idy ii+(zernModeMaxProj)*(kGs-1)];
end
idx = [];
for kGs=1:nGs
    ii = 1:zernModeMax-1;
    idx = [idx ii+(zernModeMaxProj-1)*(kGs-1)+kGs];
end

ids = [];
for kGs=1:nLayer
    ii = (1:zernModeMaxProj-1);
    ids = [ids ii+(zernWfsi.nMode)*(kGs-1)];
end

% --- matrix Palpha*SigmaPhi*Palpha' is already truncated
tmp = Palpha(idx,idy)*Sphi(ids,ids)*Palpha(idx,idy)';
Sigma_alpha = Palpha(idx,idy)*Sphi(ids,ids)*Palpha(idx,idy)';

idp = [];
for kGs=1:nSci
    ii = 1:zernModeMax-1;
    idp = [idp ii+(zernModeMax-1)*(kGs-1)+kGs];
end
Sigma_betaalpha = Pbeta(idp,idy)*Sphi*Palpha(idx,idy)';


Palpha_PR = Palpha(idx,idy);

Pbeta_PR = Pbeta(idp,idy);


PalphaSmall = Palpha_PR;
PbetaSmall  = Pbeta_PR;

fprintf('PROJECTIONS DONE \n')

%% static MMSE

Rmv = (Sigma_betaalpha*DAst')/(DAst*Sigma_alpha*DAst' + CnAst);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ---- PREDICTION ALGORITHMS ---
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if any(RecVec == 2)
    nu = logspace(-5,3,1000/2);
    ord = 1;
    [Aar1] = computeARmodel(nu,atm,telAR,zernProj,[1e-5 200],ord,10);
    Aar1 = Aar1(:);
    Aar1 = diag(Aar1);
    Sigma_betaalphaAR1 = Pbeta(idp,idy)*Aar1*Sphi*Palpha(idx,idy)';
    RmvPredAR1 = (Sigma_betaalphaAR1*DAst')/(DAst*Sigma_alpha*DAst' + CnAst);
end

if any(RecVec == 3)||any(RecVec == 5)
    nu = logspace(-5,3,1000/2);
    ord = 2;
    [Aar2 Bar2] = computeARmodel(nu,atm,telAR,zernProj,[1e-5 200],ord,10);
    
    Aar2 = Aar2(:);
    Bar2 = Bar2(:);
    S1 = diag(Aar2./(1-Bar2))*Sphi;
    S1 = (S1 + S1')/2;
    Aar2 = diag(Aar2);
    Bar2 = diag(Bar2);
    Sigma_betaalphaAar2 = Pbeta(idp,idy)*Aar2*Sphi*Palpha(idx,idy)';
    Sigma_betaalphaBar2 = Pbeta(idp,idy)*Bar2*Sphi*Palpha(idx,idy)';
    RmvApredAR2 = (Sigma_betaalphaAar2*DAst')/(DAst*Sigma_alpha*DAst' + CnAst);
    RmvBpredAR2 = (Sigma_betaalphaBar2*DAst')/(DAst*Sigma_alpha*DAst' + CnAst);
    Par2 = [PbetaSmall zeros(nMode*nSci,nModeAll*nLayer)];
end

if any(RecVec == 4)||any(RecVec == 6)
    tel = tel+atm;
    Npts = 2;
    S = cell(atm.nLayer,Npts);
    for kLayer = 1:atm.nLayer
        myAltitude = [0 1]*1e3;
        myFractionalR0 = [0 fractionalR0(kLayer)];
        myWindSpeed = [15 windSpeed(kLayer)];
        myWindDirection = [0 windDirection(kLayer)];
        myatm = atmosphere(refWavelength,r0,L0,...
            'altitude',myAltitude,...
            'fractionnalR0',myFractionalR0,...
            'windSpeed',myWindSpeed,...
            'windDirection',myWindDirection);
        phase = atm.layer(kLayer).phase;
        pupil = logical(utilities.piston(size(phase,1)));
        zern = zernike(2:zernModeMaxProj,'resolution',size(phase,1),'pupil',pupil,'D',tel.diameterAt(atm.layer(kLayer).altitude));
        
        
        for kr = 1:Npts
            alpha = (kr-1)*(frameTime)*myWindSpeed(2)/myAltitude(2);
            alpha2 = (kr-1)*(fixedLagTimeInMs*1e-3)*myWindSpeed(2)/myAltitude(2);
            alpha = 1/2*alpha*180/pi*3600;
            alpha2 = 1/2*alpha2*180/pi*3600;
            myAst = source('asterism',{[2,alpha*cougarConstants.arcsec2radian,0]});
            myAst2 = source('asterism',{[2,alpha2*cougarConstants.arcsec2radian,0]});
            S{kLayer,kr} = phaseStats.zernikeAngularCovariance(zern,myatm,myAst);
            Sfixed{kLayer,kr} = phaseStats.zernikeAngularCovariance(zern,myatm,myAst2);
        end
        A(:,:,kLayer) = S{kLayer,2}'/S{kLayer,1};
        Afixed(:,:,kLayer) = Sfixed{kLayer,2}'/Sfixed{kLayer,1};
    end
    Ammse = blkdiag(A(:,:,1), A(:,:,2), A(:,:,3));
    AmmseFixed = blkdiag(Afixed(:,:,1), Afixed(:,:,2), Afixed(:,:,3));
    
    Sigma_betaalphaAmmse = Pbeta(idp,idy)*Ammse*Sphi*Palpha(idx,idy)';
    RmvAmmse = (Sigma_betaalphaAmmse*DAst')/(DAst*Sigma_alpha*DAst' + CnAst);
    Pammse = PbetaSmall;
end

if any(RecVec == 5)
    
    SnAR2 = Sphi - Aar2*Sphi*Aar2' - Bar2*Sphi*Bar2' - Bar2*S1*Aar2' - Aar2*S1*Bar2';
    while min(eig(SnAR2))<1e-5
        SnAR2 = SnAR2 + eye(nModeAll*nLayer)*abs(min(eig(SnAR2)));
    end
    
    R = CnAst;
    Q = blkdiag(SnAR2,zeros(size(SnAR2)));
    Q = (Q+Q')./2;
    N = 0;
    
    AdAR2 = [Aar2 Bar2;eye(nModeAll*nLayer) zeros(nModeAll*nLayer)];
    
    GammaAR2 = blkdiag(eye(nLayer*nModeAll),zeros(nLayer*nModeAll));
    BdAR2 = [zeros(nModeAll*nLayer*2) GammaAR2];
    GPalphaAR2 = DAst*PalphaSmall;
    CdAR2 = [GPalphaAR2 zeros(nGs*2*wfs(1).nValidLenslet,nLayer*nModeAll)];
    
    SigmaInfAR2 = solveRiccati_2Xspeed(AdAR2,CdAR2,Q,R,10);
    MinfAR2 = (SigmaInfAR2*CdAR2')/(CdAR2*SigmaInfAR2*CdAR2'+R);
end

if any(RecVec == 6)
    SnAmmse = Sphi - Ammse*Sphi*Ammse';
    % this is to check that Sn is PSD - fix
    %(it gets stuck in a long loop when min(eig(Sn)) very small)
    
    while min(eig(SnAmmse))<1e-5
        SnAmmse = SnAmmse + eye(nModeAll*nLayer)*abs(min(eig(SnAmmse)));
    end
    
    R = CnAst;
    Q = SnAmmse;
    Q = (Q+Q')./2;
    N = 0;
    Bmmse = zeros(size(Ammse));
    AdAmmse = Ammse;
    %Ad3ms = [AmmseFixed Bmmse;eye(nModeAll*nLayer) zeros(nModeAll*nLayer)];
    
%     GammaAmmse = blkdiag(eye(nLayer*nModeAll),zeros(nLayer*nModeAll));
%     BdAmmse = [zeros(nModeAll*nLayer*2) GammaAmmse];
    GPalphaAmmse = DAst*PalphaSmall;
    CdAmmse = GPalphaAmmse;
    
    SigmaInfAmmse = solveRiccati_2Xspeed(AdAmmse,CdAmmse,Q,R,10);
    MinfAmmse = (SigmaInfAmmse*CdAmmse')/(CdAmmse*SigmaInfAmmse*CdAmmse'+R);
end

%% How To Display turbulence and residual phase
if plotPhase == true
    nRec = length(RecVec);
    if nRec == 6
    plotID1 = RecVec(1);
    plotID2 = RecVec(2);
    plotID3 = RecVec(3);
    plotID4 = RecVec(4);
    plotID5 = RecVec(5);
    plotID6 = RecVec(6);
    figure(101)
    rad2mic = 1e6/sciObjectVec{plotID1}.waveNumber;
    h = imagesc([sciObjectVec{plotID1}.meanRmPhase,sciObjectVec{plotID2}.meanRmPhase;...
        sciObjectVec{plotID3}.meanRmPhase,sciObjectVec{plotID4}.meanRmPhase;...
        sciObjectVec{plotID5}.meanRmPhase,sciObjectVec{plotID6}.meanRmPhase]*rad2mic);
    
    ax = gca;
    axis equal tight
    ylabel(colorbar,'WFE [\mum]')
    else
        disp('***Manually adjust plotIDs to display the residual phase maps****')
    end
%     figure(102)
%     plotID5 = RecVec(3);
%     plotID6 = RecVec(5);
%     
%     h2 = imagesc([sciObjectVec{plotID3}.meanRmPhase,sciObjectVec{plotID4}.meanRmPhase]*rad2mic);
%     ax = gca;
%     axis equal tight
%     ylabel(colorbar,'WFE [\mum]')
    
end

%% Open Loop

tel = tel+atm;
coefsVec = zeros(dmCal.nValidActuator,6,nIteration);
nSlope = wfs(1).nSlope;

fprintf('%5d \n',1)
for kIteration = 1:nIteration
    fprintf('\b\b\b\b\b%5d',kIteration)
    
    for iFrame = 1:frameTime/samplingTime
        +tel;  %Always advance atm by 1ms
        for iGs = 1:nGs;
            gs(iGs) = gs(iGs).*tel*wfs(iGs);
        end
        for iRec = 1:length(RecVec)
            RecID = RecVec(iRec);
            sciObjectVec{RecID} = sciObjectVec{RecID}.*tel*dmObjectVec{RecID};
        end
        
        if plotPhase == true
            set(h,'Cdata',[sciObjectVec{plotID1}.meanRmPhase,sciObjectVec{plotID2}.meanRmPhase;...
                sciObjectVec{plotID3}.meanRmPhase,sciObjectVec{plotID4}.meanRmPhase;...
                sciObjectVec{plotID5}.meanRmPhase,sciObjectVec{plotID6}.meanRmPhase]*rad2mic)
            title(ax,sprintf('#%4d/%4d',kIteration,nIteration))
            drawnow
            
%             set(h2,'Cdata',[sciObjectVec{plotID5}.meanRmPhase,sciObjectVec{plotID6}.meanRmPhase]*rad2mic)
%             title(ax,sprintf('#%4d/%4d',kIteration,nIteration))
%             drawnow
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
        
    end
    
    for iGs = 1:nGs;
        slopesStack(:,iGs,kIteration) = wfs(iGs).slopes;
    end
    
        NGS(1:nSlope,kIteration) = slopesStack(:,1,kIteration);
        NGS(nSlope+1:nSlope*2,kIteration) = slopesStack(:,3,kIteration);
        NGS(2*nSlope+1:nSlope*3,kIteration) = slopesStack(:,2,kIteration);
    
    if any(RecVec == 1)
        k1 = NGS(:,kIteration);
        coefsVec(:,1,kIteration) = Z2U*Rmv*k1(:);
    end
    
    if any(RecVec == 2)
        k1 = NGS(:,kIteration);
        coefsVec(:,2,kIteration) = Z2U*RmvPredAR1*k1(:);
    end
    
    if any(RecVec == 3)
        if kIteration > 1
            k1 = NGS(:,kIteration);
            k2 = NGS(:,kIteration-1);
            coefsVec(:,3,kIteration) = Z2U*(RmvApredAR2*k1(:) + RmvBpredAR2*k2(:));
        end
    end
    
    if any(RecVec == 4)
        k1 = NGS(:,kIteration);
        coefsVec(:,4,kIteration) = Z2U*RmvAmmse*k1(:);
    end
    
    if any(RecVec == 5)
        S_NGS = NGS(:,kIteration);
        if kIteration == 1;
            kStateAR2 = MinfAR2*S_NGS(:);
        else
            sHatAR2 = CdAR2*kStateAR2;
            kStateAR2 = kStateAR2+MinfAR2*(S_NGS(:)-sHatAR2);
        end
        
        kStateAR2 = AdAR2*kStateAR2;
        coefsVec(:,5,kIteration) = Z2U*Par2*kStateAR2;
    end
    
    if any(RecVec == 6)
        S_NGS = NGS(:,kIteration);
        if kIteration == 1;
            kStateAmmse = MinfAmmse*S_NGS(:);
        else
            sHatAmmse = CdAmmse*kStateAmmse;
            kStateAmmse = kStateAmmse+MinfAmmse*(S_NGS(:)-sHatAmmse);
        end
        
        kStateAmmse = AdAmmse*kStateAmmse;
        coefsVec(:,6,kIteration) = Z2U*Pammse*kStateAmmse;
    end
end
    for iRec = 1:length(RecVec)
        RecID = RecVec(iRec);
        strehl(RecID) = imgrObjectVec{RecID}.strehl;
        EE(RecID) = imgrObjectVec{RecID}.ee;
    end