%**************************************************************************
% Zonal Reconstructor Comparison
%
% 1 = Static SA
% 2 = Static Explicit
% 3 = Predictive SA
% 4 = Predictive Explicit
% 5 = LQG SA
% 6 = LQG Explicit
% 7 = Yoshito method
%**************************************************************************


strehl =

    0.2559    0.2506    0.2646    0.2587    0.2882    0.2990    0.2539

>> EE

EE =

    0.4852    0.4833    0.4938    0.4917    0.5143    0.5541    0.4850

strehl

%% Ensure the Raven-tailored OOMAO lib is used - **WITH MATLAB 2013 VERSION -- NOT 2015**
%rmpath(genpath('/Users/ccorreia/Documents/LAM/Simulations/_lib_matlab/OOMAO-Feb2015/'))
%rmpath(genpath('/Users/ccorreia/Documents/LAM/Simulations/OOMAO-Raven/'))

addpath(genpath('~/HarmoniSimuls/oomao_raven/'))
%addpath(genpath('/Users/ccorreia/Documents/LAM/Simulations/oomao_lam/'))


%addpath(genpath('/Users/ccorreia/Documents/LAM/Simulations/OOMAO-Raven-OLD/'))
%% Reconstructors to Compute

RecVec = [1 2 3 4 5 6 7];
%RecVec = [3 5 6 7];
%% Parameter definition

samplingTime     = 0.001;% seconds
frameTime        = 0.01;%[0.01,0.015,0.020,0.025]; seconds
fixedLagTimeInMs = 3;

nIteration              = 1000;
exposureTime            = nIteration*frameTime;
ircsSlitWidth           = 0.14;
randn('state', 25);     % sets the global random state
plotPhase = true;       % draw 2-D phase during OL execution

% System Atmosphere
refWavelength   = 500e-9;
wavelength      = photometry.R;
altitude        = [0,5.5,11.5]*1e3;
fractionalR0    = [0.596,0.223816,0.180184];
windSpeed       = [2.84*2,6,17];
windDirection   = [pi/2, pi, pi];

% altitude        = 0;
% fractionalR0    = 1;
% windSpeed       = 40;
% windDirection   = 0;
% 
% 
% altitude        = [0 5];
% fractionalR0    = [0.5 0.5];
% windSpeed       = [40 20];
% windDirection   = [0 0];


L0 = 40;
r0 = 0.19;
windSpeedNeg = -windSpeed;

atm = atmosphere(refWavelength,r0,L0,...
    'altitude',altitude,...
    'fractionnalR0',fractionalR0,...
    'windSpeed',windSpeed,...
    'windDirection',windDirection);

% atmNeg = atmosphere(refWavelength,r0,L0,...
%     'altitude',altitude,...
%     'fractionnalR0',fractionalR0,...
%     'windSpeed',windSpeedNeg,...
%     'windDirection',windDirection);

%cc added instead of the previous case
atmNeg = atmosphere(refWavelength,r0,L0,...
    'altitude',altitude,...
    'fractionnalR0',fractionalR0,...
    'windSpeed',windSpeed,...
    'windDirection',pi+windDirection);


atm.wavelength = wavelength;
nLayer = atm.nLayer;

% Telescope
nPx = 120;
tel = telescope(8,...
    'fieldOfViewInArcMin',3.5,...
    'resolution',nPx,...
    'samplingTime',samplingTime);

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
    sciStatSA = source('zenith',[arcsec(7)],'azimuth',[-0.7854],...
        'wavelength',scienceObjectWavelength,'magnitude',scienceObjectMag);
    sciObjectVec{1} = sciStatSA;
end

if any(RecVec == 2)
    sciStatExp = source('zenith',[arcsec(7)],'azimuth',[-0.7854],...
        'wavelength',scienceObjectWavelength,'magnitude',scienceObjectMag);
    sciObjectVec{2} = sciStatExp;
end

if any(RecVec == 3)
    sciPredSA = source('zenith',[arcsec(7)],'azimuth',[-0.7854],...
        'wavelength',scienceObjectWavelength,'magnitude',scienceObjectMag);
    sciObjectVec{3} = sciPredSA;
end

if any(RecVec == 4)
    sciPredExp = source('zenith',[arcsec(7)],'azimuth',[-0.7854],...
        'wavelength',scienceObjectWavelength,'magnitude',scienceObjectMag);
    sciObjectVec{4} = sciPredExp;
end

if any(RecVec == 5)
    sciLQGsa = source('zenith',[arcsec(7)],'azimuth',[-0.7854],...
        'wavelength',scienceObjectWavelength,'magnitude',scienceObjectMag);
    sciObjectVec{5} = sciLQGsa;
end

if any(RecVec == 6)
    sciLQGexp = source('zenith',[arcsec(7)],'azimuth',[-0.7854],...
        'wavelength',scienceObjectWavelength,'magnitude',scienceObjectMag);
    sciObjectVec{6} = sciLQGexp;
end

if any(RecVec == 7)
    sciYoshito = source('zenith',[arcsec(7)],'azimuth',[-0.7854],...
        'wavelength',scienceObjectWavelength,'magnitude',scienceObjectMag);
    sciObjectVec{7} = sciYoshito;
end

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
        
        dmStatSA(kPd) = deformableMirror(nActuator,...
            'modes',bif,...
            'resolution',nPx,...
            'validActuator',logical(utilities.piston(nActuator)));
    end
    dmObjectVec{1} = dmStatSA;
end

if any(RecVec ==2)
    for kPd = 1:nSci
        bif = influenceFunction('monotonic',40/100);
        
        dmStatExp(kPd) = deformableMirror(nActuator,...
            'modes',bif,...
            'resolution',nPx,...
            'validActuator',logical(utilities.piston(nActuator)));
    end
    dmObjectVec{2} = dmStatExp;
end

if any(RecVec ==3)
    for kPd = 1:nSci
        bif = influenceFunction('monotonic',40/100);
        
        dmPredSA(kPd) = deformableMirror(nActuator,...
            'modes',bif,...
            'resolution',nPx,...
            'validActuator',logical(utilities.piston(nActuator)));
    end
    dmObjectVec{3} = dmPredSA;
end

if any(RecVec ==4)
    for kPd = 1:nSci
        bif = influenceFunction('monotonic',40/100);
        
        dmPredExp(kPd) = deformableMirror(nActuator,...
            'modes',bif,...
            'resolution',nPx,...
            'validActuator',logical(utilities.piston(nActuator)));
    end
    dmObjectVec{4} = dmPredExp;
end

if any(RecVec ==5)
    for kPd = 1:nSci
        bif = influenceFunction('monotonic',40/100);
        
        dmLQGsa(kPd) = deformableMirror(nActuator,...
            'modes',bif,...
            'resolution',nPx,...
            'validActuator',logical(utilities.piston(nActuator)));
    end
    dmObjectVec{5} = dmLQGsa;
end

if any(RecVec ==6)
    for kPd = 1:nSci
        bif = influenceFunction('monotonic',40/100);
        
        dmLQGexp(kPd) = deformableMirror(nActuator,...
            'modes',bif,...
            'resolution',nPx,...
            'validActuator',logical(utilities.piston(nActuator)));
    end
    dmObjectVec{6} = dmLQGexp;
end

if any(RecVec ==7)
    for kPd = 1:nSci
        bif = influenceFunction('monotonic',40/100);
        
        dmYoshito(kPd) = deformableMirror(nActuator,...
            'modes',bif,...
            'resolution',nPx,...
            'validActuator',logical(utilities.piston(nActuator)));
    end
    dmObjectVec{7} = dmYoshito;
end


%% Command Matrix
ngs = source;
ngs.magnitude = 0;
ngs = ngs.*tel;
dmWfsCalib = calibration(dmCal,wfsCal,ngs,ngs.wavelength/2);
dmWfsCalib.threshold = 5e5;
for iPd = 1:nSci
    %   for iGs = 1:nGs
    commandMatrix(:,:,iPd) = dmWfsCalib.M;
    %  end
end

ngs.magnitude = guideStarMagnitude;

dmCal.coefs = 0;

%% Noise Covariance Zonal

nMeas = 250;
for iGs = 1:nGs
    
    wfs(iGs).camera.photonNoise = true;
    wfs(iGs).camera.readOutNoise = .2;
    wfs(iGs).framePixelThreshold = 0;
    gs(iGs) = gs(iGs).*tel*wfs(iGs);
    slopes = zeros(wfs(iGs).nSlope,nMeas);
    for kMeas=1:nMeas
        +wfs(iGs)
        slopes(:,kMeas) = wfs(iGs).slopes;
    end
    Cn(:,:,iGs) = slopes*slopes'/nMeas;
    wfs(iGs).camera.clockRate = frameTime/samplingTime;
    wfs(iGs).framePixelThreshold = 0.2;
end
CnZ = blkdiag(Cn(:,:,1), Cn(:,:,2), Cn(:,:,3));

%% Imagers

imgrObjectVec = cell(1,7);

if any(RecVec == 1)
    imgrStatSA = imager(tel);
    imgrStatSA.exposureTime = exposureTime;
    imgrStatSA.eeWidth = ircsSlitWidth; % arcsec
    imgrStatSA.imgLens.nyquistSampling = 1;
    % Get a reference frame
    sciStatSA = sciStatSA.*tel*dmStatSA*imgrStatSA;
    imgrStatSA.flush
    imgrStatSA.referenceFrame = imgrStatSA.frame;
    imgrObjectVec{1} = imgrStatSA;
    
end

if any(RecVec == 2)
    imgrStatExp = imager(tel);
    imgrStatExp.exposureTime = exposureTime;
    imgrStatExp.eeWidth = ircsSlitWidth; % arcsec
    imgrStatExp.imgLens.nyquistSampling = 1;
    % Get a reference frame
    sciStatExp = sciStatExp.*tel*...
        dmStatExp*imgrStatExp;
    imgrStatExp.flush
    imgrStatExp.referenceFrame = imgrStatExp.frame;
    imgrObjectVec{2} = imgrStatExp;
end

if any(RecVec == 3)
    imgrPredSA = imager(tel);
    imgrPredSA.exposureTime = exposureTime;
    imgrPredSA.eeWidth = ircsSlitWidth; % arcsec
    imgrPredSA.imgLens.nyquistSampling = 1;
    % Get a reference frame
    sciPredSA = sciPredSA.*tel*...
        dmPredSA*imgrPredSA;
    imgrPredSA.flush
    imgrPredSA.referenceFrame = imgrPredSA.frame;
    imgrObjectVec{3} = imgrPredSA;
end

if any(RecVec == 4)
    imgrPredExp = imager(tel);
    imgrPredExp.exposureTime = exposureTime;
    imgrPredExp.eeWidth = ircsSlitWidth; % arcsec
    imgrPredExp.imgLens.nyquistSampling = 1;
    % Get a reference frame
    sciPredExp = sciPredExp.*tel*...
        dmPredExp*imgrPredExp;
    imgrPredExp.flush
    imgrPredExp.referenceFrame = imgrPredExp.frame;
    imgrObjectVec{4} = imgrPredExp;
end

if any(RecVec == 5)
    imgrLQGsa = imager(tel);
    imgrLQGsa.exposureTime = exposureTime;
    imgrLQGsa.eeWidth = ircsSlitWidth; % arcsec
    imgrLQGsa.imgLens.nyquistSampling = 1;
    % Get a reference frame
    sciLQGsa = sciLQGsa.*tel*...
        dmLQGsa*imgrLQGsa;
    imgrLQGsa.flush
    imgrLQGsa.referenceFrame = imgrLQGsa.frame;
    imgrObjectVec{5} = imgrLQGsa;
end

if any(RecVec == 6)
    imgrLQGexp = imager(tel);
    imgrLQGexp.exposureTime = exposureTime;
    imgrLQGexp.eeWidth = ircsSlitWidth; % arcsec
    imgrLQGexp.imgLens.nyquistSampling = 1;
    % Get a reference frame
    sciLQGexp = sciLQGexp.*tel*...
        dmLQGexp*imgrLQGexp;
    imgrLQGexp.flush
    imgrLQGexp.referenceFrame = imgrLQGexp.frame;
    imgrObjectVec{6} = imgrLQGexp;
end


if any(RecVec == 7)
    imgrYoshito = imager(tel);
    imgrYoshito.exposureTime = exposureTime;
    imgrYoshito.eeWidth = ircsSlitWidth; % arcsec
    imgrYoshito.imgLens.nyquistSampling = 1;
    % Get a reference frame
    sciYoshito = sciYoshito.*tel*...
        dmYoshito*imgrYoshito;
    imgrYoshito.flush
    imgrYoshito.referenceFrame = imgrYoshito.frame;
    imgrObjectVec{7} = imgrYoshito;
end


%% Zonal Gradient Operator

%[Gamma,gridMask] = sparseGradientMatrix(wfs(1));
% CC change
[Gamma,gridMask] = sparseGradientMatrix(wfs(1),utilities.piston(nLenslet*2+1));

d = tel.D/nLenslet;
GammaBeta = Gamma/(2*pi);

Gamma = repmat({Gamma},nGs,1);
Gamma = blkdiag(Gamma{:});
Gamma = Gamma/(2*pi);


zonalMmseHR = linearMMSE(nLenslet*2+1,tel.D,atm,gs,sciCal,'pupil',gridMask,'unit',-9);
zonalMmseHRpred = linearMMSE(nLenslet*2+1,tel.D,atmNeg,gs,sciCal,'pupil',gridMask,'unit',-9,'lag',frameTime + fixedLagTimeInMs/1000);

RecStatSA = commandMatrix*GammaBeta*zonalMmseHR.Cox{1}*Gamma'/(Gamma*zonalMmseHR.Cxx*Gamma'+CnZ);
RecPredSA = commandMatrix*GammaBeta*zonalMmseHRpred.CoxLag{1}*Gamma'/(Gamma*zonalMmseHR.Cxx*Gamma'+CnZ);

%% Bilinear interpolation operator

if any(RecVec == 2)||any(RecVec == 4)||any(RecVec == 6)
    tel = tel+atm;
    tel = tel-atm;
    [Htau,maskTau] = bilinearSplineInterpMat([gs,sciCal],atm,tel,gridMask,1/100*5);
  
    Htau(nGs+1:end,:) = [];
    Htau = cell2mat(Htau);
    Htau(:,~cell2mat(maskTau)) = [];
     
    
    [H,mask] = bilinearSplineInterpMat([gs,sciCal],atm,tel,gridMask,0);
    mask = maskTau;
    Hss = H(nGs+1:end,:);
    Hss = cell2mat(Hss);
    Hss(:,~cell2mat(mask)) = [];
    H(nGs+1:end,:) = [];
    H = cell2mat(H);
    H(:,~cell2mat(mask)) = [];

end

%% Zonal layered phase covariance matrix

if any(RecVec == 2)||any(RecVec == 4)||any(RecVec == 6) any(RecVec == 7)
    
    fprintf(' Computing the layered, point-wise phase covariance matrix ...\n');
    
    CovPhi = cell(atm.nLayer);
    for kLayer=1:atm.nLayer
        for mLayer=1:atm.nLayer
            
            if kLayer == mLayer
                kLayer
                nPixLayer = sqrt(length(mask{kLayer}));
                
                [x,y] = meshgrid((0:nPixLayer-1)*d/2);
                rho = complex(x,y);
                FullCov = phaseStats.covarianceMatrix(rho,slab(atm, kLayer));
                CovPhi{kLayer, kLayer} = FullCov(mask{kLayer}, mask{kLayer});%*atm.layer(kLayer).fractionnalR0;
                
                
                ws = atm.layer(kLayer).windSpeed;
                theta = atm.layer(kLayer).windDirection;
                deltax = ws*cos(theta)*(frameTime+fixedLagTimeInMs/1000);
                deltay = ws*sin(theta)*(frameTime+fixedLagTimeInMs/1000);
                
                rho1 = complex(x+deltax,y+deltay);
                
                FullTempCov = phaseStats.covarianceMatrix(rho,rho1,slab(atm,kLayer));
                LayeredSpatioTemporalCov{kLayer, kLayer} = FullTempCov(mask{kLayer}, mask{kLayer});
                
                % --- compute the MMSE predicted phase one step ahead ---
                Az{kLayer, kLayer} = LayeredSpatioTemporalCov{kLayer,kLayer}*pinv(CovPhi{kLayer,kLayer});%FullAz(mask{kLayer}, mask{kLayer});
            end
        end
    end
    % --- populate remainder cells with zeros ---
    for kLayer=1:atm.nLayer
        for mLayer=1:atm.nLayer
            n = size(CovPhi{kLayer, kLayer},2);
            m = size(CovPhi{mLayer, mLayer},1);
            if kLayer ~= mLayer
                CovPhi{kLayer, mLayer} = zeros(n,m);
                Az{kLayer, mLayer} = zeros(n,m);
            end
        end
        
    end
    
    CPhi = cell2mat(CovPhi);
    AzPhi = cell2mat(Az);
    
    CxxExp = H*CPhi*H';
    CoxExp = Hss*CPhi*H';
    CoxExpPred = Hss*AzPhi*CPhi*H';
    
    E0 = (CoxExp*Gamma')/(Gamma*CxxExp*Gamma'+CnZ);
    E0pred = (CoxExpPred*Gamma')/(Gamma*CxxExp*Gamma'+CnZ);
    
    iF = commandMatrix(:,:,1)*GammaBeta;
    RecStatExp = iF*E0;
    RecPredExp = iF*E0pred;
end
%% CC tests to check ray-tracing and SA cov matrices are the same
%  nPts = nLenslet*2+1;
%  Cxx = phaseStats.spatioAngularCovarianceMatrix(nPts,tel.D,slab(atm,1), gs(1), 'mask', gridMask);
%  C_alphaLag = phaseStats.spatioAngularCovarianceMatrix(nLenslet*2+1,tel.D,slab(atmNeg,1), gs(1), gs(1), 'mask', gridMask,'lag', frameTime);
%  C_alphaLag = cell2mat(C_alphaLag);
% 
% AzSA = C_alphaLag/Cxx;
% imagesc(AzSA - Az{1})
% %  
%  Cxx = phaseStats.spatioAngularCovarianceMatrix(nLenslet*2+1,tel.D,atm, gs, gs, 'mask', gridMask);
%  imagesc(cell2mat(Cxx) - CxxExp)


%% Yoshito method
if any(RecVec == 7)
    
    Cxx11 = Gamma*CxxExp*Gamma';
    Cxx12 = Gamma*H*CPhi*Htau'*Gamma';
    Cxx21 = Gamma*Htau*CPhi*H'*Gamma';
    Cxx22 = Gamma*Htau*CPhi*Htau'*Gamma';
    Ry = iF*Hss*CPhi*[Gamma*H;Gamma*Htau]'/([Cxx11 Cxx12;Cxx21 Cxx22] + blkdiag(CnZ, CnZ));
end
%% Explicit Layered LQG (one state)

if any(RecVec == 6)
    Sn = CPhi - AzPhi*CPhi*AzPhi';
    
    Q = Sn;
    Q = (Q+Q')./2;
    R = blkdiag(Cn(:,:,1),Cn(:,:,2),Cn(:,:,3));
    % R = blkdiag(Cn(:,:,1),Cn(:,:,2),Cn(:,:,3),Cn(:,:,4));
    
    Ad = AzPhi;
    Cd = Gamma*H;
    
    SigmaInf = solveRiccati_2xspeed(Ad,Cd,Q,R, 10);
    Minf = (SigmaInf*Cd')/(Cd*SigmaInf*Cd'+R);
    
    Mlqg = commandMatrix(:,:,1)*GammaBeta*Hss;
end

%% Spatio-Angular LQG (one state)

if any(RecVec == 5)
    C_alphaLag = phaseStats.spatioAngularCovarianceMatrix(nLenslet*2+1,tel.D,atmNeg, gs, gs, 'mask', gridMask,'lag', frameTime + fixedLagTimeInMs/1000);
    %C_alphaLagF = phaseStats.spatioAngularCovarianceMatrix(nLenslet*2+1,tel.D,atmNeg, gs, gs, 'mask', gridMask,'lag', frameTime );
    %C_alphaLagL = phaseStats.spatioAngularCovarianceMatrix(nLenslet*2+1,tel.D,atmNeg, gs, gs, 'mask', gridMask,'lag',  fixedLagTimeInMs/1000);

    
    %C_alphaLag = phaseStats.spatioAngularCovarianceMatrix(nLenslet*2+1,tel.D,atm, gs, gs, 'mask', gridMask,'lag', frameTime);
    C_alphaLag = cell2mat(C_alphaLag);
    Cxx = phaseStats.spatioAngularCovarianceMatrix(nLenslet*2+1,tel.D,atm, gs, 'mask', gridMask);
    
    A = C_alphaLag*pinv(Cxx);
    %AF = cell2mat(C_alphaLagF)*pinv(Cxx);
    %AL = cell2mat(C_alphaLagL)*pinv(Cxx);

    %Sn = Cxx - C_alphaLag/Cxx*C_alphaLag'; = Cxx - A*Cxx*A';
    Sn = Cxx - A*Cxx*A';
    %Sn = H*Sn*H'; % equivalent to using Explicit-LQG matrices

    Q = Sn;
    Q = (Q+Q')./2;
    R = blkdiag(Cn(:,:,1),Cn(:,:,2),Cn(:,:,3));
    % R = blkdiag(Cn(:,:,1),Cn(:,:,2),Cn(:,:,3),Cn(:,:,4));
    
    AdSA = A;
    %AdSA = H*Ad/H; % equivalent to using Explicit-LQG matrices
    
    CdSA = [Gamma];
    
    SigmaInfSA = solveRiccati_2xspeed(AdSA,CdSA,Q,R, 10);
    MinfSA = (SigmaInfSA*CdSA')*pinv(CdSA*SigmaInfSA*CdSA'+R);
    %MinfSA = H*Minf; % equivalent from using Explicit-LQG matrices
    
    MlqgSA = [iF*zonalMmseHR.Cox{1}/(zonalMmseHR.Cxx)];
end


%% How To Display turbulence and residual phase
if plotPhase == true
    plotID1 = RecVec(1);
    plotID2 = RecVec(2);
    figure(101)
    rad2mic = 1e6/sciObjectVec{plotID1}.waveNumber;
    h = imagesc([sciObjectVec{plotID1}.meanRmPhase,sciObjectVec{plotID2}.meanRmPhase]*rad2mic);
    
    ax = gca;
    axis equal tight
    ylabel(colorbar,'WFE [\mum]')
    
    figure(102)
    plotID3 = RecVec(end-1);
    plotID4 = RecVec(end);
    
    h2 = imagesc([sciObjectVec{plotID3}.meanRmPhase,sciObjectVec{plotID4}.meanRmPhase]*rad2mic);
    ax = gca;
    axis equal tight
    ylabel(colorbar,'WFE [\mum]')

end
% figure(102)
% h2 = imagesc([sciSA(1).meanRmPhase,sciSA(1).meanRmPhase,sciSA(1).meanRmPhase]*rad2mic);
%
% ax2 = gca;
% axis equal tight
% ylabel(colorbar,'WFE [\mum]')

%% Open Loop

nSlopes = wfs(1).nSlope;
tel = tel+atm;
fprintf('\n %5d',1)
S_NGS = zeros(nSlopes,nGs);

coefsVec = zeros(dmCal.nValidActuator,7);

for kIteration = 1:nIteration
    fprintf('\b\b\b\b\b%5d',kIteration)
    for iFrame = 1:fixedLagTimeInMs %Frame readout, computation, commands sent to DM (3ms)
        
        +tel;
        
        for iGs = 1:nGs;
            gs(iGs) = gs(iGs).*tel*wfs(iGs);
        end
        
        for iRec = 1:length(RecVec)
            RecID = RecVec(iRec);
            sciObjectVec{RecID} = sciObjectVec{RecID}.*tel*dmObjectVec{RecID};
        end
        
        if plotPhase == true
            set(h,'Cdata',[sciObjectVec{plotID1}.meanRmPhase,sciObjectVec{plotID2}.meanRmPhase]*rad2mic)
            title(ax,sprintf('#%4d/%4d',kIteration,nIteration))
            drawnow
            
            set(h2,'Cdata',[sciObjectVec{plotID3}.meanRmPhase,sciObjectVec{plotID4}.meanRmPhase]*rad2mic)
            title(ax,sprintf('#%4d/%4d',kIteration,nIteration))
            drawnow
        end
        
        for iRec = 1:length(RecVec)
            RecID = RecVec(iRec);
            sciObjectVec{RecID} = sciObjectVec{RecID}*imgrObjectVec{RecID};
        end
        
    end
    
    for iRec = 1:length(RecVec)
        RecID = RecVec(iRec);
        dmObjectVec{RecID}.coefs = -coefsVec(:,RecID);
    end
    
    for iFrame = 1:(frameTime/samplingTime-fixedLagTimeInMs) % the rest of the exposure
        
        +tel;
        
        for iGs = 1:nGs;
            gs(iGs) = gs(iGs).*tel*wfs(iGs);
        end
        
        for iRec = 1:length(RecVec)
            RecID = RecVec(iRec);
            sciObjectVec{RecID} = sciObjectVec{RecID}.*tel*dmObjectVec{RecID};
        end
        
        if plotPhase == true
            set(h,'Cdata',[sciObjectVec{plotID1}.meanRmPhase,sciObjectVec{plotID2}.meanRmPhase]*rad2mic)
            title(ax,sprintf('#%4d/%4d',kIteration,nIteration))
            drawnow
            
            set(h2,'Cdata',[sciObjectVec{plotID3}.meanRmPhase,sciObjectVec{plotID4}.meanRmPhase]*rad2mic)
            title(ax,sprintf('#%4d/%4d',kIteration,nIteration))
            drawnow
        end
        
        for iRec = 1:length(RecVec)
            RecID = RecVec(iRec);
            sciObjectVec{RecID} = sciObjectVec{RecID}*imgrObjectVec{RecID};
        end
        
    end
    
    for iGs = 1:nGs;
        slopeStack(:,iGs, kIteration) = wfs(iGs).slopes(:);
    end
    
    S_NGS = slopeStack(:,:,kIteration);
    
    if any(RecVec == 1)
        C_staticSA = RecStatSA*S_NGS(:);
        coefsVec(:,1) = C_staticSA;
    end
    
    if any(RecVec == 2)
        C_staticExp = RecStatExp*S_NGS(:);
        coefsVec(:,2) = C_staticExp;
    end
    
    if any(RecVec == 3)
        if kIteration > 1
            C_predSA = RecPredSA*S_NGS(:);
            coefsVec(:,3) = C_predSA + 0*coefsVec(:,3);
        else
            C_predSA = RecPredSA*S_NGS(:);
            coefsVec(:,3) = C_predSA;
        end

    end
    
    if any(RecVec == 4)
        C_predExp = RecPredExp*S_NGS(:);
        coefsVec(:,4) = C_predExp;
    end
    
    if any(RecVec == 5)
        if kIteration == 1
            kStateSA = MinfSA*S_NGS(:);
        else
            sHatSA = CdSA*kStateSA;
            kStateSA = kStateSA+MinfSA*(S_NGS(:)-sHatSA);
        end
        kStateSA = AdSA*kStateSA;        
        C_LQG_SA = MlqgSA*kStateSA;
        coefsVec(:,5) = C_LQG_SA;
        
    end
    
    if any(RecVec == 6)
        if kIteration ==1
            kStateExp = Minf*S_NGS(:);
        else
            sHatExp = Cd*kStateExp;
            kStateExp = kStateExp+Minf*(S_NGS(:)-sHatExp);
        end
        kStateExp = Ad*kStateExp;
        C_LQG_Exp = Mlqg*kStateExp;
        coefsVec(:,6) = C_LQG_Exp;
    end
    
    if any(RecVec == 7)
        if kIteration > 1
            S_NGS_minus1 = slopeStack(:,:,kIteration-1);
        else
            S_NGS_minus1= slopeStack(:,:,kIteration);
        end
        coefsVec(:,7) = Ry*[S_NGS(:); S_NGS_minus1(:)];
    end
    
end

for iRec = 1:length(RecVec)
    RecID = RecVec(iRec);
    strehl(RecID) = imgrObjectVec{RecID}.strehl;
    EE(RecID) = imgrObjectVec{RecID}.ee;
end