%% DEFS
exposureTime     = 100;
startDelay       = 20;

%% SOURCE
ngs =source;

%% ATMOSPHERE

atm = atmosphere(photometry.V,15.87e-2,50,...
    'fractionnalR0',[0.5,0.3,0.2],'altitude',[0e3,5e3,12e3],...
    'windSpeed',[10,5,20],'windDirection',[0,pi/2,pi]);
%% TELESCOPE
nL   = 60;          % number of lenslets
nPx  = 10;          % number of pixels per lenslet
nRes = nL*nPx;      % resolution on the pupil plane (no of pixels)
D    = 25;          % telescope primary mirror diameter
d    = D/nL;        % lenslet pitch
samplingFreq = 500; % WFS sampling time
obstructionRatio= 0;% central obscuration ratio
fieldOfViewInArcsec = 30; %fieldOfViewInArcsec
tel = telescope(D,'resolution',nRes,...
    'obstructionRatio',obstructionRatio,'fieldOfViewInArcsec',fieldOfViewInArcsec,'samplingTime',1/samplingFreq);


%% HIGH-ORDER SH WFS
wfs = shackHartmann(nL,nRes,0.85);
%wfs.lenslets.nyquistSampling = 0.5;
%ngs.wavelength = photometry.R;

ngs = ngs.*tel*wfs;

wfs.INIT

+wfs;
figure
imagesc(wfs.camera,'parent',subplot(3,2,[1,4]))
slopesDisplay(wfs,'parent',subplot(3,2,[5,6]))

wfs.camera.frameListener.Enabled = true;
wfs.slopesListener.Enabled = true;
% WFS gain calibration such as for 1rd of tip--tilt wavefront , it will measured a slopes of 1rd.

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

% reset source on-axis
ngs.zenith = 0;
wfs.pointingDirection = [];


%% DEFORMABLE MIRROR
dmCrossCouplingCoef = 0.8;
bifa = influenceFunction('monotonic',dmCrossCouplingCoef);
figure,show(bifa,'parent',subplot(1,2,1))
title('Monototic influence function')
bifb = influenceFunction('overshoot',dmCrossCouplingCoef);
show(bifb,'parent',subplot(1,2,2))
title('Overshooted influence function')

dm = deformableMirror(nL+1,'modes',bifa,...
    'resolution',tel.resolution,...
    'validActuator',wfs.validActuator);


wfs.camera.frameListener.Enabled = false;
wfs.slopesListener.Enabled = false;


bifaLowRes = influenceFunction('monotonic',dmCrossCouplingCoef);
dmLowRes = deformableMirror(nL+1,'modes',bifaLowRes,'resolution',nL+1,...
    'validActuator',wfs.validActuator);


%% INTERACTION MATRIX
ngs = ngs.*tel;     %setup the optical path before the DM/WFS subsystem

calibDm = calibration(dm,wfs,ngs,ngs.wavelength,nL+1,'cond',1e2);

%% GENERATE ATMOSPHERE
tel = tel + atm;
figure
imagesc(tel)
ngs = ngs.*tel*wfs;


%% LOW-RES TO RECONSTRUCT PHASE AT THE TIPS OF validActuator
telLowRes = telescope(tel.D,'resolution',nL+1,...
    'obstructionRatio',obstructionRatio,'fieldOfViewInArcsec',fieldOfViewInArcsec,'samplingTime',1/samplingFreq);
telLowRes.pupil = wfs.validActuator;

%% GENERATE LOW-RES ATM
telLowRes= telLowRes + atm;
ngs = ngs.*telLowRes;
phase = ngs.meanRmOpd;

%% LGS SOURCES AND MMSE RECONSTUCTOR
lgsAst = source('asterism',{[6,arcsec(10),0]},'height',90e3);
%lgsAst = source('asterism',{[1,arcsec(15),0]},'height',90e3);

% figure, imagesc(tel, [ngs,lgsAst])
lgsAst_slmmse = slopesLinearMMSE(wfs,tel,atm,lgsAst,'mmseStar',ngs,'NF',1024);

%%
lgsAst = lgsAst.*tel*wfs;
lgsAst_ps_e = tools.meanSub( lgsAst_slmmse*wfs , wfs.validActuator );
ngs = ngs.*telLowRes*{wfs.validActuator,-lgsAst_ps_e*ngs.waveNumber};
lgsAst_ps_eRes = ngs.meanRmOpd;
lgsAst_ps_eResRms = ngs.opdRms*1e9;
 
%% TT SOURCES
%TTAst = source('asterism',{[3,arcsec(20),0]},'wavelength',photometry.H);
TTAst = source('asterism',{[1,arcsec(0),0]},'wavelength',photometry.H);
%TTAst = TTAst .* tel*OIWFS;

%% SCIENCE SOURCE
science = source('wavelength',photometry.J);
%% SCIENCE CAMERA

cam = imager(tel);
cam.clockRate    = 1;

tel = tel - atm;
science = science.*tel*cam;
figure(31416)
imagesc(cam,'parent',subplot(2,1,1))
%cam.frameListener.Enabled = true;

cam.referenceFrame = cam.frame;
+science;
fprintf('Strehl ratio: %4.1f\n',cam.strehl)


tel = tel + atm;
+science;
fprintf('Strehl ratio: %4.1f\n',cam.strehl)

%% LOW-RES DM INFLUENCE MATRIX
F = 2*bifaLowRes.modes(wfs.validActuator,:);
iF = pinv(full(F),1e-1);

%%
dm.coefs = zeros(dm.nValidActuator,1);

science = science.*tel*dm*cam;
lgsAst = lgsAst.*tel*dm*wfs;
figure(31416)
imagesc(cam,'parent',subplot(2,1,1))
subplot(2,1,2)
h = imagesc(catMeanRmPhase(science));
axis xy equal tight
colorbar
wfs.camera.frameListener.Enabled = false;
wfs.slopesListener.Enabled = false;

%TTAst = TTAst.*tel*dm*OIWFS;
%% LOOP SETUP

flush(cam)
cam.frame = cam.frame*0;
cam.exposureTime = exposureTime;
cam.startDelay = startDelay;

set(science,'logging',true)
set(science,'phaseVar',[])
lgsAst_slmmse.wavefrontSize = [dm.nValidActuator,1];
lgsAst_slmmse.warmStart = true;
cam.frameListener.Enabled = false;

gain_pol = 0.7;

gain_cl = 0.4;
u_ngs = zeros(dmLowRes.nValidActuator,1);

%OIWFS.camera.photonNoise = true;
%OIWFS.camera.readOutNoise = 0;

wfs.camera.photonNoise = 0;
%wfs.camera.readOutNoise = 0;

dm.coefs = zeros(dm.nValidActuator,1);
%% piston-tip-tilt removal 
t = ones(wfs.nValidLenslet,1);
TT = [t 0*t; 0*t,t];
SlopeTTRem = eye(2*wfs.nValidLenslet) - TT*pinv(TT);

zern = zernike(2:3,'resolution',nL+1, 'pupil',dmLowRes.validActuator);
TT = iF*zern.modes(dmLowRes.validActuator,:);
DMTTRem = eye(dmLowRes.nValidActuator) - TT*pinv(TT);

%% TT-INCLUDED HO LOOP 
% <latex>
% The loop is closed for one full exposure of the science camera.
% </latex>
nIteration = startDelay + exposureTime;

for k=1:cam.startDelay + cam.exposureTime
    tic
    k
    % Objects update
    +tel;
    +lgsAst;
    %+TTAst;
    +science;
    % Pseudo-open-loop controller  -- TTR HO LOOP   
    dm.coefs = (1-gain_pol)*dm.coefs + ...
        gain_pol*iF*( lgsAst_slmmse*( bsxfun( @minus, wfs.slopes, calibDm.D*dm.coefs ) ) );

    
    %dm.coefs = iF*( lgsAst_slmmse*( bsxfun( @minus, SlopeTTRem*wfs.slopes, calibDm.D*dm.coefs ) ) );


    % Display
     set(h,'Cdata',catMeanRmPhase(science))
     drawnow
     toc
end
imagesc(cam)
set(h,'Cdata',catMeanRmPhase(science))

%% TTR HO LOOP WITH TT FROM TT LO LOOP
% % <latex>
% % The loop is closed for one full exposure of the science camera.
% % </latex>
% nIteration = startDelay + exposureTime;
% 
% for k=1:cam.startDelay + cam.exposureTime
%     tic
%     % Objects update
%     +tel;
%     +lgsAst;
%     +TTAst;
%     +science;
%     % Pseudo-open-loop controller  -- TTR HO LOOP   
%     dm.coefs = (1-gain_pol)*dm.coefs + ...
%         gain_pol*iF*( lgsAst_slmmse*( bsxfun( @minus, SlopeTTRem*wfs.slopes, calibDm.D*dm.coefs ) ) );
% 
% 
% %dm.coefs = iF*( lgsAst_slmmse*( bsxfun( @minus, SlopeTTRem*wfs.slopes, calibDm.D*dm.coefs ) ) );
% 
% % Remove TT from DM
%     dm.coefs = DMTTRem*dm.coefs;
%     % Add TT from TT stars, % Closed-loop controller
%     u_ngs = u_ngs - TTAst(1).wavelength/8*gain_cl*TT*mean(QuadCell.slopes,2); 
%     
%     % add ngs controls back to dm commands
%     dm.coefs = dm.coefs - u_ngs;
%     % Display
%      set(h,'Cdata',catMeanRmPhase(science))
%      drawnow
%      toc
% end
% imagesc(cam)
% set(h,'Cdata',catMeanRmPhase(science))

%%
% <latex>
% The time series of wavefront variance is read from the
% \oop{stochasticWave}{phaseVar} property of the \matcall{science}
% object. 
% The Strehl ratio is estimated from the residual phase variance using the
% Marechal approximation and it is compared to the Strehl ratio derived
% from the long exposure image.
% </latex>
var_wfe_ltao = reshape(science.phaseVar(1:nIteration*2),2,[])';
wfe_ltao = sqrt(var_wfe_ltao)/science.waveNumber*1e6;
marechalStrehl_ltao = 1e2*exp(-mean(var_wfe_ltao(startDelay:end,2)))
psfStrehl_ltao =1e2*cam.strehl


