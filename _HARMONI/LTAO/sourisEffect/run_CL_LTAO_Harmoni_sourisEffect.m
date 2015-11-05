%% INTRODUCTION
%% SOURCES
% A simple on--axis \ngs object is created by calling the \oo{source} constructor without parameters:
ngs =source;
%%
% <latex>
% The first time a object from an \oomao class is created a message let the user aware he is using the library.
% A summary of the object main parameters is also displayed.
% In this example, zenith and azimuth angle are both set to zero, the star height is infinite, the wavelength is 550nm and the magnitude is set to 0.
% </latex>
%% ATMOSPHERE

% Next, we create the atmosphere the source will propagate through.

atm = atmosphere(photometry.V,16e-2,30,...
    'fractionnalR0',[0.5,0.3,0.2],'altitude',[0e3,5e3,12e3],...
    'windSpeed',[10,5,20],'windDirection',[0,pi/2,pi]);

%% TELESCOPE
nL   = 74;          % number of lenslets
nPx  = 10;          % number of pixels per lenslet
nRes = nL*nPx;      % resolution on the pupil plane (no of pixels)
D    = 37;          % telescope primary mirror diameter
d    = D/nL;        % lenslet pitch
samplingFreq = 500; % WFS sampling time
obstructionRatio=0.3;% central obscuration ratio
fieldOfViewInArcsec = 120; %fieldOfViewInArcsec
tel = telescope(D,'resolution',nRes,...
    'obstructionRatio',obstructionRatio,'fieldOfViewInArcsec',fieldOfViewInArcsec,'samplingTime',1/samplingFreq);

% cut off a portion of the pupil that is vignetted
%tel.pupil(35*nPx:42*nPx,1:7*nPx) = 0;

telFull = telescope(D,'resolution',nRes,...
    'obstructionRatio',0.3,'fieldOfViewInArcsec',fieldOfViewInArcsec,'samplingTime',1/samplingFreq);


telMasked = telescope(D,'resolution',nRes,...
    'obstructionRatio',0.3,'fieldOfViewInArcsec',fieldOfViewInArcsec,'samplingTime',1/samplingFreq);

% cut off a portion of the pupil that is vignetted
telMasked.pupil(37*nPx:40*nPx,1:3*nPx) = 0;
%telMasked.pupil = tel.pupil;
%% WAVE-FRONT SENSOR
% <latex>
% \subsubsection{The wavefront sensor}
% To be complete an \aos needs a wavefront sensor (WFS). \oomao uses the
% Shack--Hartmann \wfs. A Shack--Hartmann \wfs is made of a lenslet array
% and a detector. The numerical model follows the physical model by
% embedding a lenslet array (\oo{lensletArray}) class and a detector
% (\oo{detector}) class in the \oo{shackHartmann} class. The
% \oo{lensletArray} class perform the numerical Fraunhoffer propagation of
% opethe wavefront to the detector. The \oo{detector} class implements a CCD
% camera detection process including Poisson and read--out noise. The
% lenslet images are Nyquist sampled per default.
% \newline
% A \oo{shackHartmann} object with a $\texttt{nL}\times \texttt{nL}$
% lenslet array and a $\texttt{nRes}\times \texttt{nRes}$ resolution camera
% is created.
% The \oop{lensletArray}{minLightRatio} property of the \oo{lensletArray}
% object set the minimum ratio of light intensity between a partially and
% fully illuminated lenslet.
% In the following, \oop{lensletArray}{minLightRatio} is set to 85\%:
% </latex>
wfs = shackHartmann(nL,nRes,0.85);

% WFS initialisation
ngs = ngs.*tel*wfs;

wfs.INIT

+wfs;
figure
imagesc(wfs.camera,'parent',subplot(3,2,[1,4]))
slopesDisplay(wfs,'parent',subplot(3,2,[5,6]))

% The next 2 commands allow the displays of the frame and of the slopes to
% be updated when a new frame and new slopes are computed
wfs.camera.frameListener.Enabled = true;
wfs.slopesListener.Enabled = true;
%% WFS must be calibrated such as for 1rd of tip--tilt wavefront , it will measured a slopes of 1rd.

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

ngs.zenith = 0;
wfs.pointingDirection = [];
%% DEFORMABLE MIRROR
dmCrossCouplingCoeff = 0.4;
bifa = influenceFunction('monotonic',dmCrossCouplingCoeff);
figure,show(bifa,'parent',subplot(1,2,1))
title('Monototic influence function')
bifb = influenceFunction('overshoot',dmCrossCouplingCoeff);
show(bifb,'parent',subplot(1,2,2))
title('Overshooted influence function')

% A deformable mirror with $(\texttt{nL}+1)\times (\texttt{nL}+1)$ actuators
% and the previously defined influence function sampled at the same
% resolution than the telescope is created with
% </latex>
dm = deformableMirror(nL+1,'modes',bifa,...
    'resolution',tel.resolution,...
    'validActuator',wfs.validActuator);

bifaLowRes = influenceFunction('monotonic',dmCrossCouplingCoeff);
dmLowRes = deformableMirror(nL+1,'modes',bifaLowRes,'resolution',nL+1,...
    'validActuator',wfs.validActuator);
%% INTERACTION MATRIX
wfs.camera.frameListener.Enabled = false;
wfs.slopesListener.Enabled = false;

ngs = ngs.*tel;

calibDm = calibration(dm,wfs,ngs,ngs.wavelength,nL+1,'cond',1e2);

%% GENERATE ATMOSPHERE
tel = tel + atm;
telFull = telFull + atm;
figure
imagesc(tel)
ngs = ngs.*tel*wfs;
ngs = ngs.*tel*{telMasked.pupil, zeros(nRes)}*wfs;
%% LOW-RES TELESCOPE TO ESTIMATE PHASE 

telLowRes = telescope(tel.D,'resolution',nL+1,...
    'obstructionRatio',obstructionRatio,'fieldOfViewInArcsec',fieldOfViewInArcsec,'samplingTime',1/samplingFreq);
telLowRes.pupil = wfs.validActuator;

telLowRes= telLowRes + atm;
ngs = ngs.*telLowRes;
phase = ngs.meanRmOpd;
%% LGS SOURCES
lgsAst = source('asterism',{[6,arcsec(45),0]},'height',90e3,'magnitude',10);
%lgsAst.magnitude = 12;
% figure, imagesc(tel, [ngs,lgsAst])
lgsAst_slmmse = slopesLinearMMSE(wfs,tel,atm,lgsAst,'mmseStar',ngs,'NF',1024);


%% SCIENCE CAMERA

science = source('wavelength',photometry.K);
cam = imager(telFull);
telFull = telFull - atm;
science = science.*telFull*cam;
figure(31416)
imagesc(cam,'parent',subplot(2,1,1))

cam.referenceFrame = cam.frame;
+science;
fprintf('Strehl ratio: %4.1f\n',cam.strehl)

tel = tel + atm;
telFull = telFull + atm;
+science;
fprintf('Strehl ratio: %4.1f\n',cam.strehl)

%% INTEGRATOR

gain_cl  = 0.5;
%% PSEUDO OPEN-LOOP CONTROLLER

gain_pol = 0.6;
F = 2*bifaLowRes.modes(wfs.validActuator,:);
iF = pinv(full(F),1e-1);


%% Laser Tomography Adaptive Optics

reset(tel)
dm.coefs = zeros(dm.nValidActuator,1);
%%

science = science.*telFull*dm*cam;

% add noisy measurements
wfs.camera.photonNoise = 1;
wfs.camera.readOutNoise = 1;


%lgsAst = lgsAst.*tel*{telMasked.pupil, zeros(nRes)}*wfs;
lgsAst = lgsAst.*tel*{telMasked.pupil, zeros(nRes)}*dm*wfs;
%lgsAst = lgsAst.*tel*dm*wfs;
figure(31416)
imagesc(cam,'parent',subplot(2,1,1))
subplot(2,1,2)
h = imagesc(catMeanRmPhase(science));
axis xy equal tight
colorbar
wfs.camera.frameListener.Enabled = false;
wfs.slopesListener.Enabled = false;


%%
flush(cam)
cam.clockRate    = 1;
exposureTime     = 100;
cam.exposureTime = exposureTime;
startDelay       = 20
flush(cam)
cam.startDelay = startDelay;
set(science,'logging',true)
set(science,'phaseVar',[])
lgsAst_slmmse.wavefrontSize = [dm.nValidActuator,1];
lgsAst_slmmse.warmStart = true;
cam.frameListener.Enabled = true;
%% The loop is closed for one full exposure of the science camera.
nIteration = startDelay + exposureTime;
for k=1:cam.startDelay + cam.exposureTime
    tic
    k
    % Objects update
    +tel;
    +lgsAst;
    +science;
    % Pseudo-open-loop controller
    dm.coefs = (1-gain_pol)*dm.coefs + ...
        gain_pol*iF*( lgsAst_slmmse*( bsxfun( @minus, wfs.slopes, calibDm.D*dm.coefs ) ) );
    % Display
    %     set(h,'Cdata',catMeanRmPhase(science))
    %     drawnow
    toc
end
imagesc(cam)
set(h,'Cdata',catMeanRmPhase(science))
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
marechalStrehl_ltao = 1e2*exp(-mean(var_wfe_ltao(startDelay:end,2)));
psfStrehl_ltao =1e2*cam.strehl;
figure(fix(exp(1)*1e2))
u = 1:nIteration;
semilogy(u,wfe_lsq,u,wfe_lmmse(:,2),u,wfe_ltao(:,2))
line([1,nIteration],ones(1,2)*atm_wfe_rms,...
    'color','k','LineStyle','--','linewidth',2)
grid
xlabel('Time [2ms]')
ylabel('Wavefront rms [micron]')
legend('Full',...
    sprintf('LSR Residue  : SR: %2.0f(%2.0f)%%',psfStrel(1),marechalStrehl_lsq),...
    sprintf('LMMSE Residue: SR: %2.0f(%2.0f)%%',psfStrel(2),marechalStrehl_lmmse),...
    sprintf('LTAO Residue:  SR: %2.0f(%2.0f)%%',psfStrehl_ltao, ...
            marechalStrehl_ltao),0)
