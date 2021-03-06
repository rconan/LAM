%% INTRODUCTION
% Multiple vignetted pupils
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
telMasked.pupil(35*nPx:41*nPx,1:7*nPx) = 0;


pupilsVignetted{1} = tel.pupil;
pupilsVignetted{2} = pupilsVignetted{1}(:,end:-1:1);
pupilsVignetted{3} = pupilsVignetted{1}';
pupilsVignetted{4} = pupilsVignetted{2}';

%pupilsVignetted{1} = tel.pupil;
%pupilsVignetted{2} = pupilsVignetted{1};
%pupilsVignetted{3} = pupilsVignetted{1};
%pupilsVignetted{4} = pupilsVignetted{1};



%% WAVE-FRONT SENSOR
% <latex>
% \subsubsection{The wavefront sensor}
% To be complete an \aos needs a wavefront sensor (WFS). \oomao uses the
% Shack--Hartmann \wfs. A Shack--Hartmann \wfs is made of a lenslet array
% and a detector. The numerical model follows the physical model by
% embedding a lenslet array (\oo{lensletArray}) class and a detector
% (\oo{detector}) class in the \oo{shackHartmann} class. The
% \oo{lensletArray} class perform the numerical Fraunhoffer propagation of
% the wavefront to the detector. The \oo{detector} class implements a CCD
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

% calibration DM
wfsCal = shackHartmann(nL,nRes,0.85);
ngs = ngs.*tel*wfsCal;

wfsCal.INIT

+wfsCal;

%
nGs = 4;
for iGs = 1:nGs
    wfs(iGs) = shackHartmann(nL,nRes,0.85);
    
    % WFS initialisation
    ngs = ngs.*tel*wfs(iGs);
    
    wfs(iGs).INIT
    
    +wfs(iGs);
end

figure
imagesc(wfs(1).camera,'parent',subplot(3,2,[1,4]))
slopesDisplay(wfs(1),'parent',subplot(3,2,[5,6]))

% The next 2 commands allow the displays of the frame and of the slopes to
% be updated when a new frame and new slopes are computed
%wfs.camera.frameListener.Enabled = true;
%wfs.slopesListener.Enabled = true;
%% WFS must be calibrated such as for 1rd of tip--tilt wavefront , it will measured a slopes of 1rd.
wfsCal =  shackHartmann(nL,nRes,0.85);
ngs = ngs.*tel*wfsCal;
wfsCal.INIT;
+wfsCal;


wfsCal.pointingDirection = zeros(2,1);

pixelScale = ngs.wavelength/...
    (2*d*wfsCal.lenslets.nyquistSampling);
tipStep = pixelScale/2;
nStep   = floor(nPx/3)*2;
sx      = zeros(1,nStep+1);
u       = 0:nStep;
wfsCal.camera.frameListener.Enabled = false;
wfsCal.slopesListener.Enabled = false;
warning('off','oomao:shackHartmann:relay')
for kStep=u
    ngs.zenith = -tipStep*kStep;
    +ngs;
    drawnow
    sx(kStep+1) = median(wfsCal.slopes(1:end/2));
end
warning('on','oomao:shackHartmann:relay')
Ox_in  = u*tipStep*constants.radian2arcsec;
Ox_out = sx*ngs.wavelength/d/2*constants.radian2arcsec;
%figure
%plot(Ox_in,Ox_out)
%grid
slopesLinCoef = polyfit(Ox_in,Ox_out,1);
wfsCal.slopesUnits = 1/slopesLinCoef(1);

for iGs = 1:nGs
    wfs(iGs).slopesUnits = 1/slopesLinCoef(1);
    wfs(iGs).camera.photonNoise = 1;
    wfs(iGs).camera.readOutNoise = 1;
    wfs(iGs).framePixelThreshold = wfs(iGs).camera.readOutNoise;
end
ngs.zenith = 0;
wfsCal.pointingDirection = [];
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
    'validActuator',wfsCal.validActuator);

bifaLowRes = influenceFunction('monotonic',dmCrossCouplingCoeff);
dmLowRes = deformableMirror(nL+1,'modes',bifaLowRes,'resolution',nL+1,...
    'validActuator',wfsCal.validActuator);
%% INTERACTION MATRIX
wfsCal.camera.frameListener.Enabled = false;
wfsCal.slopesListener.Enabled = false;

ngs = ngs.*tel;

calibDm = calibration(dm,wfsCal,ngs,ngs.wavelength,nL+1,'cond',1e2);

%% GENERATE ATMOSPHERE
tel = tel + atm;
telFull = telFull + atm;
figure
imagesc(tel)
ngs = ngs.*tel*wfsCal;
%% LOW-RES TELESCOPE TO ESTIMATE PHASE 

telLowRes = telescope(tel.D,'resolution',nL+1,...
    'obstructionRatio',obstructionRatio,'fieldOfViewInArcsec',fieldOfViewInArcsec,'samplingTime',1/samplingFreq);
telLowRes.pupil = wfsCal.validActuator;

telLowRes= telLowRes + atm;
ngs = ngs.*telLowRes;
phase = ngs.meanRmOpd;
%% LGS SOURCES
%lgsAst = source('asterism',{[nGs,arcsec(45),0]},'height',90e3,'magnitude',10);

% single altitude object
lgsAst = laserGuideStar(tel.D/nL,tel.D, 90e3, [], 2e6, [],'asterism',{[nGs,arcsec(5),0]}, 'wavelength', photometry.Na,'height',9e4);
theta = linspace(0,2*pi-2*pi/nGs,nGs);
for iGS = 1:nGs
    lgsAst(iGS).viewPoint = D/2*[cos(theta(iGS)), sin(theta(iGS))];
end
% ngsAst = source('asterism',{[nGs,arcsec(45),0]},'magnitude',10);
% ngsAst = ngsAst.*tel;
% 
% for iGs=1:nGs
%     s(:,:,iGs) = pupilsVignetted{iGs};
% end
% S = {s, zeros(nRes, nRes,nGs)};
% ngsAst = ngsAst.*tel*S;



%lgsAst.magnitude = 12;
% figure, imagesc(tel, [ngs,lgsAst])

wfsCal.camera.photonNoise = 1;
wfsCal.camera.readOutNoise = 3; % just for regularisation prurposes
wfsCal.framePixelThreshold = wfsCal.camera.readOutNoise;
lgsAst_slmmse = slopesLinearMMSE(wfsCal,tel,atm,lgsAst,'mmseStar',ngs,'NF',1024);


%% SCIENCE CAMERA

science = source('wavelength',photometry.K);
cam = imager();
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

gain_pol = 0.7;
F = 2*bifaLowRes.modes(wfsCal.validActuator,:);
iF = pinv(full(F),1e-1);


%% Laser Tomography Adaptive Optics

reset(tel)
dm.coefs = zeros(dm.nValidActuator,1);
%%

science = science.*telFull*dm*cam;
%lgsAst = lgsAst.*tel*dm*wfs;
for iGS=1:nGs
    lgsAst(iGs) = lgsAst(iGs).*tel*{pupilsVignetted{iGs},zeros(nRes)}*dm*wfs(iGs);
end
figure(31416)
imagesc(cam,'parent',subplot(2,1,1))
subplot(2,1,2)
h = imagesc(catMeanRmPhase(science));
axis xy equal tight
colorbar
wfsCal.camera.frameListener.Enabled = false;
wfsCal.slopesListener.Enabled = false;
%%
flush(cam)
cam.clockRate    = 1;
exposureTime     = 220;
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
    %+lgsAst;
    for iGs=1:nGs
        lgsAst(iGs) = lgsAst(iGs).*tel*{pupilsVignetted{iGs},zeros(nRes)}*dm*wfs(iGs);
        %lgsAst(iGs) = lgsAst(iGs).*tel*dm*wfs(iGs);
        slopesk(:,iGs) = wfs(iGs).slopes;
    end
    +science;
    % Pseudo-open-loop controller
    
    dm.coefs = (1-gain_pol)*dm.coefs + ...
        gain_pol*iF*( lgsAst_slmmse*( bsxfun( @minus, slopesk, calibDm.D*dm.coefs ) ) );
    % Display
         set(h,'Cdata',catMeanRmPhase(science))
         drawnow
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
