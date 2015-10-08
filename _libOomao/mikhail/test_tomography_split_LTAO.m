clc
close all
clear all
tool_inside_path_adder('matlab_code/mikhail', 'matlab_code');

% Undefined variable "logBook" or class "logBook.checkIn".
% 
% Error in source (line 202)
%                 set(obj,'log',logBook.checkIn(obj));
% 
% Error in test_tomography_split_LTAO (line 96)
% lgsAst = source('asterism',{[3,arcsec(10),0]},'height',90e3, 'magnitude', params.guidestar_magnitude);
 

% addpath('/home/mkonnik/matlab/OOMAO_LAM/oomao_raven/OOMAO-Raven/ravenOOMAOlib')
% Undefined function or variable 'arcsec'.
% 
% Error in test_tomography_split_LTAO (line 91)
% lgsAst = source('asterism',{[3,arcsec(10),0]},'height',90e3, 'magnitude', params.guidestar_magnitude);

%  @(source)> Terminated!
% Error using repmat
% Replication factors must be a row vector of integers or integer scalars.
% 
% Error in shackHartmann (line 144)
%             obj.p_referenceSlopes = ...
% 
% Error in tool_calibrate_wfs_oomao (line 8)
% wfs = shackHartmann(nLenslet,nRes,minLightRatio );
% 
% Error in test_tomography_split_LTAO (line 166)
% QuadCell = tool_calibrate_wfs_oomao(1,8,1, tel, D, 'H');  %% since TTast is in H-band, I'm going to calibrate Quadcell in
% H-band, too


% addpath('/home/mkonnik/matlab/OOMAO_LAM/oomao_raven/OOMAO-Raven/OOMAOlibUpdated')
%  @(source)> Terminated!
% Error using repmat
% Replication factors must be a row vector of integers or integer scalars.
% 
% Error in shackHartmann (line 125)
%             obj.p_referenceSlopes = ...
% 
% Error in tool_calibrate_wfs_oomao (line 8)
% wfs = shackHartmann(nLenslet,nRes,minLightRatio );
% 
% Error in test_tomography_split_LTAO (line 139)
% QuadCell = tool_calibrate_wfs_oomao(1,8,1, tel, D, 'H');  %% since TTast is in H-band, I'm going to calibrate Quadcell in
% H-band, too

samplingTime     = 0.002;% seconds
frameTime          = 0.002;%[0.01,0.015,0.020,0.025]; seconds
exposureTime     = 100; %% [msek], equals to the Number of Iterations
startDelay           = 0; %% skip this number of frames

params = struct; %%% creating the Structure with additional parameters for the computeARmodel

% params.guidestar_type = 'NGS';
params.guidestar_type = 'LGS';

        params.guidestar_magnitude = 0;

        
params.noise_in_wfs  = 1; %% adding the noise into WFS and Quadcell (same noise everywhere)
% params.noise_in_wfs  = 0; %% noiseless WFS and Quadcell 

params.account4noise_in_TTProjections = 1; 
% params.account4noise_in_TTProjections = 0; 



%% Noisy WFS case
switch params.guidestar_type

    case 'LGS'
    params.noise.photonNoise = true;
%     params.noise.wfs_readOutNoise = 1e-6;
    params.noise.wfs_readOutNoise_SNR = 80;  %% use the mean (or median) Signal value to determine the readout noise level

    params.noise.framePixelThreshold = 0;

    case 'NGS'
    params.noise.photonNoise = true;
% % %     params.noise.wfs_readOutNoise = 1.3e5;
    params.noise.wfs_readOutNoise_SNR = 80;  %% use the mean (or median) Signal value to determine the readout noise level
    params.noise.framePixelThreshold = 0;
    
    otherwise
        error('Wrong type of GuideStar or the type is not defined')
end %% for switch

    
%% No Noise  (noiseless) WFS case:
if  (params.noise_in_wfs  == 0); %% display the diagnostic plots during the simulation.
    params.noise.wfs_readOutNoise = 0;
    params.noise.photonNoise = false;
    params.noise.framePixelThreshold = 0;
end

params.noise_account_nMeasurements =   500;%% number of measurements for the covariance matrix estimation.



if (params.account4noise_in_TTProjections == 0)
        params.noise_account_in_SlopesTT  = 0; 
        params.noise_account_in_DMTT  = 0;
else 
        params.noise_account_in_SlopesTT  = 1;     
        params.noise_account_in_DMTT  = 1;
end



switch params.guidestar_type

    case 'LGS'
%% LGS asterism

lgsAst = source('asterism',{[3,arcsec(10),0]},'height',90e3, 'magnitude', params.guidestar_magnitude);
% lgsAst = source('asterism',{[1,arcsec(10),0]},'height',90e3, 'magnitude', params.guidestar_magnitude);


    case 'NGS'
%% Temprorary use NGS asterism as LGS does not work!
lgsAst = source('asterism',{[1,arcsec(10),0]}); %%% this is NGS star!
% lgsAst = source('asterism',{[3,arcsec(20),0]},'height',inf); %%% this is NGS star! ____DOES NOT WORK!___ See Issue #14
end


%% Tip-tilt asterism (for TT measurements)
% TTAst = source('asterism',{[1,arcsec(10),0]},'wavelength',photometry.H);

TTAst = source('asterism',{[1,arcsec(20),0]},'wavelength',photometry.H);



%% Atmospheric parameters
wavelength      = photometry.R;
altitudes       = [0,    5]*1e3;
fractionalR0    = [0.7, 0.3];
windSpeed       = [2,  5];
windDirection   = [pi, pi];

L0 = 30; %% outer scale, [m]
r0 = 0.20; %% Fried parameter, [m]

% Creating the Atmosphere object
atm = atmosphere(photometry.V, r0, L0,...
    'fractionnalR0',fractionalR0, ...
    'altitude',     altitudes,...
    'windSpeed',    windSpeed,...
    'windDirection',windDirection);

nLayers = atm.nLayer;

    params.show_figures = 0; %% display the diagnostic plots during the simulation.
    params.ircsSlitWidth           = 0.14; %%% % arcsec  EEwidth for the EE estimation in Imager
    params.dm_crosscoupling_coef = 0.40; %% cross-coupling degree for the deformable mirror.


%% TELESCOPE
nLenslet   = 10;
nPx  = 12;
nRes = nLenslet*nPx;
minLightRatio = 0.65; %% large values (more than 0.8  can lead to poor reconstruction - you are truncating the phase!).

D    = 8;
d    = D/nLenslet; % lenslet pitch

tel = telescope(D,...
    'resolution',nRes,...
    'fieldOfViewInArcsec',120,...
    'samplingTime', samplingTime);

%% Quadcell calibration
% QuadCell = tool_calibrate_wfs_oomao(1,32,minLightRatio, tel, D, 'H');  %% I don't know why do we calibrate QuadCell in K band.....
% QuadCell = tool_calibrate_wfs_oomao(1,8,minLightRatio, tel, D, 'R');  %% I don't know why do we calibrate QuadCell in K band.....

% QuadCell = tool_calibrate_wfs_oomao(1,16,1, tel, D, 'H');  %% since TTast is in H-band, I'm going to calibrate Quadcell in H-band, too
QuadCell = tool_calibrate_wfs_oomao(1,8,1, tel, D, 'H');  %% since TTast is in H-band, I'm going to calibrate Quadcell in H-band, too

% QuadCell = tool_calibrate_wfs_oomao(1,4,1, tel, D, 'R');  %% I don't know why do we calibrate QuadCell in K band.....
% QuadCell = tool_calibrate_wfs_oomao(1,8,minLightRatio, tel, D, 'K');  %% I don't know why do we calibrate QuadCell in K band.....

%% SH WFS calibration
wfs = tool_calibrate_wfs_oomao(nLenslet,nPx,minLightRatio , tel, d, 'R');


%% Deformable mirror
nActuator = nLenslet + 1;
bif = influenceFunction('monotonic',params.dm_crosscoupling_coef);

dm = deformableMirror(nActuator,...
    'modes',bif,...
    'resolution',nRes,...
    'validActuator',wfs.validActuator );

wfs.camera.frameListener.Enabled = false;
wfs.slopesListener.Enabled = false;


bifaLowRes = influenceFunction('monotonic',params.dm_crosscoupling_coef);

dmLowRes = deformableMirror(nLenslet+1,...  %%% Why do we need a DM with surface of 11x11?? Is this for the TT-mirror?
    'modes',bifaLowRes,...
    'resolution',nLenslet+1,... 
    'validActuator',wfs.validActuator);


%% DM Calibration, I suppose? 
ngs = source;  %% Why do we need this Guidestar??? To remove  average tilt from wavefront? Or just to calibrate DM?
ngs.magnitude = 0;
ngs = ngs.*tel;

calibDm = calibration(dm,wfs,ngs,ngs.wavelength,nLenslet+1,'cond',1e1);

tel = tel + atm;



 %% Linear Minimum Mean Square Error (LMMSE) wavefront reconstructor
lgsAst_slmmse = slopesLinearMMSE(wfs,tel,atm,lgsAst,'mmseStar',ngs,'NF',512); 
lgsAst_slmmse.wavefrontSize = [dm.nValidActuator,1];
lgsAst_slmmse.warmStart = true;

TTAst = TTAst .* tel*QuadCell;


%% Science targets
sciStar = source('wavelength',photometry.H);
cam = imager(tel);

% The \oo{atmosphere} object is detached from the telescope and the
% science star is propagated through the telescope to the science camera
% producing a perfect diffraction limited image.
tel = tel - atm;
sciStar = sciStar.*tel*cam;

cam.referenceFrame = cam.frame;
+sciStar;

tel = tel + atm;
+sciStar;


%% Laser Tomography Adaptive Optics
sciStar = sciStar.*tel*dm*cam;
lgsAst = lgsAst.*tel*dm*wfs;


%% BEGIN::  Changing the way how we add the read-out NOISE - see Issue 14
wfs.camera.photonNoise = false;
wfs.camera.readOutNoise = 0;

+lgsAst

z_wfsCameraFrame_Noiseless = wfs.camera.frame;
z_wfsCameraFrame_Noiseless_median = median(z_wfsCameraFrame_Noiseless(:));
z_wfsCameraFrame_Noiseless_mean = mean(z_wfsCameraFrame_Noiseless(:));

z_wfsLensletIntensities_median = median(wfs.lensletIntensity);
%% END::  Changing the way how we add the read-out NOISE - see Issue 14



%% Propagating the Asterism (star constellation) for Tipt-Tilt measurements
TTAst = TTAst.*tel*dm*QuadCell;

flush(cam)
cam.clockRate    = 1;
cam.exposureTime = exposureTime;
cam.startDelay   = startDelay;

figure(31416)
imagesc(cam,'parent',subplot(2,1,1))
cam.frameListener.Enabled = true;
subplot(2,1,2)
h = imagesc(catMeanRmPhase(sciStar));
axis xy equal tight
colorbar
 


%% Pseudo-Open Loop gain:
% gain_pol = 0.7; %%% peeudo-open loop controller

%% Closed-Loop gain:
gain_cl = 0.3; 


%% Noisy or Noiseless case:
                if (params.noise_in_wfs == 0) 
                    
                        QuadCell.camera.photonNoise = false;
                        QuadCell.camera.readOutNoise = 0;
                        +QuadCell
                        
                        wfs.camera.photonNoise = false;
                        wfs.camera.readOutNoise = 0;
                        +wfs
                        
                else
                    
                        QuadCell.camera.photonNoise   = params.noise.photonNoise;
% % % %                         QuadCell.camera.readOutNoise = params.noise.wfs_readOutNoise; %% turning this noise to zero may fail the phaseRecosntruction!!!
                       QuadCell.camera.readOutNoise = z_wfsLensletIntensities_median / params.noise.wfs_readOutNoise_SNR; 
                         QuadCell.framePixelThreshold   = params.noise.framePixelThreshold; %% increase to reduce the WF error
                       +QuadCell 

                        wfs.camera.photonNoise   = params.noise.photonNoise;
% % % %                         wfs.camera.readOutNoise = params.noise.wfs_readOutNoise; %% turning this noise to zero may fail the phaseRecosntruction!!!
                       wfs.camera.readOutNoise = z_wfsLensletIntensities_median / params.noise.wfs_readOutNoise_SNR; 
                        wfs.framePixelThreshold   = params.noise.framePixelThreshold; %% increase to reduce the WF error
                        +wfs
                end %%  if (params.noise_in_wfs == 0)
                
fprintf('\n\n Readout noise value %5.5e \n\n', wfs.camera.readOutNoise)


%% Piston-tip-tilt removal from slopes
M_TT = [ones(wfs.nValidLenslet,1)   zeros(wfs.nValidLenslet,1); ...
              zeros(wfs.nValidLenslet,1),  ones(wfs.nValidLenslet,1)];  %% zonal (wfs) matrix for tip-tilt removal (a wfs.nValidLenslet X 2 matrix )

[P_w, Sigma_eta] = tool_ltao_projection4slopes(M_TT, wfs, tel, atm, params); %%% computing the P_w projection matrix, which can also account for noise
SlopeTT_Remove = eye(2*wfs.nValidLenslet) - P_w;%%%% removing tip-tilt:  Slopestt = I - P_w = I - M_TT*M_TT^T



%% Generate the Zernike with only modes 2 and 3 - that is, Tip and Tilt
zern = zernike(2:3,'resolution',nLenslet+1, 'pupil',dmLowRes.validActuator);
zern_modes = zern.modes(dmLowRes.validActuator,:);

%% zernikeStats.angularCovariance [TBD]


F = bifaLowRes.modes(wfs.validActuator,:); 
iF = 0.5*pinv(full(F));

M_TT_dm = iF*zern_modes;
P_a = tool_ltao_projection4dmtt(M_TT_dm, Sigma_eta, calibDm, params); %%% computing the P_w projection matrix, which can also account for noise
DMTT_Remove = eye(dmLowRes.nValidActuator) - P_a;    %%%% removing tip-tilt:  Dmtt = I - P_a = I - M_TT_dm*M_TT_dm^T


TTstar_wavenumber = TTAst(1).wavelength/(2*pi);
% % % % % phase_MMSE = lgsAst_slmmse*wfs.slopes;  %%% Is an MMSE phase reconstruction in micrometers
phase_MMSE_LGS = zeros(dm.nValidActuator,1);

u_ngs = zeros(dmLowRes.nValidActuator,1); % NGS control loop

%% The loop is closed for one full exposure of the science camera.
nIteration = startDelay + exposureTime;
warning off
      


%% Begin Closed Loop:
for k=1:nIteration

    % Objects update
    +tel;
    +lgsAst;
    +TTAst;
    +sciStar;
    
    
% Remove Tip-Tilt from each column of slopes SlopeTT_Remove*wfs.slopes, remove the vector of calibDm.D*dm.coefs 
%% Here   SlopeTT_Remove*wfs.slopes   removes the Tip-Tilt from the WFS slopes, and
% %             calibDm.D*DMTT_Remove*dm.coefs   removes the Tip-Tilt added by QuadCell on the previous step from the previous DM commands    
%     slopes_sans_TT  = bsxfun( @minus,  SlopeTT_Remove*wfs.slopes,    calibDm.D*DMTT_Remove*dm.coefs);
%     phase_MMSE_LGS = lgsAst_slmmse*(slopes_sans_TT); 
% 
%     %% Now add TT from TT stars, and use PI closed-loop controller for TT-loop:
% %     dmLowRes.coefs  = dmLowRes.coefs  - gain_cl * M_TT_dm * mean(QuadCell.slopes,2) * TTstar_wavenumber; 
%     dmLowRes.coefs  = dmLowRes.coefs  - gain_cl * M_TT_dm * mean(QuadCell.slopes,2) * TTAst(1).wavelength/8;  %% this lambda/8 is not just 2pi, but worked out for the tilt.
%     
% 
%     %% Add TT control signals from the TT-loop back to the DM commands
%     dm.coefs = iF*phase_MMSE_LGS  -  dmLowRes.coefs ;    


% 
% %% The code below (straight from the run_cl_Harmony  gives worse results in terms of compensation:
% %% Strehl ratio =  0.2282 compared to the 0.58 of the code above.  %% See the Issue 16 
% 
%     %% From each column of TT-removed slopes SlopeTT_Remove*wfs.slopes, remove the vector of calibDm.D*dm.coefs 
%     phase_MMSE_LGS = lgsAst_slmmse*( bsxfun( @minus,  SlopeTT_Remove*wfs.slopes,    calibDm.D*dm.coefs )); %%% this ____lgsAst_slmmse*(  ____ is an object....  so we need that extra bracket..... 
% 
%     %% Remove Tip-Tilt from the DM commands added on the previous step by the QuadCell
%     dm.coefs = DMTT_Remove*iF*phase_MMSE_LGS;
% 
%     %% Now add TT from TT stars, and use PI closed-loop controller for TT-loop:
%     u_ngs = u_ngs - gain_cl * M_TT_dm * mean(QuadCell.slopes,2) * TTstar_wavenumber; 
% 
%     %% Add TT control signals from the TT-loop back to the DM commands
%     dm.coefs = dm.coefs - u_ngs;    
% % %        

    %% Display
    set(h,'Cdata',catMeanRmPhase(sciStar))
    drawnow

end

cam.strehl


%% 3. \phi_z^{res} = Z^\dagger * \phi_{pix}^{res}
zPhase = reshape(sciStar.meanRmPhase, nRes*nRes,1);
zphi_res_modes = pinv(full(dm.modes.modes))*zPhase;
figure, bar(zphi_res_modes);

figure, imagesc(QuadCell.camera.frame); %%% this is what Quadcell camera sees - it must see the tilt!