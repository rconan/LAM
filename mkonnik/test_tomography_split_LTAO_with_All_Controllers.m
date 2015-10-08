clc
close all
clear all
tool_inside_path_adder('matlab_code/mikhail', 'matlab_code');

params = struct; %%% creating the Structure with additional parameters for the computeARmodel
    params.OOMAO_hot_start = 0;  %%% OOMAO is a catastrophic mess that cannot be efficienty hot-started because many objects modify other objects, which turns out to be interdependent mess. It is not just zeroing the camera....
 
samplingTime     = 0.002;% seconds
frameTime          = 0.002;%[0.01,0.015,0.020,0.025]; seconds

% samplingTime     = 0.01;% seconds
% frameTime          = 0.01;%[0.01,0.015,0.020,0.025]; seconds

exposureTime     = 500; %% [ms], equals to the Number of Iterations
startDelay           = 100; %% skip this number of frames

% exposureTime     = 500; %% [ms], equals to the Number of Iterations
% startDelay           = 1; %% skip this number of frames

params.asterism_variability = 'TT_angle_variable';
min_TTangle = 0;
max_TTangle = 80;
step_TTangle = 10;

% params.asterism_variability = 'TT_mag_variable';
% min_magnitude = 18;
% max_magnitude = 10;
% step_magnitude = -1;

% params.asterism_variability = 'TT_mag_variable';
% min_magnitude = 16;
% max_magnitude = 16;
% step_magnitude = -1;


for TTstarAngle=min_TTangle:step_TTangle:max_TTangle
% for TTstarmagnitude=min_magnitude:step_magnitude:max_magnitude

params.TTstars= 3;
    params.TTstars_angle = TTstarAngle; %% asterism angle in arcseconds
%     params.TTstars_angle = 20; %% asterism angle in arcseconds
    
    params.TTstar_magnitude = 10;
% params.TTstar_magnitude = TTstarmagnitude;

    
params.LGSstars = 3;
    params.LGSstars_angle = 10;  %% asterism angle in arcseconds
    params.LGSguidestar_magnitude = 10;
    

% params.controller_type_low_order = 'LQG';
% params.controller_type_low_order = 'GLAO';
% params.controller_type_low_order = 'VirtualDM';
params.controller_type_low_order = 'SAMMSE';

    
params.noise_in_wfs  = 1; %% adding the noise into WFS and OIwfs (same noise everywhere)
% params.noise_in_wfs  = 0; %% noiseless WFS and OIwfs 

params.account4noise_in_TTProjections = 1; 
% params.account4noise_in_TTProjections = 0; 

    params.show_figures = 0; %% display the diagnostic plots during the simulation
    params.dm_crosscoupling_coef = 0.40; %% cross-coupling degree for the deformable mirror.

params.noise_account_nMeasurements =   500;%% number of measurements for the covariance matrix estimation.

%% No Noise  (noiseless) WFS case:
if  (params.noise_in_wfs  == 0); %% display the diagnostic plots during the simulation.
    params.noise.wfs_readOutNoise = 0;
    params.noise.photonNoise = false;
    params.noise.framePixelThreshold = 0;
else
    params.noise.photonNoise = true;
    params.noise.wfs_readOutNoise = 2;
    params.noise.framePixelThreshold = 0;    
end


if (params.account4noise_in_TTProjections == 0)
        params.noise_account_in_SlopesTT  = 0; 
        params.noise_account_in_DMTT  = 0;
else 
        params.noise_account_in_SlopesTT  = 1;     
        params.noise_account_in_DMTT  = 1;
end


%% Gudestars: LGS and TipTilt
lgsAst = source('asterism',{[params.LGSstars,arcsec(params.LGSstars_angle),0]},...
    'height',90e3, 'magnitude', params.LGSguidestar_magnitude);

TTAst = source('asterism',{[params.TTstars,   arcsec(params.TTstars_angle),0]},...
    'magnitude', params.TTstar_magnitude, 'wavelength',photometry.H); %% Tip-tilt asterism (for TT measurements)

nTTstars = size(TTAst,2);


%% Atmospheric parameters
wavelength      = photometry.R;
altitudes       = [0,    5]*1e3;
fractionalR0    = [0.7, 0.3];
windSpeed       = [2,  5];
windDirection   = [pi, pi];

L0 = 30; %% outer scale, [m]
r0 = 0.20; %% Fried parameter, [m]


%% Creating the Atmosphere object
atm = atmosphere(photometry.V, r0, L0,...
    'fractionnalR0',fractionalR0, ...
    'altitude',     altitudes,...
    'windSpeed',    windSpeed,...
    'windDirection',windDirection);
nLayers = atm.nLayer;

atmNeg = atmosphere(photometry.V, r0, L0,...
    'fractionnalR0',fractionalR0, ...
    'altitude',     altitudes,...
    'windSpeed',    -windSpeed,...   %%% mind the negative windspeed - this is done for the SA algorithms!
    'windDirection',windDirection);


%% Wavefront Sensor
nLenslet   = 10;
nPx  = 12;
nRes = nLenslet*nPx;
minLightRatio = 0.65; %% large values (more than 0.8  can lead to poor reconstruction - you are truncating the phase!).


%% TELESCOPE
D    = 8;
d    = D/nLenslet; % lenslet pitch
tel = telescope(D,...
    'resolution',nRes,...
    'fieldOfViewInArcsec',180,...
    'samplingTime', samplingTime);

%% OIwfs calibration
OIwfs = tool_calibrate_wfs_oomao(1,8,minLightRatio, tel, D, 'H');  %% since TTast is in H-band, I'm going to calibrate OIwfs in H-band, too
% OIwfs = tool_calibrate_wfs_oomao(1,8,minLightRatio, tel, D, 'K');  %% I don't know why do we calibrate OIwfs in K band.....

%% SH WFS calibration
wfs         = tool_calibrate_wfs_oomao(nLenslet,nPx,minLightRatio , tel, d, 'R');


%% Deformable mirror
nActuator = nLenslet + 1;
bif = influenceFunction('monotonic',params.dm_crosscoupling_coef);

dm = deformableMirror(nActuator,...
    'modes',bif,...
    'resolution',nRes,...
    'validActuator',wfs.validActuator );

bifaLowRes = influenceFunction('monotonic',params.dm_crosscoupling_coef);
dmLowRes = deformableMirror(nLenslet+1,...  %%% Why do we need a DM with surface of 11x11?? Is this for the TT-mirror?
    'modes',bifaLowRes,...
    'resolution',nLenslet+1,... 
    'validActuator',wfs.validActuator);


%% DM Calibration 
ngs = source;
ngs.magnitude = 0;
ngs = ngs.*tel;
        calibDm = calibration(dm,wfs,ngs,ngs.wavelength,nLenslet+1,'cond',1e1);
tel = tel + atm;


 %% Linear Minimum Mean Square Error (LMMSE) wavefront reconstructor
lgsAst_slmmse = slopesLinearMMSE(wfs,tel,atm,lgsAst,'mmseStar',ngs,'NF',512); 
lgsAst_slmmse.wavefrontSize = [dm.nValidActuator,1];
lgsAst_slmmse.warmStart = true;


%% Science targets
sciStar = source('wavelength',photometry.H);

%% Science camera initialisation
cam = imager(tel);

tel = tel - atm;
sciStar = sciStar.*tel*cam;
    cam.referenceFrame = cam.frame;
tel = tel + atm;
+sciStar;


%% Propagating the Asterism (star constellation) for Tipt-Tilt measurements
flush(cam)
cam.frame = cam.frame*0;
cam.clockRate    = 1;
cam.exposureTime = exposureTime;
cam.startDelay   = startDelay;

        if (params.show_figures ==1)
            figure(31416)
            imagesc(cam,'parent',subplot(2,1,1))
            cam.frameListener.Enabled = true;
            subplot(2,1,2)
            h = imagesc(catMeanRmPhase(sciStar));
            axis xy equal tight
            colorbar
        end
        

%% Generate the Zernike with only modes 2 and 3 - that is, Tip and Tilt
zern = zernike(2:3,'resolution',nLenslet+1, 'pupil',dmLowRes.validActuator);
zern_modes = zern.modes(dmLowRes.validActuator,:);

[OIwfs, wfs] = tool_enforce_noise_in_wfs(OIwfs,wfs,params);%% Noisy or Noiseless case:


if  (params.OOMAO_hot_start == 0) %% HOT START!!!
    %% Piston-tip-tilt removal from slopes
    M_TT = [ones(wfs.nValidLenslet,1)   zeros(wfs.nValidLenslet,1); ...
                  zeros(wfs.nValidLenslet,1),  ones(wfs.nValidLenslet,1)];  %% zonal (wfs) matrix for tip-tilt removal (a wfs.nValidLenslet X 2 matrix )

    [P_w, Sigma_eta] = tool_ltao_projection4slopes(M_TT, wfs, tel, atm, params); %%% computing the P_w projection matrix, which can also account for noise
    SlopeTT_Remove = eye(2*wfs.nValidLenslet) - P_w;%%%% removing tip-tilt:  Slopestt = I - P_w = I - M_TT*M_TT^T


    %% Piston-tip-tilt removal from DM commands
    F = bifaLowRes.modes(wfs.validActuator,:); %% DM influence matrix
    iF = 0.5*pinv(full(F));

    M_TT_dm = iF*zern_modes;
    P_a = tool_ltao_projection4dmtt(M_TT_dm, Sigma_eta, calibDm, params); %%% computing the P_w projection matrix, which can also account for noise
    DMTT_Remove = eye(dmLowRes.nValidActuator) - P_a;    %%%% removing tip-tilt:  Dmtt = I - P_a = I - M_TT_dm*M_TT_dm^T
end  %% if  (params.OOMAO_hot_start == 0) %% HOT START!!!


%% Low-Order Loop matrices:
if  (params.OOMAO_hot_start == 0)
    Sigma_eta_OIwfs = tool_compute_wfs_noisecov(tel, atm, TTAst, OIwfs, params); %% Calculating the Noise covariance matrix for the Tip-Tilt measurements and OIWFS
    Sigma_eta_OIwfs_stacked = tool_make_blockdiag_matrix ('Sigma_eta_OIwfs', params.TTstars , Sigma_eta_OIwfs);           
end


%% Laser Tomography Adaptive Optics
sciStar = sciStar.*tel*dm*cam;
lgsAst = lgsAst.*tel*dm*wfs;
TTAst = TTAst.*tel*dm*OIwfs;

+lgsAst


%% The loop is closed for one full exposure of the science camera.
phase_MMSE_LGS = zeros(dm.nValidActuator,1);
nIteration = startDelay + exposureTime;

switch params.controller_type_low_order
    
    case 'VirtualDM'
    %%% Create extended Zernike object for Tomographic projections
    zernProj = zernike(2:5,'resolution',nLenslet+1, 'pupil',dmLowRes.validActuator); %% to calculate H_alpha and H_beta

    nLayers = 2; %% this is for the VirtualDM, where we project for 2 atmospheric layer - 0 (ground) and some other altitude (say, 11 km)
    nGs = size(TTAst,2);
    nSciobj = size(sciStar,2);

        G = eye(2)* 4/pi;
        G_large  = tool_make_blockdiag_matrix ('G', params.TTstars, G); %% % G_large = blkdiag(G,G,G);

    altitudes = [0, 5000];
    H_alpha_Cell = tool_analytical_small_footprint_expansion_with_piston(zernProj, altitudes, TTAst, tel, nLayers, nGs);
    H_alpha_Cell_TTonly = H_alpha_Cell; %% selecting only Rows 2:3 for tip-tilt

    for ii=1:nGs
        for jj=1:nLayers            
            H_alpha_Cell_TTonly{ii,jj} = H_alpha_Cell{ii,jj}(2:3,:);
        end
    end %% for ii
    
    H_alpha = cell2mat(H_alpha_Cell_TTonly);

    H_beta_Cell = tool_analytical_small_footprint_expansion_with_piston(zernProj, altitudes, sciStar, tel, nLayers, nSciobj);
    H_beta_Cell_TTonly = H_beta_Cell; %% selecting only Rows 2:3 for tip-tilt
    
    for ii=1:nLayers
            H_beta_Cell_TTonly{1,ii} = H_beta_Cell{1,ii}(2:3,:);
    end %% for ii
    H_beta = cell2mat(H_beta_Cell_TTonly);
        
    Sigma_eta_OIwfs_stacked = Sigma_eta_OIwfs_stacked + eye(2*nGs)*1e-8; %% regularization

    ReconstructorHVirtDM  =   H_beta * pinv(H_alpha' * G_large' * inv(Sigma_eta_OIwfs_stacked) * G_large * H_alpha) *...
    H_alpha' * G_large' * inv(Sigma_eta_OIwfs_stacked);


    case 'LQG'
%% Computing the cross-correlation matrix of Zernike coefficients from the modes of the atmosphere object for EACH ATM LAYER
        if (params.TTstars == 1)    %% This is for the 1 TT star
            Sigma_beta  = zernikeStats.angularCovariance(zern,atm, [TTAst, sciStar]);
            Sigma_alphaalpha = zernikeStats.angularCovariance(zern,atm, [TTAst, TTAst]);
        else
            Sigma_alphaalpha1 = zernikeStats.angularCovariance(zern, atm, [TTAst, sciStar]);
            Sigma_alphaalpha = cell2mat(Sigma_alphaalpha1(1:nTTstars, 1:nTTstars));
            Sigma_beta = cell2mat(Sigma_alphaalpha1(end, 1:nTTstars));
        end
        
TT_beta = Sigma_beta/Sigma_alphaalpha; %% Tip-Tilt reconstructor


%% Adding the estimation of LQG state space parameters:
CovarianceLag = zernikeStats.tiltsAngularCovariance(zern, atm, TTAst, 'lag', samplingTime);
% CovarianceLag = zernikeStats.tiltsAngularCovariance(zern, atmNeg, TTAst, 'lag', samplingTime); %%% use the atmNeg with -180 degr.
AmmseMatrix = CovarianceLag/Sigma_alphaalpha;

    Sigma_epsilon = Sigma_alphaalpha - AmmseMatrix*Sigma_alphaalpha*AmmseMatrix';  %% excitation noise covariance
    Sigma_epsilon = Sigma_epsilon + 1e-7*eye(size(Sigma_epsilon));  %% prevent numerical errors - use crude diagonal loading
    

%% Solve DARE for Tip-Tilt control:
Q = 100*(Sigma_epsilon+Sigma_epsilon')./2; %% we can actually boost the Q matrix to try and improve the performance.
R = Sigma_eta_OIwfs_stacked;

        G = eye(2)* 4/pi;
        G_large  = tool_make_blockdiag_matrix ('G', params.TTstars, G); %% % G_large = blkdiag(G,G,G);

Sigma_infty_Ammse = dare (AmmseMatrix, G_large', Q, R);

M_infty_Ammse = (Sigma_infty_Ammse*G_large')/(G_large*Sigma_infty_Ammse*G_large' + R);

C = G_large;
D = pinv(M_TT_dm) * 8/TTAst(1).wavelength;  %%we don't need TT_beta here

        
        
    case 'SAMMSE'
        G = eye(2)* 4/pi;
        G_large  = tool_make_blockdiag_matrix ('G', params.TTstars, G); %% % G_large = blkdiag(G,G,G);

        if (params.TTstars == 1)    %% This is for the 1 TT star
            Sigma_betaalpha   = zernikeStats.angularCovariance(zern,atm, [TTAst, sciStar]);
            Sigma_alphaalpha = zernikeStats.angularCovariance(zern,atm, [TTAst, TTAst]);
        else
            Sigma_alphaalpha1 = zernikeStats.angularCovariance(zern, atm, [TTAst, sciStar]);
            Sigma_alphaalpha   = cell2mat(Sigma_alphaalpha1(1:params.TTstars, 1:params.TTstars));
            Sigma_betaalpha     = cell2mat(Sigma_alphaalpha1(end, 1:params.TTstars));
        end

        ReconstructorMMSE =  (Sigma_betaalpha*G_large')/(G_large*Sigma_alphaalpha*G_large' + Sigma_eta_OIwfs_stacked);

end  %% for switch cases of controllers


%% Controller gains
gain_pol = 0.8;
gain_TT = 0.5;


%% Control Loop 
u_k_HO  = zeros(dmLowRes.nValidActuator,1); % High-order (LGS) control loop
u_k_LO =  zeros(dmLowRes.nValidActuator,1); %  TT (NGS) control loop

eval_command_slopes_stack = tool_make_stacked_vector ('OIwfs.slopes', 'Slopes_Stack', params.TTstars); %% Make the comand for evaluation of the SlopesStack for arbitrary number of TT guidestars

eval_command_tmp = '';    matrix_name = strcat('D * u_k_LO(: ');
for ii = 1:nTTstars
    eval_command_tmp = strcat(eval_command_tmp , matrix_name, ');' );
end
eval_command_D_times_u_k_LO = strcat(  'D_times_u_k_LO', ' = [', eval_command_tmp, '];')



%% Auto-generate the names for the masive simulations
if (params.OOMAO_hot_start == 0)
    
    switch params.asterism_variability
        case 'TT_mag_variable'
            
%             name_to_write =  strcat('SplitTomoTest_QmatrixX100_',params.controller_type_low_order, '_',num2str(params.TTstars),'TTstars_mag_variable',...
%     '_',num2str(params.LGSstars),'LGS_mag',num2str(params.LGSguidestar_magnitude),'Sampling_',num2str(1/samplingTime),...
%     'Hz_Exposure',num2str(exposureTime),'_StartDelay',num2str(startDelay),'_PhShotNoise',num2str(params.noise.photonNoise),'_ReadNoise',num2str(params.noise.wfs_readOutNoise),  '.data'); %% the name must contain the vaiable parameters of the simulation.

            name_to_write =  strcat('SplitTomoTest_',params.controller_type_low_order, '_',num2str(params.TTstars),'TTstars_mag_variable',...
    '_',num2str(params.LGSstars),'LGS_mag',num2str(params.LGSguidestar_magnitude),'Sampling_',num2str(1/samplingTime),...
    'Hz_Exposure',num2str(exposureTime),'_StartDelay',num2str(startDelay),'_PhShotNoise',num2str(params.noise.photonNoise),'_ReadNoise',num2str(params.noise.wfs_readOutNoise),  '.data'); %% the name must contain the vaiable parameters of the simulation.

        otherwise 
            
            name_to_write =  strcat('SplitTomoTest_',params.controller_type_low_order, '_',num2str(params.TTstars),'TTstars_mag_',num2str(params.TTstar_magnitude),...
    '_',num2str(params.LGSstars),'LGS_mag',num2str(params.LGSguidestar_magnitude),'Sampling_',num2str(1/samplingTime),...
    'Hz_Exposure',num2str(exposureTime),'_StartDelay',num2str(startDelay),'_PhShotNoise',num2str(params.noise.photonNoise),'_ReadNoise',num2str(params.noise.wfs_readOutNoise),  '.data'); %% the name must contain the vaiable parameters of the simulation.

    end %% for switch 
fid = fopen(name_to_write, 'w');

fprintf(fid,'# %s \n', 'This file is for study of the Strehl ratio versus statistics of algorithm convergency');
fprintf(fid,'# Strehl \t nTT stars \t nLGS stars \t TTmagnitude \t LGSmagnitude \t TTstars_angle [arcsec] \t LGSstars_angle [arcsec] \t gain_pol \t gain_TT \t samplingTime \t  exposureTime \t startDelay \n');

fclose('all');  %% close all the open data files.
end

fid = fopen(name_to_write, 'a'); %% open the file with the simulation results for Addition


accumulate_slopes = [];

%% ##### BEGIN::  Closed Loop
for k=1:nIteration

%     disp(k)

    % Objects update
    +tel;
    +lgsAst;
    +TTAst;
    +sciStar;
    
% fprintf('\n\n Iteration %d \n\n', k);

%%% BEGIN:: High-order Loop
    slopes_sans_TT  = bsxfun( @minus,  SlopeTT_Remove*wfs.slopes,    calibDm.D*DMTT_Remove*u_k_HO);
    phase_MMSE_LGS = lgsAst_slmmse*(slopes_sans_TT);     
    u_k_HO = (1-gain_pol)*u_k_HO + gain_pol*iF*phase_MMSE_LGS;
%%% END::: High-order Loop
    
            
%%% BEGIN:: Low-order (TT) Loop            
    eval(eval_command_slopes_stack); %% evaluate the command % Slopes_Stack = [ OIwfs.slopes(:,1);  OIwfs.slopes(:,2);  OIwfs.slopes(:,3) ];    as a string

    accumulate_slopes = [accumulate_slopes, Slopes_Stack(:)];
switch params.controller_type_low_order
    
    case 'LQG'
        
            if k < 2;
                    state_k = M_infty_Ammse*Slopes_Stack(:);
            else
                    state_k_minus_1 = state_k;
                    eval(eval_command_D_times_u_k_LO); %% evaluates the command D_times_u_k_LO = [D * u_k_LO(:);D * u_k_LO(:);D * u_k_LO(:);];
                    state_k = state_k_minus_1 + M_infty_Ammse*(Slopes_Stack(:) - C * state_k_minus_1 +   D_times_u_k_LO);
            end

            state_k_plus_1 = AmmseMatrix*state_k;
            u_k_LO =  M_TT_dm * TT_beta * state_k_plus_1 * TTAst(1).wavelength/(2*pi); %% converting from radians to METERS
   
            dm.coefs = u_k_HO + u_k_LO; %%% here is PLUS because of LQG


    case 'GLAO'        
                G = eye(2)* 4/pi;
                
                Slopes_Average = zeros(2,1);
                
                for kTTStar = 1:nTTstars
                    Slopes_Average = Slopes_Average + Slopes_Stack( (2*kTTStar-1) :(2*kTTStar));
                end
                
                u_k_LO = u_k_LO - gain_TT * M_TT_dm *  TTAst(1).wavelength/(2*pi)*pinv(G)*( Slopes_Average   )./ nTTstars;


    dm.coefs = u_k_HO - u_k_LO;

    
    case 'VirtualDM'
    phase_TT_reconstr = ReconstructorHVirtDM*Slopes_Stack;
    phase_TT       =  M_TT_dm *phase_TT_reconstr*TTAst(1).wavelength/(2*pi);

    u_k_LO = u_k_LO - gain_TT*phase_TT;   

    dm.coefs = u_k_HO -  u_k_LO;    

    
    
    case 'SAMMSE'
    phase_TT_reconstr = ReconstructorMMSE*Slopes_Stack;
    phase_TT       =  M_TT_dm *phase_TT_reconstr*TTAst(1).wavelength/(2*pi);

    u_k_LO = u_k_LO - gain_TT*phase_TT;   

    dm.coefs = u_k_HO -  u_k_LO;    
        
end  %% switch case of Controller types


                if (params.show_figures ==1)
                            set(h,'Cdata',catMeanRmPhase(sciStar))
                            drawnow
                end
end

%% ##### END::  Closed Loop

fprintf(fid,'%1.4f\t  %d\t  %d\t  %1.2f\t  %1.2f\t  %1.2f\t  %1.2f\t  %1.2f\t  %1.2f\t  %1.5f\t  %1.5f\t %d\t \n', ...
                 cam.strehl, params.TTstars, params.LGSstars, params.TTstar_magnitude, params.LGSguidestar_magnitude, params.TTstars_angle, ...    
                 params.LGSstars_angle, gain_pol, gain_TT,  samplingTime , exposureTime, startDelay );

fclose('all');  %% close all the open data files.

params.OOMAO_hot_start = k;

end %%% for TTstarAngle=min_TTangle:step_TTangle:max_TTangle

% cam.strehl %%% gives the Strehl ratio measurement