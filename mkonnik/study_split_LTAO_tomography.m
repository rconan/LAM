
z_NdmLowresValidActuators = dmLowRes.nValidActuator;
z_NwfsValidLenslets = wfs.nValidLenslet;

z_wfsSlopes = wfs.slopes;

z_Sigma_eta_diag = Sigma_eta_diag;

z_calibDm_D = calibDm.D;  %%% interaction matrix WFS2DM, 152x97
z_calibDm_M = calibDm.M; %% pseudoinverse of interaction matrix   97x152

z_M_TT_modal = iF*zern_modes;

z_Gm = 2*pinv(zern_modes)*F;


% z_Sigma_eta_modal = eye(z_NdmLowresValidActuators);
% z_inv_Sigma_eta_modal = inv(eye(z_NdmLowresValidActuators));

% z_M_TT_sigma_pinv = ( inv(M_TT_modal' * z_inv_Sigma_eta_modal * M_TT_modal) )*(M_TT_modal' * z_inv_Sigma_eta_modal);


slopes_dm = zeros(z_NdmLowresValidActuators,nMeas);

for kMeas=1:nMeas
    slopes_dm(:,kMeas) = calibDm.M*slopes(:,kMeas)*ngs.waveNumber;
end

% z_Sigma_eta_modal2 = slopes_dm*slopes_dm'/nMeas;
Sigma_eta_modal = slopes_dm*slopes_dm'/nMeas;

Sigma_eta_modal_diag= diag(diag(Sigma_eta_modal));
inv_Sigma_eta_modal_diag = inv(Sigma_eta_modal_diag);

M_TT_sigma_modal_pinv = ( inv(M_TT_modal' * inv_Sigma_eta_modal_diag * M_TT_modal) )*(M_TT_modal' * inv_Sigma_eta_modal_diag);
            

            
break

%% piston-tip-tilt removal 
M_TT = [ones(wfs.nValidLenslet,1)   zeros(wfs.nValidLenslet,1); ...
              zeros(wfs.nValidLenslet,1),  ones(wfs.nValidLenslet,1)];  %% zonal (wfs) matrix for tip-tilt removal (a wfs.nValidLenslet X 2 matrix )

%       if     (params.noise_in_wfs  == 0); 
P_w = M_TT*pinv(M_TT);  %%% Projection matrix for Tip-Tilt (TT) removal (works when there is no noise)
SlopeTT_Remove = eye(2*wfs.nValidLenslet) - P_w;




%%%% Generate the Zernike with only modes 2 and 3 - that is, Tip and Tilt
zern = zernike(2:3,'resolution',nLenslet+1, 'pupil',dmLowRes.validActuator);
zern_modes = zern.modes(dmLowRes.validActuator,:);

F = bifaLowRes.modes(wfs.validActuator,:); %%% WTF? Should we multiply it by lambda/2pi?, 97x97
iF = 0.5*pinv(full(F));

M_TT_modal = iF*zern_modes;

%      if     (params.noise_in_wfs  == 0); 
P_a = M_TT_modal*pinv(M_TT_modal);   % Projection matrix for Tipt-Tilt removal _in the actuator space_
DMTT_Remove = eye(dmLowRes.nValidActuator) - P_a;

            
break
%           
%       else
%              %%% Getting the Noise Covariance matrix - With noise
%             +wfs
%             ngs = source; %create a simple on-axis guide star object for WFS calibration
%                 tel = tel - atm;
% 
% %             ngs.wavelength = photometry.R;          %%% This is to match the Guidestar wavelength
%                     wfs.camera.readOutNoise = 0;
%                     wfs.camera.photonNoise = true;
%                     wfs.framePixelThreshold = 0; %% increase to reduce the WF error
%  
%                     ngs=ngs.*tel*wfs;
% 
%            %%% noise convariance matrix
%             fprintf('\n\n Initialising the noisy WFS....\n')
% 
%             nMeas = 200;
%             slopes = zeros(wfs.nSlope,nMeas);
% 
%             for kMeas=1:nMeas
%                 +wfs
%                 slopes(:,kMeas) = wfs.slopes;
%             end
%             Sigma_eta = slopes*slopes'/nMeas;
% 
%                 tel = tel + atm;
%                 +wfs
%             fprintf('... DONE! \n\n ')
% 
%             Sigma_eta_diag = diag(diag(Sigma_eta));
%             inv_Sigma_eta_diag = inv(Sigma_eta_diag);
% 
%             TT_sigma_inv = ( inv(TT' * inv_Sigma_eta_diag * TT) )*(TT'*inv_Sigma_eta_diag);
%             M_sigma = TT*TT_sigma_inv;
%             
%             SlopeTT_Remove = eye(2*wfs.nValidLenslet) - M_sigma;
%       end

      

            
            
            
break

% %%% Getting the Noise Covariance matrix - With noise
% +wfs
% ngs = source; %create a simple on-axis guide star object for WFS calibration
%     tel = tel - atm;
% ngs.wavelength = photometry.H;          %%% This is to match the Guidestar wavelength
% 
% wfs.camera.readOutNoise = 0.2;
% wfs.camera.photonNoise = true;
% wfs.framePixelThreshold = 0; %% increase to reduce the WF error
% 
% ngs=ngs.*tel*wfs;
% 
% %% noise convariance matrix
% fprintf('\n\n Initialising the noisy WFS....\n')
% 
% nMeas = 200;
% slopes = zeros(wfs.nSlope,nMeas);
% 
% for kMeas=1:nMeas
%     +wfs
%     slopes(:,kMeas) = wfs.slopes;
% end
% Sigma_eta = slopes*slopes'/nMeas;
% 
%     tel = tel + atm;
% 
%     +wfs
% 
%     fprintf('... DONE! \n\n ')
% 
% 
% Sigma_eta_diag = diag(diag(Sigma_eta));
% inv_Sigma_eta_diag = inv(Sigma_eta_diag);
% 
% TT_sigma_inv = ( inv(TT' * inv_Sigma_eta_diag * TT) )*(TT'*inv_Sigma_eta_diag);
% M_sigma = TT*TT_sigma_inv;
% 



% 
% %%%% Generate the Zernike with only modes 2 and 3 - that is, Tip and Tilt
% zern = zernike(2:3,'resolution',nLenslet+1, 'pupil',dmLowRes.validActuator);
% 
% F = bifaLowRes.modes(wfs.validActuator,:); %%% WTF? Should we multiply it by lambda/2pi?
% iF = 0.5*pinv(full(F));
% 
% iF_zern = iF*zern.modes(dmLowRes.validActuator,:);
% 
% 
% %      if     (params.noise_in_wfs  == 0); 
% 
%             M_TT = iF_zern*pinv(iF_zern);
%             
%             DMTT_Remove = eye(dmLowRes.nValidActuator) - M_TT;



% Take existing slopes and decompose them in Zernike.
% Then ???