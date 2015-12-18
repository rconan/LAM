function [P_w, Sigma_eta] = tool_ltao_projection4slopes(M_TT, wfs, tel, atm, params)

      if     (params.noise_account_in_SlopesTT  == 0); 

            fprintf('\n\n No accounting for noise in Slopes tip-tilt removal!\n\n');
            P_w = M_TT*pinv(M_TT);  %%% Projection matrix for Tip-Tilt (TT) removal (works when there is no noise)
            Sigma_eta = 0; %% this is needed for the DMTT function further in the code
            
      else
             %%% Getting the Noise Covariance matrix - With noise
            fprintf('\n\n Using noise covariane matrix for Slopes tip-tilt removal!\n\n');

            +wfs
            ngs = source; %create a simple on-axis guide star object for WFS calibration
                
            tel = tel - atm;
            
%%%%% Use the same noise settings from the above!
            ngs.wavelength = photometry.H;          %%% This is to match the Guidestar wavelength  (H-band)
            
% %                     wfs.camera.readOutNoise = 0.2;  %% turning this noise to zero may fail the phaseRecosntruction!!!
% %                     wfs.camera.photonNoise = true;
% %                     wfs.framePixelThreshold = 0.1; %% increase to reduce the WF error

                    ngs=ngs.*tel*wfs;

           %%% Begin to estimate the noise convariance matrix
            fprintf('\n\n Initialising the noisy WFS....\n')

            slopes = zeros(wfs.nSlope, params.noise_account_nMeasurements);
            for kMeas=1:params.noise_account_nMeasurements
                +wfs
                slopes(:,kMeas) = wfs.slopes;
            end
            Sigma_eta = slopes*slopes'/params.noise_account_nMeasurements; %%% estimation of the covariance matrix of WFS Noise

             tel = tel + atm;
             +wfs
             
            fprintf('... DONE! \n\n ')
           %%% Done estimation of the noise convariance matrix

            Sigma_eta_diag = diag(diag(Sigma_eta));
            inv_Sigma_eta_diag = inv(Sigma_eta_diag);

            M_TT_sigma_pinv = ( inv(M_TT' * inv_Sigma_eta_diag * M_TT) )*(M_TT' * inv_Sigma_eta_diag);
            P_w = M_TT*M_TT_sigma_pinv;            
            
      end