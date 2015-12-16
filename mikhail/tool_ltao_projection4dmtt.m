function P_a = tool_ltao_projection4dmtt (M_TT_dm, Sigma_eta, calibDm, params)

      if     (params.noise_account_in_DMTT  == 0); 

            fprintf('\n\n No accounting for noise in DM_tip-tilt removal!\n\n');
            P_a = M_TT_dm*pinv(M_TT_dm);   % Projection matrix for Tipt-Tilt removal _in the actuator space_

      else %%% if we want to account for the noise 

            fprintf('\n\n Using noise covariane matrix for Slopes tip-tilt removal!\n\n');
            
% % % %             slopes_dm = zeros(dmLowRes.nValidActuator, params.noise_account_nMeasurements);
% % % %             for kMeas=1:params.noise_account_nMeasurements
% % % %                 slopes_dm(:,kMeas) = calibDm.M*slopes(:,kMeas)*ngs.waveNumber;
% % % %             end
% % % %             Sigma_eta_dm = slopes_dm*slopes_dm'/params.noise_account_nMeasurements;

            Sigma_eta_dm =   calibDm.M * Sigma_eta * (calibDm.M') ;  %%%%    * (ngs.waveNumber)^2;
            
            Sigma_eta_dm_diag= diag( diag( Sigma_eta_dm ) );
            inv_Sigma_eta_dm_diag = inv(Sigma_eta_dm_diag);

            M_TT_sigma_dm_pinv = ( inv(M_TT_dm' * inv_Sigma_eta_dm_diag * M_TT_dm) )*(M_TT_dm' * inv_Sigma_eta_dm_diag);
            P_a = M_TT_dm*M_TT_sigma_dm_pinv ;   % Projection matrix for Tipt-Tilt removal _in the actuator space_
                 
      end