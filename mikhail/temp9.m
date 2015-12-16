% % TT = mean(QuadCell.slopes,2) * TTstar_wavenumber; 
% % 
% % TT_slopes = QuadCell.slopes;
% % 
% % % M_TT_dm 
% % 
% % TT_phase = TTAst_slmmse*(QuadCell.slopes);
% % TT_phase = M_TT_dm *TT_phase 
% % 
% % % figure, mesh(QuadCell.slopes);


Sigma_betaalpha = zernikeStats.angularCovariance(zern,atm, [TTAst, sciStar]); % Angular Cov between beta (guidestar) and alpha


break
%% We should use the zernikestat.AngularCovar
Sigma_betaalpha = zernikeStats.angularCovariance(zern,atm, [TTAst, sciStar]); %|sigma_Z
Sigma_alphaalpha= zernikeStats.angularCovariance(zern,atm, [TTAst, TTAst]); %|sigma_Z

% angularCovariance(zern,atm,src,optSrc)

gamma_T = 2;                           %% the average slope produced by TipTilt
gamma_F = 16*sqrt(3)/(3*pi); %% the average slope produced by Focus mode
gamma_A = 8*sqrt(6)/(3*pi);   %% the average slope produced by Astigmatism

G_TT1 = [gamma_T, 0; 0, gamma_T];
G_TT2 = [gamma_T, 0; 0, gamma_T];

G = blkdiag(G_TT1);  %%% this is modal matrix that translates modal coeffs from TT, TA, and TTFA modes into average slopes over the illuminatedsubregion of each subaperture.




           %%% Begin to estimate the noise convariance matrix
            fprintf('\n\n Initialising the noisy WFS....\n')
             tel = tel - atm;

            QuadCellslopes = zeros(QuadCell.nSlope, params.noise_account_nMeasurements);
            for kMeas=1:params.noise_account_nMeasurements
                +QuadCell
                QuadCellslopes(:,kMeas) = QuadCell.slopes;
            end
            Sigma_eta_quadcell = QuadCellslopes*QuadCellslopes'/params.noise_account_nMeasurements; %%% estimation of the covariance matrix of WFS Noise

             tel = tel + atm;
             
            fprintf('... DONE! \n\n ')
           %%% Done estimation of the noise convariance matrix
           
ReconstructorMMSE =  (Sigma_betaalpha*G')/(G*Sigma_alphaalpha*G' + Sigma_eta_quadcell);


phase_TT =  ReconstructorMMSE*QuadCell.slopes
