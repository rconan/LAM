close all
clc

% matlabpool % - start this if you need parallel computing
% clear
% matlabpool close - if you need to stop the parallel computing
% break
% 
% %%% This is to check the covariance of the \Sigma_z and P_\alpha * \Sigma_{\phi_z} \P_alpha^T
% % zern.D = 8;
% 
Sigma_z = zernikeStats.covariance(zern,atm); %|sigma_Z
fprintf('\n Atm zernike maxRadialDegree=%d, modes=%d\n ', maxRadialDegree, zern.nMode);

Sigma_alpha;
fprintf('\n Projected zernike maxRadialDegree=%d, modes=%d\n ', maxRadialDegreeProj, zernProj.nMode);

diff_Covar = trace(Sigma_z) - trace(Sigma_alpha) %% checking the covariance difference between the initial Zernike modes and the projected (Sigma_\alpha == P_\alpha * Sigma_phi_z P_alpha^T)


diff_Covar_relative = (trace(Sigma_z) - trace(Sigma_alpha))/trace(Sigma_z) %% checking the covariance difference between the initial Zernike modes and the projected (Sigma_\alpha == P_\alpha * Sigma_phi_z P_alpha^T)

%% Uncommet in case of Ammse
% zzz = AmmseLag{1,2};
% zzz(abs(zzz(:))<1e-4) = 0;
% 
% figure, spy(zzz)
% 

%% Illustrations for the matrices of projection
% Projection_alpha_sparse = sparse(Projection_alpha);
% figure, spy(Projection_alpha_sparse)
% 
% Projection_beta_sparse = sparse(Projection_beta);
% figure, spy(Projection_beta_sparse)

slopes_show =  wfs(1).slopes; %%% slopes for the WFS number 1

%% Minimum Variance reconstructor (static MMSE):
% % Rmv = (Sigma_betaalpha*DAst')/(DAst*Sigma_alpha*DAst' + CnAst);

% z2u_times_Rmv = Z2U*Rmv;


%% How To Display turbulence and residual phase
%     figure(101)
%     rad2mic = 1e6/sciObjectVec{1}.waveNumber;
%     h = imagesc([sciObjectVec{1}.meanRmPhase]*rad2mic);
    
%     figure(102)
%     rad2mic = 1e6/imgrObjectVec{1}.waveNumber;
%     h = imagesc([imgrObjectVec{1}.meanRmPhase]*rad2mic);


%% Show the imager's frame (science imager):
figure, mesh(dmStatic.surface); %% microemeters

% figure, mesh(imgrStatic.frame); %% microemeters

% figure, mesh(imgrObjectVec{1,1}.referenceFrame), title('No Atmosphere');
% 
% figure, mesh(imgrStatic.frame), title('After Open-Loop compensation, static MMSE');
% figure, mesh(imgrNoCorrection.frame), title('No compensation');
% figure, mesh(sciStatic.opd);


figure, mesh(imgrObjectVec{1,1}.frame);
% figure, mesh(imgrObjectVec{1,2}.frame);

figure, mesh(imgrObjectVec{1,1}.referenceFrame), title('No Atmosphere');

imgrObjectVec{1}.ee     %% encircled energy
imgrObjectVec{1}.strehl %% strehl ration

%  @(detector:relay)> starting imager integration!
% >> dm
% ___ DEFORMABLE MIRROR ___
%  11X11 actuators deformable mirror: 
%   . 97 controlled actuators
