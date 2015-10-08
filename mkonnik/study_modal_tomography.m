close all

Z2U_times_Rmv =  Z2U*Rmv;
ZCalibDM_M = calibDm.M;

norm(ZCalibDM_M - Z2U_times_Rmv)/norm(ZCalibDM_M)

% break;

% figure, imagesc(cam.referenceFrame)
% 
% figure, imagesc(cam.frame)
% figure, mesh(cam.frame)


% %%% This is to check the covariance of the \Sigma_z and P_\alpha * \Sigma_{\phi_z} \P_alpha^T
Sigma_z = zernikeStats.covariance(zern,atm); %|sigma_Z
fprintf('\n Atm zernike maxRadialDegree=%d, modes=%d\n ', maxRadialDegree, zern.nMode);

Sigma_alpha;
fprintf('\n Projected zernike maxRadialDegree=%d, modes=%d\n ', maxRadialDegreeProj, zernProj.nMode);

diff_Covar = trace(Sigma_z) - trace(Sigma_alpha) %% checking the covariance difference between the initial Zernike modes and the projected (Sigma_\alpha == P_\alpha * Sigma_phi_z P_alpha^T)

diff_Covar_relative = (trace(Sigma_z) - trace(Sigma_alpha))/trace(Sigma_z) %% checking the covariance difference between the initial Zernike modes and the projected (Sigma_\alpha == P_\alpha * Sigma_phi_z P_alpha^T)

% diff_norm_SigmaZ_Sigmaalpha = norm(Sigma_z  - Sigma_alpha) %% checking the covariance difference between the initial Zernike modes and the projected (Sigma_\alpha == P_\alpha * Sigma_phi_z P_alpha^T)


%% Illustrations for the matrices of projection
% Projection_alpha_sparse = sparse(Projection_alpha);
% figure, spy(Projection_alpha_sparse)
% 
% Projection_beta_sparse = sparse(Projection_beta);
% figure, spy(Projection_beta_sparse)

% slopes_show =  wfs(1).slopes; %%% slopes for the WFS number 1


%% this works
ZCalibDM_M = calibDm.M;

    if nGs ==1
    %%% This is for one GuideStar    
    StackSlopeLarge(1:nSlope,1) = slopesStack(:,1);
    end


if nGs ==3
 % This is for 3 NGSes
        StackSlopeLarge(1:nSlope,1) = slopesStack(:,1);
        StackSlopeLarge(nSlope+1:nSlope*2,1) = slopesStack(:,3);
        StackSlopeLarge(2*nSlope+1:nSlope*3,1) = slopesStack(:,2);    
end

U_works = calibDm.M*StackSlopeLarge(:);
[min(U_works)  mean(U_works) max(U_works)]




%% this (LS reconstruction) does not
Zphi_LS_modal = pinv(DAst*Projection_alpha(idx,idy)) * StackSlopeLarge(:);
Zphi_LS_modal_pbeta    = Projection_beta(idp,idy) * Zphi_LS_modal ;

% % % %     Z2U = calibDm.M*Dzcc;% Using zernike derivatives and slopes to DM command matrix
U_LS = calibDm.M * Dzcc * Projection_beta(idp,idy) * Zphi_LS_modal ;
[min(U_LS)  mean(U_LS) max(U_LS)]


%% this (Tomo) does not work as well
Zphi_tomo_modal =  ( (Sigma_phi*Projection_alpha(idx,idy)'*DAst'  ) / (DAst*Sigma_alpha*DAst' + CnAst) ) * StackSlopeLarge(:);
Zphi_tomo_modal_pbeta =  Projection_beta(idp,idy) * Zphi_tomo_modal;

U_Tomo = calibDm.M * Dzcc * Projection_beta(idp,idy) * Zphi_tomo_modal;
[min(U_Tomo)  mean(U_Tomo) max(U_Tomo)]


%% Comparison of the DM phase and Zernke phase:
fprintf('\n Atm zernike maxRadialDegree=%d, modes=%d\n ', maxRadialDegree, zern.nMode);
fprintf('\n Projected zernike maxRadialDegree=%d, modes=%d\n ', maxRadialDegreeProj, zernProj.nMode);

F = dm.modes.modes; %%% modal matrix

Zphi_pixels_tomo = F*U_Tomo;
Zphi_pixels_LS      = F*U_LS;


Zphi_fromZernike = zern.p * Zphi_LS_modal_pbeta;

figure, mesh(  reshape(Zphi_pixels_tomo, nRes,nRes)   ), title('Phase from Tomography reconstruction')
figure, mesh(  reshape(Zphi_pixels_LS, nRes,nRes)   ), title('Phase from Least-Square reconstruction')

figure, mesh(  reshape(Zphi_fromZernike, nRes,nRes)   ), title('Phase from Z_p * phi_beta^z')


%% Telescope's phase
tel = tel - atm;
ngs=source;
ngs = ngs.*tel*zern*wfs;

size_row      = size(ngs.phase,1);
size_col        = size(ngs.phase,2);
size_modes  = size(ngs.phase,3);
Zphi_pixels    = reshape(ngs.phase,  size_row*size_col,  size_modes);

U_dm_mod = pinv( full(F)  ) * Zphi_pixels;






% break 
% Rmv = (Sigma_betaalpha*(DAst'))/(DAst*Sigma_alpha*DAst' + CnAst);
% 
%     Z2U = calibDm.M*Dzcc;% Using zernike derivatives and slopes to DM command matrix
%     Z_slopes = StackSlopeLarge(:);
%     U_doesnotwork = Z2U*Rmv*StackSlopeLarge(:);     %%% The problem appears to be here: the slopes are a lot larger than the coefficientsVectors
% 
% [min(U_doesnotwork)  max(U_doesnotwork)]

% Zphi_beta = Rmv*StackSlopeLarge(:); %%% seems like the phase in the Beta direction - in ZERNIKE, modal!
% Zphi_beta_modal = DAst*Zphi_beta; %%% these modes are very negative - worng account for piston?!?!?!?

% figure, imagesc(Zphi_beta_modal)


% U2 = calibDm.M*Dzcc*Rmv*StackSlopeLarge(:); %%% seems like the phase in the Beta direction - in ZERNIKE, modal!
% 
% figure, imagesc(Dzcc)

% 
% 
% break
% CnAst = blkdiag( Cn , Cn , Cn );
% DAst = blkdiag( Dzcc , Dzcc , Dzcc  );
% 
% Rmv = (Sigma_betaalpha*DAst')/(DAst*Sigma_alpha*DAst' + CnAst);
% 
% Z2U = ZCalibDM_M*Dzcc;% Using zernike derivatives and slopes to DM command matrix
% 
% coefsVec(:,1,kIteration) = Z2U*Rmv*k1(:);




close all
zDz = pinv(DAst);
figure, imagesc(zDz);

zRmv_zero = (Sigma_betaalpha*DAst')/(DAst*Sigma_alpha*DAst' + CnAst);
figure, imagesc(zRmv_zero);


break
% (0.3*(tel.D/(10 *r0) )^(5/3) )*(0.5/1.65)^2
% exp(-0.3)

close all
zernTilt = zernike(2:2,'resolution',nRes,'pupil',tel.pupil,'D',tel.D);
zernTilt.c = ngs.wavelength/(2*pi);

tel=tel-atm;
ngs = ngs.*tel*zernTilt*wfs;
figure, plot(wfs.slopes)

F = dm.modes.modes;
Fz =  0.5*pinv(full(F))*zern.modes*ngs.wavelength/(2*pi);   


dm.coefs =  Fz*pinv(Dzcc)*wfs.slopes;   

sciStars = sciStars.*tel*zernTilt*dm;
phase_3 = sciStars.meanRmPhase;

figure, imagesc(phase_3), colorbar;

break

dm.coefs =  -Z2U*pinv(Dzcc)*wfs.slopes;   

sciStars = sciStars.*tel*zernTilt*dm;
phase_1 = sciStars.meanRmPhase;
figure, imagesc(phase_1); colorbar

dm.coefs =  -calibDm.M*wfs.slopes;   
sciStars = sciStars.*tel*zernTilt*dm;
phase_2 = sciStars.meanRmPhase;


figure, imagesc(phase_2);



break

Sigma_z = zernikeStats.covariance(zernProj,atm); %|sigma_Z
fprintf('\n Atm zernike maxRadialDegree=%d, modes=%d\n ', maxRadialDegree, zern.nMode);

fprintf('\n Projected zernike maxRadialDegree=%d, modes=%d\n ', maxRadialDegreeProj, zernProj.nMode);

% Rmmse = Sigma_z*Dzcc';

% Sangcov_Cell = phaseStats.zernikeAngularCovariance(zernProj,atm,astNGS); %%% Computing the Angular Covariance
% Sangcov_matrix = cell2mat(Sangcov_Cell);



% %% wavefront reconstruction minimum variance
% M = Sangcov_matrix*DzAst'/(DzAst*Sangcov_matrix*DzAst'+CznAst);
% % z = z - 1; % piston removed %% it is not there! It has been overwritten by WFS
% zern.c = reshape(M*z.c(:),z.nMode,[]);
% phaseMV = reshape(zern.p*zern.c,nRes,nRes*length(ast));


% % % Zphi_LS_modal = pinv(DAst*Projection_alpha(idx,idy)) * StackSlopeLarge(:);
% % % Zphi_LS_modal_pbeta    = Projection_beta(idp,idy) * Zphi_LS_modal ;
% % % 
% % % % % % %     Z2U = calibDm.M*Dzcc;% Using zernike derivatives and slopes to DM command matrix
% % % U_LS = calibDm.M * Dzcc * Projection_beta(idp,idy) * Zphi_LS_modal ;

% 
% wfs.slopesUnits
% temp4
% dm.modes.modes
% F = dm.modes.modes;
% zern.modes
% Zmodes = zern.modes;
% u = pinv(F)*Zmodes;
% u = pinv(full(F))*Zmodes;
% figure,imagesc(Z2U);
% figure,imagesc(u);
% figure,imagesc([u(:,1),Z2U(:,1)]);
% figure,imagesc([(ngs.wavelength/2*pi)*u(:,1),Z2U(:,1)]);
% figure,imagesc([(ngs.wavelength/2*pi)*u(:,1),-Z2U(:,1)]);
% figure,plot([(ngs.wavelength/2*pi)*u(:,1),-Z2U(:,1)]);
% figure,plot([(ngs.wavelength/pi)*u(:,1),-Z2U(:,1)]);
% figure,plot([(ngs.wavelength/2*pi)*u(:,1),-2*Z2U(:,1)]);
% test_modal_tomography_openloop_newfs
% temp4
% zernTilt.c
% temp4
% ngs.wavelength
% temp4
% test_modal_tomography_openloop_newfs
% test_simplest_closed_and_open_loops
% test_modal_tomography_openloop_newfs
% frameRef = cam.referenceFrame;
% figure, imagesc(frameRef)
% frameActual = cam.reference;
% frameActual = cam.frame
% figure, imagesc(frameActual)
% test_modal_tomography_openloop_newfs
% temp4
% test_modal_tomography_openloop_newfs
% wfs.slopes
% test_modal_tomography_openloop_newfs