%%%% this is used to be the temp4.m file;
clc


%% 3. \phi_z^{res} = Z^\dagger * \phi_{pix}^{res}
zsciStars_rms_phase = sciStars.meanRmPhase;
% figure, mesh(zsciStars_rms_phase);
zphi_res_modes = pinv(zern.p)*reshape(sciStars.meanRmPhase, nRes*nRes,1);
figure, bar(zphi_res_modes);
% break

Sigma_z = zernikeStats.covariance(zern,atm); %|sigma_Z
fprintf('\n Atm zernike maxRadialDegree=%d, modes=%d\n ', maxRadialDegree, zern.nMode);
fprintf('\n Projected zernike maxRadialDegree=%d, modes=%d\n ', maxRadialDegreeProj, zernProj.nMode);

diff_Covar = trace(Sigma_z) - trace(Sigma_alpha) %% checking the covariance difference between the initial Zernike modes and the projected (Sigma_\alpha == P_\alpha * Sigma_phi_z P_alpha^T)
diff_Covar_relative = (trace(Sigma_z) - trace(Sigma_alpha))/trace(Sigma_z) %% checking the covariance difference between the initial Zernike modes and the projected (Sigma_\alpha == P_\alpha * Sigma_phi_z P_alpha^T)


if nGs ==1      %%% This is for one GuideStar    
    StackSlopeLarge(1:nSlope,1) = wfs.slopes;
end

zPhi_LS = pinv(DAst)*StackSlopeLarge;
zdm_LS = Z2U*zPhi_LS;


% Rmv = (Sigma_betaalpha*DAst')/(DAst*Sigma_alpha*DAst' + CnAst);
% Sigma_alpha = Projection_alpha(idx,idy)*Sigma_phi(ids,ids)*Projection_alpha(idx,idy)'; %% Computing the projection of the covariance in the direction of Guidestars
% Sigma_betaalpha = Projection_beta(idp,idy)*Sigma_phi*Projection_alpha(idx,idy)'; %% Computing the Projection for the sci object:

zPalpha = Projection_alpha;
zPalpha( abs(zPalpha) < 1e-5 ) = 0;
zPalpha = sparse(zPalpha);
% figure, spy(zPalpha);

Rmv_notomo = (Sigma_z*DAst')/(DAst*Sigma_z*DAst' + CnAst);

norm (Rmv_notomo - Rmv)
norm (Rmv_notomo - Rmv) / norm (Rmv_notomo )

break 

% break
F = full(dm.modes.modes);
cond(full(F))

Fz =  0.5*pinv(full(F))*zern.modes*ngs.wavelength/(2*pi);   
cond(Fz)

Fz_info = tool_matrix_info(Fz)


%% Computing the Pseudoinverse via SVD
svd_threshold = 1e1;
[U,S,V] = svd( full(F),0); %% min singular value 0.7 max 40, economy size SVD
S( abs(S) < svd_threshold ) = 0;
S_new = S;
% F_back = U*S*V';
%%%%% Now building an inverse matrix
Fz_svd = (V*pinv(S_new)*U' ) * zern.modes*ngs.wavelength/(2*pi);

% \end{equation}
% $g_{pol}$ is the low--pass filter gain, $F$ is the matrix of the DM
% influence functions and $D$ is the poke matrix.
% </latex>

F_spie = 2*full( bif.modes(wfs.validActuator,:)) ;
Fz_spie = pinv(full(F_spie),1e-1);  %% iF in SPIE

%     % Pseudo-open-loop controller
%     dm.coefs(:,2) = (1-gain_pol)*dm.coefs(:,2) + ...
%         gain_pol*iF*( slmmse*( wfs.slopes(:,2) - calibDm.D*dm.coefs(:,2) ) );

break 
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



%% wavefront reconstruction minimum variance
M = Sangcov_matrix*DzAst'/(DzAst*Sangcov_matrix*DzAst'+CznAst);
% z = z - 1; % piston removed %% it is not there! It has been overwritten by WFS
zern.c = reshape(M*z.c(:),z.nMode,[]);
phaseMV = reshape(zern.p*zern.c,nRes,nRes*length(ast));


Zphi_LS_modal = pinv(DAst*Projection_alpha(idx,idy)) * StackSlopeLarge(:);
Zphi_LS_modal_pbeta    = Projection_beta(idp,idy) * Zphi_LS_modal ;

% % % %     Z2U = calibDm.M*Dzcc;% Using zernike derivatives and slopes to DM command matrix
U_LS = calibDm.M * Dzcc * Projection_beta(idp,idy) * Zphi_LS_modal ;


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