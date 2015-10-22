close all

z_NdmLowresValidActuators = dmLowRes.nValidActuator;
z_NwfsValidLenslets = wfs.nValidLenslet;

z_wfsSlopes = wfs.slopes;

z_calibDm_D = calibDm.D;  %%% interaction matrix WFS2DM, 152x97
z_calibDm_M = calibDm.M; %% pseudoinverse of interaction matrix   97x152

z_wfsCameraFrame= wfs.camera.frame;

figure, imagesc(wfs.camera.frame) %%% the WFS camera frame
figure, mesh(wfs.camera.frame) %%% the WFS camera frame

figure, imagesc(wfs.finiteDifferenceWavefront ) %%% the WFS camera frame


figure, imagesc(sciStar.meanRmPhase) %%% the WFS camera frame
figure, imagesc(dm.surface) %%% the WFS camera frame

break
%%%% Generate the Zernike with only modes 2 and 3 - that is, Tip and Tilt
zern = zernike(2:3,'resolution',nLenslet+1, 'pupil',dmLowRes.validActuator);
zern_modes = zern.modes(dmLowRes.validActuator,:);

F = bifaLowRes.modes(wfs.validActuator,:); %%% WTF? Should we multiply it by lambda/2pi?, 97x97
iF = 0.5*pinv(full(F));

M_TT_modal = iF*zern_modes;

%      if     (params.noise_in_wfs  == 0); 
P_a = M_TT_modal*pinv(M_TT_modal);   % Projection matrix for Tipt-Tilt removal _in the actuator space_
DMTT_Remove = eye(dmLowRes.nValidActuator) - P_a;
