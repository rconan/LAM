%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         RECONSTRUCTION  en  BOUCLE  FERMEE           %
%           COMMANDE :   INTÉGRATEUR + LS              %
%                     M_Com                            %
%                 ASO  DIFFRACTIF                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; clc; close all;
REP = 'C:\Users\mgray\AO_HARMONI\RESULTATS\';
% REP = '/data/RESULT_HARMONI/';
nb_exposureTrame = 1000; startDelay = 100; 
nb_boucle = startDelay + nb_exposureTrame; i0varM = startDelay+1;
% LARGEUR de FENTE : en Arcsec !
EE_SlitWidth = 0.02; 
% NB de TRAMES de RETARD :
Trame = 1;
% ASO  DIFFRACTIF :
BruitLec = 1;% Variance en [ (Photo-electron/pixel/trame)² ] 
SeuilBruitLec = 1;% en [ Photo-electron/pixel ]
BruitPho = true;
% MIROIR M4 : Répartition des actionnneurs sur une grille [0] carrée  ; [1] triangulaire
Typ_M4 = 0; 
% ESTIMATEUR : LS_MCom
Gain = 0.5;
% FILTRAGE des VALEURS SINGULIERES :
nbVS_MInt = 1;
% BANDES : Science / Etoile Guide Naturelle
Photo_Sci = photometry.H; Photo_Ngs = photometry.R;
% AFFICHAGE des PHASES : [1] ou [2]
AFF_PH = 0;

%% - TELESCOPE

DiamExtMin = 36.903; DiamIntMax = 11.208;
diam = DiamExtMin;
ObsCen = DiamIntMax/DiamExtMin;
nLens = 74;
nb_pxssp = 10;
nRes = nLens*nb_pxssp; d_pitch = diam/nLens;
sampFreq = 500; sampTime = 1/sampFreq;
tel = telescope(diam,'obstructionRatio',ObsCen,'resolution',nRes,'fieldOfViewInArcsec',30,'samplingTime',sampTime);
Masque = tel.pupil; MasqueLog = tel.pupilLogical;
nb_px = tel.pixelArea;

%% - PERTURBATIONS sur le M1

% load([REP 'JEFF_HSFreq_740.mat']);
% Pert_M1 = repmat(JEFF_HSFPess,1,1,nb_boucle);

% load([REP 'ESO_WL2msBS_740_S.mat']);
% load([REP 'ESO_WL2msHP_740_S.mat']);
% load([REP 'ESO_WL2msSP_740_S.mat']);
% load([REP 'MG_WL8msTT_740_S.mat']);
% load([REP 'MG_WL8msFC_740_S.mat']);
% load([REP 'MG_Transient_740_S.mat']);
% Pert_M1 = MG_WL8msTT; clear MG_WL8msTT

% SIMUL = 'NOTUR_NOBRUIT_TT_1T';

%% - ATMOSPHERE générant la TURBULENCE

seeing = 0.44; lambdaV0 = 0.5e-6;% à 0.5 um !
AngZenith = pi/6; r0V0 = lambdaV0/(seeing/3600/180*pi)*(cos(AngZenith))^(3/5);
L0 = 50; 
Alti = 0;
Frac_r0 = 1;
Wind_S = 12.5;
Wind_D = 0;
atm = atmosphere(photometry.V0,r0V0,L0,'altitude',Alti,'fractionnalR0',Frac_r0,'windSpeed',Wind_S,'windDirection',Wind_D);

%% - OBJETS SCIENTIFIQUES

Sci = source('zenith',0,'azimuth',0,'wavelength',Photo_Sci,'magnitude',5); 
lambdaSci = Photo_Sci.wavelength;
r0Sci = (lambdaSci/lambdaV0)^(6/5)*r0V0;
RadScivMicr = lambdaSci/2/pi*1e6;

%% - ETOILES GUIDES NATURELLES

Ngs = source('zenith',0,'azimuth',0,'wavelength',Photo_Ngs,'magnitude',12);
lambdaNgs = Photo_Ngs.wavelength;

%% - SOURCES DE CALIBRATION

NgsCal = source('wavelength',Photo_Ngs);% pour ASO
SciCal = source('wavelength',Photo_Sci);% pour Camera

%% - CALIBRATION  ASO DIFFRACTIF & CAMERA

seuil_ssp = 0.5;% Minimum light ratio
wfs = shackHartmann(nLens,nRes,seuil_ssp);
NgsCal = NgsCal.*tel*wfs;% Propagation of the calibration source to the WFS through the telescope
wfs.INIT
+wfs;% A new frame read-out and slopes computing
masq_ssp = wfs.validLenslet;% Map of valid sub-apertures
nb_ssp = wfs.nValidLenslet;% nb de sous-pupilles effectives
% figure(1)
% subplot(1,2,1); imagesc(wfs.camera);% WFS camera display
% subplot(1,2,2); slopesDisplay(wfs);% WFS slopes display
% Allow the displays of the frame and of the slopes to be updated when a new frame and new slopes are computed
wfs.camera.frameListener.Enabled = true;
wfs.slopesListener.Enabled = true;
% WFS must be calibrated such as for 1rd of tip--tilt wavefront , it will measured a slopes of 1rd.
% whereas the source is progressively moved off-axis
wfs.pointingDirection = zeros(2,1);
pixelScale = NgsCal.wavelength/(2*d_pitch*wfs.lenslets.nyquistSampling);
tipStep = pixelScale/2;
nStep   = floor(nb_pxssp/3)*2;
sx      = zeros(1,nStep+1);
u       = 0:nStep;
wfs.camera.frameListener.Enabled = false;
wfs.slopesListener.Enabled = false;
warning('off','oomao:shackHartmann:relay')
for kStep=u
    NgsCal.zenith = -tipStep*kStep;
    +NgsCal;
    drawnow
    sx(kStep+1) = median(wfs.slopes(1:end/2));
end
warning('on','oomao:shackHartmann:relay')
Ox_in  = u*tipStep*constants.radian2arcsec;
Ox_out = sx*NgsCal.wavelength/d_pitch/2*constants.radian2arcsec;
% figure(2)
% plot(Ox_in,Ox_out); grid
slopesLinCoef = polyfit(Ox_in,Ox_out,1);
wfs.slopesUnits = 1/slopesLinCoef(1);
% The source is reset on-axis and the WFS is set to always be aligned to
% the source by setting wfs.pointingDirection to empty.
NgsCal.zenith = 0;
wfs.pointingDirection = [];

%% - MIROIR DEFORMABLE

if Typ_M4 == 0;
    couplage = 0.4; nAct = nLens + 1;
    wfs2 = shackHartmann(nLens,nRes,0);
    NgsCal = NgsCal.*tel*wfs2; wfs2.INIT; +wfs2;
    masq_act = wfs2.validActuator;
    B_IF = influenceFunction('monotonic',couplage);
    dm = deformableMirror(nAct,'modes',B_IF,'resolution',nRes,'validActuator',masq_act);
    % figure(3); show(B_IF)
    nb_act = dm.nValidActuator;
elseif Typ_M4 == 1
    
end

%% - MATRICE d' INTERACTION ( nb_p ; nb_act ) & MATRICE de COMMANDE ( nb_act ; nb_p ) 

wfs.camera.frameListener.Enabled = false;% Switch off the display automatic update
wfs.slopesListener.Enabled = false;
NgsCal = NgsCal.*tel;% setup the optical path before the DM/WFS subsystem
calibWfsDm = calibration(dm,wfs,NgsCal,NgsCal.wavelength,nAct,'cond',1e2);
% calibWfsDm = calibration(dm,wfs,NgsCal,NgsCal.wavelength,nAct)
% calibWfsDm.nThresholded = nbVS_MInt;
VS_Int = calibWfsDm.eigenValues; Cdt_Int = calibWfsDm.cond;
M_Int = calibWfsDm.D;% MATRICE d' INTERACTION ( nb_p ; nb_act )
M_Com = calibWfsDm.M;% MATRICE de COMMANDE ( nb_act ; nb_p )
display(['Matrice INTERACTION :  Conditionnement = ',num2str(Cdt_Int)]);

%% - ESTIMATION  ZONALE

tel = tel + atm;
% figure(5); imagesc(tel)
NgsCal = NgsCal.*tel*wfs;
% The wavefront reconstruction is done by estimating the wavefront values at the corner of the lenslets.
% A new telescope is defined identical to the previous one but with a lower resolution and a pupil defined by the map of "valid actuator".
telLowRes = telescope(tel.D,'resolution',nAct,'fieldOfViewInArcsec',30,'samplingTime',sampTime);
telLowRes.pupil = wfs.validActuator;
telLowRes = telLowRes + atm;
NgsCal = NgsCal.*telLowRes;
B_IF_LowRes = influenceFunction('monotonic',couplage);
dm_LowRes = deformableMirror(nAct,'modes',B_IF_LowRes,'resolution',nAct,'validActuator',masq_act);

%% - IMAGEUR pour l'OBJET SCIENTIFIQUE 

cam = imager(tel);
tel = tel - atm;
% The science object is propagated through the telescope to the science camera
% producing a perfect diffraction limited image
SciCal = SciCal.*tel*cam;
figure(10)
imagesc(cam,'parent',subplot(2,1,1))
% cam.frameListener.Enabled = true;
% The diffraction limited image is set as the reference frame allowing to
% compute the image Strehl ratio on the subsequent frame captures.
cam.referenceFrame = cam.frame;
+SciCal;
fprintf('Strehl ratio: %4.1f\n',cam.strehl)
% The atmosphere is attached to the telescope, the science object is propagated again,
% through the atmosphere and the telescope: the atmosphere limited Strehl is obtained.
tel = tel + atm;
+SciCal;
fprintf('Strehl ratio: %4.1f\n',cam.strehl)

%% - INITIALISATION  de  l' IMAGEUR

flush(cam)
cam.frame = cam.frame*0;
cam.clockRate = 1;% means it is slaved to the sampling rate of the telescope ( sampFreq )
cam.exposureTime = nb_exposureTrame;
cam.eeWidth = EE_SlitWidth;
figure(10)
imagesc(cam,'parent',subplot(2,1,1))
% cam.frameListener.Enabled = true;
subplot(2,1,2)
Ima_Sci = imagesc(catMeanRmPhase(Sci));
axis xy equal tight; colorbar
% the DM coefficients are set to 0 and the science detector is reset
dm.coefs = 0;
flush(cam)
cam.startDelay = startDelay;
cam.frameListener.Enabled = false;
tel = tel - atm;
NgsCal.magnitude = Ngs.magnitude;
wfs.camera.photonNoise = BruitPho;
wfs.camera.readOutNoise = sqrt(BruitLec);
wfs.framePixelThreshold = SeuilBruitLec;
NgsCal = NgsCal.*tel*wfs;% re-propagate the source
% figure(4)
% intensityDisplay(wfs);% display the subaperture intensity
nMeas = 1000; slopes = zeros(wfs.nSlope,nMeas);
for kMeas = 1:nMeas
    +wfs; slopes(:,kMeas) = wfs.slopes;
end
VarBruit = mean(diag(slopes*slopes'/nMeas));% en [ Pix² ]
% tel = tel + atm;

%% SANS PERTURBATIONS sur le M1

% set(Sci,'logging',true)% The logging of the wavefront variance of the science object is turned on.
% set(Sci,'phaseVar',[])% The phase variance will be stored into the phaseVar property.
% Ngs = Ngs.*tel*dm*wfs;
% Sci = Sci.*tel*dm*cam;

                %% - FERMETURE de la BOUCLE OA

display(Gain);
% CP_TUR = zeros(nRes,nRes,nb_boucle);
% CP_RES = zeros(nRes,nRes,nb_boucle);
Ures = zeros(nb_act,1); U_k = Ures; U_kp1 = U_k; U_km1 = U_k;
% if AFF_PH == 1 % en MICRONS
%     figure(11)
%     MinMax = [-5 5];
%     Phases = imagesc([zeros(nRes),zeros(nRes),zeros(nRes)]);
%     ax = gca; axis xy equal tight; caxis(MinMax); ylabel(colorbar,'WFE  [ \mum ]'); snapnow
% end            
tic
for k_bou = 1:nb_boucle
%% PHASE TURBULENTE & PHASE ASO
    +tel;
    % SANS Perturbations sur le M1
%     +Ngs; +Sci;
    % AVEC Perturbations sur le M1
    Ngs = Ngs.*tel*{Masque;Pert_M1(:,:,k_bou)*Ngs.waveNumber}*dm*wfs;
    Sci = Sci.*tel;% CP_TUR(:,:,k_bou) = Sci.meanRmOpd;
    Sci = Sci*{Masque;Pert_M1(:,:,k_bou)*Sci.waveNumber}*dm;
    % CP_RES(:,:,k_bou) = Sci.meanRmOpd;
    Sci = Sci*cam;
    % AFFICHAGE des Phases
%     if AFF_PH == 1 % en MICRONS
%         CPCOR = CPTUR - CPRES;
%         set(Phases,'cdata',[CPTUR,CPCOR,CPRES]*RadSci1vMicr);
%         title(ax,sprintf('n° %4d / %4d',k_bou,nb_boucle)); drawnow
%     elseif AFF_PH == 2
%         h = imagesc(catMeanRmPhase(scienceCombo));
%         axis xy equal tight; colorbar
%     end
%% TENSIONS 
    Ures = - M_Com*wfs.slopes;
    if Trame == 1
        dm.coefs = dm.coefs + Gain*Ures;
    elseif Trame == 2
        U_kp1 = U_k + Gain*Ures;
        dm.coefs = U_k; U_k = U_kp1;
    elseif Trame == 3
        U_kp1 = U_k + Gain*Ures;
        dm.coefs = U_km1; U_km1 = U_k; U_k = U_kp1;
    end
end % FIN d' UNE simulation
Duree = uint16(toc);
imagesc(cam)

%% PERFORMANCES  SIMULATION

% Var_PhRes = reshape(Sci.phaseVar(1:nb_boucle*2),2,[])';
% E_Coh = exp(-mean(Var_PhRes(i0varM:end,2)))*100
VarBruit
Strehl_PSF = 100*cam.strehl
EE = cam.ee

%% PERFORMANCES  THEORIQUES

% [Sig_Fit,Sig_Alias,Sig_Bruit,Sig_Temp,Strehl_TH] = SCAO_PERFTh(Typ_M4,d_pitch,lambdaSci,r0Sci,lambdaNgs,VarBruit);

% CP_TUR = single(CP_TUR);
% CP_RES = single(CP_RES);
% save([REP 'PHASTUR_' SIMUL '.mat'],'CP_TUR','-v7.3');
% save([REP 'PHASRES_' SIMUL '.mat'],'CP_RES','-v7.3');
% save([REP 'PERF_' SIMUL '.mat'],'seeing','r0V0','lambdaSci','r0Sci','EE_SlitWidth',...
%     'Strehl_PSF','EE','Sig_Fit','Sig_Alias','Sig_Bruit','Sig_Temp','Strehl_TH');
