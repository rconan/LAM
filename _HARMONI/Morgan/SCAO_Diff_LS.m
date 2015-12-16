%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         RECONSTRUCTION  en  BOUCLE  FERMEE           %
%           COMMANDE :   INT?GRATEUR + LS              %
%                     M_Com                            %
%                 ASO  DIFFRACTIF                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; clc; close all; delete(gcp);
REP = '/data/RESULT_HARMONI/';
nb_exposureTrame = 1000; startDelay = 100; 
nb_boucle = startDelay + nb_exposureTrame;
% LARGEUR de FENTE : en Arcsec !
EE_SlitWidth = 0.02;
% NB de TRAMES de RETARD :
Trame = 3;
% ASO  DIFFRACTIF :
BruitLec = 1;% Variance en [ (Photo-electron/pixel/trame)? ] 
SeuilBruitLec = 1;% en [ Photo-electron/pixel ]
BruitPho = true;
% MIROIR M4 : R?partition des actionnneurs sur une grille [0] carr?e  ; [1] triangulaire
Typ_M4 = 0; 
% ESTIMATEUR : LS_MCom
Gain = 0.5;
% FILTRAGE des VALEURS SINGULIERES :
nbVS_MInt = 1;
% Etoile Guide Naturelle :
Photo_Ngs = photometry.R;
% Etoile Science :
Photo_Sci = photometry.H;
% AJOUT de PERTURBATIONS sur le M1 :
PERTUB = 0;
% AJOUT de l'ARAIGNEE :
ARAIG = 1;
largBras = 0.5;% en m?tres
alpha = 0/6;% 0 =< Alpha < pi/3 (sens des aiguilles d'une montre)
% TESTS sur les PETALES :
PETALE = 0;
% AFFICHAGE des PHASES :
AFF_PH = 0;

%% - TELESCOPE

% Diam?tres du M1 "UnVignetted"
diamExtUV = 38.542; diamIntUV = 11.067;
% Diam?tres du M1 "All-Glass"
diamExtAG = 36.903; diamIntAG = 11.208;
% Diam?tres utilis?s & def du tel :
diamE = diamExtAG; diamI = diamIntAG;
ObsCen = diamI/diamE;
nLens = 74; nb_pxssp = 10;
nPix = nLens*nb_pxssp; d_pitch = diamE/nLens; Res_1pix = diamE/nPix;
sampFreq = 500; sampTime = 1/sampFreq;
tel = telescope(diamE,'obstructionRatio',ObsCen,'resolution',nPix,'fieldOfViewInArcsec',30,'samplingTime',sampTime);
nb_px = tel.pixelArea;
MasqPupLog = tel.pupilLogical; MasqPup = single(tel.pupil);
 
%% - PERTURBATIONS sur le M1

if PERTUB == 0
    Pert_M1 = zeros(nPix,nPix,nb_boucle);
elseif PERTUB == 1
    %load([REP 'JEFF_HSFreq_740.mat']);
    %Pert_M1 = repmat(JEFF_HSFPess,1,1,nb_boucle);
    %load([REP 'ESO_WL2msBS_740_S.mat']);
    %load([REP 'ESO_WL2msHP_740_S.mat']);
    %load([REP 'ESO_WL2msSP_740_S.mat']);
    load([REP 'MG_WL8msTT_740_S2.mat']);
    %load([REP 'MG_WL8msFC_740_S.mat']);
    %load([REP 'MG_Transient_740_S.mat']);
    Pert_M1 = MG_WL8msTT; clear MG_WL8msTT
end

%% - PETALES d?finis par l'araign?e

if PETALE == 1
	MasqPetales = harmoniTools.EELT_PETALES(Res_1pix,alpha,largBras,MasqPupLog);
%     for kP = 1:6
%         figure(1), spy(MasqPetales(:,:,kP)); pause
%     end
%     MasqSpid = sum(MasqPetales,3);
%     figure(1), spy(MasqSpid)
    PistTur = zeros(nb_boucle,6);
    PistRes = zeros(nb_boucle,6);
end

%% - ATMOSPHERE g?n?rant la TURBULENCE

seeing = 0.65; SIMUL = 'TURs065_NoNo';
lambdaV0 = photometry.V0.wavelength;% ? 0.5 um !
AngZenith = pi/6; r0V0 = lambdaV0/(seeing/3600/180*pi)*(cos(AngZenith))^(3/5);
L0 = 50; 
Alti = 0;
Frac_r0 = 1;
Wind_S = 12.5;
Wind_D = 0;
atm = atmosphere(photometry.V0,r0V0,L0,'altitude',Alti,'fractionnalR0',Frac_r0,'windSpeed',Wind_S,'windDirection',Wind_D);

%% - OBJET SCIENTIFIQUE

Sci = source('zenith',0,'azimuth',0,'wavelength',Photo_Sci,'magnitude',5); 
lambdaSci = Sci.wavelength;
r0Sci = (lambdaSci/lambdaV0)^(6/5)*r0V0;
RadScivMicr = lambdaSci/2/pi*1e6;
k_Sci = Sci.waveNumber;

%% - ETOILE GUIDE NATURELLE

Ngs = source('zenith',0,'azimuth',0,'wavelength',Photo_Ngs,'magnitude',12);
lambdaNgs = Ngs.wavelength;
k_Ngs = Ngs.waveNumber;

%% - SOURCES DE CALIBRATION

NgsCal = source('wavelength',Photo_Ngs);% pour ASO
SciCal = source('wavelength',Photo_Sci);% pour Camera

%% - CALIBRATION  ASO DIFFRACTIF & CAMERA

seuil_ssp = 0.5;% Minimum light ratio
wfs = shackHartmann(nLens,nPix,seuil_ssp);
NgsCal = NgsCal.*tel*wfs;% Propagation of the calibration source to the WFS through the telescope
wfs.INIT
+wfs;% A new frame read-out and slopes computing
masq_ssp = wfs.validLenslet;% Map of valid sub-apertures
nb_ssp = wfs.nValidLenslet;% nb of valid sub-apertures
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

if Typ_M4 == 0;% r?partition carr?e
    couplage = 0.4; nAct = nLens + 1;
    masq_act = wfs.validActuator;
    B_IF = influenceFunction('monotonic',couplage);
    dm = deformableMirror(nAct,'modes',B_IF,'resolution',nPix,'validActuator',masq_act);
    nb_act = dm.nValidActuator;
    % figure(3); show(B_IF)
elseif Typ_M4 == 1% r?partition hexagonale
    
end
% Matrice des TENSIONS -> PHASE PIXELS @ lambdaSci !
M_UvPixSci = (4*pi/lambdaSci)*dm.modes.modes(MasqPupLog,:); 

%% - MATRICE d' INTERACTION ( nb_p ; nb_act ) & MATRICE de COMMANDE ( nb_act ; nb_p )

if ARAIG == 1
    [MasqLegs,MasqSpider] = harmoniTools.EELT_SPIDER(nPix,Res_1pix,alpha,largBras,diamI,diamE);
%     for k = 1:6
%         figure(1), spy(MasqLegs(:,:,k)); pause
%     end
%     figure(1), spy(MasqSpider)
    tel.pupil = tel.pupil.*MasqSpider;
    MasqPup = single(tel.pupil);
end
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

%% - 

tel = tel + atm;

%% - IMAGEUR pour l'OBJET SCIENTIFIQUE 

cam = imager();
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
%cam.eeWidth = EE_SlitWidth;
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
VarBruit = mean(diag(slopes*slopes'/nMeas));% en [ Pix? ]
% toto = theoreticalNoise(wfs,tel,atm,Ngs,Sci);
tel = tel + atm;

%% SANS PERTURBATIONS sur le M1

% if PERTUB == 0
%     set(Sci,'logging',true)% The logging of the wavefront variance of the science object is turned on.
%     set(Sci,'phaseVar',[])% The phase variance will be stored into the phaseVar property.
%     Ngs = Ngs.*tel*{MasqPup;zeros(nRes)}*dm*wfs;
%     Sci = Sci.*tel*{MasqPup;zeros(nRes)}*dm*cam;
% end

%%  AFFICHAGE des Phases

if AFF_PH == 1
    figure(100)% en MICRONS
    MinMax = [-0.5 0.5];
    phases = imagesc([zeros(nRes),zeros(nRes),zeros(nPix)]);%
    ax = gca; axis xy equal tight; caxis(MinMax); ylabel(colorbar,'WFE  [ \mum ]'); snapnow
end

                %% - FERMETURE de la BOUCLE OA

display(Gain);
Ures = zeros(nb_act,1); U_k = Ures; U_kp1 = U_k; U_km1 = U_k;
% CP_TUR = zeros(nPix,nPix,nb_boucle);
% CP_RES = zeros(nPix,nPix,nb_boucle);
tic
for k_bou = 1:nb_boucle
%%  PHASE TURBULENTE & PHASE ASO
    +tel;
%     if PERTUB == 0% SANS Perturbations sur le M1
%         +Ngs; +Sci;
%     elseif PERTUB == 1% AVEC ou SANS de Perturbations sur le M1
    Ngs = Ngs.*tel*{MasqPup;Pert_M1(:,:,k_bou)*k_Ngs}*dm*wfs;
    Sci = Sci.*tel; CPTUR = Sci.meanRmPhase;
    Sci = Sci*{MasqPup;Pert_M1(:,:,k_bou)*k_Sci}*dm*cam;
%     end
    CPCOR = harmoniTools.CARTES_2D(M_UvPixSci*dm.coefs,MasqPupLog);
    CPRES = CPTUR - CPCOR;
%     CP_TUR(:,:,k_bou) = CPTUR;
%     CP_RES(:,:,k_bou) = CPRES;     
%%  AFFICHAGE des PHASES en MICRONS @ lambdaSci
    if AFF_PH == 1
        figure(100)
        set(phases,'cdata',[CPTUR,CPCOR,CPRES]*RadScivMicr);% 
        title(ax,sprintf('n? %4d / %4d',k_bou,nb_boucle)); drawnow
        %saveas(gcf,['PhaseRES_' num2str(k_bou)],'jpg');
    end
%%  PISTONS & TIP / TILT sur les PETALES
    if PETALE == 1
        for kP = 1:6
            PistRes(k_bou,kP) = RadScivMicr*mean(CPRES(MasqPetales(:,:,kP)));
            PistTur(k_bou,kP) = RadScivMicr*mean(CPTUR(MasqPetales(:,:,kP)));
        end
    end
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
Duree = uint16(toc)
figure(10)
imagesc(cam)
%saveas(gcf,['HarmoniPSF_' SIMUL],'jpg');

%% PERFORMANCES  SIMULATION

% if PERTUB == 0
%     Var_PhRes = reshape(Sci.phaseVar(1:nb_boucle*2),2,[])';
%     E_Coh = exp(-mean(Var_PhRes(startDelay+1:end,2)))*100
% end
Strehl_PSF = 100*cam.strehl
EE = cam.ee;

%% SAUVEGARDES

% CP_TUR = single(CP_TUR);
% CP_RES = single(CP_RES);
% save([REP 'PHASTUR_' SIMUL '.mat'],'CP_TUR','-v7.3');
% save([REP 'PHASRES_' SIMUL '.mat'],'CP_RES','-v7.3');
% save([REP 'PERF_' SIMUL '.mat'],'seeing','r0V0','lambdaSci','r0Sci','VarBruit','EE_SlitWidth','Strehl_PSF','EE','Duree');
% save([REP 'PISTONS_' SIMUL '.mat'],'seeing','r0V0','lambdaSci','r0Sci','VarBruit','Strehl_PSF','Duree','PistTur','PistRes');
