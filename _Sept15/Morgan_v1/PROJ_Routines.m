clear; clc; close all;
REP_Result = 'C:\Users\mgray\AO_HARMONI\RESULTATS\';
REP_Pert = 'C:\Users\mgray\AO_HARMONI\PERT_DYNAMIQUES\VERSION_1_3\';
REP_M4 = 'C:\Users\mgray\AO_HARMONI\MIROIR_M4\CARRE\';
ROUT = 2;

%% - TELESCOPE

DiamExtMin = 36.903; DiamIntMax = 11.208;
diam = DiamExtMin;
ObsCen = DiamIntMax/DiamExtMin;
nLens = 74;
nb_pxssp = 10;
nRes = nLens*nb_pxssp;
tel = telescope(diam,'obstructionRatio',ObsCen,'resolution',nRes,'fieldOfViewInArcsec',30,'samplingTime',1/500);
masque = tel.pupilLogical; nb_px = tel.pixelArea;

%% PROJECTION des PHASES sur TOUS les ZERNIKE

if ROUT == 1
    % PHASES
    load([REP_Result 'MG_WL8msTT_740_S.mat']);
    PHASES = MG_WL8msTT; clear MG_WL8msTT
    nb_boucle = size(PHASES,3);
    load([REP_Pert 'MSwind8mps_tiptilt.mat']);
    Coef1Zer = MSwindtiptiltData(1:nb_boucle,2:3)';
    % MATRICE de PROJECTION sur les ZERNIKE
    ORMaxPU = 3;
    n_zPU = (ORMaxPU+1)*(ORMaxPU+2)/2;% nb de Zernike AVEC le PISTON !
    Zernik = zernike(1:n_zPU,diam,'resolution',nRes,'pupil',masque);% AVEC le PISTON !
    ZER = zeros(nb_px,n_zPU);% cartes de phase des zernike
    Ord_rad = zeros(n_zPU,1);% liste des ordres radiaux
    for k = 1:n_zPU
        Zern_k = Zernik.p(:,k); ZER(:,k) = Zern_k(masque);
        Ord_rad(k)= Zernik.n(k);
    end
    clear Zern_k ZernikPix
    ZZ = ZER'*ZER;
    [UZ,SZ]= eig(ZZ);
    UZ = fliplr(UZ);
    VPZ = diag(rot90(SZ,2));
    condZ = VPZ(1)/VPZ(n_zPU)
    iSZ = diag(1./VPZ);
    PROJ_ZER = UZ*iSZ*UZ'*ZER';
    clear ZZ UZ SZ iSZ ZER
    % PROJECTION des PHASES
    CoefZer = zeros(n_zPU,nb_boucle);
    for k_bou = 1:nb_boucle
        CPPERT = PHASES(:,:,k_bou); CP1 = CPPERT(masque);
        CoefZer(:,k_bou) = PROJ_ZER*CP1;
    end
    % COMPARAISON
    Coef2Zer = CoefZer(2:3,:)/46.7e-6;
    Err_Tip = (Coef1Zer(1,:) - Coef2Zer(1,:))./Coef1Zer(1,:);
    Err_Tilt = (Coef1Zer(2,:) + Coef2Zer(2,:))./Coef1Zer(2,:);
    
%% PROJECTION sur les FONCTIONS d'INFLUENCE

elseif ROUT == 2
    MatUvPix = 'MATUvPix_1.mat';
    if exist([REP_M4 MatUvPix],'file'), load([REP_M4 MatUvPix]);
    else
        wfs = shackHartmann(nLens,nRes,0.5);
        nAct = nLens + 1;
        masq_act = wfs.validActuator;
        B_IF = influenceFunction('monotonic',0.4);
        dm = deformableMirror(nAct,'modes',B_IF,'resolution',nRes,'validActuator',masq_act);
        M_UvPix = dm.modes.modes(masque,:);
        % INVERSION Type 1
        M_PixvU = pinv(full(M_UvPix));
        save('MATUvPix_1.mat','M_UvPix','M_PixvU','-v7.3');
        % INVERSION Type 2
%         [M_PixvU,VS,Cdt_MUvPix] = MAT_INVGEN(M_UvPix,0,1);% (M,nb_VS,Type)
%         display(['Conditionnement = ',num2str(Cdt_MUvPix)]);
%         save('MATUvPix_2.mat','M_UvPix','M_PixvU','nb_VS','Cdt_MUvPix','-v7.3');
    end
    M_UvPix = full(M_UvPix);
    load([REP_Result 'MG_WL8msTT_740_S.mat']);
    CPPERT = MG_WL8msTT(:,:,100);% en mètres
    clear MG_WL8msTT
    load([REP_Pert 'MSwind8mps_tiptilt.mat']);
    CoefTilt = -MSwindtiptiltData(100,3);% en [Arcsec] & 1 Arcsec <-> 46.7e-6 m RMS 
    CP1 = 260*(CPPERT(masque)/46.7e-6)/CoefTilt;% RMS = 1
    CP2a = M_PixvU*CP1; clear M_PixvU
    CP2 = M_UvPix*CP2a;
    CPERT1 = CARTES_2D(CP1,masque);
    CPERT2 = CARTES_2D(CP2,masque);
    CPRES = CPERT1 - CPERT2;
    figure(1), imagesc(CPERT1); axis xy equal tight; colorbar
    figure(2), imagesc(CPERT2); axis xy equal tight; colorbar
    figure(3), imagesc(CPRES); axis xy equal tight; caxis([-50 50]); colorbar
    figure(4), imagesc([CPERT1 CPERT2 CPRES]); axis xy equal tight; colorbar
    PtoV = max(CPRES(:)) - min(CPRES(:));
    RMS = sqrt(mean((CP1-CP2).^2));
end
