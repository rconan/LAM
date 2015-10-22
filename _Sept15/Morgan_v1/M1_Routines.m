clear all; clc; close all;
REP = '/data/HARMONI/SCAO/SIMUL_MORGAN/PERTURB_M1/';
nRes = 740; nb_boucle = 1100;
ROUT = 4;

%% VISUALISATION de Fichiers .FITS
if ROUT == 1
    PhaseFits = fitsread([REP '' num2str(nRes) '_' num2str(1) '.fits']);
    figure(1)
    phase = imagesc(PhaseFits*1e9);
    ax = gca; axis square; ylabel(colorbar,'nanometres'); snapnow
    for k_bou = 2:nb_boucle
        PhaseFits = fitsread([REP '' num2str(nRes) '_' num2str(k_bou) '.fits']);
        set(phase,'cdata',PhaseFits*1e9);
        title(ax,sprintf('n %4d',k_bou)); drawnow;
        %pause
    end

%% CONVERSION de Fichiers .FITS -> .MAT
elseif ROUT == 2
    ESO_WL2msSP = zeros(nRes,nRes,nb_boucle);
    for k_bou = 1:nb_boucle
        k_bou
        ESO_WL2msSP(:,:,k_bou) = fitsread([REP 'ESO_WL2msSP_' num2str(nRes) '_' num2str(k_bou) '.fits']);
    end
    ESO_WL2msSP = single(ESO_WL2msSP);
    save([REP 'ESO_WL2msSP_' num2str(nRes) Typ],'ESO_WL2msSP','-v7.3');

%% INTERPOLATION de Fichiers .MAT
elseif ROUT == 3
    nRes1 = 972; nRes2 = 740;
    load([REP 'MG_Transient_' num2str(nRes) '_S'])
    load([REP 'MG_M1AllGlass_' num2str(nRes2)])
    load([REP 'ESO_M1Pupil_' num2str(nRes)])
    [X1,Y1] = meshgrid(1:nRes1,1:nRes1);
    echant = linspace(1,nRes1,nRes2);
    [X2,Y2] = meshgrid(echant,echant);
    Transient = zeros(nRes2,nRes2,nb_boucle);
    for k_bou = 1:nb_boucle;
        k_bou
        Carte0 = MG_Transient(:,:,k_bou).*Masq_AllGlass;%1014*1014
        Carte1 = Carte0(22:993,22:993);
        Transient(:,:,k_bou) = interp2(X1,Y1,Carte1,X2,Y2,'spline');
        Transient(:,:,k_bou) = Transient(:,:,k_bou).*masqAllGlass;% 740*740
    end
    MG_Transient = single(Transient);
    save(['MG_Transient_' num2str(nRes2) '_S'],'MG_Transient','-v7.3');

%% CONVERSION de Fichiers .MAT -> .FITS
elseif ROUT == 4
    REP2 = [REP 'Transient/'];
    nom = 'MG_Transient';
    load([REP nom '.mat']);
    for k_bou = 1:nb_boucle
        k_bou
        fitswrite(MG_Transient(:,:,k_bou),[REP2 nom '_' num2str(k_bou) '.fits']);
    end
    
%% VISUALISATION de Fichiers .MAT
elseif ROUT == 5
    % load([REP 'MG_WL8msTT_' num2str(nRes)]);
    % load([REP 'MG_WL8msFC_' num2str(nRes)]);
    load([REP 'MG_Transient_' num2str(nRes) '_S.mat']);
    figure(1)
    phase = imagesc([MG_Transient(:,:,1)*1e6]);
    %phase = imagesc(CP_RES(:,:,1));% ESO_WL2msHP(:,:,1)
    ax = gca; axis equal tight; caxis([-40 40]); colorbar; snapnow
    for k_bou = 2:nb_boucle
        %set(phase,'cdata',[MG_WL8msTT(:,:,k_bou) MG_WL8msFC(:,:,k_bou)]*1e6);
        set(phase,'cdata',[MG_Transient(:,:,k_bou)*1e6]);% ESO_WL2msHP(:,:,k_bou)
        title(ax,sprintf('n %4d',k_bou)); drawnow;
        %pause
    end
end