clear all; clc; close all;
REP = '/home/mgray/OA_HARMONI/PERT_DYNAMIQUES/VERSION_1_3/';
nRes = 1014; nb_boucle = 1100; Typ = 'S';
ROUT = 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%      WIND LOAD 8m/s       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ROUT == 1 || ROUT == 2
    diamExtUV = 38.542;
    % ZERNIKES sur le M1
    telM1 = telescope(diamExtUV,'obstructionRatio',0,'resolution',nRes,'samplingTime',1/500);
    masq_telM1 = telM1.pupilLogical;
    ZernM1 = zernike(1:8,diamExtUV,'resolution',nRes,'pupil',masq_telM1);
    TipM1 = harmoniTools.CARTES_2DZER(ZernM1,2,nRes,masq_telM1);
    TiltM1 = harmoniTools.CARTES_2DZER(ZernM1,3,nRes,masq_telM1);
    FocusM1 = harmoniTools.CARTES_2DZER(ZernM1,4,nRes,masq_telM1);
    Ast0M1 = harmoniTools.CARTES_2DZER(ZernM1,5,nRes,masq_telM1);
    Ast45M1 = harmoniTools.CARTES_2DZER(ZernM1,6,nRes,masq_telM1);
    ComaXM1 = harmoniTools.CARTES_2DZER(ZernM1,7,nRes,masq_telM1);
    ComaYM1 = harmoniTools.CARTES_2DZER(ZernM1,8,nRes,masq_telM1);
end
    % PERTURBATIONS ( Tip + Tilt )
if ROUT == 1
    FICH = 'MSwind8mps_tiptilt.mat';% en [Arcsec] & 1 mas <-> 46.7 nm RMS
    load([REP FICH]);
    MG_WL8msTT = zeros(nRes,nRes,nb_boucle);
    for k_bou = 1:nb_boucle
        k_bou
        MG_WL8msTT(:,:,k_bou) = 46.7e-6*(MSwindtiptiltData(k_bou,2)*TipM1 - MSwindtiptiltData(k_bou,3)*TiltM1);
    end
    MG_WL8msTT = single(MG_WL8msTT);
    save(['MG_WL8msTT_' num2str(nRes) Typ],'MG_WL8msTT','-v7.3');
    % PERTURBATIONS ( Focus -> ComaY )
elseif ROUT == 2
    FICH = 'MSwind8mps_ZernikeNoll.mat';% en [m]
    load([REP FICH]);
    MG_WL8msFC = zeros(nRes,nRes,nb_boucle);
    for k_bou = 1:nb_boucle
        k_bou
        MG_WL8msFC(:,:,k_bou) = MSwindDataNoll(k_bou,2)*FocusM1 - MSwindDataNoll(k_bou,3)*Ast0M1 + MSwindDataNoll(k_bou,4)*Ast45M1 - MSwindDataNoll(k_bou,5)*ComaXM1 + MSwindDataNoll(k_bou,6)*ComaYM1;
    end
    MG_WL8msFC = single(MG_WL8msFC);
    save(['MG_WL8msFC_' num2str(nRes) Typ],'MG_WL8msFC','-v7.3');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%      WIND LOAD 2m/s       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ROUT == 3
    REP_M1 = '/home/mgray/OA_HARMONI/MIROIR_M1/VERSION_1_3/';
    FICH_M1 = ['MG_PupilM1_' num2str(nRes)]; load([REP_M1 FICH_M1]);
    % ZERNIKES sur 1 SEGMENT
    telSeg = telescope(D_hexaM,'obstructionRatio',0,'resolution',lg_DhexaPix,'samplingTime',1/500);
    masq_hexa = telSeg.pupilLogical;
    ZernSeg = zernike(1:8,D_hexaM,'resolution',lg_DhexaPix,'pupil',masq_hexa);
    PistonSeg = harmoniTools.CARTES_2DZER(ZernSeg,1,lg_DhexaPix,masq_hexa);
    TipSeg = harmoniTools.CARTES_2DZER(ZernSeg,2,lg_DhexaPix,masq_hexa);
    TiltSeg = harmoniTools.CARTES_2DZER(ZernSeg,3,lg_DhexaPix,masq_hexa);
    npix = nnz(PistonSeg);
    CenZernPix = [lg_DhexaPix/2+0.5 ; lg_DhexaPix/2+0.5];
    % BACK STRUCTURE / HARDPACT / SOFTPACT sur le M1 (Piston ; Tip ; Tilt)
    fn_Trans = @(pt,u) (u+pt);% translation du point pt suivant la direction du vecteur u
    FICH = 'M1wind2mps_backstructure_Optical.mat';% en [m]
    % FICH = 'M1wind2mps_hardpact_Optical.mat';% en [m]
    % FICH = 'M1wind2mps_softpact_Optical.mat';% en [m]
    load([REP FICH]);
    MG_WL2msBS = zeros(nRes,nRes,nb_boucle);
    for k_bou = 1:nb_boucle
        k_bou
        Pert_Pup = zeros(nRes);
        for k_seg = 1:nb_hexa
            Pert_Seg = zeros(nRes);
            Pert_PTT = PTTtimeDataOptical(k_bou,3*k_seg-1)*PistonSeg + PTTtimeDataOptical(k_bou,3*k_seg)*TipSeg +...
                - PTTtimeDataOptical(k_bou,3*k_seg+1)*TiltSeg;
            [iP,jP,sP] = find(Pert_PTT);
            xy_ptsPert = [jP iP]';% Calculs vectoriels en coordonnées (x,y) verticales !
            Vect_u = CentresPix(:,k_seg) - CenZernPix;
            indMasq = round(fn_Trans(xy_ptsPert,repmat(Vect_u,1,npix)));
            for k_ind = 1:npix
                Pert_Seg(indMasq(2,k_ind),indMasq(1,k_ind)) = sP(k_ind);% Indexation en (i,j) !
            end
            Pert_Seg = Pert_Seg.*Masq_SEGM(:,:,k_seg);
            Pert_Pup = Pert_Pup + Pert_Seg;
        end
        MG_WL2msBS(:,:,k_bou) = Pert_Pup;
    end
    MG_WL2msBS = single(MG_WL2msBS);
    save(['MG_WL2msBS_' num2str(nRes) Typ],'MG_WL2msBS','-v7.3');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%         TRANSIENT         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ROUT == 4
    diamExtUV = 38.542;
    % ZERNIKES sur le M1
    telM1 = telescope(diamExtUV,'obstructionRatio',0,'resolution',nRes,'samplingTime',1/500);
    masq_telM1 = telM1.pupilLogical;
    ZernM1 = zernike(1:3,diamExtUV,'resolution',nRes,'pupil',masq_telM1);
    TiltM1 = harmoniTools.CARTES_2DZER(ZernM1,3,nRes,masq_telM1);
    % PERTURBATION ( Tilt )
    FICH = 'TransientTilt.mat';% en [Arcsec] & 1 mas <-> 46.7 nm RMS
    load([REP FICH]);
    MG_Transient = zeros(nRes,nRes,nb_boucle);
    for k_bou = 1:nb_boucle
        k_bou
        MG_Transient(:,:,k_bou) = - 46.7e-6*TransientData(k_bou+1100,2)*TiltM1;
    end
    MG_Transient = single(MG_Transient);
    save(['MG_Transient_' num2str(nRes) Typ],'MG_Transient','-v7.3');
end

