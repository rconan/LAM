clear all; clc; close all;
REP = 'C:\Users\mgray\AO_HARMONI\MIROIR_M4\VERSION_1_0\';
% REP = %'/home/mgray/HARMONI/MIROIR_M4/VERSION_1_0/';%
ROUT = 3;

%% VISUALISATION  des  FONCTIONS D' INFLUENCE

if ROUT == 1
    nRes = 572;
    FICH = ['MAT_FI_M4_' num2str(nRes)]; load([REP FICH])
    for k_act = 20:50
        CarteFI2D = reshape(MatFI_M4(:,k_act),nRes,nRes);
        figure(1)
        imagesc(CarteFI2D)
        axis square; colorbar; pause
        saveas(gcf,strcat('ACT_',num2str(k_act)),'jpg');
    end
   
%% INTERPOLATION
    
elseif ROUT == 2
    nRes1 = 572; nRes2 = 740;
    FICH1 = ['MAT_FI_M4_' num2str(nRes1)]; load([REP FICH1])
    nb_act = size(MatFI_M4,2);
    [X1,Y1] = meshgrid(1:nRes1,1:nRes1);
    echant = linspace(1,nRes1,nRes2);
    [X2,Y2] = meshgrid(echant,echant);
    MatFIM4 = zeros(nRes2*nRes2,nb_act);
    for k_act = 1:nb_act;
        CarteFIres1 = reshape(MatFI_M4(:,k_act),nRes1,nRes1);
        CarteFIres2 = interp2(X1,Y1,CarteFIres1,X2,Y2,'spline');
        MatFIM4(:,k_act) = CarteFIres2(:);
        display(k_act)
    end
    MatFI_M4 = single(MatFIM4); clear MatFIM4
    save(['MAT_FI_M4_' num2str(nRes2)],'MatFI_M4','-v7.3');

%% POSITIONS des ACTIONNEURS ( considérés comme les centres d'hexagones )
    
% NE DEPEND QUE du RAYON du cercle circonscrit à l'hexagone
% CLASSEMENT des ACTIONNEURS par lignes (haut -> bas) puis par colonnes (gauche -> droite)
% REPERE : origine -> centre du Miroir ; unités -> mètres
    
%       D'après les docs fournies par l'ESO :
% - le miroir M4 a un rayon "moyen" (2.347+2.386)/2/2 = 1.18325 mètres
% - la distance (en mètres) entre 2 actionneurs dans la répartition triangulaire est de 31.5e-3 = 2*(sqrt(3)/2)*R_Hexa 
%               => R_Hexa (en mètres) = 31.5e-3/sqrt(3) = 18.187e-3   pour le M4

elseif ROUT == 3
    % Paramètres & Initialisation
    R_hexa = 31.5e-3/sqrt(3);
    R_MirExt = (2.347+2.386)/4;
    %R_MirInt = 0.54/2;% -> obstruction centrale du M4
    % ou :
    R_MirInt = 0.3037*R_MirExt;% -> obstruction centrale du M1 AllGlass
    nb_cour = 60;
    fn_NbHexa = @(x) (3*x*(x+1)+1);% nb total d'hexagones pour ( x couronnes + 1 hexagone central )
    nb_act = fn_NbHexa(nb_cour);
    CentresM = zeros(nb_act,3);% 2 coordonnées (x,y) des centres + indice de couronne
    CentresM(1,:) = [0 0 0];
    % Symétries axiales
    Alpha = zeros(11,1);
    for k = 1:11, Alpha(k) = exp(1i*k*pi/3);
    end
    % Calcul des coord des centres d'hexagones
    for k_cour = 1:nb_cour
        % Sur la partie du plan (XOX1)
        if mod(k_cour,2) == 0
            p_cou = k_cour/2;
            x = zeros(p_cou+1,1); y = zeros(p_cou+1,1);
            for k = 1:p_cou+1
                x(k) = 3*p_cou*R_hexa;
                y(k) = (k-1)*sqrt(3)*R_hexa;
            end
        elseif mod(k_cour,2) == 1
            p_cou = (k_cour-1)/2;
            x = zeros(p_cou+1,1); y = zeros(p_cou+1,1);
            for k = 1:p_cou+1
                x(k) = 3*k_cour*R_hexa/2;
                y(k) = (2*k-1)*sqrt(3)/2*R_hexa;
            end
        end
        Z0 = x + 1i*y;
        % Sur la totalité du plan
        Z = zeros(length(Z0),12);
        Z(:,1) = Z0;
        for k = 1:11
            Z(:,k+1) = Alpha(k)*conj(Z(:,k));
        end
        Z_C = Z(:);
        % Suppression des doublons
        Centres = round([real(Z_C),imag(Z_C),k_cour*ones(length(Z_C),1)],8);
        CentresM(fn_NbHexa(k_cour-1)+1:fn_NbHexa(k_cour),:) = unique(Centres,'rows');
    end
    % Classement
    CentresM = sortrows(CentresM);% coordonnées (x,y) horizontales
    IndCour = CentresM(:,3);
    CentresM(:,3) = [];
    % Suppression des Hexagones Extérieurs au miroir & Intérieurs à l'obstruction centrale
    ind = [];
    for k_hexa = 1:nb_act
        if sqrt(CentresM(k_hexa,1)^2+CentresM(k_hexa,2)^2) > R_MirExt
            ind = [ind k_hexa];
        end
        if sqrt(CentresM(k_hexa,1)^2+CentresM(k_hexa,2)^2) < R_MirInt
            ind = [ind k_hexa];
        end
    end
    CentresM(ind,:) = [];
    IndCour(ind) = [];
    nb_act = length(IndCour);
    save(['MG_M4Act_obs=' num2str(round(R_MirInt/R_MirExt*100)) '.mat'],'CentresM','IndCour','nb_act');
end

