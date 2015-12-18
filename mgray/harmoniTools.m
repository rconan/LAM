classdef harmoniTools
    
    methods (Static)
        
%% Anneau sur une carte en pixel de cot? (nRes) et ayant pour :
%   r?solution d'un pixel en m?tres (Res_1pix)
%   diam?tre int?rieur (diamInt) ; diam?tre ext?rieur (diamExt)
        function masq_Anneau = ANNEAU(nRes,Res_1pix,diamInt,diamExt)
            x = ((1:nRes) - (1+nRes)/2)*Res_1pix;
            [X,Y] = meshgrid(x);
            R = sqrt(X.^2+Y.^2);
            masq_Anneau = false(size(X));
            masq_Anneau(R(:) > diamExt/2 | R(:) < diamInt/2) = true;
            masq_Anneau = ~masq_Anneau;
        end
        
%% Carte 2D ? partir :
%   d'un vecteur des valeurs non nulles (Carte_1D)
%   du masque logique (masque_2D)
        function Carte_2D = CARTES_2D(Carte_1D,masque_2D)
            Carte_2D = zeros(size(masque_2D,1));
            Carte_2D(masque_2D == 1) = Carte_1D;
        end

%% Affichage de la carte de phase d'un Zernike
        function Zern_2D = CARTES_2DZER(Zernik,k,DiamPix,masque)
            Zern_2D = zeros(DiamPix);
            zer = Zernik.p(:,k);
            Zern_2D(masque == 1) = zer(masque);
%             figure(k)
%             imagesc(Zern_2D);
%             axis([1 DiamPix 1 DiamPix],'square'); caxis([-3 3]); ylabel(colorbar,'');
%             display(['Zernike n?',num2str(k),': Max = ',num2str(max(zer(masque))),' ; Min = ',num2str(min(zer(masque)))]);
%             saveas(gcf,strcat('ZER_n?',num2str(k)),'jpg');
        end
  
%% Inverse G?n?ralis?e avec suppression de Valeurs Singuli?res      
        function [M_IG,VS,Cdt] = MAT_INVGEN(M,nb_VS,Type)
            [UInt,SInt,VInt] = svd(full(M));
            VS = diag(SInt);
            iSInt = diag(1./VS(1:end-nb_VS));
            Cdt = VS(1)/VS(end-nb_VS);
            if Type == 1
                [nL,nC] = size(M);
                iSInt(nC,nL) = 0;
            elseif Type == 0
                if nb_VS > 0
                    [nL,nC] = size(M);
                    iSInt(nC,nL) = 0;
                end
            end
            M_IG = VInt*iSInt*UInt';
        end

%% Araign?e  de l' E-ELT sur une carte de cot? (nRes) ayant pour :
%   r?solution d'un pixel en m?tres (Res_1pix) 
%   diam?tre int?rieur (diamInt) ; diam?tre ext?rieur (diamExt) ;
% 6 Bras index?s suivant le sens des aiguilles d'une montre
        function [Masq_EELTLegs,Masq_EELTSpider] = EELT_SPIDER(nRes,Res_1pix,Alpha,largBras,diamInt,diamExt)% LargBras = 0.5;
            Theta = pi/3;
            Theta0 = pi/6+Alpha;
            x = ((1:nRes) - (1+nRes)/2) * Res_1pix;
            [X,Y] = meshgrid(x);
            Spider = false(numel(X,1),1);
            Legs = zeros(nRes*nRes,6);
            for leg = 1:6
                phi = (leg-1)*Theta + Theta0;
                P = diamInt/2 * [cos(phi), sin(phi)];
                par = [cos(phi); sin(phi)];
                per = [-sin(phi); cos(phi)];
                in = [X(:),Y(:)] * par - P * par > 0;
                thickness = abs([X(:),Y(:)] * per - P * per) - largBras/2 < 0;
                Legs(:,leg) = in.*thickness;
                Spider = Spider | (in.*thickness);
            end
            R = sqrt(X.^2+Y.^2);
            Masq_EELTSpider = false(size(X));
            Masq_EELTSpider(logical(Spider)) = true;
            Masq_EELTSpider(R(:) > diamExt/2 | R(:) < diamInt/2) = true;
            Masq_EELTSpider = ~Masq_EELTSpider;
            Masq_EELTLegs = zeros(nRes,nRes,6);
            Masq_Pup = ~reshape((R(:) > diamExt/2 | R(:) < diamInt/2),nRes,nRes);
            for leg = 1:6
                Masq_EELTLegs(:,:,leg) = reshape(Legs(:,leg),nRes,nRes).*Masq_Pup;
            end
            Masq_EELTLegs = logical(Masq_EELTLegs);
        end

%% 6 P?tales  du miroir de l'E-ELT sur une carte de cot? (nRes)
%   ayant pour r?solution d'un pixel en m?tres (Res_1pix) 
%   index?s suivant le sens des aiguilles d'une montre
%   la concat?nation des 6 p?tales -> le masque de l'araign?e !
        function Masq_EELTPetales = EELT_PETALES(Res_1pix,Alpha,largBras,MasqPupLog)
            nRes = size(MasqPupLog,1);
            Masq_EELTPetales = zeros(nRes,nRes,6);
            m1 = tan(-pi/6-Alpha); p1 = nRes/2 + largBras/2/Res_1pix/cos(-pi/6-Alpha);
            m2 = tan(pi/6-Alpha); p2 = nRes/2 - largBras/2/Res_1pix/cos(pi/6-Alpha);
            p3 = nRes/2 + largBras/2/Res_1pix/cos(pi/6-Alpha);
            u = 1:nRes;
            X = repmat(u,nRes,1);
            Y = repmat(flipud(u'),1,nRes);
            fn_Test = @(x,y,m,p) (y-m.*(x-nRes/2)-p);
            Test1 = fn_Test(nRes/2,nRes/2,m1,p1);
            Test2 = fn_Test(nRes/2,nRes/2,m2,p2);
            Test3 = fn_Test(nRes/2,nRes/2,m2,p3);
            Petale = logical(fn_Test(X,Y,m1,p1).*(Test1*ones(nRes))<0) & logical(fn_Test(X,Y,m2,p2).*(Test2*ones(nRes))<0) & MasqPupLog;
            Masq_EELTPetales(:,:,1) = Petale;
            Masq_EELTPetales(:,:,4) = rot90(Petale,2);
            if Alpha == 0
                Petale = logical(X>nRes/2+largBras/2/Res_1pix) & logical(fn_Test(X,Y,m2,p3).*(Test3*ones(nRes))<0) & MasqPupLog;
                Masq_EELTPetales(:,:,6) = Petale;
                Masq_EELTPetales(:,:,3) = rot90(Petale,2);
                Masq_EELTPetales(:,:,5) = fliplr(Petale);
                Masq_EELTPetales(:,:,2) = flipud(Petale);
            elseif Alpha > 0
                m4 = tan(pi/2-Alpha); p4 = nRes/2 - largBras/2/Res_1pix/cos(pi/2-Alpha);
                Test4 = fn_Test(nRes/2,nRes/2,m4,p4);
                Petale = logical(fn_Test(X,Y,m4,p4).*(Test4*ones(nRes))<0) & logical(fn_Test(X,Y,m2,p3).*(Test3*ones(nRes))<0) & MasqPupLog;
                Masq_EELTPetales(:,:,6) = Petale;
                Masq_EELTPetales(:,:,3) = rot90(Petale,2);
                p5 = nRes/2 - largBras/2/Res_1pix/cos(-pi/6-Alpha);
                Test5 = fn_Test(nRes/2,nRes/2,m1,p5);
                Petale = logical(fn_Test(X,Y,m4,p4).*(Test4*ones(nRes))<0) & logical(fn_Test(X,Y,m1,p5).*(Test5*ones(nRes))<0) & MasqPupLog;
                Masq_EELTPetales(:,:,2) = Petale;
                Masq_EELTPetales(:,:,5) = rot90(Petale,2);
            end
            Masq_EELTPetales = logical(Masq_EELTPetales);
        end

%% Matrices des Gradients suivant (OX) & (OY) sur les bras de l'araign?e
        function LegsGradX = EELT_LEGSGradX(MasqLegs,MasqPupLog)
            nRes = size(MasqPupLog,1);
            D = zeros(nRes); nb_lig = 0;
            for ik = 1:nRes
                if nnz(MasqLegs(ik,:)) > 1
                    nb_lig = nb_lig +1;
                    jg = find(MasqLegs(ik,:),1,'first');
                    jd = find(MasqLegs(ik,:),1,'last');
                    D(ik,jg) = -1/(jd-jg);
                    D(ik,jd) = 1/(jd-jg);
                end
            end
            D = D/nb_lig;% normalisation ?
            LegsGradX = (D(MasqPupLog))';
        end
        function LegsGradY = EELT_LEGSGradY(MasqLegs,MasqPupLog)
            nRes = size(MasqPupLog,1);
            D = zeros(nRes); nb_col = 0;
            for jk = 1:nRes
                if nnz(MasqLegs(:,jk)) > 1
                    nb_col = nb_col +1;
                    ih = find(MasqLegs(:,jk),1,'first');
                    ib = find(MasqLegs(:,jk),1,'last');
                    D(ih,jk) = -1/(ib-ih);
                    D(ib,jk) = 1/(ib-ih);
                end
            end
            D = D/nb_col;% normalisation ?
            LegsGradY = (D(MasqPupLog))';
        end

%% Miroir Segment?  de l' E-ELT :
%   sur une carte de cot? nResUV correspondant au M1 UnVignetted
%   ayant pour r?solution d'un pixel en m?tres (Res_1pix) 
        function Masq_EELTSegment = EELT_SEGMENTS(nResM1UV,Res_1pix)
            D_hexa = 1.41849;
            R_hexa = D_hexa/2;
            CentresM = importdata('C:\Users\mgray\AO_HARMONI\MIROIR_M1\VERSION_1_3\SegmentVertices.txt');
            CentresM(:,3) = [];
            x = ((1:nResM1UV) - (1+nResM1UV)/2) * Res_1pix;
            [X,Y] = meshgrid(x);
            PixelPositions = [X(:),Y(:)];
            R2_Pixels = sum(PixelPositions.^2,2);
            R2_Vertices = sum(CentresM.^2,2);
            R2_max = (sqrt(max(R2_Vertices))+(R_hexa))^2;
            in = find(R2_Pixels <= R2_max);
            xy = (PixelPositions(in,:));
            Angles =  pi/3*(0:2)';
            XY2UV = 4/3*[cos(Angles),sin(Angles)]/D_hexa;
            uvs = (XY2UV * CentresM')';
            uvsr = round(uvs);
            ReconstructedRadius = uvsr*pinv(XY2UV)';
            ReconstructedRadius = sqrt(sum(ReconstructedRadius.^2,2));
            order = 3;
            radius = sqrt(sum(CentresM.^2,2));
            m = zeros(numel(radius),order+1);
            for u = 0:order
                m(:,u+1) = radius.^u;
            end
            im = (m'*m)\m';
            c = im*ReconstructedRadius;
            radius = sqrt(sum(xy.^2,2));
            m = zeros(numel(radius),order+1);
            for u = 0:order
                m(:,u+1) = radius.^u;
            end
            xy = xy.*repmat((m*c)./radius,1,size(xy,2));
            uv = xy*XY2UV';
            uvr = round(uv);
            frac = abs(uvr - uv);
            [~,rejected] = max(frac,[],2);
            loc = zeros(size(uv,1),1);
            for k = 1:3
                rej = find(rejected==k);
                range = [mod(k,3)+1,mod(k+1,3)+1];
                [~,loc(rej)] = ismember(uvr(rej,range),uvsr(:,range),'rows');
            end
            SegmentNumber = zeros(size(PixelPositions,1),1);
            SegmentNumber(in) = loc;
            Masq_EELTSegment = logical(reshape(SegmentNumber,nResM1UV,nResM1UV));
        end

%%

    end   
        
end
