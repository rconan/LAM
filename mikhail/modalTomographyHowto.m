%% ADAPTIVE OPTICS TOMOGRAPHY HOWTO
% Demonstrate how to build a tomographic adaptive optics system

%% Atmosphere 
% atm = atmosphere(photometry.V,0.15,30,...
%     'altitude',[0,4,6,10]*1e3,...
%     'fractionnalR0',[0.6,0.25,0.1,0.05],...
%     'windSpeed',[5,10,5,20],...
%     'windDirection',[0,pi/4,0,pi]);

atm = atmosphere(photometry.V,0.15,1000,...
    'altitude',[0,4]*1e3,...
    'fractionnalR0',[1e-10 1-1e-10],...
    'windSpeed',[5,10],...
    'windDirection',[0,pi/4]);
atm.wavelength = photometry.R;

%% Telescope
nPx = 60;%60;
tel = telescope(3.6,...
    'fieldOfViewInArcMin',4,...
    'resolution',nPx,...
    'samplingTime',1/100);

%% Sources
ngs = source;
% ast = source('asterism',{[3,30*cougarConstants.arcsec2radian,0]});

%% Wavefront sensor
nLenslet = 10;
wfs = shackHartmann(nLenslet,nPx,0.75);
%setValidLenslet(wfs,utilities.piston(nPx))
ngs = ngs.*tel*wfs;
    wfs.INIT;
    wfs.camera.photonNoise = false;
    +wfs;

%% Deformable mirror
nActuator = nLenslet + 1;
bif = influenceFunction('monotonic',25/100);
dm = deformableMirror(nActuator,...
    'modes',bif,...
    'resolution',nPx,...
    'validActuator',wfs.validActuator,...
    'zLocation',0e3);

%% Building the system
ngs=ngs.*tel*dm*wfs;
%wfs.referenceSlopes = wfs.slopes;
+ngs;
figure
subplot(1,2,1)
imagesc(wfs.camera)
subplot(1,2,2)
slopesDisplay(wfs)

  %% DM/WFS calibration
% dm.coefsDefault = 0;
% stroke = ngs.wavelength/2;
% dm.coefs = eye(dm.nValidActuator)*stroke;
% +ngs;
% dm.coefs = 0;
% calibrationMatrix = wfs.slopes./stroke;
% figure(10)
% subplot(1,2,1)
% imagesc(calibrationMatrix)
% xlabel('DM actuators')
% ylabel('WFS slopes [px]')
% ylabel(colorbar,'slopes/actuator stroke')
% 
% %% Command matrix derivation
% [nS,nC] = size(calibrationMatrix);
% [U,S,V] = svd(calibrationMatrix);
% nThresholded = 4;
% eigenValues = diag(S);
% subplot(1,2,2)
% semilogy(eigenValues,'.')
% xlabel('Eigen modes')
% ylabel('Eigen values')
% iS = diag(1./eigenValues(1:end-nThresholded));
% iS(nC,nS) = 0;
% commandMatrix = V*iS*U';

%% Zernike measurement
maxRadialDegree = 9;
zern = zernike(2:zernike.nModeFromRadialOrder(maxRadialDegree),'resolution',nPx,'pupil',tel.pupil);
zern.lex = false;
% figure(10)
% imagesc(zern.phase)
zern.c = eye(zern.nMode);
ngs=ngs.*zern*wfs;
z = zernike(1:zernike.nModeFromRadialOrder(maxRadialDegree))\wfs;
Dz = z.c;%(2:end,:);

% %% With noise
% ngs.wavelength = photometry.R;
% ngs.magnitude = 0;
% wfs.camera.readOutNoise = 5;
% wfs.camera.photonNoise = true;
% wfs.framePixelThreshold = 0;
% ngs=ngs.*tel*wfs;
% 
% %% noise convariance matrix
% nMeas = 250;
% slopes = zeros(wfs.nSlope,nMeas);
% parfor kMeas=1:nMeas
%     +wfs
%     slopes(:,kMeas) = wfs.slopes;
% end
% Cn = slopes*slopes'/nMeas;
% figure(5)
% subplot(1,2,1)
% imagesc(Cn)
% axis equal tight
% colorbar
% wfs.slopes = slopes;
% z = z\wfs;
% Czn = z.c*z.c'/nMeas;
% subplot(1,2,2)
% imagesc(Czn)
% axis equal tight
% colorbar


%% TOMOGRAPHY
tel = tel+atm;
 %% Tomography
ast = source('asterism',{[3,30*cougarConstants.arcsec2radian,0]},...
    'wavelength',photometry.R,'magnitude',0)
ast = source('asterism',{[2,30*cougarConstants.arcsec2radian,0]},...
    'wavelength',photometry.R,'magnitude',0)
rr = [0 10 20 30];
th = [0 0 0 0]
rr = [120 120 120];
th = [0 0 pi/2];
ast = source('zenith',rr*constants.arcsec2radian,'azimuth',th); 
% ast = source('asterism',{[3,4*30*cougarConstants.arcsec2radian,0]},...
%     'wavelength',photometry.R,'magnitude',0)
src = source('asterism',{[0,0],[5,60*constants.arcsec2radian,0]},'wavelength',photometry.H,'magnitude',[8 10 12 9 11 14])
pd = source('zenith',0*cougarConstants.arcsec2radian,'wavelength',photometry.R);
figure(3)
imagesc(tel,[ast,pd])
wpd = ones(size(pd));
gs = ast;
nGs = length(ast);
nPd = length(pd);
zernModeMax = zernike.nModeFromRadialOrder(maxRadialDegree);
wfsMaxRadialDegree = zernModeMax;
zernWfs = zernike(2:zernModeMax,'resolution',nPx,'pupil',tel.pupil);
%tel = tel+atm;
% Zernike section expansion
%% Projection
nDm = length(dm);
P = cell(nPd,nDm);
for kDm = 1:nDm
    fprintf('@(Projection)> ')
    for kPd = 1:nPd
        fprintf('pd#%d/dm#%d - ',kPd,kDm)
        src = pd(kPd);
        delta = dm.zLocation.*tan(src.zenith).*...
            [cos(src.azimuth),sin(src.azimuth)];
        delta = delta*2/tel.diameterAt(dm.zLocation);
        alpha = tel.diameterAt(dm.zLocation)./tel.D;
        P{kPd,kDm} = smallFootprintExpansion(zernWfs,delta,alpha);
    end
    fprintf('\n')
end
%% Turbulence covariance matrix
tic
SFile = sprintf('S%d__18-Jan-2010.mat',maxRadialDegree);
if exist(SFile,'file')==2
    load(SFile)
else
    S = phaseStats.zernikeAngularCovariance(zernWfs,atm,ast);
    SFile = sprintf(['S%d__',date,'.mat'],maxRadialDegree);
    save(SFile,'S','maxRadialDegree')
    fprintf(' S matrix saved to %s\n',SFile)
end
toc

%--------------------------------------------------------------------------
%% Tomographic phase estimation using explicit reconstruction
tic;
nL = atm.nLayer;
Ptheta = cell(nL,nGs);
altitudes = [atm.layer.altitude];
for kLayer = 1:nL
    fprintf('@(Projection)> ')
    parfor kPd = 1:nGs
        fprintf('pd#%d/dm#%d - ',kPd,kLayer)
        src = gs(kPd);
        delta = altitudes(kLayer).*tan(src.zenith).*...
            [cos(src.azimuth),sin(src.azimuth)];
        delta = delta*2/tel.diameterAt(altitudes(kLayer));
        alpha = tel.diameterAt(altitudes(kLayer))./tel.D;
        Ptheta{kPd,kLayer} = smallFootprintExpansion(zernWfs,delta,alpha);
    end
    fprintf('\n')
end
Ptheta = cell2mat(Ptheta);

%AnalyticalSmallFootprintExpansion
[Palpha Palphacell] = AnalyticalSmallFootprintExpansion(zernWfs,tel,ast,atm);


Slayered = cell(nL,nL);
fr0 = [atm.layer.fractionnalR0];
altitudes = [atm.layer.altitude];
ws = [atm.layer.windSpeed];
wd = [atm.layer.windDirection];

%myatm = atm;
for kLayer = 1:nL
    for kLayer1 = 1:nL
        if kLayer == kLayer1
%             myatm = atmosphere(photometry.V,0.15,30,...
%                 'altitude',altitudes(kLayer),...
%                 'fractionnalR0',1,...
%                 'windSpeed',ws(kLayer),...
%                 'windDirection',wd(kLayer));
%             atm.wavelength = photometry.R;

            zernWfsi = zernike(2:zernModeMax,'resolution',nPx,'pupil',tel.pupil);
            zernWfsi.D = tel.diameterAt(altitudes(kLayer))/tel.D*2;
            Slayered{kLayer,kLayer1} = phaseStats.zernikeCovariance(zernWfsi,atm)*fr0(kLayer);
        else
            Slayered{kLayer,kLayer1} = zeros(zernWfs.nMode);
        end
    end
end
Slayered = cell2mat(Slayered);
idy = [];
for kGs=1:nL
    ii = (2:zernWfs.nMode+1)-1;
    idy = [idy ii+zernWfs.nMode*(kGs-1)+kGs];
end
idx = [];
for kGs=1:nGs
    ii = (2:zernWfs.nMode+1)-1;
    idx = [idx ii+zernWfs.nMode*(kGs-1)+kGs];
end
%Sigma_theta = Ptheta(idx,:)*Slayered*Ptheta(idx,:)';
Sigma_theta = Palpha(idx,idx)*Slayered*Palpha(idx,idx)';
toc
Smat = cell2mat(S);

%% figures for the theta directions
for kGs=1:nGs
    subplot(1,nGs,kGs)
    semilogy(diag(Sigma_theta,zernWfs.nMode*(kGs-1)),'k');
    hold
    semilogy(diag(Smat,zernWfs.nMode*(kGs-1)),'r--');
    legend('MV','L&A')
    title(['Diagonals of #1 with #' num2str(kGs)])
end

figure
subplot(1,2,1)
imagesc(Sigma_theta)
title('MV')
subplot(1,2,2)
imagesc(Smat)
title('L&A')

mymat = Sigma_theta(1:nMode,nMode+1:end);
figure
loglog(svd(Sigma_theta),'k')
hold
loglog(svd(mymat(1:44,1:44)),'k')
loglog(svd(Smat(1:44,1:44)),'r');
title('Singular value decomposition')
    legend('MV','L&A')
axis tight
    
Palpha = cell(nPd,nL);
parfor kLayer = 1:nL
    fprintf('@(Projection)> ')
    for kPd = 1:nPd
        fprintf('pd#%d/dm#%d - ',kPd,kLayer)
        src = pd(kPd);
        delta = altitudes(kLayer).*tan(src.zenith).*...
            [cos(src.azimuth),sin(src.azimuth)];
        delta = delta*2/tel.diameterAt(altitudes(kLayer));
        alpha = tel.diameterAt(altitudes(kLayer))./tel.D;
        Palpha{kPd,kLayer} = smallFootprintExpansion(zernWfs,delta,alpha);
    end
    fprintf('\n')
end

Palpha = cell2mat(Palpha);
Palpha = Palpha(2:end,:);


Sigma_thetaalpha = Palpha*Slayered*Ptheta(idx,:)';

Rmv = (Sigma_thetaalpha* DzAst')/(DzAst*Sigma_theta*DzAst' + CznAst);
fprintf('PROJECTIONS DONE')
%--------------------------------------------------------------------------


%% wavefront reconstruction least square fit
ast=ast.*tel;
ps = [ast.meanRmPhase];
ast=ast*wfs
z = z\wfs;
zern.c = Dz\z.c;%(2:end,:);
zern.lex = true;
phaseLS = reshape(zern.p*zern.c,nPx,nPx*length(ast));

%% wavefront reconstruction minimum variance
%S = cell2mat(S);
CznAst = blkdiag( Czn , Czn , Czn );
DzAst = blkdiag( Dz , Dz , Dz  );
%M = S/(S+CznAst);
M = S*DzAst'/(DzAst*S*DzAst'+CznAst);
z = z - 1; % piston removed
zern.c = reshape(M*z.c(:),z.nMode,[]);
phaseMV = reshape(zern.p*zern.c,nPx,nPx*length(ast));
figure(12)
imagesc([ps;phaseLS;phaseMV])
axis equal tight xy

%% Data/Target covariance
CFile = sprintf('C%d__31-Mar-2010.mat',maxRadialDegree);
if exist(CFile,'file')==2
    load(CFile)
else
    C = phaseStats.zernikeAngularCovariance(zernWfs,atm,gs,pd);
    CFile = sprintf(['C%d__',date,'.mat'],maxRadialDegree);
    save(CFile,'C','maxRadialDegree')
    fprintf(' C matrix saved to %s\n',CFile)
end
%% Target matrix
T = 0;
R = 0;
for kPd = 1:nPd
    PnDM = cell2mat(P(kPd,:));
    T = T + wpd(kPd)*DzAst*cell2mat(C(:,kPd))*PnDM(2:end,:);
    R = R + wpd(kPd)*(PnDM'*PnDM);
end

%% Command matrix
M = R\T'/(DzAst*S*DzAst'+CznAst);
ast=ast.*tel*wfs;
z = zernike(1:zernike.nModeFromRadialOrder(maxRadialDegree))\wfs;
z = z - 1;
zern.c = M*z.c(:);
pd = pd.*tel;
figure
imagesc([reshape(zern.p*zern.c,nPx,nPx),pd.meanRmPhase])
axis equal tight xy
colorbar
