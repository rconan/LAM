clear
ngs = source;
nPx = 40;
tel = telescope(1,'resolution',nPx);
zern = zernike(tel,3);
zern.c = ngs.wavelength/4/16;
pyr = pyramid(nPx,'modulation',0,'binning',1);
%pyr.modulation = 1;
%pyr.binning = 2;
%%
ngs.*tel*pyr;
pyr.INIT
ngs.*tel*zern*pyr;

figure(101)
subplot(2,1,1)
imagesc(pyr.camera.frame)
axis equal tight
subplot(2,1,2)
imagesc(pyr.slopesMap.*pyr.validSlopes)
axis equal tight

%%
ast = source('asterism',{[3,arcsec(30),0]});
ast.*tel*zern*pyr

figure(102)
subplot(2,1,1)
imagesc(pyr.camera.frame)
axis equal tight
subplot(2,1,2)
imagesc(pyr.slopesMap)
axis equal tight

%%
bif = influenceFunction('monotonic',0.5);
dm = deformableMirror(20,'modes',bif,'resolution',nPx,...
    'validActuator',tools.piston(20,'type','logical'));
dm.coefs = eye(dm.nValidActuator)*ngs.wavelength/10;

ngs.*tel*dm;
calib = calibration(dm,pyr,ngs,ngs.wavelength/10,100);

