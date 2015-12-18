function wfs = tool_calibrate_wfs_oomao(nLenslet,nPx,minLightRatio , tel, d, photometry_lambda)

ngs = source; %% make a calibration source
eval(strcat('ngs.wavelength = photometry.', photometry_lambda, ';'));

nRes = nLenslet*nPx; %% total number of pixels.

wfs = shackHartmann(nLenslet,nRes,minLightRatio );
%%%wfs.lenslets.nyquistSampling = 0.5;
ngs = ngs.*tel*wfs;

wfs.INIT
+wfs;

% The next 2 commands allow the displays of the frame and of the slopes to
% be updated when a new frame and new slopes are computed
wfs.camera.frameListener.Enabled = true;
wfs.slopesListener.Enabled = true;

% The WFS must be calibrated such as for 1rd of tip--tilt wavefront , it will
% measured a slopes of 1rd.
wfs.pointingDirection = zeros(2,1);  %%% DON'T COMMENT - important for calibration

% whereas the source is progressively moved off-axis
pixelScale = ngs.wavelength/(2*d*wfs.lenslets.nyquistSampling);
tipStep = pixelScale/2;
nStep   = floor(nPx/3)*2;
sx      = zeros(1,nStep+1);
u       = 0:nStep;

wfs.camera.frameListener.Enabled = false;
wfs.slopesListener.Enabled = false;

warning('off','oomao:shackHartmann:relay')
%%% Begin calibration on the tilt angle 
            for kStep=u
                ngs.zenith = -tipStep*kStep;  %%% changing zenith - check the NGS name!
                +ngs;
                drawnow
                sx(kStep+1) = median(wfs.slopes(1:end/2));
            end
%%% Begin calibration on the tilt angle 
warning('on','oomao:shackHartmann:relay')

Ox_in  = u*tipStep*constants.radian2arcsec;
Ox_out = sx*ngs.wavelength/d/2*constants.radian2arcsec;

% figure, plot(Ox_in,Ox_out), grid

slopesLinCoef = polyfit(Ox_in,Ox_out,1);
wfs.slopesUnits = 1/slopesLinCoef(1);

% The source is reset on--axis and the WFS is set to always be aligned to
% the source by setting \oop{shackHartmann}{pointingDirection} to empty.
wfs.pointingDirection = []; %%% DO NOT UNCOMMENT IT!