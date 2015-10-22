
% Script for wfs calibration
% Calibrates slope units and plots response over and beyond linear range

% Set wfs on-axis
wfs.pointingDirection = zeros(2,1);
% Start with ngs on axis
ngs.zenith = 0;

% The source is progressively moved off-axis
% Rather than measure dx or dy in pixels (or m) the angle is measured:
% 1 pixel is equivalent to an anglular change of source direction.
pixelScale = ngs.wavelength/...
    (2*d*wfs.lenslets.nyquistSampling);
% Move source across in steps of 1/2 pixels
tipStep = pixelScale/2;
% Chose number of steps so we stay within the linear regime, i.e.
% nStep < nPxLens (we haven't moved across more pixels than there are per 
% lenslet 
kStep = 0:floor(nLPx/3)*2;
sx      = zeros(1,length(kStep));
% Don't replot camera/slopes
wfs.camera.frameListener.Enabled = false;
wfs.slopesListener.Enabled = false;
warning('off','oomao:shackHartmann:relay')
% Move source and record wfs measurement
for k=1:length(kStep)
    % Move source
    ngs.zenith = -tipStep*kStep(k);
    +ngs;
    drawnow
    % Get median wfs signal (for tip-tilt should be constant across
    % lenslets when within linear regime)
    sx(k) = median(wfs.slopes(1:end/2));
end
warning('on','oomao:shackHartmann:relay')
% Convert input and output to arcsecs
Ox_in  = kStep*tipStep*constants.radian2arcsec;
Ox_out = sx*ngs.wavelength/d/2*constants.radian2arcsec;
% Plot result
figure
    plot(Ox_in,Ox_out)
    xlabel('WF tilt [arcsecs]')
    ylabel('WFS signal [arcsecs]')
    legend('Diffractive SH response (linear range)')
    grid

% Calibrate units
slopesLinCoef = polyfit(Ox_in,Ox_out,1);
wfs.slopesUnits = 1/slopesLinCoef(1) * wfs.slopesUnits;

% Then test plot over greater range: to non-linear part
kStep = linspace(0,nPx/nLenslet/2+1,100);
warning('off','oomao:shackHartmann:relay')
for k=1:length(kStep)
    % Move source (where 1 kStep corresponds to 1 pixel)
    ngs.zenith = -tipStep*kStep(k)*2;
    +ngs;
    % Get median wfs signal (for tip-tilt should be constant across
    % lenslets when within linear regime)
    sx(k) = median(wfs.slopes(1:end/2));
end
warning('on','oomao:shackHartmann:relay')

% Plot result
figure
    plot(kStep,sx)
    xlabel('WF tilt [pixels]')
    ylabel('WFS signal [pixels]')
    legend('Diffractive SH response')
    grid
    
% Then test plot over small range
% Want to check limits of wfs to sense distortion
kStep = linspace(0,nPx/nLenslet/100,100);
warning('off','oomao:shackHartmann:relay')
for k=1:length(kStep)
    % Move source (where 1 kStep corresponds to 1 pixel)
    ngs.zenith = -tipStep*kStep(k)*2;
    +ngs;
    % Get median wfs signal (for tip-tilt should be constant across
    % lenslets when within linear regime)
    sx(k) = median(wfs.slopes(1:end/2));
end
warning('on','oomao:shackHartmann:relay')

% Plot result
figure
    plot(kStep,sx)
    xlabel('WF tilt [pixels]')
    ylabel('WFS signal [pixels]')
    legend('Diffractive SH response')
    grid

% The source is reset on--axis and the WFS is set to always be aligned to
% the source by setting pointing direction [].
ngs.zenith = 0;
+ngs
wfs.pointingDirection = [];

