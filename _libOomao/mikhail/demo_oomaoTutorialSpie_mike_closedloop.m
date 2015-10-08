clc
close all
clear all
warning off


%% \subsubsection{The source}
% A simple on--axis \ngs object is created by calling the \oo{source} constructor without parameters:
ngs = source; %create a simple on-axis guide star object

% dAsterism               = 0.5; % 0.5, 0.75, 1.0, 1.5, 2 diameter in arcmin
% guideStarWavelength     = photometry.R;
% guideStarMagnitude      = 14;%[14,14.5,15,15.5,16,16.5,17]; % lower later
% 
% %% Guide Stars with Regular asterism
% ngs = source('asterism',{[3,arcmin(dAsterism/2),0]},...
%             'wavelength',guideStarWavelength,... 
%             'magnitude',guideStarMagnitude);


% In the following, an on-axis natural guide star in V band is defined.
% ngsCal = source('wavelength',photometry.J);


%% \subsubsection{The atmosphere}
altitude = [0, 5, 12]*1e3;
fractionalR0 = [0.5,0.3,0.2];
windSpeed = [10,5,20];
windDirection = [0,0,0];

L0 = 40; %% outer scale, [m]
r0 = 0.19; %% Fried parameter, [m]

atm = atmosphere(photometry.V, r0, L0,...
    'fractionnalR0',fractionalR0, ...
    'altitude',     altitude,...
    'windSpeed',    windSpeed,...
    'windDirection',windDirection);


%% \subsubsection{The telescope}
% Lest first make the NGSAO model scalable by defining the number of
% lenslet \texttt{nL} of the wavefront sensor, the number of pixel per lenslet \texttt{nPx}, the
% telescope diameter \texttt{D} in meter and the sampling frequency \texttt{samplingFreq} in Hz.

nL   = 16; % number of lenslets
nPx  = 8; % pixels per lenslet
nRes = nL*nPx; % this Resolution variable means that the "phase" var will be of the same size
D    = 8;  % telescope diameter [m]
d    = D/nL; % lenslet pitch
samplingFreq = 500;  % sampling frequency [Hz]

% Now, we create a telescope of diameter \texttt{D} with a pupil sampled with \texttt{nRes} pixel.
% The telescope field--of--view is 30arcsec and the temporal evolution of
% wavefront in the telescope pupil will be sampled at \texttt{samplingFreq}.

tel = telescope(D, ...
    'resolution',nRes,...
    'fieldOfViewInArcsec',30, ...
    'samplingTime',1/samplingFreq);

%% \subsubsection{The wavefront sensor}
minLightRatio = 0.85; % The \oop{lensletArray}{minLightRatio} property of the \oo{lensletArray}
% object set the minimum ratio of light intensity between a partially and
% fully illuminated lenslet. In the following, \oop{lensletArray}{minLightRatio} is set to 85\%:
wfs = shackHartmann(nL, nRes, minLightRatio);

% Propagation of the calibration source to the WFS through the telescope
ngs = ngs.*tel*wfs;

%% \paragraph{Wavefront sensor initialization}
% The WFS must be calibrated such as for 1rd of tip--tilt wavefront , it will
% measured a slopes of 1rd.
% To do so, the \oop{shackHartmann}{pointingDirection} property is set on-axis
    wfs.INIT;
    wfs.camera.photonNoise = false;
    wfs.pointingDirection = zeros(2,1);
% % % A new frame read-out and slopes computing:
    +wfs;
break



break
%% The closed loop
% Combining the atmosphere and the telescope
tel = tel+atm;

% Propagation throught the atmosphere to the telescope
ngs=ngs.*tel;

% Saving the turbulence aberrated phase
turbPhase = ngs.meanRmPhase; %% this is a turbulent phase

% Propagation to the WFS
ngs=ngs*wfs;


% Display of turbulence and residual phase
figure(11), h = imagesc([turbPhase]); %%%  turbPhase = ngs.meanRmPhase; %% this is a turbulent phase
axis equal tight

% closing the loop
nIteration = 100;
kmv_pausetime = 0.1;

total  = zeros(1,nIteration);
residue = zeros(1,nIteration);

%% This is actuall closed-loop simualtions with kIteration steps.
for kIteration=1:nIteration
    
    % Propagation throught the atmosphere to the telescope, +tel means that
    % all the layers move of one step based on the sampling time and the
    % wind vectors of the layers
    ngs=ngs.*+tel; 

    % Saving the turbulence aberrated phase
    turbPhase = ngs.meanRmPhase;

    % Finally, we get the phase of the atmospheric turbulence.s
    phase = ngs.meanRmOpd; %% atmospheric phase, in \mu m

    % Variance of the atmospheric wavefront
    total(kIteration) = var(ngs);

    % Propagation to the WFS
    ngs=ngs*wfs; 
    
    % The first method to estimate the wavefront is the finite difference (FD)
    % method using the inverse of the sparse gradient matrix computed with
    % the method \oom{shackHartmann}{sparseGradientMatrix}.
    phaseEst = tools.meanSub(wfs.finiteDifferenceWavefront*ngs.wavelength, wfs.validActuator);

    % Display of turbulence and residual phase
    set(h,'Cdata',[turbPhase])
    
    figure(111),imagesc(phaseEst)
    pause(kmv_pausetime)
    drawnow
end