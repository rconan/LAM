close all
clc
warning off

kmv_clean_cache = 0; %% put 1 here if you want to clean all and start over again (otherwise leave 0)s
if ( kmv_clean_cache == 1)
    clear all
end

if (exist('nL') == 0 ) %% if the simulator is started first time - do all the preparations, otherwise skip.
% \subsubsection{The source}
% A simple on--axis \ngs object is created by calling the \oo{source} constructor without parameters:
ngs = source; %create a simple on-axis guide star object

% \subsubsection{The atmosphere}
% A 3--layer atmosphere with a 20cm Fried parameter in the
% visible and a 30m outer scale is created with:

atm = atmosphere(photometry.V,20e-2,30,... % r0 = 20cm; L0 = 30m
    'fractionnalR0',[0.5,0.3,0.2], ...
    'altitude',     [0e3,5e3,12e3],...
    'windSpeed',    [10,5,20],...
    'windDirection',[0,0,0]);


% \subsubsection{The telescope}
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

% \subsubsection{The wavefront sensor}
minLightRatio = 0.85; % The \oop{lensletArray}{minLightRatio} property of the \oo{lensletArray}
% object set the minimum ratio of light intensity between a partially and
% fully illuminated lenslet. In the following, \oop{lensletArray}{minLightRatio} is set to 85\%:
wfs = shackHartmann(nL, nRes, minLightRatio);

% \paragraph{Wavefront sensor initialization}
% The next component to define is a deformable mirror (DM).
% The WFS and the DM are optically conjugated to the telescope pupil and
% the DM actuators are aligned to the WFS lenslet array following the Fried
% geometry i.e. they are located at the lenslet corners.

ngs = ngs.*tel*wfs;


% A warning is issued because the lenslets at the corner of the array are in the dark. 
% This issue is solved by selecting the lenslets that are receiving enough light to be useful to the wavefront detection.
wfs.INIT

% The WFS must be calibrated such as for 1rd of tip--tilt wavefront , it will
% measured a slopes of 1rd.
% To do so, the \oop{shackHartmann}{pointingDirection} property is set on-axis

wfs.pointingDirection = zeros(2,1);

% Lets switch off the display automatic update:
% wfs.camera.frameListener.Enabled = false;
% wfs.slopesListener.Enabled = false;


%% The closed loop
% Combining the atmosphere and the telescope
tel = tel+atm;
figure
imagesc(tel)

% Propagation throught the atmosphere to the telescope
ngs=ngs.*tel;

% Saving the turbulence aberrated phase
turbPhase = ngs.meanRmPhase; %% this is a turbulent phase

% Propagation to the WFS
ngs=ngs*wfs;

end %% if (exist('nL') == 0) %% if the simulator is started first time - do all the preparations, otherwise skip.


% Display of turbulence and residual phase
figure(11), h = imagesc([turbPhase]);
axis equal tight

% closing the loop
nIteration = 100;
kmv_pausetime = 0.1;

total  = zeros(1,nIteration);
residue = zeros(1,nIteration);

%% This is actuall closed-loop simualtions with kIteration steps.
tic
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
toc