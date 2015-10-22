
% Simulation parameters for OOMAO comparisons of different reconstructors

%% Set parameters
%--------------------------------------------------------------------------
% Atmosphere parameters
r0 = 0.1549;
L0 = 25;
altitude = 0;
fractionalR0 = [0.5 0.2 0.3];
windSpeed = [10 5 15];
windDirection = [0 pi/2 pi];

% Source for wave front sensing
lambda_wfs = photometry.V;
magnitude_wfs = 8;

% Science source
lambda_sci = photometry.K;
magnitude_sci = 8;

% Telescope parameters
samplingFreq = 1000;
samplingTime = 1/samplingFreq;
D = 8;
fieldOfView = 2.5;          % In arcmins

% SH parameters
nLenslet = 40;
nPx = 6*nLenslet;           % 6 pixels per lenslet
illRatio = 0.75;

% DM parameters
% Coupling coeficient for dm influence function
DMcoef = 0.35;

% Derived parameters
% Lenslet dimensions
d = D/nLenslet;
% Pixels per lenslet
nLPx = nPx/nLenslet;
% Pixel size (m)
dx = D/nPx;
% No. of actuators
nActuators = nLenslet+1;
%--------------------------------------------------------------------------

%% Load cube of phase (or regenerate) if variable filename exists
% If filename doesn't exist we have no cube of phase, must use phase
% generated in simulation loop (i.e. for closed loop)
if exist('filename','var')
    if ~exist(filename,'file')
        % Generate atmosphere
        atm = atmosphere(lambda_wfs,r0,L0,'altitude',altitude,...
            'fractionnalR0',fractionalR0,'windSpeed',windSpeed,'windDir',windDirection);
        phi = fourierPhaseScreen(atm,D,nPx,nLenslet,nScreens);
        save(filename,'phi')
    else
        load(filename)
    end
end



