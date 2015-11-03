

% Telescope, atmosphere and wfs parameters
D = 8;
r0 = 20e-2;
L0 = 30;
nLenslet = 20;
nPxLowRes = nLenslet+1;
nPx = 200;
sFreq = 500;
% 85% illumination cut-off
illRat = 0.85;

% Derived parameters
% Size of one lenslet
d = D/nLenslet;
% Pixels per lenslet
nPxLens = nPx/nLenslet;

% 1 layer atmosphere
atm = atmosphere(photometry.V,r0,L0,'windSpeed',15,'windDir',0);

phi = fourierPhaseScreen(atm,D,nPx,nScreens);
%phiLowRes = fourierPhaseScreen(atm,D,nPxLowRes);

%save(filename,'-v7.3','D','r0','L0','nLenslet','nPx','sFreq','illRat','d','nPxLens','atm','nScreens','phi','nPxLowRes','phiLowRes')
save(filename,'-v7.3','D','r0','L0','nLenslet','nPx','sFreq','illRat','d','nPxLens','atm','nScreens','phi')
