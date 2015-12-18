close all
samplingTime     = 0.05;% seconds

%% Atmospheric parameters
wavelength      = photometry.R;
altitudes       = [0]*1e3;
fractionalR0    = [1];
windSpeed       = [ 5];
windDirection   = [ pi];

L0 = 30; %% outer scale, [m]
r0 = 0.20; %% Fried parameter, [m]


%% Creating the Atmosphere object
atm = atmosphere(photometry.V, r0, L0,...
    'fractionnalR0',fractionalR0, ...
    'altitude',     altitudes,...
    'windSpeed',    windSpeed,...
    'windDirection',windDirection);
nLayers = atm.nLayer;


%% TELESCOPE
D    = 8;
nLenslet = 10;
d    = D/nLenslet; % lenslet pitch
nRes = 4*nLenslet;
tel = telescope(D,...
    'resolution',nRes,...
    'fieldOfViewInArcsec',180,...
    'samplingTime', samplingTime);


%% Generate the Zernike with only modes 2 and 3 - that is, Tip and Tilt
zern = zernike(2:3,'resolution',nRes, 'pupil',tel.pupil);


TTAst = source('asterism',{[3,   arcsec(20),0], [1,   0, 0]},...
    'magnitude', 10, 'wavelength',photometry.H); %% Tip-tilt asterism (for TT measurements)

TTprojection =[];
tel = tel+atm;

nPts=100;
for k = 1:nPts
    
    +tel
    TTAst = TTAst.*tel;
    zern \ TTAst(1).phase(:);
    TTprojection = [TTprojection,  zern.c ];

end

%% Plots 
figure, plot(TTprojection(2,:))
figure, plot(TTprojection(1,:))

%%

figure 

Acorr_TTproj = xcorr(TTprojection(2,:));

plot(Acorr_TTproj./max(Acorr_TTproj(:)));

hold on
Acorr_TTproj = xcorr(TTprojection(2,:));
plot(Acorr_TTproj./max(Acorr_TTproj(:)));

%% Here
tel = tel-atm;
S = compute_spatio_angular_temporal_correlation(atm, tel, zern, 50);


%% Loop
for x = 1:50
    Acorr_SA(x) = S{x,1}(1);
end

time_acorr = [0:samplingTime:samplingTime*49]
time_SA = [0:samplingTime:samplingTime*(nPts-1)]

figure, plot(time_SA,  Acorr_TTproj(nPts:end)./max(Acorr_TTproj(:) ))
hold on
plot(time_acorr,  Acorr_SA./max(Acorr_SA(:) ))

%%
nu = logspace(-3,log10(0.5/samplingTime),201);

specZern = lamStats.temporalSpectrum(nu,atm,tel,3);

zern = zernike(4,'resolution',nRes, 'pupil',tel.pupil);
zern.D = 8;
specZern2 = zernikeStats.temporalSpectrum(nu,atm,zern);

specZern = lamStats.multipleSpectra(nu, atm, zern, tel);

nuLin = 0.005:0.1:0.5/samplingTime;

specZernLin = interp1(nu, specZern, nuLin);

specZernLin = [specZernLin specZernLin(end:-1:2)];


temporalCorr = real(fft(specZernLin));

dt  = 1/(2*(0.5/samplingTime))
time = 0:dt:5;
len = numel(time);
plot(time, temporalCorr(1:len)/max(temporalCorr))


