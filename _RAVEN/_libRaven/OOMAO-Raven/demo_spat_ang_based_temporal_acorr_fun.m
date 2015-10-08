% Hi Mikhail, here's a script that computes the Zernike temporal correlation functions using the analytical angular covariance routines.
% 
% Carlos

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ---- SPATIO-ANGULAR_BASED TEMPORAL AUTO-CORRELATION FUNCTIONS ---
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------------------------------------------------------------------------
% MULTI-LAYERED CASE
% -------------------------------------------------------------------------
% THEORY: Ammse = < || phi_{k+1}, phi_{k+1} ||^2 > * inv { < || phi_k, phi_k ||^2>}

%% Atmosphere
% --->>> alternative parameters ---
 altitude       = [0]*1e3;
 fractionalR0   = [1];
 windSpeed      = [15];
 windDirection  = [0]*pi/180;
 nLayer = 1;
 L0 = 30;
 refWavelength = 0.5e-6;
 r0 = 0.156;
 

%----------------------------------
atm = atmosphere(refWavelength,r0,L0,...
    'altitude',altitude,...
    'fractionnalR0',fractionalR0,...
    'windSpeed',windSpeed,...
    'windDirection',windDirection);

%atm.wavelength = wavelength;

%% Telescope

nPx_lenslet = 15;
nPx = 128; %% pixels across the telescope's pupil size

samplingTime = 2e-3*5;

dTel = 8;
telFieldOfView = 3.5;
tel = telescope(dTel,'fieldOfViewInArcMin',telFieldOfView,...
    'resolution',nPx,...
    'samplingTime',samplingTime);
%% Zernike measurement - Explicit Tomography
maxRadialDegree = 9;
zernModeMax = zernike.nModeFromRadialOrder(maxRadialDegree);
zern = zernike(2:zernModeMax,'resolution',nPx,'pupil',tel.pupil,'D',tel.D);
zern.lex = false;
zern.c = eye(zern.nMode);


Npts     = 11;                  % for the phase covar and 1-step phase covar
S        = cell(Npts, nLayer); % cell with the phase covar matrices
AmmseLag = cell(1,nLayer);     % cell with the 1-step MMSE temporal predictor

for kLayer = 1:nLayer
    % create a temp atmosphere with one layer at ground and another at 1km
    % (user option) to compute spatio-angular covariance matrices. The
    % option for a 2-layer temp atmosphere is to comply with the existing
    % code. One-layer atm object would create an integration error (Nov'12)
    refH = 1e3;  %% reference altitude Heighte
    altitude      = [0 refH];
    fractionalR0  = [0 atm.layer(kLayer).fractionnalR0];
    windSpeed     = [0 atm.layer(kLayer).windSpeed];
    windDirection = [0 atm.layer(kLayer).windDirection];

    myatm = atmosphere(refWavelength,r0,L0,...
        'altitude',altitude,...
        'fractionnalR0',fractionalR0,...
        'windSpeed',windSpeed,...
        'windDirection',windDirection);
    %myatm.wavelength = wavelength;

% Npts - number of points, time domain!

    for kr = 1:Npts
        alpha        = (kr-1)*tel.samplingTime*windSpeed(2)/altitude(2);
        alpha        = 1/2*alpha*180/pi*3600;
        ast          = source('asterism',{[2,alpha*cougarConstants.arcsec2radian,0]});
        S{kr,kLayer} = phaseStats.zernikeAngularCovariance(zern,myatm,ast);%*(dTel/tel.diameterAt(refH))^(5/3);
    end
    AmmseLag{kLayer} = S{2,kLayer}*pinv(S{1,kLayer});
end
% takes about 25 sec to compute with 54 Zernike modes


% check results: 1) sum(S{1}(:)) should be the full atm covariance matrix
FullAtmCov = phaseStats.zernikeCovariance(zern,atm);

FullAtmCov_fromLayers = zeros(zernModeMax-1);
for kLayer = 1:nLayer
    FullAtmCov_fromLayers = FullAtmCov_fromLayers + S{1,kLayer};
end

plot(diag(FullAtmCov))
hold on
plot(diag(FullAtmCov_fromLayers),'r--')


%S1_ang = S{2};
%Sn_ang = Sphi - A*Sphi*A' - B*Sphi*B' - B*S1_ang*A' - A*S1_ang*B';

t = 0:tel.samplingTime: Npts * tel.samplingTime-samplingTime;
for kr = 1:Npts
    myS(:,:,kr) = S{kr};
end
%myS = cell2mat(S);
%myS = reshape(myS,zernModeMax-1,zernModeMax-1,Npts);




% break

figure
m = {'--+','--o','--*','--.','--x','--s','--d','--^','--v','-->','--<','--p','--h'};
set(gca(), 'LineStyleOrder',m, 'ColorOrder',[0 0 0], 'NextPlot','replacechildren','LineWidth',1)

idx = [1 2 3 4 6 7 8];

tmp = myS(3,idx,:)./S{1}(3,3);
tmp = squeeze(tmp);
%plot(t,tmp,'--','LineWidth',1)
plot(t, tmp)
box on
legend('Z_2','Z_3','Z_4','Z_5','Z_7','Z_8','Z_9','Z_{10}')
title('Temporal auto- and cross-correlation functions with focus','fontsize',14,'fontweight','bold')
ylabel('Temporal auto-correlation','fontsize',14)
xlabel('Time [s]','fontsize',14)
%print('-depsc2','-tiff','Temporal_correlations_focus.eps')

figure
m = {'--+','--o','--*','--.','--x','--s','--d','--^','--v','-->','--<','--p','--h'};
set(gca(), 'LineStyleOrder',m, 'ColorOrder',[0 0 0], 'NextPlot','replacechildren','LineWidth',1)

tmp = myS(1,idx,:)./S{1}(1,1);
tmp = squeeze(tmp);
%plot(t,tmp,'--','LineWidth',1)
plot(t, tmp)
box on
legend('Z_2','Z_3','Z_4','Z_5','Z_7','Z_8','Z_9','Z_{10}')
title('Temporal auto- and cross-correlation functions with tilt','fontsize',14,'fontweight','bold')
ylabel('Temporal auto-correlation','fontsize',14)
xlabel('Time [s]','fontsize',14)
axis([0, (Npts-1) * tel.samplingTime, -0.2,1.01])
%print('-depsc2','-tiff','Temporal_correlations_tilt.eps')


