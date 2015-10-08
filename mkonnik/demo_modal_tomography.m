clc
close all
clear all
% warning off
%%% using the modalTomographyHowto.m example for setting up 3NGSs

tool_inside_path_adder('matlab_code/mikhail', 'matlab_code');

% randn('state', 25);     % sets the global random state


params = struct; %%% creating the Structure with additional parameters for the computeARmodel
params.show_figures = 1; %% display the diagnostic plots during the simulation.

%% \subsubsection{The source}
% A simple on--axis \ngs object is created by calling the \oo{source} constructor without parameters:
ngs = source; %create a simple on-axis guide star object
%%% ngs.magnitude = 15;


%% \subsubsection{The atmosphere}
altitudes        = [0, 10]*1e3;
fractionalR0    = [0.7,0.3];
windSpeed       = [10,5];
windDirection   = [0,0];

L0 = 40; %% outer scale, [m]
r0 = 0.19; %% Fried parameter, [m]


if (exist('atm') == 0)
atm = atmosphere(photometry.V, r0, L0,...
    'fractionnalR0',fractionalR0, ...
    'altitude',     altitudes,...
    'windSpeed',    windSpeed,...
    'windDirection',windDirection);

end %% if (exist('atm') 

%% \subsubsection{The telescope}
nLenslet    = 16; % number of lenslets
nPx         = 8; % pixels per lenslet
nRes        = nLenslet*nPx; % this Resolution variable means that the "phase" var will be of the same size
D           = 8;  % telescope diameter [m]
samplingTime= 0.001;% seconds

% Now, we create a telescope of diameter \texttt{D} with a pupil sampled with \texttt{nRes} pixel.
tel = telescope(D, ...
    'resolution',nRes,...
    'fieldOfViewInArcMin',2, ...  %%% field of view, can be in ArcMins
    'samplingTime',samplingTime);
%     'fieldOfViewInArcsec',90, ...


%% \subsubsection{The wavefront sensor}
minLightRatio = 0.6; %% large values (more than 0.8  can lead to poor reconstruction - you are truncating the phase!).

% The \oop{lensletArray}{minLightRatio} property of the \oo{lensletArray}
% object set the minimum ratio of light intensity between a partially and
% fully illuminated lenslet. In the following, \oop{lensletArray}{minLightRatio} is set to 85\%:
wfs = shackHartmann(nLenslet, nRes, minLightRatio);

% Propagation of the calibration source to the WFS through the telescope
ngs = ngs.*tel*wfs;

%% \paragraph{Wavefront sensor initialization}
% The WFS must be calibrated such as for 1rd of tip--tilt wavefront , it will
% measured a slopes of 1rd.
% To do so, the \oop{shackHartmann}{pointingDirection} property is set on-axis
    wfs.INIT;
%     wfs.camera.photonNoise = false;  %%% DON'T TOUCH IT! IT IS CALIBRAION!!!
% % % %     wfs.pointingDirection = zeros(2,1);  %% do NOT use this line as it locks the WFS into upright position
    wfs.pointingDirection = []; %%% Keep it here!
% % % A new frame read-out and slopes computing:
    +wfs;


% % Saving the turbulence aberrated phase
turbPhase = ngs.meanRmPhase; %% this is a turbulent phase


%% Deformable mirror
nActuator = nLenslet + 1;
bif = influenceFunction('monotonic',25/100);
dm = deformableMirror(nActuator,...
    'modes',bif,...
    'resolution',nRes,...
    'validActuator',wfs.validActuator,...
    'zLocation',0);  %% the DM is assumed to be on the ground layer

    
%% Building the system
ngs=ngs.*tel*dm*wfs;
+ngs;



%% Zernike coefs to WFS slopes calibration
maxRadialDegree = 3; %% - setting the Max Radial Degree = how many modes we will use
zernModeMax = zernike.nModeFromRadialOrder(maxRadialDegree);

zern = zernike(2:zernModeMax,'resolution',nRes,'pupil',tel.pupil); %% Without Piston!!!!!!!!!!!!
%%% zern.lex = false;  %% don't use this crapy option (prevent lexicographic ordering) at all! Default is true!!!
zern.c = eye(zern.nMode);

    ngs=ngs.*zern*wfs;  %%% propagating the ngs to WFS via Zernike

%% Caibrating the WFS - the modes from the NGS measurements
z = zernike(1:zernModeMax)\wfs;
Dz = z.c;%(2:end,:);


% %% Measurements With noise
wfs.camera.photonNoise = true;
% wfs.camera.photonNoise = false;
ngs=ngs.*tel*wfs;  %% Propagating the NGS through telescope onto WFS

% %% noise convariance matrix
nMeas = 10;
slopes = zeros(wfs.nSlope,nMeas);

for kMeas=1:nMeas
    +wfs
    slopes(:,kMeas) = wfs.slopes;
end

Cn = slopes*slopes'/nMeas;

z = z\wfs;  %% projecting the wfs measurements onto Zernike
Czn = z.c*z.c'/nMeas; %% Covariance matrix of noisy measurements

if (params.show_figures == 1)
    figure(5)
    subplot(1,2,1)
    imagesc(Cn)
    axis equal tight
    colorbar
    subplot(1,2,2)
    imagesc(Czn)
    axis equal tight
    colorbar
end


%% Constructing matrices for the Rmv reconstruction
CznAst = blkdiag( Czn , Czn , Czn ); %%% block-diag matrix, the size of N Ngs
DzAst = blkdiag( Dz , Dz , Dz  );


%%
%% BEGIN::::: TOMOGRAPHY ::::::::::::::::
%%
tel = tel+atm;
params.do_projections = 'via_analytical_small_footprint'; %%% faster
% params.do_projections = 'via_numerical_small_footprint';  %%% slower

%% Define the Asterism (star constellation) for the tomography
%%% Don't forget that the FOV of the telescope (meta-pupil) must cover all the guide stars
dAsterism               = 1.0; % 0.5, 0.75, 1.0, 1.5, 2 diameter in arcmin
guideStarWavelength     = photometry.R;
guideStarMagnitude      = 14;%[14,14.5,15,15.5,16,16.5,17]; % lower later

%% Guide Stars with Regular asterism
ast = source('asterism',{[3,arcmin(dAsterism/2),0]},...
            'wavelength',guideStarWavelength,... 
            'magnitude',guideStarMagnitude);

%% %% This is Science object
sci_object = source('zenith',0*cougarConstants.arcsec2radian,'wavelength',photometry.R); %%% this is probably Projection Disk (like the telescope pupil)

%% Plot the Tomography setup - NGS and Sci Object positions. 
if (params.show_figures == 1)
    figure(20), imagesc(tel,[ast,sci_object])
end



%% Determining the number of Guide stars, layers and sci objects.
nSciobj = length(sci_object);
nDm = length(dm);
nGs = length(ast);
nLayers = atm.nLayer;

% zernModeMax = zernike.nModeFromRadialOrder(maxRadialDegree);
% zernWfs = zernike(2:zernModeMax,'resolution',nRes,'pupil',tel.pupil);
zernWfs = zern;

%% Zernike Projection for the DM
altitudes = dm.zLocation;
Projection_dm_Cell = tool_do_modal_projection(zernWfs, altitudes, sci_object, tel, nDm, nSciobj, params);
Projection_dm = cell2mat(Projection_dm_Cell);


%% Zernike projection for atmosphere layers
altitudes = [atm.layer.altitude];
Projection_atmlayers_Cell = tool_do_modal_projection(zernWfs, altitudes, ast, tel, nLayers, nGs, params);
Projection_atmlayers = cell2mat(Projection_atmlayers_Cell);


%% Modal Projection of the science object 
altitudes = [atm.layer.altitude];
Projection_sciobj_Cell = tool_do_modal_projection(zernWfs, altitudes, sci_object, tel, nLayers, nSciobj, params);

Projection_sciobj = cell2mat(Projection_sciobj_Cell);
Projection_sciobj = Projection_sciobj(2:end,:); %%% Removing piston mode?



%% Computing Turbulence Angular covariance matrix
Sangcov_Cell = phaseStats.zernikeAngularCovariance(zernWfs,atm,ast); %%% Computing the Angular Covariance
Sangcov_matrix = cell2mat(Sangcov_Cell);


%% Computing the cross-correlation matrix of Zernike coefficients from the modes of the atmosphere object for EACH ATM LAYER
Scov_layer_Cell = cell(nLayers,nLayers);  %%% composite cell array for Atm Covariance

% zernWfsi = zernike(2:zernModeMax,'resolution',nRes,'pupil',tel.pupil);
zernWfsi = zern;

for kLayer = 1:nLayers
    for kLayer1 = 1:nLayers
        
        if kLayer == kLayer1
            zernWfsi.D = tel.diameterAt( atm.layer(kLayer).altitude )/tel.D*2;
            Scov_layer_Cell{kLayer,kLayer1} = ...
                phaseStats.zernikeCovariance(zernWfsi,atm)*atm.layer(kLayer).fractionnalR0; %% %% this is from zernike.m
        else
            Scov_layer_Cell{kLayer,kLayer1} = zeros(zernWfs.nMode);
        end
    end
end  %%% for kLayer = 1:nLayers

Scov_layer_matrix = cell2mat(Scov_layer_Cell); %% matrix for the Atmosphere Covariance



%% Projection of the Atmosphere Cross-correlation between ATM and on DM
idy = [];
for kL=1:nLayers
    ii = (2:zernWfs.nMode+1)-1;
    idy = [idy ii+zernWfs.nMode*(kL-1) + kL];
end

idx = [];
for kGs=1:nGs
    ii = (2:zernWfs.nMode+1)-1;
    idx = [idx ii+zernWfs.nMode*(kGs-1)+kGs];
end

Sigma_theta = Projection_atmlayers(idx,:)*Scov_layer_matrix*Projection_atmlayers(idx,:)'; %% 

% 
% if (params.show_figures == 1)
%     %% figures for the theta directions
%     for kGs=1:nGs
%         subplot(1,nGs,kGs)
%         semilogy(diag(Sigma_theta,zernWfs.nMode*(kGs-1)),'k');
%         hold on
%         semilogy(diag(Sangcov_matrix,zernWfs.nMode*(kGs-1)),'r--');
%         legend('MV','L&A')
%         title(['Diagonals of #1 with #' num2str(kGs)])
%         hold off
%     end
% end
% 
% if (params.show_figures == 1)
%     figure
%     subplot(1,2,1)
%     imagesc(Sigma_theta)
%     title('MV')
%     subplot(1,2,2)
%     imagesc(Sangcov_matrix)
%     title('L&A')
% end
 
% mymat = Sigma_theta(1:zernWfs.nMode,zernWfs.nMode+1:end);
% 
% if (params.show_figures == 1)
%     figure
%     loglog(svd(Sigma_theta),'k')
%     hold on
%     loglog(svd(mymat(1:zern.nMode,1:zern.nMode)),'k')
%     loglog(svd(Sangcov_matrix(1:zern.nMode,1:zern.nMode)),'r');
%     title('Singular value decomposition')
%         legend('MV','L&A') %%% L&A - learn and apply
%     axis tight
%     hold off
% end





%% Computing the Projection for the sci object:
Sigma_thetaalpha = Projection_sciobj*Scov_layer_matrix*Projection_atmlayers(idx,:)';
% % % Projection_atmlayers(idx,:) - taking all the modes except piston?

Rmv = (Sigma_thetaalpha* DzAst')/(DzAst*Sigma_theta*DzAst' + CznAst);
%%% DzAst here consists of NaNs, so does Rmv

fprintf('\n PROJECTIONS DONE\n \n')
%% End calculating the projections



ast=ast.*tel*wfs;  %%% This obviously changes the projections - WHY?!???!
% ast=ast.*tel*wfs;  %%% This obviously changes the projections - WHY?!???!
+wfs;  %%% either this, or propagate twice
% 
phaseAst = [ast.meanRmPhase];
% figure, mesh(phaseAst);

%% wavefront reconstruction least square fit
z = z\wfs;
zern.c = Dz\z.c;%(2:end,:);
%%% zern.lex = true; %%% WHY????
phaseLS = reshape(zern.p*zern.c, nRes, nRes*length(ast));  %%% first run - you have NaNs, but second - it works.
figure, mesh(phaseLS);


%% wavefront reconstruction minimum variance
M = Sangcov_matrix*DzAst'/(DzAst*Sangcov_matrix*DzAst'+CznAst);
% z = z - 1; % piston removed %% it is not there! It has been overwritten by WFS
zern.c = reshape(M*z.c(:),z.nMode,[]);
phaseMV = reshape(zern.p*zern.c,nRes,nRes*length(ast));

figure,imagesc([phaseAst;phaseLS;phaseMV])
axis equal tight xy



%% Data/Target covariance
Sangcov_target_Cell = phaseStats.zernikeAngularCovariance(zernWfs,atm,ast,sci_object);
Sangcov_target_matrix = cell2mat(Sangcov_target_Cell);

%%% The target matrix T depends on the observed target location \theta, but also on the DM parameters and atmospheric parameters.
T = 0;
R = 0;
numSciobj = ones(size(sci_object));

for kSci = 1:nSciobj
    ProjectionDM = cell2mat(Projection_dm_Cell(kSci,:));

    T = T + numSciobj(kSci)*DzAst*cell2mat( Sangcov_target_Cell(:,kSci) )*ProjectionDM(2:end,:);
    R = R + numSciobj(kSci)*(ProjectionDM'*ProjectionDM);
end


%% Command matrix
M2_nominator = R\T';
M2_denominator = (DzAst*Sangcov_matrix*DzAst'+CznAst);
M2 = M2_nominator/M2_denominator;

ast=ast.*tel*wfs;
z = zernike(1:zernike.nModeFromRadialOrder(maxRadialDegree))\wfs;
z = z - 1;
zern.c = M2*z.c(:);
sci_object = sci_object.*tel;

phaseSciobject = [ast.meanRmPhase];

figure
imagesc([reshape(zern.p*zern.c,nRes,nRes), sci_object.meanRmPhase])
axis equal tight xy
colorbar
