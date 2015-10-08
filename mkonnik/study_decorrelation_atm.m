clc
clear all
% close all
samplingTime     = 0.05;% seconds
% samplingTime     = 0.002;% seconds
tool_inside_path_adder('matlab_code/mikhail', 'matlab_code');

%% Atmospheric parameters (single layer ATM
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
zern = zernike(2:10,'resolution',nRes, 'pupil',tel.pupil);
zern.D = D;

TTAst = source('asterism',{[3,   arcsec(20),0], [1,   0, 0]},...
    'magnitude', 10, 'wavelength',photometry.H); %% Tip-tilt asterism (for TT measurements)


TTprojection =[];
tel = tel+atm;

        nPts=1001; %% the number of points of the phase (make it large)

for k = 1:nPts
    
    +tel
    TTAst = TTAst.*tel;
    zern \ TTAst(1).phase(:);
    TTprojection = [TTprojection,  zern.c ];

end

%% Plots 
% figure, plot(TTprojection(2,:))
% figure, plot(TTprojection(1,:))



%% Compute the spatio-angular temporal correlation
tel = tel-atm;
N_points_SA = 100;
S = compute_spatio_angular_temporal_correlation(atm, tel, zern, N_points_SA);






%% Compute for other modes
zernike_mode_min  =1;
zernike_mode_max = 4;

for zz = zernike_mode_min:zernike_mode_max
zernike_mode = zz;

dt  = 1/(2*(0.5/samplingTime));
time = 0:dt:(samplingTime*nPts);
len = numel(time);

%% Computing Autocorrelation
Acorr_TTproj = xcorr(TTprojection(zernike_mode,:));
figure, plot(time, Acorr_TTproj(1:len)./max(Acorr_TTproj(:))  );
%     hold on


%% Loop too collect the Spatio-Angular autocorrelation.
for x = 1:N_points_SA
    Acorr_SA(x) = S{x,1}(zernike_mode,zernike_mode); 
end


time_SA = [0:samplingTime:samplingTime*(N_points_SA - 1)];
time_acorr = [0:samplingTime:samplingTime*(nPts-1)];

plot(time_acorr,  Acorr_TTproj(nPts:end)./max(Acorr_TTproj(nPts:end) ))
hold on
plot(time_SA,  Acorr_SA./max(Acorr_SA(:) ));


%%  Specturm of Zernike
number_of_points_sampling = round(0.5*1/samplingTime);

min_freq = 1e-3; %% smallest frequency we are interested in.
max_freq = 2*number_of_points_sampling; %% largest frequency we are interested in.

nu = logspace(log10(min_freq),log10(max_freq), number_of_points_sampling); %% frequencies for Zernike temporal spectrum
%%%%% y = logspace(a,b,n) generates n points between decades 10^a and 10^b.

        params.freq_part = 0.2;  %% default is 200, or 0.2 of maximum freq.
        params.fN = 1001; % number of points in the frequency domain - a smaller number of points leads to wrong results in the decorrelation functions cc, Jul2013

        FRange =  [min_freq round(params.freq_part*max_freq)]; %Frequency Range

        f1 = linspace(nu(1), max(FRange), params.fN)'; %% nu(1) here is the smallest frequency we are interested in. 
        f  = f1(2:end)';

%% New calculation of power specturm
PSDi = lamStats.temporalSpectrum(nu,atm,tel,zernike_mode);;

% ---  Interpolation taking place ---
sp = interp1(nu, PSDi, f1)'; % returns interpolated values of a 1-D function. Vector nu contains the sample 
% points, PSDi contains the corresponding values, and vector f1 contains the coordinates of the query points.
sp_full  = [0 sp(end:-1:2)  sp]; %% we flip the sp_full vector from left to right and take all elements except the last one (so is ...:2) stands for.
Cphi_aoa = real(fftshift(fft(fftshift(sp_full))));        % autocorrelation thru Wiener-Khinchine theorem
Cphi_aoa = Cphi_aoa/max(Cphi_aoa);
correl   = Cphi_aoa(1,end/2+1:end-1); % taking only a half of the correlation function.

dt_correl = 1/max_freq;
time_end = 0.5*( params.fN*params.freq_part )*samplingTime;
time_correl = 0:dt_correl:time_end;

plot(time_correl, correl(1:length(time_correl)))

% legend ('Autocorrelation, TT projection zern/TTast(1)','Autocorrelation, Spatio-Angular', 'Temporal correlation of ZernikeLinear', 'From ARmodel covar. via powerspec')
legend ('Autocorrelation, TT projection zern/TTast(1)','Spatio-Angular Acorr', 'Acorr from ARmodel via powerspec')
hold off


end %% for zz = 1:8 zernike_mode = zz;


% figure, loglog(nu, specZern)
% hold on
% loglog(nu, nu.^(-17/3))
% hold off

figure, plot(Acorr_SA./max(Acorr_SA(:)))
figure, plot(time_correl, correl(1:length(time_correl))./max(correl)  )
break
% %% Computing via AR model
% params.flag_correlation_method = 'corr_via_PSD';
% 
% params.order  = 1;  %% order of the AR model to fit into Zernike
% params.Npts = 20; %% Npts is Number of points to consider in the fitting
% % params.tau_f = 16e-3; %% time for estimation of the AR model (instead of Npts)
% %     params.tau_f = 'auto'; %% time for estimation of the AR model (instead of Npts), use the 1/e decay model.
% 
% params.fN = 1001; % number of points in the frequency domain - a smaller number of points leads to wrong results in the decorrelation functions cc, Jul2013
% params.freq_part = 0.45;  %% default is 200, or 0.2 of maximum freq.
% 
% [ARmodel, ZPSD, correl_all] = lam_utilities.tool_compute_AR_model_via_covariance_bfgs_fitting(atm,tel,zern,params); %% fits the AR model into atm data.

% 
% 
% plot(sp)
% hold on
% plot(specZernLin,'r')
% whos sp
%   Name      Size              Bytes  Class     Attributes
% 
%   sp        1x4001            32008  double              
% 
% whos sp_full
%   Name         Size              Bytes  Class     Attributes
% 
%   sp_full      1x8002            64016  double              
% 
% plot(sp_full)
% f1(end)
% 
% ans =
% 
%      5
% 
% whos nu1
%   Name      Size              Bytes  Class     Attributes
% 
%   nu1       1x1001             8008  double              
% 
% nu1(end)
% 
% ans =
% 
%     10
% 
% whos nu
%   Name      Size              Bytes  Class     Attributes
% 
%   nu        1x1001             8008  double              
% 
% nu(end)
% 
% ans =
% 
%    100