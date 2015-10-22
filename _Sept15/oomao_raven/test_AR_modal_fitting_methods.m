%> @file test_AR_modal_fitting_methods.m
%> @brief This is the start-up script that sets the atmosphere at telescope 
%> parameters for testing the AR models fitting methods.
%> @author Mikhail Konnik (from the file "modalMethodsFeb2015.m")
%> @date   24 April 2015
%> 
%>  Modal Reconstructor Comparison
%>  1 = static MMSE
%>  2 = AR1 prediction
%>  3 = AR2 prediction
%>  4 = Ammse prediction
%>  5 = LQG with AR2 prediction
%>  6 = LQG with Ammse prediction

clear all
close all
clc
format long

addpath('./OOMAO-Raven/OOMAOlibUpdated/')  % Only use OOMAOlibUpdated!
addpath('./OOMAO-Raven/')

%% Setting up Parameters
samplingTime    = 0.001;% seconds
frameTime       = 0.01; % seconds
fixedLagTimeInMs = 3;
nIteration              = 1000;
exposureTime            = nIteration*frameTime;
randn('state', 25);     % sets the global random state


%% Atmosphere parameters
refWavelength   = 500e-9;
wavelength      = photometry.R;

altitude = [1]*1e3;
fractionalR0 = [1];
windSpeed = [10];
windDirection = [pi/2];

L0 = 40;
r0 = 0.19;

atm = atmosphere(refWavelength,r0,L0,...
    'altitude',altitude,...
    'fractionnalR0',fractionalR0,...
    'windSpeed',windSpeed,...
    'windDirection',windDirection);
atm.wavelength = wavelength;


%% Telescope
nPx = 128; %% pixels across the telescope's pupil size
D = 8; %% teelscope diameter [m]
telAR = telescope(D, ...
    'fieldOfViewInArcMin',3.5,...
    'resolution',nPx,...
    'samplingTime',samplingTime);

%%% Here we preparing the Zernike functions evaluated on a telescope pupil
%%% of size tel.D and pixels nPx, for we get zernProj.modes
maxRadialDegreeProj = 3;
zernModeMaxProj = zernike.nModeFromRadialOrder(maxRadialDegreeProj);
zernProj = zernike(2:zernModeMaxProj,'resolution',nPx,'pupil',telAR.pupil,'D',telAR.D);
% % % zernProj.lex = false;
% % % zernProj.c = eye(zernProj.nMode);



%% Switch between the methods:
flag_method = 'ARmodel';  %% use the compute AR model fitting
% flag_method = 'spat_ang_based';  %% computes the Zernike temporal correlation functions using the analytical angular covariance routines

%% Parameters for the fitting methods:
params = struct; %%% creating the Structure with additional parameters for the computeARmodel

% params.flag_correlation_method = 'corr_via_spatio_angular';
params.flag_correlation_method = 'corr_via_PSD';

params.order  = 1;  %% order of the AR model to fit into Zernike
params.Npts = 20; %% Npts is Number of points to consider in the fitting
% params.tau_f = 16e-3; %% time for estimation of the AR model (instead of Npts)
%     params.tau_f = 'auto'; %% time for estimation of the AR model (instead of Npts), use the 1/e decay model.

params.fN = 1001; % number of points in the frequency domain - a smaller number of points leads to wrong results in the decorrelation functions cc, Jul2013
params.freq_part = 0.45;  %% default is 200, or 0.2 of maximum freq.

params.flag_messages = 1; %% show or hide test messages from the function
params.flag_do_plots = 1;





switch flag_method
    
    case 'ARmodel'
%       [ARmodel, ZPSD, correl_all] = tool_compute_AR_model_via_covariance_bfgs_fitting(atm,telAR,zernProj,params); %% fits the AR model into atm data.
    [ARmodel, ZPSD, correl_all] = lam_utilities.tool_compute_AR_model_via_covariance_bfgs_fitting(atm,telAR,zernProj,params); %% fits the AR model into atm data.

        diag(ARmodel.A) %% show the identified model


    case 'spat_ang_based'
%        [AmmseLag, S]   = tool_compute_AR_model_via_spatio_angular_temporal_corr(atm,telAR,zernProj,params);
    [AmmseLag, S] = lam_utilities.tool_compute_AR_model_via_spatio_angular_temporal_corr(atm,telAR,zernProj,params);

%%%% Convert the cell array of covariances into the covariance matrix
        Layer_Number = 1;
        number_of_modes = length(zernProj.n);

        correl_all = zeros(number_of_modes, params.Npts);
        for ii=1:params.Npts
            correl_all(:,ii) = diag(S{ii,Layer_Number}); 
        end

        for ii=1:number_of_modes
            correl_all(ii,:) = correl_all(ii,:)./(max(correl_all(ii,:)));
        end

     diag(AmmseLag{1}) %% show the identified model

end %%% for switch


%% Plot the decorrelation graphs
if (params.flag_do_plots == 1)        
        % %%% BEGIN: Plotting for spat_ang_based
        plot_decorrelation_modes(correl_all, zernProj.n, params.Npts);
        % %%% END: Plotting for spat_ang_based
end %if (flag_do_plots == 1)        
