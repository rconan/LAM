%> @file tool_compute_AR_model_via_spatio_angular_temporal_corr.m
%> @brief A script that computes the Zernike temporal correlation functions using the analytical angular covariance routines.
%> @author Mikhail Konnik
%> @date   24 April 2015
%> @section spatangbased Spatio-angular-based temporal auto-correlation functions
%> This script computes the spatio-angular-based temporal auto-correlation functions for the multi-layered case.
%> 
%> THEORY: Ammse = < ||phi_k+1, phi_k||^2 > * inv { < ||phi_k, phi_k||^2>}
%> 
%> A script that computes the Zernike temporal correlation functions using the analytical angular covariance routines.
%======================================================================
%> @param atm =  atmosphere object.
%> @param tel =  telescope object.
%> @param zern = Zernike decomposition of the atmosphere.
%> @param params = the structure that contains other parameters like the number of points to consider (Npts) and order of the AR model.
%> @retval AmmseLag = the covariance matrix.
%> @retval S = the covariance cell array (can be converted into the matrix).
% ====================================================================== 
function [AmmseLag, S] = tool_compute_AR_model_via_spatio_angular_temporal_corr(atm,tel,zern,params)
Ts  = tel.samplingTime;

    if (isfield(params,'Npts') == 1);          %%% if there is an explicit parameter for number of points - use them
        Npts = params.Npts;  %% Npts is Number of points to consider in the fitting
    else
        if (isfield(params,'tau_f') == 1) && (isnumeric(params.tau_f) == 1)
            Npts = round(params.tau_f / Ts);
        else 
% coherence function decay for theta0 and tau0 computation (as defined in the atmosphere.m)
        % . Roddier (default): coherence function 1/e decay
        % . Fried: structure function equal 1rd^2 (1/\sqrt(e) decay)
           Npts = round( (0.5/exp(1))*(1/Ts) );
        end
    end

S        = cell(Npts, atm.nLayer); % cell with the phase covar matrices
AmmseLag = cell(1,atm.nLayer);     % cell with the 1-step MMSE temporal predictor

for kLayer = 1:atm.nLayer
    % create a temp atmosphere with one layer at ground and another at 1km
    % (user option) to compute spatio-angular covariance matrices. The
    % option for a 2-layer temp atmosphere is to comply with the existing
    % code. One-layer atm object would create an integration error (Nov'12)

    refH = 1e3;  %% reference altitude Height
    
    altitude      = [0 refH];
    fractionalR0  = [0 atm.layer(kLayer).fractionnalR0];
    windSpeed     = [0 atm.layer(kLayer).windSpeed];
    windDirection = [0 atm.layer(kLayer).windDirection];

    myatm = atmosphere(atm.wavelength, atm.r0, atm.L0,...
        'altitude',altitude,...
        'fractionnalR0',fractionalR0,...
        'windSpeed',windSpeed,...
        'windDirection',windDirection);

    for kr = 1:Npts  % Npts - number of points, time domain!
        alpha        = (kr-1)*tel.samplingTime*windSpeed(2)/altitude(2);
        alpha        = 1/2*alpha*180/pi*3600;
        ast          = source('asterism',{[2,alpha*cougarConstants.arcsec2radian,0]});
        S{kr,kLayer} = phaseStats.zernikeAngularCovariance(zern,myatm,ast);%*(dTel/tel.diameterAt(refH))^(5/3);
    end
    
    AmmseLag{kLayer} = S{2,kLayer}*pinv(S{1,kLayer}); %% gives the matrix A assuming AR(1),
    %%% although the matrix can have non-diagonal elements - coupling between the Zernike modes of the atmosphere.
end

%% check results: 1) sum(S{1}(:)) should be the full atm covariance matrix
% FullAtmCov = phaseStats.zernikeCovariance(zern,atm);
% 
% maxRadialDegreeProj = zern.nMode;
% 
% FullAtmCov_fromLayers = zeros(maxRadialDegreeProj);
% for kLayer = 1:atm.nLayer
%     FullAtmCov_fromLayers = FullAtmCov_fromLayers + S{1,kLayer}; %% S{1,klayer} here is a zenith
% end

% plot(diag(FullAtmCov))
% hold on
% plot(diag(FullAtmCov_fromLayers),'r--')