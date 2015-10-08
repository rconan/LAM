%> @file tool_compute_AR_model_via_covariance_bfgs_fitting.m
%> @brief A script that computes the Zernike temporal correlation functions using the analytical angular covariance routines.
%> @author Mikhail Konnik
%> @date   24 April 2015
%> @section spatangbased Spatio-angular-based temporal auto-correlation functions
%> 
%> Use Broyden-Fletcher-Goldfarb-Shanno (BFGS) Method to fit parameters of
%> AR1 or AR2 models (see Chapter 7, section 7.6 of Practical Optimization:
%> Algorithms and Engineering Applications, A. Antoniou and W.-S. Lu, 2007.
%> 
%> Note: gradient is a numerical approximation, not a closed form equation.
%> 
%> Objective         ::  Find the coefficients of a perturbation model that best fit the Npts first
%>                        steps of the (de)correlation curve correl
%> Comments:             Based on identical procedures developped at ONERA
%>  
%> INPUT VARS
%> correl            :: The correlation sequence with at least Npts length, normalised
%> t                 :: The temporal vector of points were the correlation is evaluated
%> Npts              :: Number of points to consider in the fitting
%> ord               :: Order of the model to fit
%> 
%> OUTPUT VARS
%> coeff             :: Coefficients of the fitting model
%> 
%> Created by         :: C. Correia
%> Creation date      :: 06/05/2009
%======================================================================
%> @param atm =  atmosphere object.
%> @param tel =  telescope object.
%> @param zern = Zernike decomposition of the atmosphere.
%> @param params = the structure that contains other parameters like the number of points to consider (Npts) and order of the AR model.
%> @retval ARn  = the covariance matrix.
%> @retval ZPSD = the matrix of zernike Angular Covariance.
%> @retval correl_all = correlation matrix.
% ====================================================================== 
function [ARn, ZPSD, correl_all] = tool_compute_AR_model_via_covariance_bfgs_fitting(atm,tel,zern,params)
ord = params.order; %% order of the AR model to fit into Zernike
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

number_of_points_sampling = round(0.5*1/Ts);
min_freq = 1e-5; %% smallest frequency we are interested in.
max_freq = 2*number_of_points_sampling; %% largest frequency we are interested in.
nu = logspace(log10(min_freq), log10(max_freq), number_of_points_sampling); %% frequency spectrum for sampling in the freq. domain.

FRange =  [min_freq round(params.freq_part*max_freq)]; %Frequency Range

f1 = linspace(nu(1), max(FRange), params.fN)'; %% nu(1) here is the smallest frequency we are interested in. 
f  = f1(2:end)';

ARn = struct; %%% create the structure array for the function's answer


    switch params.flag_correlation_method

        case 'corr_via_PSD'
        %%% This Zernike Miltipe spectral calculation is TERRIBLY inefficient and may have numerical problems (it uses some integration methods)
        ZPSD = zernikeStats.multipleSpectra(nu,atm,zern,tel);
        correl_all = [];

        case 'corr_via_spatio_angular'
        S = compute_spatio_angular_temporal_correlation(atm,tel,zern,Npts);
        ZPSD = [];

    end %%% switch flag_correlation_method




%% %%%%%%%%  Begin fitting the AR model for each atm layer and each mode %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for kLayer = 1:atm.nLayer

    mydisp(params.flag_messages, strcat('Processing atm.Layer#', num2str(kLayer) ));

    if (strcmp(params.flag_correlation_method, 'corr_via_spatio_angular') == 1)
        correl_all = compute_layer_correlation_spatio_angular(S, zern, kLayer, Npts);
    end

    for kmode=1:zern.nMode

        switch params.flag_correlation_method

            case 'corr_via_PSD'
                %----------------------------------------------------------------------
                % --- make fit and compute model temporal auto-correlation ---
                % 1) use params for a continuous time model in Laplace domain,
                % 2) extract magnitude,
                % 3) compute the auto-correlation by Wiener-Khinchine---
                %----------------------------------------------------------------------               
                PSDi = ZPSD(kmode,:,kLayer);
                % ---  Interpolation taking place ---
                sp = interp1(nu, PSDi, f1)'; % returns interpolated values of a 1-D function. Vector nu contains the sample 
                % points, PSDi contains the corresponding values, and vector f1 contains the coordinates of the query points.
                sp_full  = [0 sp(end:-1:2)  sp]; %% we flip the sp_full vector from left to right and take all elements except the last one (so is ...:2) stands for.
                Cphi_aoa = real(fftshift(fft(fftshift(sp_full))));        % autocorrelation thru Wiener-Khinchine theorem
                Cphi_aoa = Cphi_aoa/max(Cphi_aoa);
                correl   = Cphi_aoa(1,end/2+1:end-1); % taking only a half of the correlation function.

                correl_all = [correl_all; correl]; %% collect all covariance vectors into one for further analysis

            case 'corr_via_spatio_angular' %% use the Spatio-Angular Covariance method
                correl = correl_all(kmode,:);

        end %%% switch flag_correlation_method
           

%%% Issue 4: Convergence in BFGS depends heavily on the starting point x0
        x0 = 3*ones(1,ord); %%% set initial starting point for the BFGS algorithm. 3 seems to be OK.
%%% The AR model is computed in computeARmodel.m using BFGS algorithm, but
%%% the convergence may break if you choose wrong starting point x0.
%%% Besides, the starting point with minimal iterations seem to depend on the order of AR model  
%%% However, I changed this code into this (initially):
%%%          x0 = ones(1,ord); %%% set initial starting point for the BFGS algorithm.
%%% which worked well for AR(1), but it does not work for AR(2) and AR(3)!
%%% The BFGS takes ages to converge and it actually fails to do so (takes too many iterations).  
%%% 
%%% [Comment by Carlos Correia from Issue 4]This is a normal feature for a gradient-based non-linear minimisation
%%% method that gets trapped in local minima. We've come across that in the
%%% past and our approach was the following:  
%%% 
%%% 1) From a physical point of view, we want the model parameters to be
%%% meaningful. I.e., the decorrelation functions we get are in good
%%% agreement with theory despite the identified params not being the global minimum  
%%% 
%%% 2) We chose starting points that relate to the range of values we would
%%% normally find for the first tens of Zernike modes. We know we don't
%%% start off from afar, but have no certainty whether the convergence to a
%%% suitable point can be made. Our previous experience with the
%%% identification method tell us that it does, but there is no
%%% mathematical proof. In my opinion, the starting point ought to be a
%%% function of the Zernike radial order, for example.      


%% This calls for the BFGS optimisation algorithm - every time, for each layer!
        [coeff,resid,nIter] = bfgs(x0,1e-6,correl,f,Npts,ord);

        % minimization parameter substitution: p1 = coeff(1)^2, p2 = coeff(2)^2 to
        % enforce constraint [p1, p2] > 0.
        coeff = coeff.^2;
        res(kLayer,kmode) = resid;
               
switch ord  %%% depending on the ord (order of ARmodel), convert the coefficients into AR
    
    case 1  %% First-order autoregressive function
            z1 = exp(-coeff(1)*Ts);
                       
            Hz = zpk([], [z1], 1,Ts); % a normalization by the gain@0 may be needed here
            Hz = tf(Hz); % convert to tf
            pp = get(Hz,'den');
            pp = pp{1};
            
            A(kmode,kLayer) = -pp(2);
        
    case 2 %% Second-order autoregressive function
            z1 = exp(-coeff(1)*Ts);
            z2 = exp(-coeff(2)*Ts);
                       
            Hz = zpk([], [z1 z2], 1,Ts); % a normalization by the gain@0 may be needed here
            Hz = tf(Hz); % convert to tf
            pp = get(Hz,'den');
            pp = pp{1};
            
            A(kmode,kLayer) = -pp(2);
            B(kmode,kLayer) = -pp(3);
            

    case 3 %% Third-order autoregressive function
            z1 = exp(-coeff(1)*Ts);
            z2 = exp(-coeff(2)*Ts);
            z3 = exp(-coeff(3)*Ts);
            
            Hz = zpk([], [z1 z2 z3], 1,Ts); % a normalization by the gain@0 may be needed here
            Hz = tf(Hz); % convert to tf
            pp = get(Hz,'den');
            pp = pp{1};
            
            A(kmode,kLayer) = -pp(2);
            B(kmode,kLayer) = -pp(3);
            C(kmode,kLayer) = -pp(4);
            
    end %% for switch
        
    end %% for each zern.nMode
    
end %% for each atm.nLayer

%% Finally, assemble all the answers into the structure ARn for the output:
ARn.res = res;
    switch ord
        case 1            
            ARn.A = A;

        case 2    
            ARn.A = A;
            ARn.B = B;

        case 3
            ARn.A = A;
            ARn.B = B;
            ARn.C = C;

        otherwise
          fprintf('The ord value is wrong! \n')
    end %% for Switch 

 mydisp(params.flag_messages, strcat('The number of points Npts=', num2str(Npts) ));

end %%% for the ENTIRE function


%%
%%%%%%%%% ##### The rest of the function are auxiliary sub-functions necessary for the computations 
%%
function out = EvalCorrFunction(x,correl,f,Npts,order)  %% evaluating the value of the correlation function. 
%%% this function does not evaluate the correlation function, but compares the corr.fun of the estimated AR model
% % % with the correlation matrix computed from the theoretical PSD
    switch order
        case 1            
            syswind = tf(x^2, [1 x^2]);

        case 2    
            syswind= tf(x(1)^2*x(2)^2,  [1 x(1)^2+x(2)^2 x(1)^2*x(2)^2]);

        case 3
            syswind= zpk([],[-x(1)^2 -x(2)^2 -x(3)^2],x(1)^2*x(2)^2*x(3)^2);

        otherwise
          fprintf('The order value is wrong! \n');
          return;
    end

    [mag]   = bode(syswind,2*pi*f);
    mag     = reshape(mag, 1,length(mag));
    sp_full = [0 mag(length(mag):-1:1)  mag];
    Cphi    = fftshift(fft(fftshift(sp_full.^2)));     % autocorrelation thru Wiener-Khinchine theorem
    Cphi    = real(Cphi);
    Cphi    = Cphi/max(Cphi);
    Np      = length(Cphi);

    out = var(correl(1:Npts) - Cphi(floor(Np/2+1):floor(Np/2+1)+Npts-1));
end %% end for  function out = EvalCorrFunction(x,correl,f,Npts,order)



function out = gradient(x,delta,correl,f,Npts,order)
fnc = EvalCorrFunction(x,correl,f,Npts,order);
n = length(x);

%%% We construct the gradient of the function this way since fDelta will
%%% give you an overall delta, as a single number, instead of vector of
%%% deltas for each component.
for i = 1:n
    x(i) = x(i)+delta;
    fDelta = EvalCorrFunction(x,correl,f,Npts,order);
    out(i) = (fDelta-fnc)/delta;
    x(i) = x(i)-delta;
end
end  %% for Gradient function




%% This is where the computeARmodel actually ends and the BFGS optimization algorithm, which solves the fitting problem, starts
function [xs,fs,k] = bfgs(x0,epsi1,correl,f,Npts,order)
flag_messages = 0;
mydisp(flag_messages, 'Program bfgs.m')

n = length(x0);
I = eye(n);
k = 1;
xk = x0';
Sk = I;

fk = EvalCorrFunction(xk,correl,f,Npts,order);
gk = gradient(xk,1e-5,correl,f,Npts,order)';
dk = -Sk*gk;
ak = inex_lsearch(xk,dk,correl,f,Npts,order);
dtk = ak*dk;
xk_new = xk + dtk;
fk_new = EvalCorrFunction(xk_new,correl,f,Npts,order);
dfk = abs(fk - fk_new);
err = max(dfk,norm(dtk));

%%% Main part of the BFGS (iteration) starts below:
while err >= epsi1,
    gk_new = gradient(xk_new,1e-5,correl,f,Npts,order)';
    gmk = gk_new - gk;
    D = dtk'*gmk;
    if D <= 0,
        Sk = I;
    else
        sg = Sk*gmk;
        sw0 = (1+(gmk'*sg)/D)/D;
        sw1 = dtk*dtk';
        sw2 = sg*dtk';
        Sk = Sk + sw0*sw1 - (sw2'+sw2)/D;
    end
    fk = fk_new;
    gk = gk_new;
    xk = xk_new;
    dk = -Sk*gk;
    ak = inex_lsearch(xk,dk,correl,f,Npts,order);
    dtk = ak*dk;
    xk_new = xk + dtk;
    fk_new = EvalCorrFunction(xk_new,correl,f,Npts,order);
    dfk = abs(fk - fk_new);
    err = max(dfk,norm(dtk));
    k = k + 1;
end %% while err (main BFGS loop)

% mydisp(flag_messages, 'solution point:');
xs = xk_new;
% mydisp(flag_messages, 'objective function at the solution point:');

fs = EvalCorrFunction(xs,correl,f,Npts,order);
mydisp(flag_messages, strcat('number of iterations at convergence: ', num2str(k) ));

end


%%% Inexact line search for the BFGS optimization algorithm
function z = inex_lsearch(xk,s,correl,f,Npts,order)
k = 0;
m = 0;
tau = 0.1;
chi = 0.75;
rho = 0.1;
sigma = 0.1;
mhat = 400;
epsilon = 1e-10;
xk = xk(:);
s = s(:);
% evaluate given parameters:
% compute f0 and g0
f0 = EvalCorrFunction(xk,correl,f,Npts,order);
gk = gradient(xk,1e-5,correl,f,Npts,order)';
m = m+2;
deltaf0 = f0;

% step 2 Initialize line search
dk = s;
aL = 0;
aU = 1e99;
fL = f0;
dfL = gk'*dk;
if abs(dfL) > epsilon,
    a0 = -2*deltaf0/dfL;
else
    a0 = 1;
end
if ((a0 <= 1e-9)|(a0 > 1)),
    a0 = 1;
end

%step 3
while 1,
    deltak = a0*dk;
    f0 = EvalCorrFunction(xk+deltak,correl,f,Npts,order);    
    m = m + 1;
    %step 4
    if ((f0 > (fL + rho*(a0 - aL)*dfL)) & (abs(fL - f0) > epsilon) & (m < mhat))
        if (a0 < aU)
            aU = a0;
        end
        % compute a0hat using equation 7.65
        a0hat = aL + ((a0 - aL)^2*dfL)/(2*(fL - f0 + (a0 - aL)*dfL));
        a0Lhat = aL + tau*(aU - aL);
        if (a0hat < a0Lhat)
            a0hat = a0Lhat;
        end
        a0Uhat = aU - tau*(aU - aL);
        if (a0hat > a0Uhat)
            a0hat = a0Uhat;
        end
        a0 = a0hat;
    else
        gtemp = gradient(xk+a0*dk,1e-5,correl,f,Npts,order)';
        df0 = gtemp'*dk;
        m = m + 1;

        % step 6
        if (((df0 < sigma*dfL) & (abs(fL - f0) > epsilon) & (m < mhat) & (dfL ~= df0)))
            deltaa0 = (a0 - aL)*df0/(dfL - df0);
            if (deltaa0 <= 0)
                a0hat = 2*a0;
            else
                a0hat = a0 + deltaa0;
            end
            a0Uhat = a0 + chi*(aU - a0);
            if (a0hat > a0Uhat)
                a0hat = a0Uhat;
            end
            aL = a0;
            a0 = a0hat;
            fL = f0;
            dfL = df0;
        else
            break;
        end
    end
end % while 1
if a0 < 1e-5,
    z = 1e-5;
else
    z = a0;
end
end

% ======================================================================
%> @brief Turns diagnostic messages off or on depending on the flag
%> @param flag  = can be 0 (no messages) or 1 (messages will appear)
%> @param ss    = message string.
% ======================================================================
function mydisp(flag,ss)  %%% truns on or off the diagnostic messages via the __flag__ variable.
if (flag == 1)
    disp(ss);
end
end


% ======================================================================
%> @brief The function computes the temporal correlation matrix via spatio-angular method (like in AmmseLag).
%> @param atm =  atmosphere object.
%> @param tel =  telescope object.
%> @param zern = Zernike decomposition of the atmosphere.
%> @param Npts = Npts is Number of points to consider in the fitting.
%> @retval S = the covariance cell array (can be converted into the matrix).
% ======================================================================
function [S] = compute_spatio_angular_temporal_correlation(atm, tel, zern, Npts)
S = cell(Npts, atm.nLayer); % cell with the phase covar matrices

    for kLayer = 1:atm.nLayer
        refH = 1e3; %% reference altitude Height

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
            alpha        = 1/2*alpha*180/pi*3600; %% converting the alngle to arcseconds
            guidestar    = source('asterism',{[2,alpha*cougarConstants.arcsec2radian,0]}); %%% creating a fake guide star
            S{kr,kLayer} = phaseStats.zernikeAngularCovariance(zern,myatm,guidestar); % S contains the covariance matrices for different alpha (angles, where S{Npts,1}) is the covariance for zenith
        end %% for kr = 1:Npts 

    end %% for kLayer

end %%% for function compute_spatio_angular_temporal_correlation(atm,Npts)



% ======================================================================
%> @brief The function converts the cell array into a matrix of covariances for a specified atmospheric layer.
%> @param S =  the covariance cell array (can be converted into the matrix).
%> @param zern = Zernike decomposition of the atmosphere.
%> @param kLayer = Current atmospheric layer.
%> @param Npts = Npts is Number of points to consider in the fitting.
%> @retval correl_all = the covariance matrix for the specified (kLayer) atmospheric layer.
% ======================================================================
function [correl_all] = compute_layer_correlation_spatio_angular(S, zern, kLayer, Npts)

        number_of_modes = length(zern.n);

        correl_all = zeros(number_of_modes, Npts);
        for ii=1:Npts
            correl_all(:,ii) = diag(S{ii,kLayer}); 
        end
        for ii=1:number_of_modes
            correl_all(ii,:) = correl_all(ii,:)./(max(correl_all(ii,:)));
        end

end %% function [correl_all] = compute_layer_correlation_spatio_angular