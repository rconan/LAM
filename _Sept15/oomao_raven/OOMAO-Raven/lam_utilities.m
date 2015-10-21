%> @class lam_utilities 
%> @brief A class with the numerical methods for OOMAO and the tools developed for RAVEN.
%> @author Mikhail Konnik
%> @date   24 April 2015
classdef lam_utilities

    methods (Static)

        % ======================================================================
        %> @brief This simple demo is intended to show how to compute a simple atmosphere covariance and its spectrum.
        %> @author Mikhail Konnik
        %> @date   24 April 2015
        %> 
        %> @par Demonstation of Zernike covariance
        %> Example:% @n
        %> atm = atmosphere(photometry.V,0.15,30); % @n
        %> tel = telescope(10);% @n
        %> modes = 1:15;% @n
        %> figure% @n
        %> spy(zernikeStats.covariance(modes,atm))% @n
        % ======================================================================    
        function [] = demo_zernike_cov()            
            L0 = 40;
            r0 = 0.19;
            max_mode = 15;
            max_freq = 100;

            atm = atmosphere(photometry.V,r0,L0);
            tel = telescope(8);
            modes = 2:max_mode;
            figure, spy(zernikeStats.covariance(modes,atm))

            f = [0:0.1:max_freq];
            out = phaseStats.spectrum(f,atm);
            figure, plot(out)
        end %%% END OF FUNCTION demo_zernike_cov
    
        


%% BEGIN functions used by spat_ang_based_temporal_acorr method

        %======================================================================
        %> @brief Computes the spatio-angular-based temporal auto-correlation functions for the multi-layered case. Ammse = < ||phi_k+1, phi_k||^2 > * inv { < ||phi_k, phi_k||^2>}
        %> First, it computes the Zernike temporal correlation functions using the analytical angular covariance routines.
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
        end %%% END OF FUNCTION spat_ang_based_temporal_acorr
        
        
        



        %======================================================================
        %> @brief Find the coefficients of a perturbation model that best fit the Npts first
        %> steps of the (de)correlation curve correl. The algorithm uses Broyden-Fletcher-Goldfarb-Shanno (BFGS) Method to fit parameters of
        %> AR1 or AR2 models 
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
                S = lam_utilities.compute_spatio_angular_temporal_correlation(atm,tel,zern,Npts);
                ZPSD = [];

            end %%% switch flag_correlation_method




            %% Begin fitting the AR model for each atm layer and each mode
    for kLayer = 1:atm.nLayer
        lam_utilities.mydisp(params.flag_messages, strcat('Processing atm.Layer#', num2str(kLayer) ));

        if (strcmp(params.flag_correlation_method, 'corr_via_spatio_angular') == 1)
            correl_all = lam_utilities.compute_layer_correlation_spatio_angular(S, zern, kLayer, Npts);
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
                    [coeff,resid,nIter] = lam_utilities.bfgs(x0,1e-6,correl,f,Npts,ord);

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

            end %%% for the ENTIRE function


            %%
            %%%%%%%%% ##### The rest of the function are auxiliary sub-functions necessary for the computations 
            %%
            function out = EvalCorrFunction(x,correl,f,Npts,order)  %% evaluating the value of the correlation function. 

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
            fnc = lam_utilities.EvalCorrFunction(x,correl,f,Npts,order);
            n = length(x);

            %%% We construct the gradient of the function this way since fDelta will
            %%% give you an overall delta, as a single number, instead of vector of
            %%% deltas for each component.
            for i = 1:n
                x(i) = x(i)+delta;
                fDelta = lam_utilities.EvalCorrFunction(x,correl,f,Npts,order);
                out(i) = (fDelta-fnc)/delta;
                x(i) = x(i)-delta;
            end
            end  %% for Gradient function




            %% This is where the computeARmodel actually ends and the BFGS optimization algorithm, which solves the fitting problem, starts
            function [xs,fs,k] = bfgs(x0,epsi1,correl,f,Npts,order)
            flag_messages = 0;
            lam_utilities.mydisp(flag_messages, 'Program bfgs.m')

            n = length(x0);
            I = eye(n);
            k = 1;
            xk = x0';
            Sk = I;

            fk = lam_utilities.EvalCorrFunction(xk,correl,f,Npts,order);
            gk = lam_utilities.gradient(xk,1e-5,correl,f,Npts,order)';
            dk = -Sk*gk;
            ak = lam_utilities.inex_lsearch(xk,dk,correl,f,Npts,order);
            dtk = ak*dk;
            xk_new = xk + dtk;
            fk_new = lam_utilities.EvalCorrFunction(xk_new,correl,f,Npts,order);
            dfk = abs(fk - fk_new);
            err = max(dfk,norm(dtk));

            %%% Main part of the BFGS (iteration) starts below:
            while err >= epsi1,
                gk_new = lam_utilities.gradient(xk_new,1e-5,correl,f,Npts,order)';
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
                ak = lam_utilities.inex_lsearch(xk,dk,correl,f,Npts,order);
                dtk = ak*dk;
                xk_new = xk + dtk;
                fk_new = lam_utilities.EvalCorrFunction(xk_new,correl,f,Npts,order);
                dfk = abs(fk - fk_new);
                err = max(dfk,norm(dtk));
                k = k + 1;
            end %% while err (main BFGS loop)

            % lam_utilities.mydisp(flag_messages, 'solution point:');
            xs = xk_new;
            % lam_utilities.mydisp(flag_messages, 'objective function at the solution point:');

            fs = lam_utilities.EvalCorrFunction(xs,correl,f,Npts,order);
            lam_utilities.mydisp(flag_messages, strcat('number of iterations at convergence: ', num2str(k) ));

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
            f0 = lam_utilities.EvalCorrFunction(xk,correl,f,Npts,order);
            gk = lam_utilities.gradient(xk,1e-5,correl,f,Npts,order)';
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
                f0 = lam_utilities.EvalCorrFunction(xk+deltak,correl,f,Npts,order);    
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
                    gtemp = lam_utilities.gradient(xk+a0*dk,1e-5,correl,f,Npts,order)';
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

%% END functions used by spat_ang_based_temporal_acorr method



%% BEGIN functions used by AnalyticalSmallFootprintExpansion method

            %======================================================================
            %> @brief A script that computes the Zernike temporal correlation functions using the analytical angular covariance routines.
            %> @author Carlos Correia, improved by Mikhail Konnik
            %> @date   June 2012, revised May 2015
            %> 
            %> @par Spatio-angular-based temporal auto-correlation functions
            %> The function @b AnalyticalSmallFootprintExpansion computes the modal projection
            %> of Zernike polynomials analytically using Noll's indexing convention onto a
            %> smaller telescope pupil, displaced by \f$\Delta x\f$ and \f$\Delta y\f$ and (possibly)
            %> rotated. In our case, the pupil is not rotated at all.
            %> 
            %> The major benefit of the @b AnalyticalSmallFootprintExpansion function is the
            %> speed: since there is no numerical integration, the expansion is very fast,
            %> typically 10x faster than the @b tel.footprintProjection.
            %> 
            %> The @b AnalyticalSmallFootprintExpansion has the following major steps:
            %> 
            %> - the Zernike ordering is changed from Noll conventions to ANSI, because the
            %> @b TransformC function uses ANSI:
            %> 
            %> 	A2N   =  ANSI2Noll(nMode); %% converts the ANSI order of the Zernike modes into Noll convention 
            %> 
            %> - calculating telescope's diameter at an altitude of kLayer of the atmosphere:
            %> 
            %> 	D0 = tel.diameterAt(altitudes)*1000,  where:  out = obj.D + 2.*height.*tan(obj.fieldOfView/2);
            %> 
            %> here diameterAt is from telescopeAbstract, that takes into account the Field of
            %> View (FOV) and height of the layer. Then the D0 is converted to millimeters
            %> (because @b TransformC computes in these units);
            %> 
            %> - Then for each atmospheric layer, and for each mode in that layer, we compute
            %> the projection matrix that accounts for displacement and rotation via @b
            %> TransformC;
            %> 
            %> - The matrix of projections (Proj) is computed, but for the ANSI conventions -
            %> we have to sort the columns of the projection matrix back.
            %> 
            %> The major portion of the computation time is spent in @b TransformC (about 60\%
            %> of time), the rest is for cell2mat conversion (about 8\%)  and displacement
            %> computations (7\%).
            %> @param zernModeMaxProj   = maximum order of Zernike function [integer number]
            %> @param tel               = telescope object.
            %> @param gs                = the guide-star asterism
            %> @param atm               = atmosphere
            %> @retval Proj             = Projection matrix with [nGs x Nz] by [Nz x Nlayers] dimensions
            %> @retval Projcell         = Cell version of Proj
            % ======================================================================
            function [Proj Projcell] = tool_analytical_small_footprint_expansion(zernModeMaxProj,tel,gs,atm)
            nMode = zernModeMaxProj-1;          %% number of Zernike modes.

                A2N   =  lam_utilities.ANSI2Noll(nMode); %% converts the ANSI order of the Zernike modes into Noll convention

            altitudes = [atm.layer.altitude];    %% atmospheric layers altitudes
            nGs = length(gs);                    %% number of Guide Stars
            re_sorting_vector = zeros(1,nMode+1);%% preallocation of a re-sorting vector for the Zernike modes.

            %%% diameterAt is from telescopeAbstract.m, that takes into account the Field
            %%% of View (FOV) and height of the layer:  out = obj.D + 2.*height.*tan(obj.fieldOfView/2);
            Dh = tel.D*1000;

            Proj = zeros((nMode+1), (nMode+1),  atm.nLayer, nGs);
            Projcell = cell(atm.nLayer, nGs);  %%% Cell array of size NumberOfNGS x number of atmospheric layers.

            for kLayer = 1:atm.nLayer  %% iterating for each atmosphere layer atm.nLayer;
                for kGuidestar = 1:nGs %% iterating through Guide Stars

                    D0 = tel.diameterAt(altitudes(kLayer))*1000; %% telescope's diameter at an altitude of kLayer, changes with the altitude!

                    src = gs(kGuidestar); %% taking the properties of Guide Star for each one.
                    delta = altitudes(kLayer)*tan(src.zenith)*[cos(src.azimuth),sin(src.azimuth)]; %displacement due to altitude and FoV of a telescope

                    %% Computing the Delta_x and Delta_y displacement of the pupil due to asterism        
                    tx = delta(1)*1000; %% to convert the tx into millimeters, which is standard unit in TransformC function
                    ty = delta(2)*1000; %% to convert the ty into millimeters, which is standard unit in TransformC function

                    for kMode = 2:nMode+1
                        vec0 = zeros(1,nMode+1);
                        vec0(kMode) = 1; %% considering only current mode, thus 1 in vec0 corresponds to kMode (current Mode)

                        %% Transforming the Zernike coeffs for displaced telescope pupil
                        % % % function C2 = TransformC (C1,dia2,tx,ty,thetaR) returns transformed Zernike coefficient set, C2,
                        % % % from the original set, C1,   !!!____with the pupil diameter___ in mm as the first vector element!!!.
                        C21 = lam_utilities.TransformC( [D0 (vec0(1:kMode)) ], Dh, tx, ty, 0); %% this function is from the article 

                        C2x = zeros(nMode+1,1);
                        C2x(1:kMode) = C21(2:kMode+1); %% take all elements from C21 since TransformC(1) is a pupil diameter
                        C2x_converted = A2N*C2x;

                        re_sorting_vector(kMode) = find(C2x_converted, 1, 'last');  %% finds the last index of non-zero entry in the vector converted_tmp; 
                        Proj(:,kMode, kLayer, kGuidestar) = C2x_converted;
                    end %% for kMode

                    Proj(1,1, kLayer, kGuidestar) = 1; %%% accounting for the piston mode

                    [~, resort_indices] = sort(re_sorting_vector); %% we need a vector of resort_indices that will recover the correct order in ProjectionMatrix below.
                    Projcell{kLayer, kGuidestar} = Proj(:,resort_indices,kLayer,kGuidestar); %% Store the projected coefficients into Cell array

                end %%     for kGuidestar = 1:nGs %% iterating through Guide Stars    
            end %% for kLayer = 1:atm.nLayer  %% iterating for each atmosphere layer nL = atm.nLayer;

            Proj = cell2mat(Projcell); %% converting the Projections cell back to a matrix.

            end %%% function [Proj Projcell] = AnalyticalSmallFootprintExpansion(zernModeMaxProj,tel,gs,atm)



            %======================================================================
            %> @brief Convert ANSI ordering to Noll ordering.
            %> @author Carlos Correia
            %> @date   18 June 2012
            %> 
            %> @par ANSI2Noll - Convert ANSI ordering to Noll ordering
            %> There are many conventions on the Zernike expansions, which differ in the 
            %> ordering of modes.
            %> 
            %> The Zernike ordering scheme used in the telescope optics domain (e.g. ZEMAX 
            %> standard Zernike coefficients), based on @b Noll's @b concept, is shown below 
            %> for the first 21 terms. Here, the ordering number j is determined by ordering 
            %> polynomial with lower radial order first, and for given radial order with odd 
            %> number for sine function and even number for the cosine.
            %> 
            %> @image html zernike_noll.png
            %> 
            %> @b Noll's @b scheme also has a vertical expansion symmetry similar to Wyant's, 
            %> but since it directly uses radial order n, it unfolds somewhat differently. So 
            %> for n=6, the term for m=0 is secondary spherical aberration (j=22), for m=2 it 
            %> is tertiary astigmatism (j=23 and 24 for the sine and cosine function, 
            %> respectively), for m=4 it is secondary quadrafoil (j=25 and 26) and for m=6 the 
            %> hexafoil (j=27 and 28).
            %> 
            %> Zernike scheme commonly used in ophthalmology - the @b ANSI @b standard @b 
            %> Zernike @b expansion - is different from those used in assessing telescope 
            %> aberrations. The ANSI standard single-index scheme (Zernike pyramid) routinely 
            %> used in ophthalmology for evaluating eye aberrations. It is graphically 
            %> presented as a pyramid resulting from Zernike term expansion as a function of 
            %> radial order n and angular frequency m, with the latter being numerically 
            %> positive for cosine function of \f$\theta\f$, and negative for the sine 
            %> function 
            %> 
            %> @image html zernike_pyramide.png
            %> 
            %> Zernike terms expansion pyramid is a function of term's radial degree (or order) 
            %> n and azimuthal frequency m. It is the basis for classifying aberrations as 
            %> lower \f$(n\leq 2)\f$ and higher-order \f$(n>2)\f$ in ophthalmology. Shown are 
            %> the top 20 terms.
            %> @param nModes     = The number of modes
            %> @retval ANSIt	 = The nMode x nModes matrix 
            % ======================================================================
            function ANSIt = ANSI2Noll(nModes)
            [Nj, Mj] = lam_utilities.make_zernike_index(nModes); %Output the radial order and azimuthal frequency indexes for up to a given number of Nmodes
            ANSIt = sparse(numel(Nj), numel(Nj));

            for n = 0:max(Nj)
                Nollj = find(Nj == n);
                nn = numel(Nollj);

                if lam_utilities.iseven(nn)

                    if iseven(Nollj(min(numel(Nollj),2)))
                        coss = Nollj(2:2:end);
                        sins = Nollj(3:2:end);
                    else
                        sins = Nollj(2:2:end);
                        coss = Nollj(3:2:end);
                    end

                    ANSIj = [sins(end:-1:1) Nollj(1) coss];

                else

                    if iseven(Nollj(min(numel(Nollj),2)))
                        coss = Nollj(2:2:end);
                        sins = Nollj(3:2:end);
                    else
                        sins = Nollj(2:2:end);
                        coss = Nollj(3:2:end);
                    end

                    ANSIj = [sins(end:-1:1) Nollj(1) coss];

                end

                for rr = 1:numel(Nollj)
                    ANSIt(ANSIj(rr), Nollj(rr)) = 1;
                end

            end

            end %% function ANSIt = ANSI2Noll(nModes)



            %======================================================================
            %> @brief   Output the radial order and azimuthal frequency indexes for up to a given number of Nmodes
            %> @author Carlos Correia (Improved by Mike Konnik)
            %> @date   15 March 2004 (revised 17 May 2015)
            %> @param Nmodes     = The number of modes
            %> @retval Nj	 = The radial order
            %> @retval Mj	 = The azimuthal frequency
            %> 
            %> @par Zernike radial polynomials approximation
            %> The radial polynomials \f$R^m_n\f$ are defined as: @n
            %> \f$ R^m_n(\rho) = \! \sum_{k=0}^{(n-m)/2} \!\!\! \frac{(-1)^k\,(n-k)!}{k!\,((n+m)/2-k)!\,((n-m)/2-k)!} \;\rho^{n-2\,k} \f$@n
            %> for \f$(n - m)\f$ even. If \f$(n - m)\f$ is odd then \f$R^m_n(\rho) = 0\f$.
            %> Here:
            %> - index \f$n = 0,1,2,\dots\f$ is called the @b degree of the function or polynomial;
            %> - index \f$m = -n \dots +n\f$ is called the @b order of the function or polynomial;
            %> - \f$\rho\f$ is the radial distance \f$0 \le\rho\le 1\f$
            %> Indexes $m$ and $n$ are non-negative integers with $n\geq m$. 
            %>
            %> For \f$m = 0\f$, the even definition is used, which reduces to \f$R_n^0(\rho)\f$.
            %> 
            %> As one can see from Eq. above, Zernike radial polynomials are usually defined using their series representation as a finite sum of powers of \f$\rho^2\f$. Although the expression in Eq.~\ref{eq:zernikepoly} looks complicated, the \f$R^m_n(\rho)\f$ are just simple polynomials: @n
            %> \f$R^0_0(\rho) = 1\f$ @n
            %> \f$R^1_1(\rho) = \rho\f$ @n
            %> \f$R^0_2(\rho) = 2\rho^2 - 1\f$ @n
            %> \f$R^2_2(\rho) = \rho^2\f$ @n
            %> \f$R^1_3(\rho) = 3\rho^3 - 2\rho\f$ @n
            %> \f$R^3_3(\rho) = \rho^3\f$ @n
            %> \f$R^0_4(\rho) = 6\rho^4 - 6\rho^2 + 1\f$ @n
            %> \f$R^2_4(\rho) = 4\rho^4 - 3\rho^2\f$ @n
            %> \f$R^4_4(\rho) = \rho^4\f$ @n
            %> 
            %> @note This implementation of Zernike WORKS with any MATLAB versions and does not contain stupid stunt with unreadable vector notations. There is actually a widely used function tool_zernike_poly written by Paul Fricker. Unfortunately it is practically useless due to its implementation: Mr. Fricker uses a lot of scout boy's stunt aka ``am-a-thirteen-year-old-hacker-and-gonna-write-a-code-that-nobody-can-read''. Such a style makes the code unserviceable and potentially dangerous (there were reported some problems with norm).
            %> 
            %> @par Normalisation of Zernike polynomials
            %> The normalization of Zernike functions is given by\cite{thibos2002standards}: @n
            %> \f$N_n^m = \sqrt{\frac{2(n+1)}{1+\delta_{m0}}}\f$ @n
            %> where \f$\delta_{m0}\f$ is the Kronecker delta function (i.e., \f$\delta_{m0} =1\f$ for \f$m = 0\f$, and \f$\delta_{m0} = 0\f$ for \f$m \neq 0\f$). Note that the value of \f$n\f$ is a positive integer or zero. For a given \f$n\f$, \f$m\f$ can only take on values \f$-n, -n + 2, -n +4, \dots n\f$.
            %> @par Converting Degrees, Orders and Mode Number of Zernike polynomials.
            %> The Zernike Polynomial definitions used are taken from: Thibos, L., Applegate, R.A., Schweigerling, J.T., Webb, R., VSIA Standards Taskforce Members, "Standards for Reporting the Optical Aberrations of Eyes", OSA Trends in Optics and Photonics Vol. 35, Vision Science and its Applications, Lakshminarayanan,V. (ed) (Optical Society of America, Washington, DC 2000), pp: 232-244.
            %>  
            %> Occasionally, a single indexing scheme is useful  for describing Zernike expansion coefficients. Since the polynomials actually depend on two parameters, \f$n\f$ and \f$m\f$, ordering of a single indexing scheme is arbitrary. To avoid confusion, a standard single indexing scheme should be used, and this scheme should only be used for bar plots of expansion coefficients.
            %> 
            %> To obtain the single index, \f$j\f$, it is convenient to lay out the polynomials in a pyramid with row number n and column number \f$m\f$. The single index, j, starts at the top of the pyramid and steps down from left to right. To convert between \f$j\f$ and the values of \f$n\f$ and \f$m\f$, the following relationships can be used:@n
            %> mode number:  \f$z = \frac{n(n+2)+m}{2}\f$ @n
            %> degree: \f$n = round[\frac{-3+\sqrt{9+8\cdot z}}{2} ] \f$ @n
            %> order:  \f$m = 2\cdot z - n(n+2) \f$ @n
            %> 
            %> @par Comments in the code
            %> In the process of replacing and testing the code, I've noticed that the new code (which is more compact and neat) is noticeably slower than the old, hard-coded one.
            %> 
            %> For maxRadialDegreeProj = 4 we have: @n
            %> Old Code (hardcoded one)
            %> - AnalyticalSmallFootprintExpansion takes 0.1665 seconds;
            %> - AnalyticalSmallFootprintExpansion takes 0.1340 seconds;
            %> - AnalyticalSmallFootprintExpansion takes 0.1203 seconds;
            %> 
            %> New Code (uses zernike class)@n
            %> - AnalyticalSmallFootprintExpansion takes 0.2557 seconds;
            %> - AnalyticalSmallFootprintExpansion takes 0.2276 seconds;
            %> - AnalyticalSmallFootprintExpansion takes 0.2311 seconds;
            %> 
            %> Not much, but noticeable.
            %> I will leave both: I'm going to comment the hardcoded code and leave new (that uses Zernike class). [Mikhail Konnik]
            % ======================================================================
            function [Nj, Mj] = make_zernike_index(Nmodes)
            %%%% The code below is hard-coded, but seems to be at least 2x faster than the New Code below    
            % % % Nj = 0;
            % % % n = 1;
            % % % while size(Nj,2) < Nmodes +1
            % % %     Nj = [Nj n*ones(1,n+1)];
            % % %     n = n+1;
            % % % end
            % % % 
            % % % even = [0 2 2 4 4 6 6 8 8 10 10 12 12 14 14 16 16 18 18 20 20 22 22 24 24 26 26 28 28 30 30 32 32 34 34 36 36 38 38 40 40 42 42 44 44 46 46 48 48 50 50 52 52 54 54 56 56 58 58 60 60 ];
            % % % odd  = [1 1 3 3 5 5 7 7 9 9 11 11 13 13 15 15 17 17 19 19 21 21 23 23 25 25 27 27 29 29 31 31 33 33 35 35 37 37 39 39 41 41 43 43 45 45 47 47 49 49 51 51 53 53 55 55 57 57 59 59 61 61];
            % % % 
            % % % Mj = 0;
            % % % isodd = 1;
            % % % n = 1;
            % % % while size(Mj,2) < Nmodes +1
            % % %     if isodd
            % % %         Mj = [Mj odd(1:n+1)];
            % % %         isodd = 0;
            % % %     else
            % % %         Mj = [Mj even(1:n+1)];
            % % %         isodd = 1;
            % % %     end
            % % %     n = n+1;
            % % % end

            %%%% the New Code below
            zern = zernike(1:Nmodes);
            maxRadialDegree = max(zern.n);
            zern = zernike(1:zernike.nModeFromRadialOrder(maxRadialDegree))
            Nj = zern.n;
            Mj = zern.m;
            end %%% function [Nj, Mj] = make_zernike_index(Nmodes)



            %======================================================================
            %> @brief  Returns true if the input is an even number. Syntax: tf = iseven(x)
            %> @author Dan Kominsky, Copyright 2012  Prime Photonics, LC.
            %> @date   2012
            %> @param tf = The number of modes
            %> @retval x	 = numeric (real) input of any number of dimensions and sizes
            %> @retval Mj	 =  true for each element that is an even whole number
            % ======================================================================
            function tf = iseven(x)
              if ~isreal(x)
                error('iseven:badinput','iseven requires real inputs');
              else
                tf = mod(x,2)==0;
              end

            end %% function tf = iseven(x)

            
            
%======================================================================
%> @brief TransformC returns transformed Zernike coefficient set, C2, from the original set, C1, 
%> both in standard ANSI order, with the pupil diameter in mm as the first term. Scaling and translation is performed first and then rotation.
%> @author Linda Lundstrom and Peter Unsbo
%> @date   March 2007 (from the original article)
%> 
%> @par TransformC - transforming Zernike coeffs for shifted and rotated pupils
%> The @b TransformC function was proposed in the article ``Transformation of 
%> Zernike coefficients: scaled, translated, and rotated wavefronts with circular 
%> and elliptical pupils'' by Linda Lundstrom and Peter Unsbo, published in 
%> Vol. 24, No. 3, March 2007,  J. Opt. Soc. Am. A (pages 569-577).
%> 
%> The algorithm presents the means to transform Zernike coefficients analytically 
%> with regard to concentric scaling, translation of pupil center, and rotation. 
%> The transformations are described both for circular and elliptical pupils. 
%> 
%> TransformC returns transformed Zernike coefficient set, C2, from the original 
%> set, C1, both in standard ANSI order, with the pupil diameter in mm as the first term. 
%> Scaling and translation is performed first and then rotation.
%> @param C1     = the original set of Zernike coefficients [vector], with the pupil diameter in mm as the first term.
%> @param dia2	 = new (desired) pupil diameter [mm] for concentric scaling
%> @param tx	 = Cartesian translation coordinates [mm]
%> @param ty	 = Cartesian translation coordinates [mm]
%> @param thetaR = angle of rotation [degrees]
%> @retval C2	 = transformed Zernike coefficient set in standard ANSI order with the pupil diameter in mm as the first term.
% ======================================================================
function C2 = TransformC (C1,dia2,tx,ty,thetaR)
% ''TransformC'' returns transformed Zernike coefficient set, C2, from the original set, C1,
% both in standard ANSI order, with the pupil diameter in mm as the first term.
% dia2 - new pupil diameter [mm]
% tx, ty - Cartesian translation coordinates [mm]
% thetaR - angle of rotation [degrees]
% Scaling and translation is performed first and then rotation.

dia1=C1(1); % Original pupil diameter
C1=C1(2:end);
etaS=dia2/dia1; % Scaling factor
etaT=2*sqrt(tx^2+ty^2) / dia1; % Translation in Cartesian coordinates
thetaT=atan2(ty, tx);
thetaR=thetaR*pi/180; % Rotation in radians
jnm=length(C1)-1; 

nmax=ceil((-3+sqrt(9+8*jnm))/2); 
jmax=nmax*(nmax+3)/2;

S=zeros(jmax+1,1); S(1:length(C1))=C1; C1=S; clear S
P=zeros(jmax+1); % Matrix P transforms from standard to Campbell order
N=zeros(jmax+1); % Matrix N contains the normalization coefficients
R=zeros(jmax+1); % Matrix R is the coefficients of the radial polynomials

CC1=zeros(jmax+1,1); % CC1 is a complex representation of C1
counter=1;

for m=-nmax:nmax % Meridional indexes
    for n=abs(m):2:nmax % Radial indexes
        jnm=(m+n*(n+2))/2;
        P(counter,jnm+1)=1;
        N(counter,counter)=sqrt(n+1);
        for s=0:(n-abs(m))/2
            R(counter-s,counter)=(-1)^s*factorial(n-s)/(factorial(s)*factorial((n+m)/2-s)*factorial((n-m)/2-s));
        end
        if m<0, CC1(jnm+1)=(C1((-m+n*(n+2))/2+1)+1i*C1(jnm+1))/sqrt(2);
        elseif m==0, CC1(jnm+1)=C1(jnm+1);
        else
            CC1(jnm+1)=(C1(jnm+1)-1i*C1((-m+n*(n+2))/2+1))/sqrt(2);
        end
        counter=counter+1;
    end
end

ETA=[]; % Coordinate-transfer matrix

for m=-nmax:nmax
    for n=abs(m):2:nmax
        ETA=[ETA P*(lam_utilities.transform(n,m,jmax,etaS,etaT,thetaT,thetaR))];
    end
end

C=inv(P)*inv(N)*inv(R)*ETA*R*N*P;

CC2=C*CC1;
C2=zeros(jmax+1,1); % C2 is formed from the complex Zernike coefficients, CC2

for m=-nmax:nmax
    for n=abs(m):2:nmax
        jnm=(m+n*(n+2))/2;
        if m<0, C2(jnm+1)=imag(CC2(jnm+1)-CC2((-m+n*(n+2))/2+1))/sqrt(2);
        elseif m==0, C2(jnm+1)=real(CC2(jnm+1));
        else
            C2(jnm+1)=real(CC2(jnm+1)+CC2((-m+n*(n+2))/2+1))/sqrt(2);
        end 
    end 
end
C2=[dia2;C2];
end %% function C2 = TransformC (C1,dia2,tx,ty,thetaR)




%======================================================================
%> @brief The sub-function @b transform returns 
%> @param n      = Meridional indexes of Zernike coefficients
%> @param m      = Radial indexes of Zernike coefficients
%> @param jmax	 = maximum Meridional index
%> @param etaS	 = Scaling factor, etaS=dia2/dia1; 
%> @param etaT	 = Translation in Cartesian coordinates
%> @param thetaT = angle of rotation [degrees]
%> @param thetaR = angle of rotation [degrees]
%> @retval Eta	 =  The matrix \f$\eta\f$ is formed successively, column by column, with the function @b transform, which
%> performs the calculations of Eqs. (20), (23), and (25) from the article ``Transformation of 
%> Zernike coefficients: scaled, translated, and rotated wavefronts with circular 
%> and elliptical pupils'' by Linda Lundstrom and Peter Unsbo, published in 
%> Vol. 24, No. 3, March 2007,  J. Opt. Soc. Am. A (pages 569-577). The
%> implemented code uses complex Zernike coefficients for
%> the calculations and includes Eqs. (7) and (8) to convert
%> between complex and real coefficients. The matrix  \f$\eta \f$ will be a diagonal matrix with each element equal
%> to \f$\eta_s^n \f$, where \f$n\f$ is the exponent of the corresponding \f$\rho\f$ term in \f$ \langle \rho M|\f$.
% ======================================================================
function Eta=transform(n,m,jmax,etaS,etaT,thetaT,thetaR)
% Returns coefficients for transforming a ro^n * exp(i*m*theta)-term (from Eq. 20) into '-terms
Eta=zeros(jmax+1,1);

for p=0:(n+m)/2
    for q=0:(n-m)/2
        nnew=n-p-q; 
        mnew=m-p+q;
        jnm=(mnew+nnew*(nnew+2))/2;
        Eta(floor(jnm+1))=Eta(floor(jnm+1))+nchoosek((n+m)/2,p)*nchoosek((n-m)/2,q)...
            * etaS^(n-p-q)*etaT^(p+q)*exp(1i*((p-q)*(thetaT-thetaR)+m*thetaR));
        %%% nchoosek returns the binomial coefficient, defined as n!/((n–k)! k!). 
    end 
end
end %% function Eta=transform(n,m,jmax,etaS,etaT,thetaT,thetaR)
            %% END   functions used by AnalyticalSmallFootprintExpansion method

                end %%% for methods (Static)

            end %% for lam_utulities class