%**************************************************************************
% Use Broyden-Fletcher-Goldfarb-Shanno (BFGS) Method to fit parameters of
% AR1 or AR2 models (see Chapter 7, section 7.6 of Practical Optimization:
% Algorithms and Engineering Applications, A. Antoniou and W.-S. Lu, 2007.
%
% Note: gradient is a numerical approximation, not a closed form equation.
% *************************************************************************

function [A B C res, ZPSD] = computeARmodel(nu,atm,tel,zern,FRange,FitOrd,NptsFit)

%lambda      = atm.wavelength;
%R           = tel.D/2;
%V           = atm.layer.windSpeed;
Fmax        = max(FRange);
%Fmin        = min(FRange);
%L0          = atm.L0;
Nmodes      = zern.nMode;
%Nlayers     = atm.nLayer;
%theraW      = atm.layer.windDirection;
fN          = 1001; % a smaller number of points leads to wrong results in the decorrelation functions cc, Jul2013
f1          = linspace(0,Fmax,fN)';
f           = f1(2:end)';
dnu         = mean(diff(f1));
dt          = 1/(2*Fmax);
tmax        = 1/dnu;
t           = linspace(-tmax/2,tmax/2-dt, tmax/dt);
tvec        = t(1,end/2+1:end);
ord         = FitOrd;

Ts          = tel.samplingTime;
Tfit        = NptsFit*Ts;

fMaxAR      = 1/(2*Ts);
dtz         = 1/(2*fMaxAR);
z           = tf('z',Ts);
fz          = linspace(0,fMaxAR,fN*(fMaxAR/Fmax));
tz          = linspace(-tmax/2,tmax/2-dtz, tmax/dtz);
tz          = tz(1,end/2+1:end);

% Be careful with this, the ZPSD file is specific to a given atm
% filename = 'Zpsd_54modes_robustness.mat';
% if ~exist(filename,'file')
    ZPSD = zernikeStats.multipleSpectra(nu,atm,zern,tel);
%     save(filename, 'ZPSD');
% else
%     load(filename)
% end

tic
for kLayer = 1:atm.nLayer
    for kmode=1:Nmodes
        
        %{
------------HEADER-----------------
Objective         ::  Find the coefficients of a perturbation model that best fit the Npts first
                      steps of the (de)correlation curve correl

Comments:             Based on identical procedures developped at ONERA

INPUT VARS
 correl            :: The correlation sequence with at least Npts length, normalised
 t                 :: The temporal vector of points were the correlation is evaluated
 Npts              :: Number of points to consider in the fitting
 ord               :: Order of the model to fit

OUTPUT VARS
 coeff             :: Coefficients of the fitting model
Created by         :: C. Correia
Creation date      :: 06/05/2009
                      
Change Record:     :: 4/8/2011
                      - Bug found for identification of AR1, now corrected for
------------HEADER END----------------
        %}
        
        
             
        
        % ---  Interpolation taking place ---
        PSDi = ZPSD(kmode,:,kLayer);
        %PSDi = Fzer(kmode,:);
        f1(1) = nu(1);
        tmp = interp1(nu, PSDi, f1);
        
        PSDi = [f1, tmp];
        
        %         figure(46)
        %         loglog(f1, tmp,'r--')
        %         hold on
        
        sp = PSDi(1:end,2)';
        
        sp_full  = [0 sp(length(sp):-1:2)  sp];              %to keep symmetry. Otherwise fft(sp_full) presents imaginary components...
        Cphi_aoa = fftshift(fft(fftshift(sp_full)));        % autocorrelation thru Wiener-Khinchine theorem
        Cphi_aoa = real(Cphi_aoa);
        Cphi_aoa = Cphi_aoa/max(Cphi_aoa);
        
        
        
        %----------------------------------------------------------------------
        % --- make fit and compute model temporal auto-correlation ---
        % 1) use params for a continuous time model in Laplace domain,
        % 2) extract magnitude,
        % 3) compute the auto-correlation by Wiener-Khinchine---
        %----------------------------------------------------------------------
        correl   = Cphi_aoa(1,end/2+1:end-1);
        do_plots = 1;
            if do_plots
         figure(54), hold on
         plot(tvec, correl,'k--')
            end
        if ord == 1
            x0 = 1;
        elseif ord == 2
            %if ~exist('x0','var')
                x0 = [5,5];
            %else
            %    x0 = coeff';
            %end
        elseif ord == 3
            x0 = [5,5,5];
        end
       
        [coeff,resid,nIter] = bfgs(x0,1e-6,correl,f,NptsFit,ord);
        % minimization parameter substitution: p1 = coeff(1)^2, p2 = coeff(2)^2 to
        % enforce constraint [p1, p2] > 0.
        coeff = coeff.^2;
        res(kLayer,kmode) = resid;
        
               
        if ord == 1
            p1 = coeff(1);
            a1 = -exp(-Ts*p1);
            A(kmode,kLayer) = -a1;
        elseif ord == 2
            
            p1 = coeff(1);
            p2 = coeff(2);
            
            %----------------------------------------------------------------------
            % --- Compute equivalent (same indicial response) AR model ---
            %----------------------------------------------------------------------
            z1 = exp(-p1*Ts);
            z2 = exp(-p2*Ts);
           
            
            Hz = zpk([], [z1 z2], 1,Ts); % a normalization by the gain@0 may be needed here
            Hz = tf(Hz); % convert to tf
            pp = get(Hz,'den');
            pp = pp{1};
            A(kmode,kLayer) = -pp(2);
            B(kmode,kLayer) = -pp(3);
            
            
            % --- See Petit08: First laboratory validation of vibration filtering
            % with LQG control law for Adaptive Optics ---
            %             a1       = -2*exp(-Ts*1/2*(p1+p2))*cos(sqrt(p1*p2)*Ts*sqrt(1-1/4*(p1^2+2*p1*p2+p2^2)/(p1*p2)));
            %             b1       = exp(-(p1+p2)*Ts);
            %
            %             A(kmode,kLayer) = -a1;
            %             B(kmode,kLayer) = -b1;
            do_plots = 0;
            if do_plots
                % --- See Petit08: First laboratory validation of vibration filtering
                        % with LQG control law for Adaptive Optics ---
                        a1       = -2*exp(-Ts*1/2*(p1+p2))*cos(sqrt(p1*p2)*Ts*sqrt(1-1/4*(p1^2+2*p1*p2+p2^2)/(p1*p2)));
                        b1       = exp(-(p1+p2)*Ts);
                        H        = (p1*p2*Ts)/(1 + a1/z + b1/(z^2));
                        [mag,~,~]= bode(H,2*pi*fz);
                        mag      = reshape(mag, 1,length(mag));
                        corr     = utilities.correl_fcn_WienerKhinchin(mag.^2);
                        Np       = length(corr);
                        
                        figure(54), hold on
                plot(tz, corr(floor(Np/2+1):end),'k.');
                title(sprintf('Temporal auto-correlation functions mode %4d',kmode+1),'fontsize',14,'fontweight','bold')
                ylabel('Temporal auto-correlation','fontsize',14)
                legend('Disturbance auto-correlation','Identified model auto-correlation (20pts=0.1s)','Discrete ARn model')
                xlabel('Time [s]','fontsize',14)
                
            end
        elseif ord == 3
            p1 = coeff(1);
            p2 = coeff(2);
            p3 = coeff(3);
            
            z1 = exp(-p1*Ts);
            z2 = exp(-p2*Ts);
            z3 = exp(-p3*Ts);
            
            Hz = zpk([], [z1 z2 z3], 1,Ts); % a normalization by the gain@0 may be needed here
            Hz = tf(Hz); % convert to tf
            pp = get(Hz,'den');
            pp = pp{1};
            A(kmode,kLayer) = -pp(2);
            B(kmode,kLayer) = -pp(3);
            C(kmode,kLayer) = -pp(4);
        end
        
    end
    toc
end
end
function out = evaluateFncOrder1(x,correl,f,Npts)

syswind = tf(x^2,[1 x^2]);
[mag]   = bode(syswind,2*pi*f);
mag     = reshape(mag, 1,length(mag));
sp_full = [0 mag(length(mag):-1:1)  mag];
Cphi    = fftshift(fft(fftshift(sp_full.^2)));     % autocorrelation thru Wiener-Khinchine theorem
Cphi    = real(Cphi);
Cphi    = Cphi/max(Cphi);
Np      = length(Cphi);

out = var(correl(1:Npts) - Cphi(floor(Np/2+1):floor(Np/2+1)+Npts-1));
end
function out = evaluateFnc(x,correl,f,Npts)

syswind= tf(x(1)^2*x(2)^2,[1 x(1)^2+x(2)^2 x(1)^2*x(2)^2]);
[mag]  = bode(syswind,2*pi*f);
mag    = reshape(mag, 1,length(mag));
sp_full= [0 mag(length(mag):-1:1)  mag];
Cphi   = fftshift(fft(fftshift(sp_full.^2)));     % autocorrelation thru Wiener-Khinchine theorem
Cphi   = real(Cphi);
Cphi   = Cphi/max(Cphi);
Np     = length(Cphi);

out = var(correl(1:Npts) - Cphi(floor(Np/2+1):floor(Np/2+1)+Npts-1));
end

function out = evaluateFncOrder3(x,correl,f,Npts)
            syswind= zpk([],[-x(1)^2 -x(2)^2 -x(3)^2],x(1)^2*x(2)^2*x(3)^2);
            [mag]  = bode(syswind,2*pi*f);            
            mag    = reshape(mag, 1,length(mag));
            sp_full= [0 mag(length(mag):-1:1)  mag];
            Cphi   = fftshift(fft(fftshift(sp_full.^2)));     % autocorrelation thru Wiener-Khinchine theorem
            Cphi   = real(Cphi);
            Cphi   = Cphi/max(Cphi);
            Np     = length(Cphi);
           
            out = var(correl(1:Npts) - Cphi(floor(Np/2+1):floor(Np/2+1)+Npts-1));
end
        
function out = gradient(x,delta,correl,f,Npts,order)

if order == 1
    fnc = evaluateFncOrder1(x,correl,f,Npts);
elseif order == 2
    fnc = evaluateFnc(x,correl,f,Npts);
elseif order == 3
    fnc = evaluateFncOrder3(x,correl,f,Npts);
end
n = length(x);
for i = 1:n
    x(i) = x(i)+delta;
    if order == 1
        fDelta = evaluateFncOrder1(x,correl,f,Npts);
    elseif order == 2
        fDelta = evaluateFnc(x,correl,f,Npts);
    elseif order == 3
        fDelta = evaluateFncOrder3(x,correl,f,Npts);
    end
    out(i) = (fDelta-fnc)/delta;
    x(i) = x(i)-delta;
end
end

function [xs,fs,k] = bfgs(x0,epsi1,correl,f,Npts,order)
disp(' ')
disp('Program bfgs.m')
n = length(x0);
I = eye(n);
k = 1;
xk = x0';
Sk = I;
if order == 1
    fk = evaluateFncOrder1(xk,correl,f,Npts);
elseif order == 2
    fk = evaluateFnc(xk,correl,f,Npts);
elseif order == 3
    fk = evaluateFncOrder3(xk,correl,f,Npts);
end
gk = gradient(xk,1e-5,correl,f,Npts,order)';
dk = -Sk*gk;
ak = inex_lsearch(xk,dk,correl,f,Npts,order);
dtk = ak*dk;
xk_new = xk + dtk;
if order ==1
    fk_new = evaluateFncOrder1(xk_new,correl,f,Npts);
elseif order ==2
    fk_new = evaluateFnc(xk_new,correl,f,Npts);
elseif order == 3
    fk_new = evaluateFncOrder3(xk_new,correl,f,Npts);
end
dfk = abs(fk - fk_new);
err = max(dfk,norm(dtk));
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
    if order == 1
        fk_new = evaluateFncOrder1(xk_new,correl,f,Npts);
    elseif order == 2
        fk_new = evaluateFnc(xk_new,correl,f,Npts);
    elseif order == 3
        fk_new = evaluateFncOrder3(xk_new,correl,f,Npts);
    end
    dfk = abs(fk - fk_new);
    err = max(dfk,norm(dtk));
    k = k + 1;
end
format long
disp('solution point:')
xs = xk_new
disp('objective function at the solution point:')
if order ==1
    fs = evaluateFncOrder1(xs,correl,f,Npts);
elseif order == 2
    fs = evaluateFnc(xs,correl,f,Npts);
elseif order == 3
    fs = evaluateFncOrder3(xs,correl,f,Npts);
end
format short
disp('number of iterations at convergence:')
k
end

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
if order == 1
    f0 = evaluateFncOrder1(xk,correl,f,Npts);
elseif order == 2
    f0 = evaluateFnc(xk,correl,f,Npts);
elseif order == 3
    f0 = evaluateFncOrder3(xk,correl,f,Npts);
end
gk = gradient(xk,1e-5,correl,f,Npts,order)';
%                 eval(['f0 = ' F '(xk' parameterstring ');']);
%                 eval(['gk = ' G '(xk' parameterstring ');']);
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
    if order == 1
        f0 = evaluateFncOrder1(xk+deltak,correl,f,Npts);
    elseif order == 2
        f0 = evaluateFnc(xk+deltak,correl,f,Npts);
    elseif order == 3
        f0 = evaluateFncOrder3(xk+deltak,correl,f,Npts);
    end
    %eval(['f0 = ' F '(xk+deltak' parameterstring ');']);
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
        %eval(['gtemp =' G '(xk+a0*dk' parameterstring ');']);
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
