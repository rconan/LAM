%--------------------------------------------------------------------------
% Run1_open_loop.m
%
% Runs an open loop simulation for a range of wavefront reconstructors and
% their associated parameters.
%
%
% Charlotte Bond        14.09.15
%--------------------------------------------------------------------------
%

clear all;

filename = '~/Documents/MATLAB/FFT_reconstruction/testing/e2e/phaseCubes/Allphasecube.mat';
basename = 'Run_211015';

% Variable parameters
%--------------------------------------------------------------------------
% TIMES
% Camera exposure time (averaging time for psf, strehl etc.) and delay
% This is in frames, to convert to time, * sampling time.
exposureTime = 20;
startDelay = 0;         
nScreens = 2000;

% NOISE
photonNoise = true;
readOutNoise = 2;
% trying m=1!
% RECONSTRUCTORs and relevant parameters
sims = [4 7 10 13 17];%1:17;
Nsims = 18;
NFFTRsims = 15;
Recon{1} = 'Direct';
Recon{2} = 'LSQ';
Recon{3} = 'MMSE';
Recon(4:Nsims) = {'FFTR'};

Filt(1:3) = {'Hudgin'};
Exten(1:3) = {'simple' 'Gerchberg' 'Gerchberg'};
Gerch(1:3) = {'' 'Hudgin' 'Exact'};
Filt(4:6) = {'Fried'};
Exten(4:6) = {'simple' 'Gerchberg' 'Gerchberg'};
Gerch(4:6) = {'' 'Fried' 'Exact'};
Filt(7:9) = {'Southwell'};
Exten(7:9) = {'simple' 'Gerchberg' 'Gerchberg'};
Gerch(7:9) = {'' 'Southwell' 'Exact'};
Filt(10:12) = {'antialias'};
Exten(10:12) = {'simple' 'Gerchberg' 'Gerchberg'};
Gerch(10:12) = {'' 'antialias' 'Exact'};
Filt(13:14) = {'Exact'};
Exten(13:14) = {'simple' 'Gerchberg'};
Gerch(13:14) = {'' 'Exact'};
Filt(15) = {'e2eAA'};
Exten(15) = {'Gerchberg'};
Gerch(15) = {'Exact'};

Gnbiter = 7;
ext = 3;
kn = 1;
ka = 1;
%noiseVar = 0;

% 1: Global waffle removal (direct) (up strehl, ~94.2%, but gets worse with
% anti alias)
% 2: global and local waffle removal (direct) (down in strehl, ~93.6%, not
% such a big difference with anti-alias, but still no improvement)
% 3: G waffle removal fourier (same result from direct)
% 4: G and l waffle removal fourier
%
% ka +/- 0.1,0.2,0.3,0.4,0.5 (using direct G waffle removal)
% ka +/- 0.1,0.2,0.3,0.4,0.5 (using direct G and L waffle removal)
modes = [1 0 0];
%--------------------------------------------------------------------------

for h=1:length(sims)

    % Simulation
    %----------------------------------------------------------------------
    reconstructor = Recon{sims(h)};
    if sims(h)>Nsims-NFFTRsims
        method = Exten{sims(h)-(Nsims-NFFTRsims)};
        filter = Filt{sims(h)-(Nsims-NFFTRsims)};
        Gfilter = Gerch{sims(h)-(Nsims-NFFTRsims)};
        savename = sprintf('%s_%s_filter:%s_extension:%s%s.mat',basename,reconstructor,...
            filter,method,Gfilter);
    else
        savename = sprintf('%s_%s',basename,reconstructor);
    end
    
    % Load simulation parameters and cube of phase
    simulationParameters

    % Setup OOMAO, initiate objects, calibrate etc.
    % All prepared for open loop wavefront reconstruction
    setupOOMAO

    % Run open loop simulation
    open_loop_simulator
    %----------------------------------------------------------------------
    if sims(h)>Nsims-NFFTRsims
        save(savename,'reconstructor','FFTR','psd','rms','cam','psdS','wS')
    else
        save(savename,'reconstructor','wfs','atm','tel','psd','rms','cam')
    end
    
    fprintf('strehl (total)      = %g \n',cam(1).strehl)
    fprintf('strehl (low order)  = %g \n',cam(2).strehl)
    fprintf('strehl (high order) = %g \n',cam(3).strehl)
    
    close all;
    clear wfs dm tel atm psd rms cam ngs;
    
end



