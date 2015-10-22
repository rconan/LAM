%--------------------------------------------------------------------------
% open_loop_simulator.m
%
% Script which runs an open loop simulation.  Has options to specify which
% type of reconstructor is used.
%
% Input parameters = Joel's parameters (with exception of number of
% pixels, which is set to 6 per lenslet, an even number).
% Results = strehl, residual phase psd, psfs
%
% Charlotte Bond        14.09.15
%--------------------------------------------------------------------------
%


%% Set parameters for reconstructors

noiseVar = mean(theoreticalNoise(wfs,tel,atm,ngs,sci));

if strcmp(reconstructor,'MMSE')
    % MMSE 
    % create inverse of matrix F, where F is matrix of DM influence functions 
    % to go from phase to DM commands. A tolerance of 0.1 is used, 
    % values > 0.1 will be set to zero.
    F = 2*bifLowRes.modes(wfs.validActuator,:);
    iF = pinv(full(F),0.1); 
    MMSE = slopesLinearMMSE(wfs,tel,atm,ngs,'mmseStar',ngs);
elseif strcmp(reconstructor,'FFTR')
    FFTR = fourierReconstructor(wfs,tel,atm,'extensionMethod',method,...
    'extension',ext,'dm',dm,'stroke',ngs.wavelength/2,'filter',filter,...
    'GerchbergFilter',Gfilter,'GerchbergIterations',Gnbiter,...
    'modeRemoval',modes,'kAliasing',ka,'kNoise',kn,'noiseVariance',noiseVar);
elseif strcmp(reconstructor,'Direct')
    % For direct projection (i.e. high resolution phase)
    Fhr = 2*bif.modes(tel.pupilLogical,:);
    iFhr = pinv(full(Fhr),0.1);
end

%%

% Vectors for psds and the hanning window
psd = zeros(nPx,nPx,3,2);
N = nLenslet+2*ext+2;
% 2 slopes (x,y)
% 2 screens (total/low order)
% 2 methods (weighted, not weighted)
psdS = zeros(N,N,2,2,2);
wxy = mywindow('squareTukey',N,1-2*ext/N);
wS = wxy(:,:,1).*wxy(:,:,2);
w = mywindow('hann',nPx);
k = 1/D*(-nPx/2:nPx/2-1);

nIter = exposureTime + startDelay;


scrsz = get(groot,'ScreenSize');
figure('Position',[1 scrsz(4)/2 2*scrsz(3)/3 scrsz(4)])
%figure
j=0;
for i=1:nIter
    j=j+1;
    if j==round(0.01*nIter) || i==nIter
        fprintf('Loop cycle %g of %g \n',i,nIter)
        j=0;
    end
    
    for n=1:3
        % Propagate source and save turbulent phase, take psd of orginal phase
        ngs = ngs.*tel;
        ngs.phase = phi(:,:,n,i);
        phase = ngs.meanRmOpd;
        if i>startDelay
            psd(:,:,n,1) = psd(:,:,n,1) + (1/exposureTime)*abs(fftshift(fft2(w.*phase))).^2;
        end
        if n==1
            phaseT = phase;
        end

        % Propagate wavefront to wfs
        ngs = ngs*wfs;
        % Reconstruct
        if strcmp(reconstructor,'LSQ')
            % LSQ
            % Use pseudo inverse of interaction matrix to get dm coefficients
            dm.coefs = -calibDm.M*wfs.slopes;
        elseif strcmp(reconstructor,'MMSE')
            % MMSE
            % Get dm coefficients by first calculating the low resolution phase
            % using an MMSE reconstructor, then getting the commands using the dm
            % filter functions (the inverse of the influence functions).
            phaseRecon = MMSE*wfs.slopes;
            dm.coefs = iF*phaseRecon(wfs.validActuator);
        elseif strcmp(reconstructor,'FFTR')
            % FFTR
            [~,coefs] = FFTR*wfs;
            dm.coefs = coefs;
            if n~=3
                psdS(:,:,1,n,1) = psdS(:,:,1,n,1) + (1/exposureTime)*abs(fftshift(fft2(FFTR.Sx))).^2;
                psdS(:,:,2,n,1) = psdS(:,:,2,n,1) + (1/exposureTime)*abs(fftshift(fft2(FFTR.Sy))).^2;
                psdS(:,:,1,n,2) = psdS(:,:,1,n,2) + (1/exposureTime)*abs(fftshift(fft2(wS.*FFTR.Sx))).^2;
                psdS(:,:,2,n,2) = psdS(:,:,2,n,2) + (1/exposureTime)*abs(fftshift(fft2(wS.*FFTR.Sy))).^2;
            end
        elseif strcmp(reconstructor,'Direct')
            OPD = ngs.opd;
            dm.coefs = iFhr*OPD(tel.pupilLogical);
        end

        % Propagate to dm
        ngs = ngs*dm;
        % Propagate science case
        sci = sci.*tel;
        sci.phase = phi(:,:,n,i).*ngs.wavelength/sci.wavelength;
        sci = sci*dm*cam(n);

        % Phases and results
        phaseRes = ngs.meanRmOpd;
        phaseDM = 2*dm.surface.*tel.pupil;
        rms(i,n) = ngs.opdRms*1e9;
        if i>startDelay
            psd(:,:,n,2) = psd(:,:,n,2) + (1/exposureTime)*abs(fftshift(fft2(w.*phaseRes))).^2;
        end
        if n==1
            phaseResT = phaseRes;
            phaseDMT = phaseDM;
        end
    end
    %----------------------------------------------------------------------
    % Plots
    % Phase screens
    subplot(3,3,[1 2 3])
    imagesc([phaseT phaseDMT phaseResT])
    axis tight
    axis equal
    colorLim = max(max(abs(phaseT)));
    set(gca,'CLim',[-colorLim colorLim])
    drawnow;
    
    % Rms residual phase
    for n=1:3
        subplot(3,3,n+3)
        plot(1:i,rms(1:i,n))
        set(gca,'XLim',[0 nIter])
        set(gca,'YLim',[0 200])
        xlabel('iteration')
        ylabel('rms [nm]')
    end
    
    % PSD
    for n=1:3
        subplot(3,3,n+6)
        loglog(k,psd(nPx/2+1,:,n,1),k,psd(nPx/2+1,:,n,2))
        legend('uncorrected phase','residual phase')
        xlabel('spatial frequency [m^{-1}]')
        ylabel('PSD(\phi)')
    end
end