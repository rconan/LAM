
clear all;

% load('Run_181015_LSQ')
% psdLSQ = psd;
% cam(1).strehl
% camLSQ = cam(1);
% 
% load('Run_181015_MMSE')
% psdMMSE = psd;
% cam(1).strehl
% camMMSE = cam(1);

% load('Run_171015_FFTR_filter:Exact_extension:simple')
% psdFFTRExactS = psd;
% cam(1).strehl
% camExactS = cam(1);

% load('Run_171015_FFTR_filter:antialias_extension:simple')
% psdFFTRAAS = psd;
% cam(1).strehl
% camAAS = cam(1);

load('Run_211015_FFTR_filter:Exact_extension:GerchbergExact')
psdFFTRExactG = psd;
cam(1).strehl
camExactG = cam(1);
FFTRE = FFTR;

load('Run_211015_FFTR_filter:antialias_extension:Gerchbergantialias')
psdFFTRe2e = psd;
cam(1).strehl
camAAG = cam(1);
FFTRAA = FFTR;

% load('Run_211015_FFTR_filter:antialias_extension:GerchbergExact')
% psdFFTRAAGE = psd;
% cam(1).strehl
% camAAGE = cam(1);
% FFTRAAE = FFTR;

nPx = size(psd,1);
D = 8;
k = 1/D*(-nPx/2:nPx/2-1);

figure
loglog(k,psd(nPx/2+1,:,1,1),...
    k,psdFFTRExactG(nPx/2+1,:,1,2),...
    k,psdFFTRe2e(nPx/2+1,:,1,2))
legend('uncorrected phase',...
    'Fourier: Exact, Gerchberg','Fourier: e2e, Exact Gerchberg')
xlabel('spatial frequency [m^{-1}]')
ylabel('PSD(\phi)')

cmap = FT_flip_colormap();
clim = max(max(abs((psdFFTRe2e(:,:,1,2)-psdFFTRExactG(:,:,1,2))./psdFFTRExactG(:,:,1,2))));
figure
    imagesc(k,k,(psdFFTRe2e(:,:,1,2)-psdFFTRExactG(:,:,1,2))./psdFFTRExactG(:,:,1,2))
    axis equal
    axis tight
    colorbar()
    set(gca,'clim',[-clim clim])
    colormap(cmap)
    set(gca,'YDir','normal')
    xlim([-2.5 2.5])
    ylim([-2.5 2.5])
    
% figure
%     imagesc(k,k,(psdFFTRAAGE(:,:,1,2)-psdFFTRExactG(:,:,1,2))./psdFFTRExactG(:,:,1,2))
%     axis equal
%     axis tight
%     colorbar()
%     set(gca,'clim',[-clim clim])
%     colormap(cmap)
%     set(gca,'YDir','normal')
%     xlim([-2.5 2.5])
%     ylim([-2.5 2.5])
%     
%     subplot(2,2,2)
%     imagesc(k,k,psdFFTRAA_20ka(:,:,1,2)./psdFFTRExact(:,:,1,2)-1)
%     axis equal
%     axis tight
%     colorbar()
%     set(gca,'clim',[-0.7 0.7])
%     colormap(cmap)
%     set(gca,'YDir','normal')
%     
%     subplot(2,2,3)
%     imagesc(k,k,psdFFTRAA_10ka(:,:,1,2)./psdFFTRExact(:,:,1,2)-1)
%     axis equal
%     axis tight
%     colorbar()
%     set(gca,'clim',[-0.7 0.7])
%     colormap(cmap)
%     set(gca,'YDir','normal')
%     
%     subplot(2,2,4)
%     imagesc(k,k,psdFFTRAA_n10ka(:,:,1,2)./psdFFTRExact(:,:,1,2)-1)
%     axis equal
%     axis tight
%     colorbar()
%     set(gca,'clim',[-0.7 0.7])
%     colormap(cmap)
%     set(gca,'YDir','normal')
%  
% xPSF = size(cam(1).frame,1);    
%     
% figure
%     plot(1:xPSF/2,cam(1).referenceFrame(xPSF/2+1,xPSF/2+1:xPSF),...
%         1:xPSF/2,camLSQ.frame(xPSF/2+1,xPSF/2+1:xPSF),...
%         1:xPSF/2,camMMSE.frame(xPSF/2+1,xPSF/2+1:xPSF),...
%         1:xPSF/2,camExact.frame(xPSF/2+1,xPSF/2+1:xPSF),...
%         1:xPSF/2,camAA_100ka.frame(xPSF/2+1,xPSF/2+1:xPSF),...
%         1:xPSF/2,camAA_20ka.frame(xPSF/2+1,xPSF/2+1:xPSF),...
%         1:xPSF/2,camAA_10ka.frame(xPSF/2+1,xPSF/2+1:xPSF),...
%         1:xPSF/2,camAA_n10ka.frame(xPSF/2+1,xPSF/2+1:xPSF))
%     set(gca,'YScale','log')
%     legend('reference','LSQ','MMSE','Exact','AA: 100%','AA: 20%','AA: 10%','AA: -10%')
%     
%     
% figure
%     subplot(2,1,1)
%     imagesc(camExact.frame)
%     axis equal
%     axis tight
%     colorbar()
%     %set(gca,'clim',[-65 65])
%     colormap(cmap)
%     set(gca,'YDir','normal')
%     
%     subplot(2,1,2)
%     imagesc(camAA_100ka.frame)
%     axis equal
%     axis tight
%     colorbar()
%     %set(gca,'clim',[-65 65])
%     colormap(cmap)
%     set(gca,'YDir','normal')
%     