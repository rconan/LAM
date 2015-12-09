%--------------------------------------------------------------------------
% Script to read in and process data from experiment and simulations of
% Pyramid measurements: poke of central actuator

clear all;
close all;

% Actuator pitch in nanometre
%pAmp = 500;
pAmp = 1;
D = 38;
R = D/2;
nb_poke = 41;
poke = linspace(-0.8,0.8,41);
poke = pAmp*poke;

% Load data
imageEx = fitsread('PymEx_pupil_images_081015.fits');
imageSim = fitsread('PymSim_pupil_images_081015.fits');

% Plot raw images (first poke, -80%)
figure
    subplot(2,1,1)
    imagesc(imageEx(:,:,1))
    axis equal
    axis tight
    set(gca,'YDir','normal')
    
    subplot(2,1,2)
    imagesc(imageSim(:,:,1))
    axis equal
    axis tight
    set(gca,'YDir','normal')
    
%% Create pupil
x = -(D-1)/2:D-(D+1)/2;
[X,Y] = meshgrid(x,x);
r = sqrt(X.^2+Y.^2);
pupil = zeros(D,D);
pupil(r<=D/2) = 1;
smallpup = zeros(D,D);
smallpup(r<=D/2) = 1;

figure
    imagesc(pupil)
    axis equal
    axis tight


%% Cropping imageEx and extracting 4 pupils
% centre coordinates
x0 = [36 194 39 199];
y0 = [194 198 34 37];

% Selecting different pupil images for sim data and applying pupil mask
for n=1:nb_poke
    I_Ex(:,:,1,n) = pupil.*imageEx(y0(1)-R+1:y0(1)+R,x0(1)-R+1:x0(1)+R,n);
    I_Ex(:,:,2,n) = pupil.*imageEx(y0(2)-R+1:y0(2)+R,x0(2)-R+1:x0(2)+R,n);
    I_Ex(:,:,3,n) = pupil.*imageEx(y0(3)-R+1:y0(3)+R,x0(3)-R+1:x0(3)+R,n);
    I_Ex(:,:,4,n) = pupil.*imageEx(y0(4)-R+1:y0(4)+R,x0(4)-R+1:x0(4)+R,n);
    I_Sim(:,:,1,n) = smallpup.*imageSim(D+1:2*D,1:D,n);
    I_Sim(:,:,2,n) = smallpup.*imageSim(D+1:2*D,D+1:2*D,n);
    I_Sim(:,:,3,n) = smallpup.*imageSim(1:D,1:D,n);
    I_Sim(:,:,4,n) = smallpup.*imageSim(1:D,D+1:2*D,n);
end

figure
    imagesc([I_Ex(:,:,3,1) I_Ex(:,:,4,1); I_Ex(:,:,1,1) I_Ex(:,:,2,1)])
    axis equal
    axis tight
    set(gca,'YDir','normal')
    set(gca,'CLim',[0 max(max([I_Ex(:,:,3,1) I_Ex(:,:,4,1); I_Ex(:,:,1,1) I_Ex(:,:,2,1)]))])
    %colormap(jet)
    axis off
    
figure
    imagesc([I_Sim(:,:,3,1) I_Sim(:,:,4,1); I_Sim(:,:,1,1) I_Sim(:,:,2,1)])
    axis equal
    axis tight
    set(gca,'YDir','normal')
    set(gca,'CLim',[0 max(max([I_Sim(:,:,3,1) I_Sim(:,:,4,1); I_Sim(:,:,1,1) I_Sim(:,:,2,1)]))])
    %colormap(jet)
    axis off
    
%% Reference signals (i.e. flat DM at poke = 0)
n0 = find(poke==0);

% Flat reference measurements
I_Ex_ref = I_Ex(:,:,:,n0);
I_Sim_ref = I_Sim(:,:,:,n0);

%% Pupil images
% 1) raw images
% Postive and negative pokes (+/- 0, 16%, 32%, 48%, 64%, 80%)
ps = 21:-3:1;
ps = [ps 21:3:nb_poke];

tmp_I_Ex = zeros(D,D*length(ps),4);
tmp_I_Sim = zeros(D,D*length(ps),4);
for i=1:length(ps)
    tmp_I_Ex(:,(i-1)*D+1:i*D,:) = I_Ex(:,:,:,ps(i));
    tmp_I_Sim(:,(i-1)*D+1:i*D,:) = I_Sim(:,:,:,ps(i));
end

c_max_E = max(max(max(tmp_I_Ex)));
c_min_E = min(min(min(tmp_I_Ex))); 
c_max_S = max(max(max(tmp_I_Sim)));
c_min_S = min(min(min(tmp_I_Sim)));

%cmap = winter;

scrsz = get(groot,'ScreenSize');
figure('Position',[1 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2])
    for p=1:4
        subplot(4,2,2*p-1)
        imagesc([tmp_I_Ex(:,1:end/2,p); tmp_I_Ex(:,end/2+1:end,p)])
        axis equal
        axis tight
        axis off
        set(gca,'YDir','normal')
        set(gca,'clim',[0 c_max_E])
        %colormap(cmap)
        
        subplot(4,2,2*p)
        imagesc([tmp_I_Sim(:,1:end/2,p); tmp_I_Sim(:,end/2+1:end,p)])
        axis equal
        axis tight
        axis off
        set(gca,'YDir','normal')
        set(gca,'clim',[0 c_max_S])
        %colormap(cmap)
    end
    
% 2) Reference removed
tmp_I_Ex = zeros(D,D*length(ps),4);
tmp_I_Sim = zeros(D,D*length(ps),4);
for i=1:length(ps)
    tmp_I_Ex(:,(i-1)*D+1:i*D,:) = I_Ex(:,:,:,ps(i)) - I_Ex_ref;
    tmp_I_Sim(:,(i-1)*D+1:i*D,:) = I_Sim(:,:,:,ps(i)) - I_Sim_ref;
end

c_max_E = max(max(max(tmp_I_Ex)));
c_min_E = min(min(min(tmp_I_Ex))); 
c_max_S = max(max(max(tmp_I_Sim)));
c_min_S = min(min(min(tmp_I_Sim))); 

scrsz = get(groot,'ScreenSize');
figure('Position',[1 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2])
    for p=1:4
        subplot(4,2,2*p-1)
        imagesc([tmp_I_Ex(:,1:end/2,p); tmp_I_Ex(:,end/2+1:end,p)])
        axis equal
        axis tight
        axis off
        set(gca,'YDir','normal')
        set(gca,'clim',[0 c_max_E])
        %colormap(cmap)
        
        subplot(4,2,2*p)
        imagesc([tmp_I_Sim(:,1:end/2,p); tmp_I_Sim(:,end/2+1:end,p)])
        axis equal
        axis tight
        axis off
        set(gca,'YDir','normal')
        set(gca,'clim',[0 c_max_S])
        %colormap(cmap)
    end
    

%% Slopes: Old calculation, non-flat map, Rigazoni method
% Sx = I1+I3-I2-I4/(I1+I2+I3+I4)
% Experimental
Sx_Ex_ref = (I_Ex_ref(:,:,1)+I_Ex_ref(:,:,3)-I_Ex_ref(:,:,2)-I_Ex_ref(:,:,4))./...
    (I_Ex_ref(:,:,1)+I_Ex_ref(:,:,2)+I_Ex_ref(:,:,3)+I_Ex_ref(:,:,4));
Sy_Ex_ref = (I_Ex_ref(:,:,1)+I_Ex_ref(:,:,2)-I_Ex_ref(:,:,3)-I_Ex_ref(:,:,4))./...
    (I_Ex_ref(:,:,1)+I_Ex_ref(:,:,2)+I_Ex_ref(:,:,3)+I_Ex_ref(:,:,4));

for n=1:nb_poke
    Sx_Ex(:,:,n) = (I_Ex(:,:,1,n)+I_Ex(:,:,3,n)-I_Ex(:,:,2,n)-I_Ex(:,:,4,n))./...
        (I_Ex(:,:,1,n)+I_Ex(:,:,2,n)+I_Ex(:,:,3,n)+I_Ex(:,:,4,n));
    Sx_Ex(:,:,n) = Sx_Ex(:,:,n) - Sx_Ex_ref; 
    Sy_Ex(:,:,n) = (I_Ex(:,:,1,n)+I_Ex(:,:,2,n)-I_Ex(:,:,3,n)-I_Ex(:,:,4,n))./...
        (I_Ex(:,:,1,n)+I_Ex(:,:,2,n)+I_Ex(:,:,3,n)+I_Ex(:,:,4,n));
    Sy_Ex(:,:,n) = Sy_Ex(:,:,n) - Sy_Ex_ref; 
end

Sx_Ex(isnan(Sx_Ex)) = 0;
Sy_Ex(isnan(Sy_Ex)) = 0;

% Sim
Sx_Sim_ref = (I_Sim_ref(:,:,1)+I_Sim_ref(:,:,3)-I_Sim_ref(:,:,2)-I_Sim_ref(:,:,4))./...
    (I_Sim_ref(:,:,1)+I_Sim_ref(:,:,2)+I_Sim_ref(:,:,3)+I_Sim_ref(:,:,4));
Sy_Sim_ref = (I_Sim_ref(:,:,1)+I_Sim_ref(:,:,2)-I_Sim_ref(:,:,3)-I_Sim_ref(:,:,4))./...
    (I_Sim_ref(:,:,1)+I_Sim_ref(:,:,2)+I_Sim_ref(:,:,3)+I_Sim_ref(:,:,4));

for n=1:nb_poke
    Sx_Sim(:,:,n) = (I_Sim(:,:,1,n)+I_Sim(:,:,3,n)-I_Sim(:,:,2,n)-I_Sim(:,:,4,n))./...
        (I_Sim(:,:,1,n)+I_Sim(:,:,2,n)+I_Sim(:,:,3,n)+I_Sim(:,:,4,n));
    Sx_Sim(:,:,n) = Sx_Sim(:,:,n) - Sx_Sim_ref; 
    Sy_Sim(:,:,n) = (I_Sim(:,:,1,n)+I_Sim(:,:,2,n)-I_Sim(:,:,3,n)-I_Sim(:,:,4,n))./...
        (I_Sim(:,:,1,n)+I_Sim(:,:,2,n)+I_Sim(:,:,3,n)+I_Sim(:,:,4,n));
    Sy_Sim(:,:,n) = Sy_Sim(:,:,n) - Sy_Sim_ref; 
end

Sx_Sim(isnan(Sx_Sim)) = 0;
Sy_Sim(isnan(Sy_Sim)) = 0;

%% Plot slopes
cmap = 1-copper;
%cmap = FT_flip_colormap();;
%cmap = parula;
%cmap = hot;
plot_slopes
    
%% Plot maxima
maxSx_Ex = squeeze(max(max(Sx_Ex)));
maxSy_Ex = squeeze(max(max(Sy_Ex)));
minSx_Ex = squeeze(min(min(Sx_Ex)));
minSy_Ex = squeeze(min(min(Sy_Ex)));

maxSx_Sim = squeeze(max(max(Sx_Sim)));
maxSy_Sim = squeeze(max(max(Sy_Sim)));
minSx_Sim = squeeze(min(min(Sx_Sim)));
minSy_Sim = squeeze(min(min(Sy_Sim)));

c0 = 21; 
mx1 = abs(maxSx_Ex(c0+2)/poke(c0+2));
mx2 = abs(maxSx_Ex(c0-2)/poke(c0-2));
mx3 = abs(minSx_Ex(c0+2)/poke(c0+2)); 
mx4 = abs(minSx_Ex(c0-2)/poke(c0-2)); 
mx = (mx1+mx2+mx3+mx4)/4;
linX = mx*poke((c0-3):(c0+3));

my1 = abs(maxSy_Ex(c0+2)/poke(c0+2));
my2 = abs(maxSy_Ex(c0-2)/poke(c0-2));
my3 = abs(minSy_Ex(c0+2)/poke(c0+2)); 
my4 = abs(minSy_Ex(c0-2)/poke(c0-2)); 
my = (my1+my2+my3+my4)/4;
linY = my*poke((c0-3):(c0+3));


figure
    plot(poke,maxSx_Ex,'r',poke,minSx_Ex,'b','LineWidth',2)
    hold on;
    plot(poke,maxSx_Sim,'r--',poke,minSx_Sim,'b--','LineWidth',2)
    plot(poke((c0-3):(c0+3)),linX,'k')
    plot(poke((c0-3):(c0+3)),-linX,'k')
    grid on;
    legend('max. Sx: experiment','min. Sx: experimenet',...
        'max. Sx: simulation','min. Sx: simulation')
    ylim([-3 3])
    xlabel('stroke [nm]')
    ylabel('Sx')
    
figure
    plot(poke,maxSy_Ex,'r',poke,minSy_Ex,'b','LineWidth',2)
    hold on;
    plot(poke,maxSy_Sim,'r--',poke,minSy_Sim,'b--','LineWidth',2)
    plot(poke((c0-3):(c0+3)),linY,'k')
    plot(poke((c0-3):(c0+3)),-linY,'k')
    grid on;
    legend('max. Sy: experiment','min. Sy: experimenet',...
        'max. Sy: simulation','min. Sy: simulation')
    ylim([-3 3])
    xlabel('stroke [nm]')
    ylabel('Sy')


%% Slopes: Old calculation, non-flat map but constant (i.e. normalisation 
% measured from reference)
% Sx = I1+I3-I2-I4/(I1_f+I2_f+I3_f+I4_f)
% Experimental
Sx_Ex_ref = (I_Ex_ref(:,:,1)+I_Ex_ref(:,:,3)-I_Ex_ref(:,:,2)-I_Ex_ref(:,:,4))./...
    (I_Ex_ref(:,:,1)+I_Ex_ref(:,:,2)+I_Ex_ref(:,:,3)+I_Ex_ref(:,:,4));
Sy_Ex_ref = (I_Ex_ref(:,:,1)+I_Ex_ref(:,:,2)-I_Ex_ref(:,:,3)-I_Ex_ref(:,:,4))./...
    (I_Ex_ref(:,:,1)+I_Ex_ref(:,:,2)+I_Ex_ref(:,:,3)+I_Ex_ref(:,:,4));

for n=1:nb_poke
    Sx_Ex(:,:,n) = (I_Ex(:,:,1,n)+I_Ex(:,:,3,n)-I_Ex(:,:,2,n)-I_Ex(:,:,4,n))./...
        (I_Ex_ref(:,:,1)+I_Ex_ref(:,:,2)+I_Ex_ref(:,:,3)+I_Ex_ref(:,:,4));
    Sx_Ex(:,:,n) = Sx_Ex(:,:,n) - Sx_Ex_ref; 
    Sy_Ex(:,:,n) = (I_Ex(:,:,1,n)+I_Ex(:,:,2,n)-I_Ex(:,:,3,n)-I_Ex(:,:,4,n))./...
        (I_Ex_ref(:,:,1)+I_Ex_ref(:,:,2)+I_Ex_ref(:,:,3)+I_Ex_ref(:,:,4));
    Sy_Ex(:,:,n) = Sy_Ex(:,:,n) - Sy_Ex_ref; 
end

Sx_Ex(isnan(Sx_Ex)) = 0;
Sy_Ex(isnan(Sy_Ex)) = 0;

% Sim
Sx_Sim_ref = (I_Sim_ref(:,:,1)+I_Sim_ref(:,:,3)-I_Sim_ref(:,:,2)-I_Sim_ref(:,:,4))./...
    (I_Sim_ref(:,:,1)+I_Sim_ref(:,:,2)+I_Sim_ref(:,:,3)+I_Sim_ref(:,:,4));
Sy_Sim_ref = (I_Sim_ref(:,:,1)+I_Sim_ref(:,:,2)-I_Sim_ref(:,:,3)-I_Sim_ref(:,:,4))./...
    (I_Sim_ref(:,:,1)+I_Sim_ref(:,:,2)+I_Sim_ref(:,:,3)+I_Sim_ref(:,:,4));

for n=1:nb_poke
    Sx_Sim(:,:,n) = (I_Sim(:,:,1,n)+I_Sim(:,:,3,n)-I_Sim(:,:,2,n)-I_Sim(:,:,4,n))./...
        (I_Sim_ref(:,:,1)+I_Sim_ref(:,:,2)+I_Sim_ref(:,:,3)+I_Sim_ref(:,:,4));
    Sx_Sim(:,:,n) = Sx_Sim(:,:,n) - Sx_Sim_ref; 
    Sy_Sim(:,:,n) = (I_Sim(:,:,1,n)+I_Sim(:,:,2,n)-I_Sim(:,:,3,n)-I_Sim(:,:,4,n))./...
        (I_Sim_ref(:,:,1)+I_Sim_ref(:,:,2)+I_Sim_ref(:,:,3)+I_Sim_ref(:,:,4));
    Sy_Sim(:,:,n) = Sy_Sim(:,:,n) - Sy_Sim_ref; 
end

Sx_Sim(isnan(Sx_Sim)) = 0;
Sy_Sim(isnan(Sy_Sim)) = 0;

%% Plot slopes
plot_slopes
    
%% Plot maxima
maxSx_Ex = squeeze(max(max(Sx_Ex)));
maxSy_Ex = squeeze(max(max(Sy_Ex)));
minSx_Ex = squeeze(min(min(Sx_Ex)));
minSy_Ex = squeeze(min(min(Sy_Ex)));

maxSx_Sim = squeeze(max(max(Sx_Sim)));
maxSy_Sim = squeeze(max(max(Sy_Sim)));
minSx_Sim = squeeze(min(min(Sx_Sim)));
minSy_Sim = squeeze(min(min(Sy_Sim)));


figure
    plot(poke,maxSx_Ex,'b',poke,minSx_Ex,'r')
    hold on;
    plot(poke,maxSx_Sim,'b--',poke,minSx_Sim,'r--')
    grid on;
    
figure
    plot(poke,maxSy_Ex,'b',poke,minSy_Ex,'r')
    hold on;
    plot(poke,maxSy_Sim,'b--',poke,minSy_Sim,'r--')
    grid on;


%% Slopes: using flat map, non-constant normalisation (taken from each 
% measurement)
% Sx = I1+I3-I2-I4/(sum(sum(I1+I2+I3+I4))
Npup = sum(sum(pupil));
% Experimental
Sx_Ex_ref = (I_Ex_ref(:,:,1)+I_Ex_ref(:,:,3)-I_Ex_ref(:,:,2)-I_Ex_ref(:,:,4))/...
    (sum(sum(I_Ex_ref(:,:,1)+I_Ex_ref(:,:,2)+I_Ex_ref(:,:,3)+I_Ex_ref(:,:,4)))/Npup);
Sy_Ex_ref = (I_Ex_ref(:,:,1)+I_Ex_ref(:,:,2)-I_Ex_ref(:,:,3)-I_Ex_ref(:,:,4))/...
    (sum(sum(I_Ex_ref(:,:,1)+I_Ex_ref(:,:,2)+I_Ex_ref(:,:,3)+I_Ex_ref(:,:,4)))/Npup);

for n=1:nb_poke
    Sx_Ex(:,:,n) = (I_Ex(:,:,1,n)+I_Ex(:,:,3,n)-I_Ex(:,:,2,n)-I_Ex(:,:,4,n))/...
        (sum(sum(I_Ex(:,:,1,n)+I_Ex(:,:,2,n)+I_Ex(:,:,3,n)+I_Ex(:,:,4,n)))/Npup);
    Sx_Ex(:,:,n) = Sx_Ex(:,:,n) - Sx_Ex_ref; 
    Sy_Ex(:,:,n) = (I_Ex(:,:,1,n)+I_Ex(:,:,2,n)-I_Ex(:,:,3,n)-I_Ex(:,:,4,n))/...
        (sum(sum(I_Ex(:,:,1,n)+I_Ex(:,:,2,n)+I_Ex(:,:,3,n)+I_Ex(:,:,4,n)))/Npup);
    Sy_Ex(:,:,n) = Sy_Ex(:,:,n) - Sy_Ex_ref; 
end

% Sim
Sx_Sim_ref = (I_Sim_ref(:,:,1)+I_Sim_ref(:,:,3)-I_Sim_ref(:,:,2)-I_Sim_ref(:,:,4))/...
    (sum(sum(I_Sim_ref(:,:,1)+I_Sim_ref(:,:,2)+I_Sim_ref(:,:,3)+I_Sim_ref(:,:,4)))/Npup);
Sy_Sim_ref = (I_Sim_ref(:,:,1)+I_Sim_ref(:,:,2)-I_Sim_ref(:,:,3)-I_Sim_ref(:,:,4))/...
    (sum(sum(I_Sim_ref(:,:,1)+I_Sim_ref(:,:,2)+I_Sim_ref(:,:,3)+I_Sim_ref(:,:,4)))/Npup);

for n=1:nb_poke
    Sx_Sim(:,:,n) = (I_Sim(:,:,1,n)+I_Sim(:,:,3,n)-I_Sim(:,:,2,n)-I_Sim(:,:,4,n))/...
        (sum(sum(I_Sim(:,:,1,n)+I_Sim(:,:,2,n)+I_Sim(:,:,3,n)+I_Sim(:,:,4,n)))/Npup);
    Sx_Sim(:,:,n) = Sx_Sim(:,:,n) - Sx_Sim_ref; 
    Sy_Sim(:,:,n) = (I_Sim(:,:,1,n)+I_Sim(:,:,2,n)-I_Sim(:,:,3,n)-I_Sim(:,:,4,n))/...
        (sum(sum(I_Sim(:,:,1,n)+I_Sim(:,:,2,n)+I_Sim(:,:,3,n)+I_Sim(:,:,4,n)))/Npup);
    Sy_Sim(:,:,n) = Sy_Sim(:,:,n) - Sy_Sim_ref; 
end

%% Plot slopes
plot_slopes

%% Plot maxima
maxSx_Ex = squeeze(max(max(Sx_Ex)));
maxSy_Ex = squeeze(max(max(Sy_Ex)));
minSx_Ex = squeeze(min(min(Sx_Ex)));
minSy_Ex = squeeze(min(min(Sy_Ex)));

maxSx_Sim = squeeze(max(max(Sx_Sim)));
maxSy_Sim = squeeze(max(max(Sy_Sim)));
minSx_Sim = squeeze(min(min(Sx_Sim)));
minSy_Sim = squeeze(min(min(Sy_Sim)));

c0 = 21; 
mx1 = abs(maxSx_Ex(c0+3)/poke(c0+3));
mx2 = abs(maxSx_Ex(c0-3)/poke(c0-3));
mx3 = abs(minSx_Ex(c0+3)/poke(c0+3)); 
mx4 = abs(minSx_Ex(c0-3)/poke(c0-3)); 
mx = (mx1+mx2+mx3+mx4)/4;
linX = mx*poke((c0-5):(c0+5));

my1 = abs(maxSy_Ex(c0+3)/poke(c0+3));
my2 = abs(maxSy_Ex(c0-3)/poke(c0-3));
my3 = abs(minSy_Ex(c0+3)/poke(c0+3)); 
my4 = abs(minSy_Ex(c0-3)/poke(c0-3)); 
my = (my1+my2+my3+my4)/4;
linY = my*poke((c0-5):(c0+5));


figure
    plot(poke,maxSx_Ex,'r',poke,minSx_Ex,'b','LineWidth',2)
    hold on;
    plot(poke,maxSx_Sim,'r--',poke,minSx_Sim,'b--','LineWidth',2)
    plot(poke((c0-5):(c0+5)),linX,'k')
    plot(poke((c0-5):(c0+5)),-linX,'k')
    grid on;
    legend('max. Sx: experiment','min. Sx: experimenet',...
        'max. Sx: simulation','min. Sx: simulation')
    ylim([-3 3])
    xlabel('stroke [nm]')
    ylabel('Sx')
    
figure
    plot(poke,maxSy_Ex,'r',poke,minSy_Ex,'b','LineWidth',2)
    hold on;
    plot(poke,maxSy_Sim,'r--',poke,minSy_Sim,'b--','LineWidth',2)
    plot(poke((c0-5):(c0+5)),linY,'k')
    plot(poke((c0-5):(c0+5)),-linY,'k')
    grid on;
    legend('max. Sy: experiment','min. Sy: experimenet',...
        'max. Sy: simulation','min. Sy: simulation')
    ylim([-3 3])
    xlabel('stroke [nm]')
    ylabel('Sy')

%% Slopes: using flat map, constant normalisation (taken from reference)
% Sx = I1+I3-I2-I4/(sum(sum(I1+I2+I3+I4))
% Experimental
Sx_Ex_ref = (I_Ex_ref(:,:,1)+I_Ex_ref(:,:,3)-I_Ex_ref(:,:,2)-I_Ex_ref(:,:,4))/...
    (sum(sum(I_Ex_ref(:,:,1)+I_Ex_ref(:,:,2)+I_Ex_ref(:,:,3)+I_Ex_ref(:,:,4))));
Sy_Ex_ref = (I_Ex_ref(:,:,1)+I_Ex_ref(:,:,2)-I_Ex_ref(:,:,3)-I_Ex_ref(:,:,4))/...
    (sum(sum(I_Ex_ref(:,:,1)+I_Ex_ref(:,:,2)+I_Ex_ref(:,:,3)+I_Ex_ref(:,:,4))));

for n=1:nb_poke
    Sx_Ex(:,:,n) = (I_Ex(:,:,1,n)+I_Ex(:,:,3,n)-I_Ex(:,:,2,n)-I_Ex(:,:,4,n))/...
        (sum(sum(I_Ex_ref(:,:,1)+I_Ex_ref(:,:,2)+I_Ex_ref(:,:,3)+I_Ex_ref(:,:,4))));
    Sx_Ex(:,:,n) = Sx_Ex(:,:,n) - Sx_Ex_ref; 
    Sy_Ex(:,:,n) = (I_Ex(:,:,1,n)+I_Ex(:,:,2,n)-I_Ex(:,:,3,n)-I_Ex(:,:,4,n))/...
        (sum(sum(I_Ex_ref(:,:,1)+I_Ex_ref(:,:,2)+I_Ex_ref(:,:,3)+I_Ex_ref(:,:,4))));
    Sy_Ex(:,:,n) = Sy_Ex(:,:,n) - Sy_Ex_ref; 
end

% Sim
Sx_Sim_ref = (I_Sim_ref(:,:,1)+I_Sim_ref(:,:,3)-I_Sim_ref(:,:,2)-I_Sim_ref(:,:,4))/...
    (sum(sum(I_Sim_ref(:,:,1)+I_Sim_ref(:,:,2)+I_Sim_ref(:,:,3)+I_Sim_ref(:,:,4))));
Sy_Sim_ref = (I_Sim_ref(:,:,1)+I_Sim_ref(:,:,2)-I_Sim_ref(:,:,3)-I_Sim_ref(:,:,4))/...
    (sum(sum(I_Sim_ref(:,:,1)+I_Sim_ref(:,:,2)+I_Sim_ref(:,:,3)+I_Sim_ref(:,:,4))));

for n=1:nb_poke
    Sx_Sim(:,:,n) = (I_Sim(:,:,1,n)+I_Sim(:,:,3,n)-I_Sim(:,:,2,n)-I_Sim(:,:,4,n))/...
        (sum(sum(I_Sim_ref(:,:,1)+I_Sim_ref(:,:,2)+I_Sim_ref(:,:,3)+I_Sim_ref(:,:,4))));
    Sx_Sim(:,:,n) = Sx_Sim(:,:,n) - Sx_Sim_ref; 
    Sy_Sim(:,:,n) = (I_Sim(:,:,1,n)+I_Sim(:,:,2,n)-I_Sim(:,:,3,n)-I_Sim(:,:,4,n))/...
        (sum(sum(I_Sim_ref(:,:,1)+I_Sim_ref(:,:,2)+I_Sim_ref(:,:,3)+I_Sim_ref(:,:,4))));
    Sy_Sim(:,:,n) = Sy_Sim(:,:,n) - Sy_Sim_ref; 
end

%% Plot slopes
plot_slopes

%% Plot maxima
maxSx_Ex = squeeze(max(max(Sx_Ex)));
maxSy_Ex = squeeze(max(max(Sy_Ex)));
minSx_Ex = squeeze(min(min(Sx_Ex)));
minSy_Ex = squeeze(min(min(Sy_Ex)));

maxSx_Sim = squeeze(max(max(Sx_Sim)));
maxSy_Sim = squeeze(max(max(Sy_Sim)));
minSx_Sim = squeeze(min(min(Sx_Sim)));
minSy_Sim = squeeze(min(min(Sy_Sim)));

figure
    plot(poke,maxSx_Ex,'b',poke,minSx_Ex,'r')
    hold on;
    plot(poke,maxSx_Sim,'b--',poke,minSx_Sim,'r--')
    grid on;
    
figure
    plot(poke,maxSy_Ex,'b',poke,minSy_Ex,'r')
    hold on;
    plot(poke,maxSy_Sim,'b--',poke,minSy_Sim,'r--')
    grid on;
    