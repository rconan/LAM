% Plot slopes
tmp_Sx_Ex = zeros(D,D*length(ps));
tmp_Sy_Ex = zeros(D,D*length(ps));
tmp_Sx_Sim = zeros(D,D*length(ps));
tmp_Sy_Sim = zeros(D,D*length(ps));
for i=1:length(ps)
    tmp_Sx_Ex(:,(i-1)*D+1:i*D) = Sx_Ex(:,:,ps(i));
    tmp_Sy_Ex(:,(i-1)*D+1:i*D) = Sy_Ex(:,:,ps(i));
    tmp_Sx_Sim(:,(i-1)*D+1:i*D) = Sx_Sim(:,:,ps(i));
    tmp_Sy_Sim(:,(i-1)*D+1:i*D) = Sy_Sim(:,:,ps(i));
end

clim = max(max(abs([tmp_Sx_Ex tmp_Sy_Ex tmp_Sx_Sim tmp_Sy_Sim])));

%scrsz = get(groot,'ScreenSize');
%figure('Position',[1 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2])
figure
    imagesc([tmp_Sx_Sim(:,end/2+1:end); tmp_Sx_Ex(:,end/2+1:end)])
    axis equal
    axis tight
    axis off
    set(gca,'YDir','normal')
    set(gca,'clim',[-clim clim])
    colormap(cmap)

figure
    imagesc([tmp_Sx_Sim(:,1:end/2); tmp_Sx_Ex(:,1:end/2)])
    axis equal
    axis tight
    axis off
    set(gca,'YDir','normal')
    set(gca,'clim',[-clim clim])
    colormap(cmap)

figure
    imagesc([tmp_Sy_Sim(:,end/2+1:end); tmp_Sy_Ex(:,end/2+1:end)])
    axis equal
    axis tight
    axis off
    set(gca,'YDir','normal')
    set(gca,'clim',[-clim clim])
    colormap(cmap)

figure
    imagesc([tmp_Sy_Sim(:,1:end/2); tmp_Sy_Ex(:,1:end/2)])
    axis equal
    axis tight
    axis off
    set(gca,'YDir','normal')
    set(gca,'clim',[-clim clim])
    colormap(cmap)
    
