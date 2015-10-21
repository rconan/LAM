function Zern = CARTES_2DZER(Zernik,k,DiamPix,masque)

Zern = zeros(DiamPix); 
zer = Zernik.p(:,k); 
Zern(masque == 1) = zer(masque);
figure(k)
imagesc(Zern);
axis([1 DiamPix 1 DiamPix],'square'); caxis([-3 3]); ylabel(colorbar,'');
display(['Zernike n°',num2str(k),': Max = ',num2str(max(zer(masque))),' ; Min = ',num2str(min(zer(masque)))]);
%saveas(gcf,strcat('ZER_n°',num2str(k)),'jpg');
