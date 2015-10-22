function Carte_2D = CARTES_2D(Carte_1D,masque_2D)

Carte_2D = zeros(size(masque_2D,1));
Carte_2D(masque_2D == 1) = Carte_1D;