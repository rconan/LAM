function [G] = tool_compose_G_matrix(how_many_components)
%%%% Implementation and improvement of the article "Icreased sky coverage with optimal correction of tilt-anisoplanatism ..." by Correia

gamma_T = 2;                           %% the average slope produced by TipTilt
gamma_F = 16*sqrt(3)/(3*pi); %% the average slope produced by Focus mode
gamma_A = 8*sqrt(6)/(3*pi);   %% the average slope produced by Astigmatism

G_TT1 = [gamma_T, 0; 0, gamma_T];
G_TT2 = [gamma_T, 0; 0, gamma_T];

G_TTFA_Tiltpart = [  G_TT1; G_TT1; G_TT1; G_TT1 ];
G_TTFA_Astigmatismpart = [  gamma_F,    gamma_A,    gamma_A ; ...
                                                   gamma_F,   gamma_A,    -gamma_A ; ...
                                                   -gamma_F,  gamma_A,    -gamma_A ; ...
                                                   gamma_F,   -gamma_A,   -gamma_A ; ...
                                                   -gamma_F,  -gamma_A,   gamma_A ; ...
                                                   -gamma_F,  -gamma_A,   gamma_A ; ...
                                                   gamma_F,  -gamma_A,   gamma_A ; ...
                                                   -gamma_F,  gamma_A,   gamma_A ];
G_TTFA = [G_TTFA_Tiltpart, G_TTFA_Astigmatismpart];



switch how_many_components
    case 1
            G = blkdiag(G_TT1);  %%% this is modal matrix that translates modal coeffs from TT, TA, and TTFA modes into average slopes over the illuminatedsubregion of each subaperture.
        
    case 2
            G = blkdiag(G_TT1, G_TT2);  %%% this is modal matrix that translates modal coeffs from TT, TA, and TTFA modes into average slopes over the illuminatedsubregion of each subaperture.

    case 3
            G = blkdiag(G_TT1, G_TT2, G_TTFA);  %%% this is modal matrix that translates modal coeffs from TT, TA, and TTFA modes into average slopes over the illuminatedsubregion of each subaperture.
    
end %% for switch how_many_components