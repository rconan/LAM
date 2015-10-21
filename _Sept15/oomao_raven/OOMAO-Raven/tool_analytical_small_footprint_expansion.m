%> @file tool_analytical_small_footprint_expansion.m
%> @brief A script that computes the Zernike temporal correlation functions using the analytical angular covariance routines.
%> @author Carlos Correia, improved by Mikhail Konnik
%> @date   June 2012, revised May 2015 (based on AnalyticalSmallFootprintExpansion.m)
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
%======================================================================
%> @param zernModeMaxProj   = maximum order of Zernike function [integer number]
%> @param tel               = telescope object.
%> @param gs                = the guide-star asterism
%> @param atm               = atmosphere
%> @retval Proj             = Projection matrix with [nGs x Nz] by [Nz x Nlayers] dimensions
%> @retval Projcell         = Cell version of Proj
% ======================================================================
function [Proj Projcell] = tool_analytical_small_footprint_expansion(zernModeMaxProj,tel,gs,atm)
nMode = zernModeMaxProj-1;          %% number of Zernike modes.
 
    A2N   =  ANSI2Noll(nMode); %% converts the ANSI order of the Zernike modes into Noll convention

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
            C21 = TransformC( [D0 (vec0(1:kMode)) ], Dh, tx, ty, 0); %% this function is from the article 
           
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