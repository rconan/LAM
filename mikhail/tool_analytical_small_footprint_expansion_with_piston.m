function Projection_cell = tool_analytical_small_footprint_expansion_with_piston(zernWfs, altitudes, src, tel, M, N)
nMode = zernWfs.nMode;          
A2N   =  ANSI2Noll(nMode); %% converts the ANSI order of the Zernike modes into Noll convention
nMode = length(A2N)-1; %% number of Zernike modes.

re_sorting_vector = zeros(1,nMode);%% preallocation of a re-sorting vector for the Zernike modes.
Dh = tel.D*1000;

Proj = zeros((nMode+1), (nMode+1));
Projection_cell = cell(N,M); %%% Cell array of size NumberOfNGS x number of atmospheric layers.
 
for kM = 1:M  %% iterating for each atmosphere layer atm.nLayer;

    for kN = 1:N%% iterating through Guide Stars

                D0 = tel.diameterAt(altitudes(kM))*1000; %% telescope's diameter at an altitude of kM, changes with the altitude!
                delta = altitudes(kM)*tan(src(kN).zenith)*[cos(src(kN).azimuth),sin(src(kN).azimuth)]; %displacement due to altitude and FoV of a telescope

                tx = delta(1)*1000; %% to convert the tx into millimeters, which is standard unit in TransformC function
                ty = delta(2)*1000; %% to convert the ty into millimeters, which is standard unit in TransformC function

 
                for kMode = 2:nMode+1
                    vec0 = zeros(1,nMode+1);
                    vec0(kMode) = 1; %% considering only current mode, thus 1 in vec0 corresponds to kMode (current Mode)

                    C21 = TransformC( [D0 (vec0(1:kMode)) ], Dh, tx, ty, 0); %% this function is from the article 
                    C2x = zeros(nMode+1,1);
                    C2x(1:kMode) = C21(2:kMode+1); %% take all elements from C21 since TransformC(1) is a pupil diameter
                    C2x_converted = A2N*C2x;

                    re_sorting_vector(kMode) = find(C2x_converted, 1, 'last');  %% finds the last index of non-zero entry in the vector converted_tmp; 
                    Proj(:,kMode) = C2x_converted;
                end

                [~, resort_indices] = sort(re_sorting_vector); %% we need a vector of resort_indices that will recover the correct order in ProjectionMatrix below.
                Projection_cell{kN, kM} = Proj(:,resort_indices); %% Store the projected coefficients into Cell array

    end %%     for kGuidestar = 1:nGs %% iterating through Guide Stars    
end %% for kLayer = 1:atm.nLayer  %% iterating for each atmosphere layer nL = atm.nLayer;
