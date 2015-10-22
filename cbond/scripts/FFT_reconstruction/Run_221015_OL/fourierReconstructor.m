classdef fourierReconstructor < handle
    %% fastFourierTransformReconstructor 
    % Wavefront reconstructor computed in the fourier domain.
    %
    % Required inputs: wfs, tel, atm
    %
    % Optional inputs       | Defaults
    %----------------------------------------------------------------------
    % dmFunction            | []
    % stroke                | 0 (must specify a stroke if dmFunction is
    %                       | defined)
    % extension             | 3
    % extensionMethod       | 'Gerchberg'
    % GerchbergFilter       | {}        (if method =/= 'Gerchberg')
    %                       | 'filter'  (i.e. same as recon filter)
    % GerchbergIterations   | 6
    % filter                | 'exact'
    % modeRemoval           | []
    %                       | A: global waffle
    %                       | B: local waffle
    %                       | C: piston
    % interpolationFactor   | 4
    % 
    
    
    
    properties
        tag = 'fourierReconstructor';
        % Wavefront sensor (Shack-Hartmann object)
        wfs;
        % Telescope object
        tel;
        % Atmosphere object
        atm;
        % Parameters (resolution, no. of pixels etc.)
        params;
        % Method for extending data
        extensionMethod;
        % Extension of data (no. of pixels)
        extension;
        % Filter used for reconstrction
        filter;
        % Coefficient for noise within filter (i.e. gamma in literature)
        kNoise;
        % Variance of the noise.
        noiseVar;
        % Coefficients for aliasing (i.e. factors of aliased PSD included
        % in the filter)
        kAliasing;
        % Shifts to increase corespondence between model and measurement.
        % Should be a vector [xshift yshift]
        shifts;
        % Mode removal
        modeRemoval;
        doModeRemoval;
        % Shannon Interpolation factor (for the interpolation of low
        % resolution phase to high resolution phase)
        interpolationFactor
        % Valid actuators (for padded data)
        validActuator;
        % Valid lenslets (for padded data)
        validLenslet;
        % Number of data points for reconstructed phase
        N;
        % Filters (Fourier domain) for x and y slopes
        RX;
        RY;
        % wfs response to spatial frequencies
        GX;
        GY;
        % FFT filter coresponding to DM influence function
        dm;
        % Gerchberg filter used for extension (if Gerchberg is used as
        % extension method (should this just be same as filter?)
        Gerchberg;
        
        % N.B. add dependent variables at some point (ones which will
        % update)
        
        % Orginal slopes
        xSlopes;
        ySlopes;
        % Extended slopes
        Sx;
        Sy;
        
    end
    
    methods
        
        %% Constructor
        function obj = fourierReconstructor(wfs,tel,atm,varargin)
           
            inputs = inputParser;
            inputs.addRequired('wfs',@(x) isa(x,'shackHartmann') || isa(x,'geomGrad'));
            inputs.addRequired('tel',@(x) isa(x,'telescope'));
            inputs.addRequired('atm',@(x) isa(x,'atmosphere'));
            inputs.addOptional('dm',[],@(x) isa(x,'deformableMirror'));
            inputs.addOptional('stroke',0,@isnumeric);
            inputs.addOptional('extension',3,@isnumeric);
            inputs.addOptional('extensionMethod','simple',@ischar);
            inputs.addOptional('GerchbergFilter',[],@ischar);
            inputs.addOptional('GerchbergIterations',6,@isnumeric);
            inputs.addOptional('filter','exact',@ischar);
            inputs.addOptional('kNoise',[],@isnumeric);
            inputs.addOptional('noiseVariance',0,@isnumeric);
            inputs.addOptional('kAliasing',[],@isnumeric);
            inputs.addOptional('shifts',[],@isnumeric);
            inputs.addOptional('modeRemoval',[],@isnumeric);
            inputs.addOptional('interpolationFactor',4,@isnumeric);
            
            inputs.parse(wfs,tel,atm,varargin{:});
            
            % Required parameters
            obj.wfs                 = inputs.Results.wfs;
            obj.tel                 = inputs.Results.tel;
            obj.atm                 = inputs.Results.atm;
            
            % Parameters
            obj.params.nLenslet     = inputs.Results.wfs.lenslets.nLenslet;
            obj.params.vActuator    = inputs.Results.wfs.validActuator;
            obj.params.vLenslet     = inputs.Results.wfs.validLenslet;
            obj.params.nPx          = inputs.Results.tel.resolution;
            obj.params.stroke       = inputs.Results.stroke;
            obj.params.r0           = inputs.Results.atm.r0;
            obj.params.L0           = inputs.Results.atm.L0;
            obj.params.d            = inputs.Results.tel.D/obj.params.nLenslet;
            obj.extensionMethod     = inputs.Results.extensionMethod;
            obj.extension           = inputs.Results.extension;
            obj.filter              = inputs.Results.filter;
            obj.noiseVar            = inputs.Results.noiseVariance;
            obj.shifts              = inputs.Results.shifts;
            obj.modeRemoval         = inputs.Results.modeRemoval;
            obj.interpolationFactor = inputs.Results.interpolationFactor;
            
            % If we have an even number of lenslets this coresponds to an
            % odd number of actuator points.  We want to always use an even
            % number of points as it makes filter out frequencies such as
            % the waffle easier.
            nL                  = obj.params.nLenslet;
            if obj.params.nLenslet/2==round(obj.params.nLenslet/2)
                obj.N               = obj.params.nLenslet + 2 + 2*obj.extension;
                % First pad with 1 additional point
                % For actuators padd first row, column (have data centre at
                % N/2 +1, for NxN even grid)
                % For lenslets add zeros to all sides: padding to make
                % lenslet and slope arrays consistent with size of actuator
                % array.
                obj.validActuator   = logical([zeros(1,nL+2); zeros(nL+1,1) obj.params.vActuator]);
                obj.validLenslet    = logical([zeros(1,nL+2); zeros(nL,1) obj.params.vLenslet zeros(nL,1); zeros(1,nL+2)]);
                % Then pad to extension
                obj.validActuator   = logical(fourierReconstructor.padarray(obj.validActuator,obj.extension));
                obj.validLenslet    = logical(fourierReconstructor.padarray(obj.validLenslet,obj.extension));
            else
                obj.N               = obj.params.nLenslet + 1 + 2*obj.extension;
                % No additional padding for actuator array (we already have
                % an even number of points)
                % Extend lenslets to no. of actuator points
                obj.validLenslet    = logical([obj.params.vLenslet zeros(nL,1) ;zeros(1,nL+1);]);
                % Pad actuator and lenslet arrays
                obj.validActuator   = logical(fourierReconstructor.padarray(obj.params.vActuator,obj.extension));
                obj.validLenslet    = logical(fourierReconstructor.padarray(obj.validLenslet,obj.extension));
            end
            
            % If a noise variance is specified and no value is given for
            % the noise coefficient, it is set to 1.
            if isempty(inputs.Results.kNoise) && obj.noiseVar~=0
                obj.kNoise          = 1;
            else
                obj.kNoise          = inputs.Results.kNoise;
            end
                
            % If the aliasing filter is specified but no aliasing
            % coefficient is given it is set to 1.
            if isempty(inputs.Results.kAliasing) && strcmp(obj.filter,'antialias')
                obj.kAliasing       = 1;
            else
                obj.kAliasing       = inputs.Results.kAliasing;
            end
            
            [RX,RY,div,GX,GY] = fftFilters(obj);
            obj.RX                  = RX;
            obj.RY                  = RY;
            obj.GX                  = GX;
            obj.GY                  = GY;

            if strcmp(obj.extensionMethod,'Gerchberg')
                if isempty(inputs.Results.GerchbergFilter)
                    obj.Gerchberg.filter = obj.filter;
                else
                    obj.Gerchberg.filter = inputs.Results.GerchbergFilter;
                end
                obj.Gerchberg.nIter      = inputs.Results.GerchbergIterations;
            end
            
            if ~isempty(obj.modeRemoval) && (obj.modeRemoval(1)==1 || ...
                obj.modeRemoval(2)==1 || obj.modeRemoval(3)==1)
                obj.doModeRemoval = 1;
            else
                obj.doModeRemoval = 0;
            end
            
            if ~isempty(inputs.Results.dm)
                obj.dm.function         = inputs.Results.dm.modes.bezier; 
                obj.dm.modes            = inputs.Results.dm.modes.modes;
                obj.dm.stroke           = inputs.Results.stroke;
                obj.dm.filter           = fourierReconstructor.dmFilter(obj.dm.function,obj.extension,obj.params.nLenslet);
                obj.dm.compensation     = 1;
            else
                obj.dm.compensation     = 0;
            end
            
            % Set up empty matrices to store x/y slopes
            obj.xSlopes = zeros(obj.params.nLenslet,obj.params.nLenslet);
            obj.ySlopes = zeros(obj.params.nLenslet,obj.params.nLenslet);
            
        end    
        

        %% Wavefront reconstruction
        function [Phi,DMcoeffs,PhiDM,PhiHR] = mtimes(obj,wfs)
            
            % Returns:
            % Phi:      lenslet+1 sampled phase
            % DMcoeffs: dm coefficients for compensation
            % PhiDM:    nPx resolution phase coresponding to phaes of the
            %           dm.
            % PhiHR:    high resolution phase equivalent to Phi, where the
            %           resolution is determined by the interpolation
            %           factor
            
            % To Do: add mode removal
            
            % Can use either SH object or the wfs slopes themselves
            if isa(wfs,'shackHartmann')
                obj.xSlopes(obj.params.vLenslet) = wfs.slopes(1:end/2);
                obj.ySlopes(obj.params.vLenslet) = wfs.slopes(end/2+1:end);
            else
                obj.xSlopes(obj.params.vLenslet) = wfs(1:end/2);
                obj.ySlopes(obj.params.vLenslet) = wfs(end/2+1:end);
            end
           
            N = obj.N;
            nL = obj.params.nLenslet;
            mask = obj.validActuator;
            
            Sx = obj.xSlopes;
            Sy = obj.ySlopes;
            % Pad slopes to number of actuators points, with extension 
            % (always pad to even number)
            if nL/2==round(nL/2)
                % Extend slopes to the same size as the no. of actuator
                % points (i.e. recreated phase points)
                Sx = [zeros(1,nL+2); zeros(nL,1) Sx zeros(nL,1); zeros(1,nL+2)];
                Sy = [zeros(1,nL+2); zeros(nL,1) Sy zeros(nL,1); zeros(1,nL+2)];
                Sx = fourierReconstructor.padarray(Sx,obj.extension);
                Sy = fourierReconstructor.padarray(Sy,obj.extension);
            else
                % Extend slopes to the same as the number of actuator
                % points
                Sx = [Sx zeros(nL,1); zeros(1,nL+1)];
                Sy = [Sy zeros(nL,1); zeros(1,nL+1)];
                Sx = fourierReconstructor.padarray(Sx,obj.extension);
                Sy = fourierReconstructor.padarray(Sy,obj.extension);
            end
            
            
            
            % Extend slopes to edges of grid (with given method)
            [Sx,Sy] = extendSlopes(Sx,Sy,obj);

            obj.Sx = Sx;
            obj.Sy = Sy;
            
            % FFT of slopes
            fft_Sx = 1/N^2 * fftshift(fft2(Sx,N,N));
            fft_Sy = 1/N^2 * fftshift(fft2(Sy,N,N));
            
            % Use filters to create FFT of phase
            fft_Phi = obj.RX.*fft_Sx + obj.RY.*fft_Sy;
            % Remove piston
            fft_Phi(floor(N/2+1),floor(N/2+1)) = 0;
            
            % Modal removal (in direct space)
            if obj.doModeRemoval==1
                % Remove selected mdoes
                %fft_Phi = fourierReconstructor.removeModes(fft_Phi,N,obj.validLenslet,obj.modeRemoval,'fourier');
                Phi = N^2*ifft2(ifftshift(fft_Phi),N,N);
                Phi = real(Phi); 
                Phi = fourierReconstructor.removeModes(Phi,N,obj.validLenslet,obj.modeRemoval,'direct');
                fft_Phi = 1/N^2 * fftshift(fft2(Phi,N,N));
                %
            else
                % Remove piston over the pupil
                %Phi = fourierReconstructor.removeModes(Phi,N,obj.validLenslet,[0 0 1],'direct');
                fft_Phi = fourierReconstructor.removeModes(fft_Phi,N,obj.validLenslet,[0 0 1],'fourier');
                % Now use inverse fourier transform to get phase
                Phi = N^2 * ifft2(ifftshift(fft_Phi),N,N);
                % Return real value of phase
                Phi = real(Phi);
            end            
            
            % Shannon interpolation to return high resolution phase
            nSh = (N-1) * obj.interpolationFactor;
            fft_PhiHR = zeros(nSh,nSh);
            % Coordinates for central frequency (high res)
            nSh0 = floor(nSh/2+1);
            % Central FFT is equivalent to FFT of phi
            fft_PhiHR(nSh0-floor(N/2):nSh0+ceil(N/2-1),nSh0-floor(N/2):nSh0+ceil(N/2-1)) = fft_Phi;
                        
            % Use inverse FFT to get high resolution phase
            PhiHR = nSh^2 * ifft2(fftshift(fft_PhiHR),nSh,nSh);
            PhiHR = real(PhiHR);
            
            % Impose pupil (i.e. only return valid measurements)
            Phi = mask.*Phi;
            % Remove extended results
            if nL/2==round(nL/2)
                Phi = Phi(2+obj.extension:N-obj.extension,2+obj.extension:N-obj.extension);
            else
                Phi = Phi(1+obj.extension:N-obj.extension,1+obj.extension:N-obj.extension);
            end
            
            % And remove extended results for high resolution phase
            M = obj.extension*obj.interpolationFactor;
            if nL/2==round(nL/2)
                PhiHR = PhiHR(1+obj.interpolationFactor+M:nSh-M,1+obj.interpolationFactor+M:nSh-M);
            else
                PhiHR = PhiHR(1+M:nSh-M,1+M:nSh-M);
            end
            
            % Calculate dm coefficients
            DMcoeffs = [];
            PhiDM = [];
            if obj.dm.compensation==1
               
                % To do: check DM filter -> understand and check 'aliasing
                % effect' is properly taken into account.
                coeffs = N^2 * ifft2(ifftshift(fft_Phi).*obj.dm.filter,N,N);
                % Take real part and calibrate with stroke (and factor of 2
                % to take reflection from dm into account)
                coeffs = obj.params.stroke/2 * real(coeffs);
                DMcoeffs = coeffs(obj.validActuator);
             
                % Equivalent phase of the deformable mirror with the given
                % coefficients
                PhiDM = obj.dm.modes * DMcoeffs;
                PhiDM = reshape(PhiDM,obj.params.nPx,obj.params.nPx);
            end
            
        end
        
        %% Calculation of filters
        function [RX, RY, div, GX, GY] = fftFilters(obj)

            % ------------HEADER-----------------
            % N.B. Need to edit comments etc.
            %
            % Filters:  
            %
            % Objective         ::  Implement several filters to reconstruct (LSR & MVR) the
            %                       phase in Fourier domain
            % INPUT VARS
            % filter            ::  the type of filter
            %                       'exact'
            %                       'Hudgin'
            %                       'modifiedHudgin'
            %                       'Freischlad86'
            %                       'Poyneer02'
            % N                 ::  The number of computational points
            % hshift, vshift    ::  The horizontal and vertical shifts
            %
            % OUTPUT VARS
            % CX                ::  Filter for the X slopes
            % CY                ::  Filter for the Y slopes
            % div               ::  The denominator, div = abs(CX)^2+abs(CY)^2
            % CXo               ::  Filter for the X slopes without being divided by
            %                       div
            % CYo               ::  Filter for the Y slopes without being divided by
            %                       div
            %Created by         ::  Carlos Correia
            %Creation date      ::  16/02/2007
            %Change Record:     ::  CCorreia, July09
            %                       - Added two output variables, to have the x and y
            %                       filters without being divided by the denominator
            %                       CCorreia, Jun09
            %                       - horizontal and vertical shifts were incorrectly computed. Added the
            %                       frequency vectors fx and fy with the correct numbering to allow for
            %                       fractional shifts

            %                       JTeixeira, July14
            %                       - Added option to Wiener filter;
            %                       - input parameter changed:
            %                           params :: structure containing parameters:
            %                              . filter
            %                              . N
            %                              . r0 & L0 to estimate PSD's
            %                              . kn & ka to the regularization parameter
            %                                coeficients, kn for noise, ka for aliasing;
            %
            %
            %       CORRECTIONS TO DO: Extra-alignment in Hudgin & Southwell is
            %                              incorrect, each filters should be shifted in
            %                              both directions x and y.

            % ------------HEADER END----------------

            N = obj.N;
            d = obj.params.d;
            % No. of pixels per lenslet
            nPxL = (obj.params.nPx/obj.params.nLenslet);
            % Averaging length (d-size of 1 pixel)
            d_av = d - d/nPxL;
            
            if isempty(obj.shifts)
                xshift = 0;
                yshift = 0;
            else
                xshift = obj.shifts(1);
                yshift = obj.shifts(2);
            end

            % Spatial frequency over defined phase points (plus extension)
            D = N*d;
            
            
            if N/2==round(N/2)
                kx = -N/2*1/D:1/D:(N/2-1)*1/D;
            else
                kx = -floor(N/2)*1/D:1/D:floor(N/2)*1/D;
            end

            ky = kx;
            [KX,KY] = meshgrid(kx,ky);


            % PSDs for noise (i.e. for Wiener filters)
            if ~isempty(obj.kNoise)
                
                % Calculate phase PSD (in band -> no aliasing)
                m = 0;
                Wphi = fourierReconstructor.createPhasePSD(kx,m,d,nPxL,obj.params.r0,obj.params.L0);
                
                % Calculate white noise PSD
                Wnoise = fourierReconstructor.createNoisePSD(kx,obj.noiseVar);
                
                % Noise PSD normalised against phase PSD
                Rnoise = obj.kNoise * Wnoise./Wphi;
            else
                Rnoise = 0;
            end


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %--------------  filters   ----------------
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if strcmp(obj.filter,'modifiedHudgin') || strcmp(obj.filter,'Hudgin')|| strcmp(obj.filter,'hudgin')
                % -------------------------------------------
                % ---------------------------------
                %----------  modified Hudgin  -----
                % ---------------------------------

                % No noise case
                Dx = exp(-2*1i*pi*d*KX)-1;
                Dy = exp(-2*1i*pi*d*KY)-1;
                Ax = 1;
                Ay = 1;
                E = exp(-1/4*1i*pi*d*(KX+KY));

                % CZ and CY with additional horizontal and vertical shifts
                GX = Dx.*Ax.*E;
                GY = Dy.*Ax.*E;
                
                div = 1./(abs(conj(Dx).*conj(Ax)).^2 + abs(conj(Dy).*conj(Ay)).^2 + Rnoise);
                %div = 1./(abs(GX).^2+abs(GY).^2);

                
            elseif strcmp(obj.filter,'exact') || strcmp(obj.filter,'Exact')
                % --------------------------------
                %----------   Exact  --------------

                % No noise case
%                 Dx = exp(-2*1i*pi*d*KX)-1;
%                 Dy = exp(-2*1i*pi*d*KY)-1;
%                 Ax = sinc(d_av*KY).*exp(-1i*pi*d*KY);
%                 Ay = sinc(d_av*KX).*exp(-1i*pi*d*KX);
%                 E = 1;
% 
%                 GX = Dx.*Ax.*E;
%                 GY = Dy.*Ay.*E;
                Dx = -2*1i*pi*d*KX;
                Dy = -2*1i*pi*d*KY;
                Ax = sinc(d_av*KX).*sinc(d*KY).*exp(-1i*pi*d*(KX+KY));
                Ay = sinc(d*KX).*sinc(d_av*KY).*exp(-1i*pi*d*(KX+KY));

                GX = Dx.*Ax;
                GY = Dy.*Ay;
                div = 1./(abs(conj(Dx).*conj(Ax)).^2 + abs(conj(Dy).*conj(Ay)).^2 + Rnoise);



            elseif strcmp(obj.filter,'Freischlad86')
                % ---------------------------------------
                %----------  Freischlad86  --------------
                % ---------------------------------------
                
                % N.B. doesn't work for dm coeffs????
                
                GX = exp(-1i*2*pi*d*KX) - 1;
                GY = exp(-1i*2*pi*d*KY) - 1;
                div = 1./(4*(sin(pi*d*KX).^2 + sin(pi*d*KY).^2) + Rnoise);
                
                
            elseif    strcmp(obj.filter,'Poyneer02') || strcmp(obj.filter,'Fried')|| strcmp(obj.filter,'fried')
                % -------------------------------------------
                % -------------------------------------------
                %----------  Poyneer02, exp. 13  ------------
                % -------------------------------------------
                % theoretical filter for the Fried geometry (See Poyneer02, 
                % exp. 13) THE APPLICATION OF LISA'S FILTER TO THE FRIED 
                % GEOMETRY AFTER LINEAR RECOMBINATION OF GRIDS (EXP. 18 AND
                % 19) GIVE VERY ACCURATE RESULTS, CLOSE OR EVEN BETTER THAN
                % THOSE OBTAINED WITH FREISCHLAD'S FILTER...

                % TO DO: add back in LISA's addition to this filter (waffle
                % removal) done
                % TO DO: check shifts

                Dx = exp(-2*1i*pi*d*KX)-1;
                Dy = exp(-2*1i*pi*d*KY)-1;
                Ax = (1/2)*(1+exp(-2*1i*pi*d*KY));
                Ay = (1/2)*(1+exp(-2*1i*pi*d*KX));

                % CX and CY with additional horizontal and vertical shifts
                GX = Dx.*Ax;
                GY = Dy.*Ay;

                div = 1./(abs(conj(Dx).*conj(Ax)).^2 + abs(conj(Dy).*conj(Ay)).^2 + Rnoise);
                
                if N/2==round(N/2)
                    % If even number of points remove waffle mode
                    GX(1,1) = 0;
                    GY(1,1) = 0;
                end

            elseif strcmp(obj.filter,'DCT')
        %         % TO DO: what is this filter?
        %         % re-write in my format
        %         for k=1:N
        %             for l=1:N
        %                 CX(k,l) =  (exp(-imath*pi*(k-1)/N)-1); %2* (sin(pi*(k-1)/(2*N)));
        %                 CY(k,l) = (exp(-imath*pi*(l-1)/N)-1); %2* (sin(pi*(l-1)/(2*N)));%
        %                 div(k,l) = (4*((sin(pi*(k-1)/(2*N)))^2 + (sin(pi*(l-1)/(2*N)))^2))^(-1);
        %             end
        %         end


            elseif strcmp(obj.filter,'Southwell')|| strcmp(obj.filter,'southwell')
                %----------  Southwell_80  --------------

                Dx = exp(-2*1i*pi*d*KX)-1;
                Dy = exp(-2*1i*pi*d*KY)-1;
                Ax = 1;
                Ay = 1;
                Sx = (1/2)*(1 + exp(2*1i*pi*d*KX));
                Sy = (1/2)*(1 + exp(2*1i*pi*d*KY));
                E = exp(-1i*pi*d*(KX+KY));

                % CX and CY with additional horizontal and vertical shifts
                GX = Dx.*Sx.*E;
                GY = Dy.*Sy.*E;
        
                div = 1./(abs(conj(Dx).*conj(Ax)).^2 + abs(conj(Dy).*conj(Ay)).^2 + Rnoise);


            elseif strcmp(obj.filter,'antialias')
                % -------------------------------------------
                % -------------------------------------------
                % --------------------------------------
                %   ANTI-ALIASING filter
                % --------------------------------------
                
                % CX and CY filters (from exact approach)
                %CX = (-2*1i*pi*d*KX) .* sinc(KX*d) .* sinc(KY*d) .* exp(-1i*pi*d*(KX+KY));
                %CY = (-2*1i*pi*d*KY) .* sinc(KY*d) .* sinc(KX*d) .* exp(-1i*pi*d*(KX+KY));
%                 Dx = exp(-2*1i*pi*d*KX)-1;
%                 Dy = exp(-2*1i*pi*d*KY)-1;
%                 Ax = sinc(d_av*KY).*exp(-1i*pi*d*KY);
%                 Ay = sinc(d_av*KX).*exp(-1i*pi*d*KX);
%                 E = 1;
% 
%                 GX = Dx.*Ax.*E;
%                 GY = Dy.*Ay.*E;

                Dx = -2*1i*pi*d*KX;
                Dy = -2*1i*pi*d*KY;
                Ax = sinc(d_av*KX).*sinc(d*KY).*exp(-1i*pi*d*(KX+KY));
                Ay = sinc(d*KX).*sinc(d_av*KY).*exp(-1i*pi*d*(KX+KY));

                GX = Dx.*Ax;
                GY = Dy.*Ay;

                % PSD of the phase (and X and Y slopes) in band (no
                % aliasing)
                [Wphi,WSx0,WSy0] = fourierReconstructor.createPhasePSD(kx,0,d,nPxL,obj.params.r0,obj.params.L0);
                % PSDs of the slopes including aliasing components.  Using
                % m = 2 aliasing factor (two folds in) 
                m = 2;
                [WphiAA,WSx,WSy] = fourierReconstructor.createPhasePSD(kx,m,d,nPxL,obj.params.r0,obj.params.L0);
                % The slope PSDs with only the aliasing component
                WSxa = WSx - WSx0;
                WSya = WSy - WSy0;
                
                div = 1./(abs(GX).^2 + abs(GY).^2 + ...
                    (obj.kAliasing*(WSxa+WSya)./Wphi) + Rnoise);
                
            elseif strcmp(obj.filter,'e2eAA')
                % -------------------------------------------
                % -------------------------------------------
                % --------------------------------------
                %   Filter using slope PSD measured from e2e simulation
                % --------------------------------------
                
                Dx = -2*1i*pi*d*KX;
                Dy = -2*1i*pi*d*KY;
                Ax = sinc(d_av*KX).*sinc(d*KY).*exp(-1i*pi*d*(KX+KY));
                Ay = sinc(d*KX).*sinc(d_av*KY).*exp(-1i*pi*d*(KX+KY));

                GX = Dx.*Ax;
                GY = Dy.*Ay;

                % PSD of the phase (and X and Y slopes) in band (no
                % aliasing)
                [Wphi,WSx0,WSy0] = fourierReconstructor.createPhasePSD(kx,0,d,nPxL,obj.params.r0,obj.params.L0);
                % load e2e calculated slopes
                load('e2e_slopes')
                c = 0.9;
                SxAA = c*SxAA - WSx0;
                SyAA = c*SyAA - WSy0;
                
                div = 1./(abs(GX).^2 + abs(GY).^2 + ...
                    (obj.kAliasing*(SxAA+SyAA)./Wphi) + Rnoise);
                
                
            end

            
            % Normalising
            RX = GX.*div;
            RY = GY.*div;
            
            % Remove piston
            RX(floor(N/2+1),floor(N/2+1)) = 0;
            RY(floor(N/2+1),floor(N/2+1)) = 0;

        end
        
        %% Function to measure wfs filter
        function [Rx,Ry,Gx,Gy] = measureFilter(obj)
            
            % To do: check capability for extended phase
            
            
            % Number of points
            N = obj.N;
            nPx = obj.params.nPx;
            wfs = obj.wfs;
            tel = obj.tel;
            atm = obj.atm;
            
            % x (n) and y (m) vectors in indices (high resolution phase)
            [n,m] = meshgrid(0:nPx-1,0:nPx-1);
            
            % Define x (l) and y (k) frequency vectors.  Slightly different
            % if we have even or odd number of points
            if N/2==round(N/2)
                % N even
                l = -N/2:(N/2-1);
                k = -N/2:(N/2-1);
                % Central index
                i0 = N/2+1;
            else
                % N odd
                l = -(N-1)/2:(N-1)/2;
                k = -(N-1)/2:(N-1)/2;
                % Central index
                i0 = (N+1)/2;
            end
            
            % Number of points corresponding to high resolution phase for
            % actuator points (i.e. nPx is high resolution points for phase
            % over number of lenslets)
            nM = nPx*N/(N-1);
            
            % Check atmosphere is removed from telescope (we will propagate
            % phase defined by spatial frequencies)
            tel = tel-atm;
            
            % Cycle through each spatial frequency, propagating this
            % through the wfs
            for i=1:N
                for j=1:N
                    % Initiate slope vectors
                    Sx = zeros(N,N);
                    Sy = zeros(N,N);
                    
                    % Normalisation factors (depend on which frequency we
                    % are looking at).  In cases with combinations of 0
                    % frequency of waffle frequency the sin term goes to
                    % zero
                    if N/2==round(N/2)
                        % N even
                        if (i==1 && j==1) || (i==1 && j==i0) || (i==i0 && j==1) || (i==i0 && j==i0)
                            % If sine term is 0 the normalisation constant
                            % is 1 and we only have one frequency in fft
                            % (i.e. no negative/positive counterpart)
                            Dkl = 1;
                            M = 1;
                        else
                            % Otherwise normalisation is sqrt(2) and we
                            % take into account power split over
                            % positive/negative frequencies
                            Dkl = sqrt(2);
                            M = 2;
                        end
                    else
                        % N odd (only no sine term for DC term)
                        if (i==i0 && j==i0)
                            Dkl = 1;
                            M = 1;
                        else
                            Dkl = sqrt(2);
                            M = 2;
                        end
                    end
                    
                    % Eigen function (only use cos term)
                    C_kl = Dkl/N * cos(2*pi/nM*(-l(i)*n-k(j)*m));
                    
                    % Propagate ngs through telescope and add custom phase.
                    % Then propagate to wfs
                    ngs = ngs.*+tel;
                    ngs.phase = C_kl;
                    ngs = ngs*wfs;
                    
                    % Select the valid measurement of the slopes
                    Sx(wfs.validLenslet) = wfs.slopes(1:end/2);
                    Sy(wfs.validLenslet) = wfs.slopes(end/2+1:end);
                    
                    % Pad data with 1 extra data point (going from N-1
                    % measurements to N actuator points)
                    Sx = [Sx zeros(N-1,1);zeros(1,N)];
                    Sy = [Sy zeros(N-1,1);zeros(1,N)];
                    [Sx,Sy] = extendSlopes(Sx,Sy,FFTR);  
                    
                    % Take fourier transform of slopes
                    fftSx = (1/N)*fftshift(fft2(Sx,N,N));
                    fftSy = (1/N)*fftshift(fft2(Sy,N,N));
                    
                    % Record response for given frequency
                    Gx(j,i) = M*(1/Dkl)*fft(j,i);
                    GY(j,i) = M*(1/Dkl)*fftSy(j,i);
                    
                end
            end
            
            
            % Calculate filter from response functions
            Rx = conj(Gx)./(abs(Gx).^2+abs(Gy).^2);
            Ry = conj(Gy)./(abs(Gx).^2+abs(Gy).^2);
            
            % Zero the 0th frequency component
            Rx(i0,i0)=0;
            Ry(i0,i0)=0;
            
            
        end
        
        %%
        function [Sx,Sy] = extendSlopes(Sx,Sy,obj)

            % Extension methods             | Description
            %-----------------------------------------------------------
            % 'simple'/'simpleextension'    | 
            %
            % 'FriedClosedExtension'/       |
            % 'mesh2x2'                     |
            %
            % 
            mask = obj.validLenslet;

            if strcmp(obj.extensionMethod,'simple') || strcmp(obj.extensionMethod,'simpleextension')
                %----------------------------------------------------------------------
                % Extension method: simple
                % Extends the slopes along x (for x slopes) and along y (for y slopes)
                % using the last value in the valid apperture
                %----------------------------------------------------------------------
                                
                for i = 1:length(mask)
                    % Extend x slopes up and down out of aperture
                    idx = find(mask(:,i));
                    Sx(1:min(idx)-1,i) = Sx(min(idx),i);
                    Sx(max(idx)+1:end,i) = Sx(max(idx),i);
                    
                    % Extend y slopes left and right out of aperture
                    idy = find(mask(i,:));
                    Sy(i,1:min(idy)-1) = Sy(i,min(idy));
                    Sy(i,max(idy)+1:end) = Sy(i,max(idy));
                end

            elseif strcmp(obj.extensionMethod,'FriedClosedExtension') ||strcmp(obj.extensionMethod,'mesh2x2')
                %-------------------------------------
                % Extension method for the FRIED geometry
                %-------------------------------------
                N = length(mask);
                masktmp = mask;
                surrounded = 1;
                %for debug purposes
                xs = Sx;
                ys = Sy;
                while surrounded > 0
                    surrounded = 0;
                    for i = 2:length(mask)-1
                        % --- X slopes ---
                        idx = find(masktmp(i,:));
                        %solve conflicts - a point that is surronded by three other inside
                        %the aperture
                        if ~isempty(idx) %if subapertures have been found!
                            point = min(idx) - 1;
                            [Sx,Sy,nn] = solve_conflicts(i,point,N,masktmp,Sx,Sy);
                            surrounded = surrounded + nn;
                            point = max(idx) + 1;
                            [Sx,Sy,nn] = solve_conflicts(i,point,N,masktmp,Sx,Sy);
                            surrounded = surrounded + nn;
                        end
                    end
                    %enlarge the new masktmp to account for the conflict-free
                    %subapertures
                    maskidx = find(Sx+Sy);
                    masktmp(maskidx) = 1;
                end
                for i = 2:length(masktmp)-1
                    % --- Y slopes ---
                    idx = find(masktmp(:,i));
                    Sx(1:min(idx)-1,i) = Sx(min(idx),i);
                    Sx(max(idx)+1:end,i) = Sx(max(idx),i);
                    %create the waffle-like pattern in the y-slope direction
                    if ~isempty(idx)
                        pos = 1:min(idx)-1;
                        Sy(1:min(idx)-1,i) = Sy(min(idx),i)*pattern(length(pos),'inv')';
                        pos = max(idx)+1:N;
                        Sy(max(idx)+1:end,i) = Sy(max(idx),i)*pattern(length(pos))';
                    end

                    % --- Y slopes ---
                    idy = find(masktmp(i,:));
                    Sy(i,1:min(idy)-1) = Sy(i,min(idy));
                    Sy(i,max(idy)+1:end) = Sy(i,max(idy));
                    %create the waffle-like pattern in the x-slope direction
                    if ~isempty(idy)
                        pos = 1:min(idy)-1;
                        Sx(i,1:min(idy)-1) = Sx(i,min(idy))*pattern(length(pos),'inv');
                        pos = max(idy)+1:N;
                        Sx(i,max(idy)+1:end) = Sx(i,max(idy))*pattern(length(pos));
                    end
                end

                maskidx = find(Sx+Sy);
                masktmp(maskidx) = 1;
                surrounded = 1;
                while surrounded > 0
                    surrounded = 0;
                    for i = 1:length(mask)
                        % --- X slopes ---
                        idx = find(masktmp(i,:));
                        %solve conflicts - a point that is surronded by three other inside
                        %the aperture
                        if ~isempty(idx) & length(idx) < N%if subapertures have been found!
                            point = min(idx) - 1;
                            [Sx,Sy,nn] = solve_conflicts(i,point,N,masktmp,Sx,Sy);
                            surrounded = surrounded + nn;
                            point = max(idx) + 1;
                            [Sx,Sy,nn] = solve_conflicts(i,point,N,masktmp,Sx,Sy);
                            surrounded = surrounded + nn;
                        end
                    end
                    %enlarge the new masktmp to account for the conflict-free
                    %subapertures
                    maskidx = find(Sx+Sy);
                    masktmp(maskidx) = 1;
                end

            elseif strcmp(obj.extensionMethod,'HudginSimpleExtension') || strcmp(obj.extensionMethod,'Hudginsimpleextension')
                %-------------------------------------
                % Extension method for the HUDGIN geometry
                %   (takes into account three-slopes subapertures)
                %-------------------------------------
                for i = 2:length(mask)-1
                    idx = find(mask(:,i));
                    Sx(1:min(idx)-1,i) = Sx(min(idx),i);

                    idy = find(mask(i,:));
                    Sy(i,1:min(idy)-1) = Sy(i,min(idy));
                    if (mask(max(idx),i+1))%if there is a right-neighbour
                        if isempty(idx)
                            idx = 1;
                        end
                        Sx(max(idx)+1:end,i) = Sx(max(idx),i) - Sy(max(idx),i) + Sy(max(idx),i+1);
                    else
                        Sx(max(idx)+1:end,i) = Sx(max(idx),i);
                    end
                    if (mask(i+1,max(idy)))%if there is a bottom-neighbour
                        if isempty(idy)
                            idy = 1;
                        end
                        Sy(i,max(idy)+1:end) = Sy(i,max(idy)) - Sx(i,max(idy)) + Sx(i+1,max(idy));
                    else
                        Sy(i,max(idy)+1:end) = Sy(i,max(idy));
                    end
                end

            elseif strcmp(obj.extensionMethod,'FriedSkewedExtension')
                %-------------------------------------
                % Extension method for the Fried geometry...
                %-------------------------------------
                [B,d] = spdiags(Sy);
                for i=1:size(B,2)
                    idx = find(B(:,i));
                    B(1:min(idx)-1,i) = B(min(idx),i);
                    B(max(idx)+1:end,i) = B(max(idx),i);
                end
                %   --- re-shape the B matrix so it resembles a circle again and not
                %   an elipse ---
                %     shift = -27:27;
                %     for i=1:size(B,2)
                %         B(:,i) = circshift(B(:,i), floor(-shift(i)/2));
                %     end
                Sy = spdiags(B,d,Sy);
                Sy = full(Sy);
                circ = 1:2:length(mask);
                circ = [circ circ];
                for i = 1:length(mask)
                    Sx(i,:) = Sx(i,length(Sx):-1:1);
                end
                [B,d] = spdiags(Sx);
                for i=1:size(B,2)
                    idx = find(B(:,i));
                    B(1:min(idx)-1,i) = B(min(idx),i);
                    B(max(idx)+1:end,i) = B(max(idx),i);
                end
                Sx = spdiags(B,d,Sx);
                Sx = full(Sx);
                for i = 1:length(mask)
                    Sx(i,:) = Sx(i,length(Sx):-1:1);
                end

            elseif strcmp(obj.extensionMethod,'Gerchberg')
                %-------------------------------------
                % Gerchberg
                %-------------------------------------
                
                [Sx, Sy] = GerchbergMethod(Sx,Sy,obj);
                
                
            elseif strcmp(obj.extensionMethod,'edge')
                %-------------------------------------
                % Edge correction (Poyneer05)
                % Requires at least 1 point of padding to work
                % N.B. produce error message if no padding
                %-------------------------------------
                for i = 1:length(mask)
                    idx = find(mask(i,:));
                    tmp = -1/2*sum(Sx(i,idx));
                    Sx(i,min(idx)-1) = tmp;
                    Sx(i,max(idx)+1) = tmp;

                    idy = find(mask(:,i));
                    tmp = -1/2*sum(Sy(idx,i));
                    Sy(min(idy)-1,i) = tmp;
                    Sy(max(idy)+1,i) = tmp;
                end
                
            elseif strcmp(obj.extensionMethod,'new')
                % Also requires at least ext = 1
                % Same idea as simple extension, but in this case the x
                % slopes are extended in x and the y slopes are extended in
                % y
                for i = 1:length(mask)
                    % Extending x slopes left to right
                    idx = find(Sx(i,:));
                    Sx(i,1:min(idx)-1) = Sx(i,min(idx));
                    Sx(i,max(idx)+1:end) = Sx(i,max(idx));
                    % Extending y slopes top to bottom
                    idy = find(Sy(:,i));
                    Sy(1:min(idy)-1,i) = Sy(min(idy),i);
                    Sy(max(idy)+1:end,i) = Sy(max(idy),i);
                end
            else
                fprintf('No extension method chosen...\n')
            end

            % EDIT 15/10/15
            % Testing if we should use a window for cyclic conditions,
            % rather than the sum method.  Also enforce in x and y for both
            % slopes.  So use square window with 1 in pupil area.  May need
            % large extension for this???
%             dEdge = 2*obj.extension;
%             fracCentre = 1 - dEdge/obj.N;
%             % Square Tukey window
%             w = mywindow('squareTukey',obj.N,fracCentre);
%             Sx = w(:,:,2).*Sx;
%             Sy = w(:,:,1).*Sy;
            
            % Enforcing cyclic condition on slopes
            Sx(:,end) = -sum(Sx(:,1:end-1)');
            Sy(end,:) = -sum(Sy(1:end-1,:));
        end
        
        %% Gerchberg extension method
        function [Sx,Sy] = GerchbergMethod(Sx,Sy,obj)
            
            % Number of iterations
            nIter = obj.Gerchberg.nIter;
            N = obj.N;
            nL = obj.params.nLenslet;
            ext = obj.extension;

            % Defining size of slopes measurement data
            mask = obj.validLenslet;
            [mx,my] = size(mask);

            % Inverse mask (for use later)
            imask = 1-mask;

            % Valid slopes measurements
            idx = find(mask);

            % Prepare slopes: subtract average
            mSx = mean(Sx(idx));
            Sx = mask.*(Sx-mSx);
            mSy = mean(Sy(idx));
            Sy = mask.*(Sy-mSy);

            % Estimated slopes initially set to original slopes
            Sx_est = Sx;
            Sy_est = Sy;

            % Calculate filters
            RX = obj.RX;
            RY = obj.RY;
            GX = conj(obj.GX);
            GY = conj(obj.GY);
            absG2 = abs(GX).^2+abs(GY).^2;
            absG2(absG2<1e-3)=1e-3;
            
            wxy = mywindow('squareTukey',N,N-2*ext/N);
            w = wxy(:,:,1).*wxy(:,:,2);
            S2 = sqrt(sum(sum(w.^2)));
            
            for n=1:nIter

                % Add estimated slopes to measured ones
                % I.e. area within valid pupil is measured slopes, outside this
                % area we have the estimated slopes
                Sx_est = Sx.*mask + Sx_est.*imask;
                Sy_est = Sy.*mask + Sy_est.*imask;
                
                % FFT of new slopes
                fft_Sx = fftshift(fft2(w.*Sx_est))*N/S2;
                fft_Sy = fftshift(fft2(w.*Sy_est))*N/S2;
        
                % Estimate FFT slopes by filtering to produce FFT of phase (using
                % our reconstructors) and then use the reverse process (s = G*phi)
                % to get an estimate of the slopes (in Fourier domain)
                fft_Sx_est = (GX.*conj(GX).*fft_Sx + GX.*conj(GY).*fft_Sy)./absG2;
                fft_Sy_est = (GY.*conj(GX).*fft_Sx + GY.*conj(GY).*fft_Sy)./absG2;

                % Now get estimate of slopes using inverse FFT
                Sx_est = (ifft2(ifftshift(fft_Sx_est)));
                Sy_est = (ifft2(ifftshift(fft_Sy_est)));

            end

            % Final slopes are estimated slopes with mean added back on
            Sx = real(Sx.*mask + Sx_est.*imask + mSx);
            Sy = real(Sy.*mask + Sy_est.*imask + mSy);

        end
        
        
        
        
                   
            
    end

    
    methods (Static)
        
        %% Noise PSD
        function [Wnoise] = createNoisePSD(k,noiseVar)
            
            % Creates noise PSD from noise variance, assuming white noise
            % (constant noise over all frequencies).
            
            N = length(k);
            Wnoise = noiseVar/N^2 * ones(N,N);
            
        end
        
        %% Phase PSD
        function [Wphi,Sx,Sy] = createPhasePSD(k,m,d,nPxL,r0,L0)
            
            % creates PSD of phase over given frequencies.  Option to
            % include aliasing
            % k:    frequency vector
            % m:    factor dictating maximum frequency to be aliased (i.e.
            %       how many times spectrum is folded).
            % d:    sampling period, max(k) = 1/(2*d).
            % nPxL: pixels per lenslet
            % r0,L0:atmospheric parameters
            %
            % N.B. no piston removal, can be added if needed.
            
            [kx,ky] = meshgrid(k,k);
            
            Wphi = 0;
            Sx = 0;
            Sy = 0;
            % Averaging d for sinc function
            d_av = d - d/nPxL;
            
            for mi = -m:m
                for ni = -m:m
                   
                    % Frequency range
                    km = kx - mi/d;
                    kn = ky - ni/d;
                    
                    % Phase PSD 
                    % For current range of frequencies
                    Wmn = (abs(km.^2 + kn.^2) + (1/L0)^2).^(-11/6);
                    % Adding to total PSD
                    Wphi = Wphi + Wmn;
                    
                    % Slopes PSD
                    % Averaging term
                    Avx = abs(sinc(d_av*km).*sinc(d*kn).*exp(1i*pi*d*(km+kn))).^2;
                    Avy = abs(sinc(d*km).*sinc(d_av*kn).*exp(1i*pi*d*(km+kn))).^2;
                    Sx = Sx + Wmn.*abs(1i*2*pi*d*km).^2.*Avx;
                    Sy = Sy + Wmn.*abs(1i*2*pi*d*kn).^2.*Avy;
                    
                end
            end
            
            % Normalisation factors (~0.0229*r0^(-5/3)
            Wphi = (24.*gamma(6./5)./5).^(5./6).*...
                (gamma(11./6).^2./(2.*pi.^(11./3)))*r0^(-5/3)*Wphi;
            Sx = (24.*gamma(6./5)./5).^(5./6).*...
                (gamma(11./6).^2./(2.*pi.^(11./3)))*r0^(-5/3)*Sx;
            Sy = (24.*gamma(6./5)./5).^(5./6).*...
                (gamma(11./6).^2./(2.*pi.^(11./3)))*r0^(-5/3)*Sy;
            
                 
        end
        
        
        %% Remove modes
        function [Phi] = removeModes(Phi,N,pupil,modes,method)
            
            % methods = 'direct'  (in spatial domain)
            %           'fourier' (spatial frequency domain)
            % N.B. if the 'fourier' method is used then the phase (Phi)
            % must be passed to the function in fourier space (and will
            % also be returned in fourier space).  This phase should be
            % fftshifted to have the 0 frequency component in the centre of
            % the matrix/vector
            %
            % N: number of measurements
            % pupil: logical matrix describing valid elements (i.e. inside
            % the pupil)
            %
            % N.B: only compatible for waffle when Phi has an even number
            % of points.  I.e. we remove the waffle frequency at (1,1),
            % which must be an eigenmode of the measurement space. a vector
            % of 1,-1 across a grid is only an eigenmode over an even
            % number of points
            
            if strcmp(method,'direct')
                % Global waffle
                if modes(1)==1
                    % Waffle mode
                    waffle = ones(N,N);
                    for i=1:N
                        if mod(i,2) == 0
                            waffle(1:2:N,i) = -1;
                        else
                            waffle(2:2:N,i) = -1;
                        end
                    end
                    
                    % Waffle coefficient in Phi (calculate using inner
                    % product)
                    cw = sum(Phi(pupil).*waffle(pupil)) / sum(waffle(pupil).*waffle(pupil));
                    Phi = Phi - cw*waffle;
                end
                                
                % Local waffle removal
                if modes(2)==1
                    lwaffle = [-1 1; 1 -1];
                    tmp_phi = Phi;
                    for n=1:N/2
                        for m=1:N/2
                            v = zeros(N,N);
                            v(2*n-1:2*n,2*m-1:2*m) = lwaffle;
                            c(n,m) = sum(sum(tmp_phi.*v))./sum(sum(v.*v));
                            tmp_phi = tmp_phi - c(n,m)*v;
                        end
                    end
                    %Phi_waf = conv2(Phi,[-1 1; 1 -1],'same');
                    Phi = tmp_phi;
                end
                                
                % Piston removal
                if modes(3)==1
                    Phi(pupil) = Phi(pupil) - mean(Phi(pupil));
                end
                
            elseif strcmp(method,'fourier')
                % Both the global and local waffle modes can be filtered in
                % the Fourier space (no need to go to direct space); 
                % N.B.: As has been outlined in Pyoneer02 (OSA paper) the 
                % removal in Fourier space is not as much efective as in 
                % direct space due to artifacts originated by the extension 
                % method.
                
                % Global waffle
                if modes(1)==1
                    
                    % Waffle frequency (assuming an even number of points
                    % in the grid)
                    Phi(1,1) = 0;
                    %Phi(1,N/2+1) = 0;
                    %Phi(N/2+1,1) = 0;
                end
                
                % Local waffle
                if modes(2)==1
                    % the Fourier transform of the direct-space filter is used (See Poyneer
                    % 2002 (SPIE conf. paper) for the original filter)
                    % Local waffle mode in direct space
%                     Wmode = zeros(N,N);
%                     Wmode(N/2:N/2+1,N/2:N/2+1) = [1 -1; -1 1];
%                     
%                     % Apply filter (created using fft) to phase
%                     % 4 from normalisation: sum(sum(wmode.*wmode))
%                     Phi = Phi.*fft2(ifftshift(Wmode))/4;
                    [k,l] = meshgrid(-N/2:N/2-1,-N/2:N/2-1);
                    L = 1/(4)*(3+exp(-1i*2*pi*k*(N-1)/N)+exp(-1i*2*pi*l*(N-1)/N)...
                        -exp(1i*2*pi*(k+l)*(N-1)/N));
                    Phi = L.*Phi;
                end
                
                % Global piston
                if modes(3)==1
                    Phi(floor(N/2+1),floor(N/2+1)) = 0;
                end
            end
        end
        
        %% Pad data
        function mat = padarray(mat,vv)
            % Function to pad array
            if length(vv) == 1
                x = vv;
                y = vv;
            elseif length(vv) == 2
                x = vv(1);
                y = vv(2);
            else
                fprintf('Bad input')
            end

            % Padding with zeros
            mat = [zeros(size(mat,1),x) mat zeros(size(mat,1),x)];
            mat = [zeros(y,size(mat,2)); mat; zeros(y,size(mat,2))];
        end

        %% DM filters in fourier domain
        function out = dmFilter(dmFunction,ext,nLenslet)

            % DM filter
            aux = (-nLenslet/2:-max(dmFunction(:,1)))';
            aux_x = sort([ aux; flipud(-dmFunction(2:end,1)); dmFunction(:,1); -aux]);
            aux_y = [ 0*aux; flipud(dmFunction(2:end,2)); dmFunction(:,2); 0*aux];
            poly = spline( aux_x, aux_y);
            sub_h = ppval(poly,round(-nLenslet/2):nLenslet/2);
            h = sub_h'*sub_h;
            hh = h;
            
            % Pad filter, always to an even number of points
            if nLenslet/2==round(nLenslet/2)
                hh = [zeros(1,nLenslet+2); zeros(nLenslet+1,1) hh];
            end
            
            hh = fourierReconstructor.padarray(hh,ext);
            out = 1./(fft2(ifftshift(hh))); % DM filter (1/fft2(DM influence func))


        end

    end
    
end
        
function p = pattern(length,arg)
    p = ones(1,length);
    for pi = 1:2:length
        p(pi) = -1;
    end
    if nargin > 1
        p = p(end:-1:1);
    end
end

function [xslopes,yslopes,nn] = solve_conflicts(i,point,N,mask,xslopes,yslopes)
    localmask = mask(max(i-1,1):min(i+1,N),max(point-1,1):min(point+1,N));
    nneighbours = nnz(localmask);
    if nneighbours > 2 & nnz(sum(localmask,2) > 0) > 1 & nnz(sum(localmask,1) > 0) > 1% the conditions stands for: more than two neighbours,
        %with at least one at the left/right and one above/below
        nn = 1;
        if i < N/2+1 && point < N/2+1 %1st quadrant
            Q = 1;
        elseif i < N/2+1 && point > N/2+1 %2nd quadrant
            Q = 2;
        elseif i > N/2+1 && point < N/2+1%3rd quadrant
            Q = 3;
        elseif i > N/2+1 && point > N/2+1%4th quadrant
            Q = 4;
        end
        if Q == 1
            xslopes(i,point) = -xslopes(i,point+1) + xslopes(i+1,point) + xslopes(i+1,point+1);
            yslopes(i,point) = -yslopes(i+1,point) + yslopes(i,point+1) + yslopes(i+1,point+1);
        elseif Q == 2
            xslopes(i,point) = -xslopes(i,point-1) + xslopes(i+1,point) + xslopes(i+1,point-1);
            yslopes(i,point) = -yslopes(i+1,point) + yslopes(i+1,point-1) + yslopes(i,point-1);
        elseif Q == 3
            xslopes(i,point) = xslopes(i-1,point) + xslopes(i-1,point+1) - xslopes(i,point+1);
            yslopes(i,point) = -yslopes(i-1,point) + yslopes(i-1,point+1) + yslopes(i,point+1);
        elseif Q == 4
            xslopes(i,point) = -xslopes(i,point-1) + xslopes(i-1,point-1) + xslopes(i-1,point);
            yslopes(i,point) = -yslopes(i-1,point) + yslopes(i-1,point-1) + yslopes(i,point-1);
        end
    else
        nn = 0;
    end
end






