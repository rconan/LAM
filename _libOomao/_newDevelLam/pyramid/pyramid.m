classdef pyramid < handle
    %pyramid
    %A class that implements a Pyramid WFS
    %%
    properties 
        nLenslet                % # of lenslet
        resolution;             % resolution of the single quadrant detector images
        camera;                 %the detector used in the end of the chain
        %satisfy the Shannon's condition). Does not seem to work if modified
        %from its default value.
        %TODO fix c
        modulation;             %amplitude of the modulation in \lambda/D units, a decimal number
        referenceSlopes;        %the slopes vector of reference
        referenceSlopesMap;     %the slopes map of reference
        % slopes display handle
        slopesDisplayHandle;
        % slope listener
        slopesListener;
        tag = 'PYRAMID';
    end
    properties(GetAccess = 'public', SetAccess = 'private')
        validActuator
        lightMap;               %map of the light passed through the pyramid
        slopesMap;              %the slopes as read by processing the sensor
        wave;                   %incoming wave, a square complex matrix
        validSlopes;            % logical insindex indicating the validSlopes
        validIntensityPupil;    % logical index indicating the validIntensityPupil
        pyrMask;                % pyramid face transmitance and phase
        slopesUnits=1;          % normalisation factor to have the 1-to-1 gain btw input and output
    end
    
    properties(SetAccess=private,SetObservable=true)
        slopes;                 %the slopes as read by processing the sensor
    end
    
    properties (Dependent)
        alpha;                  %angle of incoming rays
        c;            %Multiplicative term of the Nyquist sampling (which is 2 by default). With  c = 2 (default) the quadrants are spread by a factor 2. End result :: c=2 -> 4pixels across \lambda/D
        binning;                % binning factor, default 1
    end
    
    properties (Dependent,SetAccess=private)
        nSlope;                 % # of slopes
    end
    
    properties (Access=private)
        p_alpha;
        p_c;
        p_binning = 1;          % binning factor, default 1
        isInitialized;
        u;
        log;
    end
    
    %%
    methods
        %% Constructor
        function pwfs = pyramid(nLenslet,resolution,varargin)
            p = inputParser;
            addRequired(p,'nLenslet',@isnumeric);
            addRequired(p,'resolution',@isnumeric);
            addParameter(p,'modulation',0,@isnumeric);
            addParameter(p,'binning',1,@isnumeric);
            parse(p,nLenslet,resolution,varargin{:})
            
            pwfs.nLenslet   = p.Results.nLenslet;
            pwfs.resolution = p.Results.resolution;
            pwfs.modulation = p.Results.modulation;
            
            pwfs.p_alpha=pi/2;
            pwfs.p_c=2;
            pwfs.camera=detector(2*pwfs.c*nLenslet);

            pwfs.binning    = p.Results.binning;
            pwfs.validSlopes = [pwfs.validIntensityPupil pwfs.validIntensityPupil];
            
            pwfs.u = 1+pwfs.resolution*(2*pwfs.c-1)/2:pwfs.resolution*(2*pwfs.c+1)/2;
            pwfs.camera.frameGrabber=pwfs; %
            makePyrMask(pwfs)
            
            pwfs.slopesListener = addlistener(pwfs,'slopes','PostSet',...
                @(src,evnt) pwfs.slopesDisplay );
            pwfs.slopesListener.Enabled = false;            
            pwfs.log = logBook.checkIn(pwfs);
        end
        
        %% Destructor
        function delete(pwfs)
            if isvalid(pwfs.camera)
                delete(pwfs.camera)
            end
            checkOut(pwfs.log,pwfs);
        end
        
        % get nSlope
        function out = get.nSlope(pwfs)
            out = sum(pwfs.validSlopes(:));
        end
        
        % setget alpha
        function out = get.alpha(pwfs)
            out = pwfs.p_alpha;
        end
        function set.alpha(pwfs,val)
            pwfs.p_alpha = val;
            makePyrMask(pwfs);
        end
        
        % setget c
        function out = get.c(pwfs)
            out = pwfs.p_c;
        end
        function set.c(pwfs,val)
            pwfs.p_c = val;
            pwfs.u = 1+pwfs.resolution*(2*pwfs.c-1)/2:pwfs.resolution*(2*pwfs.c+1)/2;

            makePyrMask(pwfs);
            
            pwfs.validIntensityPupil = utilities.piston(pwfs.nLenslet/pwfs.binning,...
                pwfs.nLenslet*pwfs.c/pwfs.binning,'type','logical');
            pwfs.validSlopes = [pwfs.validIntensityPupil pwfs.validIntensityPupil];
            
        end
        
        % setget binning
        function out = get.binning(pwfs)
            out = pwfs.p_binning;
        end
        function set.binning(pwfs,val)
            pwfs.p_binning = val;
            pwfs.isInitialized = false;
            
            pwfs.validIntensityPupil = utilities.piston(pwfs.nLenslet/pwfs.p_binning,...
                pwfs.nLenslet*pwfs.c/pwfs.p_binning,'type','logical');
            pwfs.validSlopes = [pwfs.validIntensityPupil pwfs.validIntensityPupil];
            
        end
        
        %% Get valid actuators
        function val = get.validActuator(obj)
            nElements            = 2*obj.nLenslet+1; % Linear number of lenslet+actuator
            validLensletActuator = zeros(nElements);
            index                = 2:2:nElements; % Lenslet index
            w = floor(obj.nLenslet/2)+1:floor(obj.nLenslet/2)+obj.nLenslet;
            validLensletActuator(index,index) = obj.validIntensityPupil(w,w);
            for xLenslet = index
                for yLenslet = index
                    if validLensletActuator(xLenslet,yLenslet)==1
                        xActuatorIndice = [xLenslet-1,xLenslet-1,...
                            xLenslet+1,xLenslet+1];
                        yActuatorIndice = [yLenslet-1,yLenslet+1,...
                            yLenslet+1,yLenslet-1];
                        validLensletActuator(xActuatorIndice,yActuatorIndice) = 1;
                    end
                end
            end
            index = 1:2:nElements; % Actuator index
            val   = logical(validLensletActuator(index,index));
        end
        
        %% Relay
        % Method that allows compatibility with the overloaded mtimes
        % operator, allowing things like source=(source.*tel)*wfs
        function relay(pwfs, src)
            pyramidTransform(pwfs,src.catWave);
            grab(pwfs.camera,pwfs.lightMap);
            dataProcessing(pwfs);
        end
        
        %% slopesDisplay
        function slopesDisplay(pwfs)
            if ishandle(pwfs.slopesDisplayHandle)
                set(pwfs.slopesDisplayHandle,'Cdata',pwfs.slopesMap)
            else
                pwfs.slopesDisplayHandle = imagesc(pwfs.slopesMap);
                ax = gca;
                pos = get(ax,'position');
                axis xy equal tight
                ylabel(colorbar('location','EastOutside'),'Pixel')
                set(ax,'position',pos)
            end
        end
        
        %% Propagate the wavefront transformed by the pyram WFS to the detector
        function pwfs = pyramidTransform(pwfs,wave)
            [n1,n2,n3] = size(wave);
            nWave = n2*n3/n1;
            px_side  = n1*2*pwfs.c;
            q = zeros(px_side,px_side,nWave);
            u = 1+n1*(2*pwfs.c-1)/2:n1*(2*pwfs.c+1)/2;
            q(u,u,:) = reshape(wave,n1,n1,nWave);
            
                        
            %             figure()
            %             h = imagesc(I4Q);
            %             axis square
            %             colorbar
            %             % drawnow
            
            if pwfs.modulation>0
                [u,v] = ndgrid((0:(px_side-1))./px_side);
                [o,r] = cart2pol(u,v);
                nTheta = round(4*pi*pwfs.modulation);

                I4Q = zeros(px_side,px_side,nWave,nTheta); %
                for kTheta = 1:nTheta
                    theta = (kTheta-1)*2*pi/nTheta;
                    fftPhasor = exp(-1i.*pi.*8*pwfs.modulation*r.*cos(o+theta));
                    buf = bsxfun(@times,q,fftPhasor);
                    buf = bsxfun(@times,fft2(buf),pwfs.pyrMask);
                    I4Q(:,:,:,kTheta) = abs(fft2(buf)).^2;
                    
                    %I4Q(:,:,kTheta) = fft2(fft2(q.*fftPhasor).*fpym);
                    %imagesc(fftshift(abs(fft2(q.*fftPhasor))))
                    %drawnow
                    %pause
                    %         set(h,'Cdata',I4Q(:,:,kTheta))
                    %         drawnow
                end
                I4Q = sum(I4Q,4);
                
                %fftPhasor = 2./(pwfs.modulation*r).*(cos(pi*pwfs.modulation*r) + sin(pi*pwfs.modulation*r));
                %fftPhasor(isinf(fftPhasor)) = pwfs.modulation;
                %I4Q =  abs(fft2(fft2(q).*fftPhasor.*fpym)).^2;
                %                set(h,'Cdata',I4Q)
                %                drawnow
            else
                I4Q = bsxfun(@times,fft2(q),pwfs.pyrMask);
                I4Q = abs(fft2(I4Q)).^2;
            end
            pwfs.lightMap = reshape(I4Q,px_side,px_side*nWave);
        end
        
        %% make the pyramid phase mask
        function makePyrMask(pwfs)
            px_side  = pwfs.resolution*2*pwfs.c;
            
            
            [fx,fy] = freqspace(px_side,'meshgrid');
            fx = fx.*floor(px_side/2);
            fy = fy.*floor(px_side/2);
            %pym = zeros(px_side);
            
            % pyramid face transmitance and phase for fx>=0 & fy>=0
            mask  = heaviside(fx).*heaviside(fy);
            phase = -pwfs.alpha.*(fx+fy);
            pym   = mask.*exp(1i.*phase);
            % pyramid face transmitance and phase for fx>=0 & fy<=0
            mask  = heaviside(fx).*heaviside(-fy);
            phase = -pwfs.alpha.*(fx-fy);
            pym   = pym + mask.*exp(1i.*phase);
            % pyramid face transmitance and phase for fx<=0 & fy<=0
            mask  = heaviside(-fx).*heaviside(-fy);
            phase = pwfs.alpha.*(fx+fy);
            pym   = pym + mask.*exp(1i.*phase);
            % pyramid face transmitance and phase for fx<=0 & fy>=0
            mask  = heaviside(-fx).*heaviside(fy);
            phase = -pwfs.alpha.*(-fx+fy);
            pym   = pym + mask.*exp(1i.*phase);
            pwfs.pyrMask   = fftshift(pym./sum(abs(pym(:))));
        end
        
        %%
        function INIT(pwfs)
            
            %makePyrMask(pwfs)
            %pyramidTransform(pwfs)
            %dataProcessing(pwfs)
            
            pwfs.referenceSlopesMap = pwfs.slopesMap;
            pwfs.referenceSlopes = pwfs.slopes;
            
            gainCalibration(pwfs)
            
            pwfs.isInitialized = true;
            
        end
        
        %%
        function uplus(pwfs)
            %This function updates the detector with whatever light there
            %is in the pyramid, and then procedes to process the output of
            %the detector in order to obtain the slopes
            %pwfs.camera.frameGrabber=pwfs;
            %grab(pwfs.camera,pwfs.lightMap);
            %dataProcessing(pwfs);
            
            grab(pwfs.camera,pwfs.lightMap);
            dataProcessing(pwfs);
            
        end
        
        %% Grab a frame and crop signal areas
        function dataProcessing(pwfs)
            %This function computes the slopes corresponding to a flat
            %incoming light wave.
            
            px_side  = size(pwfs.camera.frame,1);
            
            I4Q=utilities.binning(pwfs.camera.frame,size(pwfs.camera.frame)/pwfs.binning);             % binning

            n = px_side/pwfs.binning;
            I4Q = reshape(I4Q,n,n,[]);
            
            half = floor(px_side/2/pwfs.binning);
            xc = half+1;
            
            is = xc-half+1;
            ie = xc;
            js = xc-half+1;
            je = xc;
            I1 = I4Q(is:ie,js:je,:);
            is = xc;
            ie = xc+half-1;
            I2 = I4Q(is:ie,js:je,:);
            is = xc;
            ie = xc+half-1;
            js = xc;
            je = xc+half-1;
            I3 = I4Q(is:ie,js:je,:);
            is = xc-half+1;
            ie = xc;
            I4 = I4Q(is:ie,js:je,:);
            
            computeSlopes(pwfs, I1, I2, I3, I4);% modiffff 
        end
        
        %% Compute slopes from 4 intensity maps
        function computeSlopes(pwfs, I1, I2, I3, I4)
            
            nPhase = size(I1,3);
            % normalisation options
            %    1) normalisation pixel-wise by the intensity
            I = (I1+I2+I3+I4);      %
            %    2) normalisation by the integrated flux
            I2D = utilities.toggleFrame(I,2);
            I = sum(I2D(pwfs.validIntensityPupil))*ones(size(I));
            
            SyMap = (I1-I2+I4-I3)./I;
            SxMap = (I1-I4+I2-I3)./I;

            if pwfs.isInitialized == true
                pwfs.slopesMap = bsxfun(@minus,[SxMap,SyMap],pwfs.referenceSlopesMap);
            else
                pwfs.slopesMap=[SxMap,SyMap];
            end
            [n1,n2,n3] = size(pwfs.slopesMap);
            slopesMap_ = reshape(pwfs.slopesMap,n1*n2,n3);
            pwfs.slopes = slopesMap_(pwfs.validSlopes(:),:)*pwfs.slopesUnits;
            pwfs.slopesMap = reshape( slopesMap_ , n1, n2*n3);
        end
        
        %% gain calibration
        function gainCalibration(pwfs)
            nPx = pwfs.resolution;
            ngs = source('wavelength',photometry.R);
            tel = telescope(1,'resolution',nPx);
            zer = zernike(3,'resolution', nPx);
            for i = 1:5
                zer.c = (i-3)*0.1/(ngs.waveNumber);
                ngs = ngs.*tel*zer*pwfs;
                sx(i) = mean(pwfs.slopes(1:end/2));
                sy(i) = mean(pwfs.slopes(end/2+1:end));
            end
            Ox_in = 4*[-2:2]*0.1;%/ngs.waveNumber;
            Ox_out = sy;
            slopesLinCoef = polyfit(Ox_in,Ox_out,1);
            pwfs.slopesUnits = 1/slopesLinCoef(1);
            
            ngs = ngs.*tel*pwfs;
        end
        %% HIlbert transform approximation for the slopes
        function hilbertSlopes(pwfs,telescope,source)
            P=telescope.pupil;
            phi=angle(source.wave);
            %            Sx = ???PHx (P??) + P??Hx (P) ??? Hxy (P)Hy (P??) + Hxy (P??)Hy (P),
            %            Sy = ???PHy (P??) + P??Hy (P) ??? Hxy (P)Hx (P??) + Hxy (P??)Hx (P).
            %            Sx = -P*imag(hilbert(P*phi)) + P*phi*imag(hilbert(P)) - imag(hilbert(imag(hilbert(transpose(P)))))*imag(hilbert(transpose(P*phi))) + imag(hilbert(imag(hilbert(transpose(P*phi)))))*imag(hilbert(transpose(P))) ;
            %            Sy = -P*imag(hilbert(transpose(P*phi))) + P*phi*imag(hilbert(transpose(P))) - imag(hilbert(imag(hilbert(transpose(P)))))*imag(hilbert(P*phi)) + imag(hilbert(imag(hilbert(transpose(P*phi)))))*imag(hilbert(P)) ;
            %Sx = -P.*imag(hilbert(P.*phi)) + P.*phi.*imag(hilbert(P)) - imag(hilbert(transpose(imag(hilbert(P))))).*imag(hilbert(transpose(P.*phi))) + imag(hilbert(transpose(imag(hilbert(P.*phi))))).*imag(hilbert(transpose(P))) ;
            %Sy = -P*imag(hilbert(transpose(P*phi))) + P*phi*imag(hilbert(transpose(P))) - imag(hilbert(imag(hilbert(transpose(P)))))*imag(hilbert(P*phi)) + imag(hilbert(imag(hilbert(transpose(P*phi)))))*imag(hilbert(P)) ;
            myHilbertX = @(x) imag(hilbert(x));
            myHilbertY = @(x) imag(hilbert(x')');
            myHilbertXY = @(x) imag(hilbert(imag(hilbert(x')')));
            Sx = -P.*myHilbertX(P.*phi) + P.*phi.*myHilbertX(P) - myHilbertXY(P).*myHilbertY(P.*phi) + myHilbertXY(P.*phi).*myHilbertY(P);
            Sy = -P.*myHilbertY(P.*phi) + P.*phi.*myHilbertY(P) - myHilbertXY(P).*myHilbertX(P.*phi) + myHilbertXY(P.*phi).*myHilbertX(P) ;
            pwfs.slopes=[Sx,Sy];
            
        end
        
    end
end




















