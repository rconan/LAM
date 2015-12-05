classdef pyramid < handle
    %pyramid
    %A class that implements a Pyramid WFS
    %%
    properties
        resolution;             % resolution of the single quadrant detector images
        camera;                 %the detector used in the end of the chain
        %satisfy the Shannon's condition). Does not seem to work if modified
        %from its default value.
        %TODO fix c
        modulation;             %amplitude of the modulation in \lambda/D units, a decimal number
        referenceSlopes;        %the slopes vector of reference
        referenceSlopesMap;     %the slopes map of reference
    end
    properties(GetAccess = 'public', SetAccess = 'private')
        lightMap;               %map of the light passed through the pyramid
        slopes;                 %the slopes as read by processing the sensor
        slopesMap;              %the slopes as read by processing the sensor
        wave;                   %incoming wave, a square complex matrix
        validSlopes;            % logical index indicating the validSlopes
        validIntensityPupil;    % logical index indicating the validIntensityPupil
        pyrMask;                % pyramid face transmitance and phase
        slopesUnits=1;          % normalisation factor to have the 1-to-1 gain btw input and output
    end
    
    properties (Dependent)
        alpha;                  %angle of incoming rays
        multNyqSamp;            %Multiplicative term of the Nyquist sampling (which is 2 by default). With  multNyqSamp = 2 (default) the quadrants are spread by a factor 2. End result :: multNyqSamp=2 -> 4pixels across \lambda/D
        binning;                % binning factor, default 1
    end
    
    properties (Access=private)
        p_alpha;
        p_multNyqSamp;
        p_binning = 1;          % binning factor, default 1
        isInitialized;
        u;
    end
    
    %%
    methods
        %% Constructor
        function pwfs = pyramid(resolution)
            pwfs.resolution = resolution;
            pwfs.modulation = 0;
            pwfs.wave=ones(floor(resolution));
            pwfs.p_alpha=pi/2;
            pwfs.p_multNyqSamp=2;
            pwfs.camera=detector(2*pwfs.multNyqSamp*floor(resolution));
            pwfs.isInitialized=false;
            pwfs.validIntensityPupil = utilities.piston(resolution,resolution*pwfs.multNyqSamp,'type','logical');
            pwfs.validSlopes = [pwfs.validIntensityPupil pwfs.validIntensityPupil];
            
            pwfs.u = 1+pwfs.resolution*(2*pwfs.multNyqSamp-1)/2:pwfs.resolution*(2*pwfs.multNyqSamp+1)/2;
            pwfs.camera.frameGrabber=pwfs; %
            makePyrMask(pwfs)
        end
        
        %% Destructor
        function delete(pwfs)
            %             if isvalid(obj.slopesListener)
            %                 delete(obj.slopesListener)
            %             end
            %             if isvalid(obj.intensityListener)
            %                 delete(obj.intensityListener)
            %             end
            %             if isvalid(obj.paceMaker)
            %                 if strcmp(obj.paceMaker.Running,'on')
            %                     stop(obj.paceMaker)
            %                 end
            %                 delete(obj.paceMaker)
            %             end
            %             if ishandle(pwfs.slopesDisplayHandle)
            %                 delete(pwfs.slopesDisplayHandle)
            %             end
            %             if ishandle(pwfs.intensityDisplayHandle)
            %                 delete(pwfs.intensityDisplayHandle)
            %             end
            
            if isvalid(pwfs.camera)
                delete(pwfs.camera)
            end
            if ~isempty(pwfs.log)
                checkOut(pwfs.log,pwfs);
            end
        end
        %% Gets and sets
        % The get functions are not necessary here, because all properties
        % are defined as public with respect to GetAccess
        
        % setget alpha
        function out = get.alpha(pwfs)
            out = pwfs.p_alpha;
        end
        function set.alpha(pwfs,val)
            pwfs.p_alpha = val;
            makePyrMask(pwfs);
        end
        
        % setget multNyqSamp
        function out = get.multNyqSamp(pwfs)
            out = pwfs.p_multNyqSamp;
        end
        function set.multNyqSamp(pwfs,val)
            pwfs.p_multNyqSamp = val;
            pwfs.u = 1+pwfs.resolution*(2*pwfs.multNyqSamp-1)/2:pwfs.resolution*(2*pwfs.multNyqSamp+1)/2;

            makePyrMask(pwfs);
            
            pwfs.validIntensityPupil = utilities.piston(pwfs.resolution/pwfs.binning,pwfs.resolution*pwfs.multNyqSamp/pwfs.binning,'type','logical');
            pwfs.validSlopes = [pwfs.validIntensityPupil pwfs.validIntensityPupil];
            
        end
        
        % setget binning
        function out = get.binning(pwfs)
            out = pwfs.p_binning;
        end
        function set.binning(pwfs,val)
            pwfs.p_binning = val;
            pwfs.isInitialized = false;
            
            pwfs.validIntensityPupil = utilities.piston(pwfs.resolution/pwfs.binning,pwfs.resolution*pwfs.multNyqSamp/pwfs.binning,'type','logical');
            pwfs.validSlopes = [pwfs.validIntensityPupil pwfs.validIntensityPupil];
            
        end
        
        
        %         %%Get and set modulation
        %         function modulation = getmodulation(pwfs)
        %             modulation = pwfs.modulation;
        %         end
        %         function pwfs = setmodulation(pwfs,modulation)
        %             pwfs.modulation = modulation;
        %         end
        
        %         %%Get and set wave
        %         function wave = getwave(pwfs)
        %             wave = pwfs.wave;
        %         end
        %         function pwfs = setwave(pwfs,wave)
        %             if mod(length(wave),2)==0
        %                 pwfs.wave = wave;
        %             else
        %                 display('WARNING The wave must be a square matrix of even dimension. The new wave setting has been disregarded')
        %             end
        %         end
        
        
        %         %%Get and set binning
        %         function pwfs = setBinning(pwfs, binning)
        %             pwfs.binning = binning;
        %         end
        %         function binning = getBinning(pwfs)
        %             binning = pwfs.binning;
        %         end
        
        %% Relay
        % Method that allows compatibility with the overloaded mtimes
        % operator, allowing things like source=(source.*tel)*wfs
        function relay(pwfs, src)
            pwfs.wave=src.wave;
            pyramidTransform(pwfs);
            grab(pwfs.camera,pwfs.lightMap);
            dataProcessing(pwfs);
        end
        
        %% slopesDisplay
        function slopesDisplay(pwfs)
            imagesc(pwfs.slopesMap.*pwfs.validSlopes);
        end
        
        %% Propagate the wavefront transformed by the pyramid WFS to the detector
        function pwfs = pyramidTransform(pwfs)
            if ~pwfs.isInitialized
                makePyrMask(pwfs)
            end
            %n = length(pwfs.wave);
            px_side  = pwfs.resolution*2*pwfs.multNyqSamp;
            q = zeros(px_side);
            %u = 1+n*(2*pwfs.multNyqSamp-1)/2:n*(2*pwfs.multNyqSamp+1)/2;
            q(pwfs.u,pwfs.u) = pwfs.wave;
            
            
            I4Q = abs(fft2(fft2(q).*fftshift(pwfs.pyrMask))).^2;
            
            %             figure()
            %             h = imagesc(I4Q);
            %             axis square
            %             colorbar
            %             % drawnow
            
            if pwfs.modulation>0
                [u,v] = ndgrid((0:(px_side-1))./px_side);
                [o,r] = cart2pol(u,v);
                nTheta = round(4*pi*pwfs.modulation);
                I4Q = zeros(px_side,px_side,nTheta); %
                fpym = fftshift(pwfs.pyrMask);
                for kTheta = 1:nTheta
                    theta = (kTheta-1)*2*pi/nTheta;
                    fftPhasor = exp(-1i.*pi.*8*pwfs.modulation*r.*cos(o+theta));
                    I4Q(:,:,kTheta) = abs(fft2(fft2(q.*fftPhasor).*fpym)).^2;
                    
                    %I4Q(:,:,kTheta) = fft2(fft2(q.*fftPhasor).*fpym);
                    %imagesc(fftshift(abs(fft2(q.*fftPhasor))))
                    %drawnow
                    %pause
                    %         set(h,'Cdata',I4Q(:,:,kTheta))
                    %         drawnow
                end
                I4Q = sum(I4Q,3);
                
                %fftPhasor = 2./(pwfs.modulation*r).*(cos(pi*pwfs.modulation*r) + sin(pi*pwfs.modulation*r));
                %fftPhasor(isinf(fftPhasor)) = pwfs.modulation;
                %I4Q =  abs(fft2(fft2(q).*fftPhasor.*fpym)).^2;
                %                set(h,'Cdata',I4Q)
                %                drawnow
            end
            pwfs.lightMap=I4Q;
        end
        
        %% make the pyramid phase mask
        function makePyrMask(pwfs)
            px_side  = length(pwfs.wave)*2*pwfs.multNyqSamp;
            
            
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
            pwfs.pyrMask   = pym./sum(abs(pym(:)));
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
            
            
            %pwfs.wave=src.wave;
            pyramidTransform(pwfs);
            grab(pwfs.camera,pwfs.lightMap);
            dataProcessing(pwfs);
            
        end
        
        %% Grab a frame and crop signal areas
        function dataProcessing(pwfs)
            %This function computes the slopes corresponding to a flat
            %incoming light wave.
            
            px_side  = length(pwfs.camera.frame);
            
            I4Q=utilities.binning(pwfs.camera.frame,[px_side px_side]/pwfs.binning);             % binning
            
            half = floor(px_side/pwfs.multNyqSamp/pwfs.binning);
            xc = half+1;
            
            is = xc-half+1;
            ie = xc;
            js = xc-half+1;
            je = xc;
            I1 = I4Q(is:ie,js:je);
            is = xc;
            ie = xc+half-1;
            I2 = I4Q(is:ie,js:je);
            is = xc;
            ie = xc+half-1;
            js = xc;
            je = xc+half-1;
            I3 = I4Q(is:ie,js:je);
            is = xc-half+1;
            ie = xc;
            I4 = I4Q(is:ie,js:je);
            
            
            computeSlopes(pwfs, I1, I2, I3, I4);
        end
        
        %% Compute slopes from 4 intensity maps
        function computeSlopes(pwfs, I1, I2, I3, I4)
            
            % normalisation options
            %    1) normalisation pixxel-wise by the intensity
            I = (I1+I2+I3+I4);      %
            %    2) normalisation by the integrated flux
            I = sum(I(pwfs.validIntensityPupil))*ones(size(I));
            
            SyMap = (I1-I2+I4-I3)./I;
            SxMap = (I1-I4+I2-I3)./I;
            if pwfs.isInitialized == true
                pwfs.slopesMap=[SxMap,SyMap]-pwfs.referenceSlopesMap;
            else
                pwfs.slopesMap=[SxMap,SyMap];
            end
            pwfs.slopes = pwfs.slopesMap(pwfs.validSlopes)*pwfs.slopesUnits;
        end
        
        %% gain calibration
        function gainCalibration(pwfs)
            nPx = length(pwfs.wave);
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




















