classdef pyramid < handle
    %pyramid
    %A class that implements a Pyramid WFS
    %% 
    properties
        camera;%the detector used in the end of the chain
    end
    properties(GetAccess = 'public', SetAccess = 'private')
        alpha;%angle of incoming rays
        c;%The Nyquist sampling (which is 2 by default, in order to
          %satisfy the Shanon???s condition). Does not semm to work if modified
          %from its default value.
        %TODO fix c
        isInitialized;
        lightmap;%map of the light passed through the pyramid
        modulation;%amplitude of the modulation, a decimal number
        referenceSlopes;%the slopes of reference
        slopes;%the slopes as read by processing the sensor
        wave;%incoming wave, a square complex matrix
        validSlopes; % logical index indicating the validSlopes
    end
    %%
    methods
        %% Constructor
        function pwfs = pyramid(resolution)
            pwfs.modulation = 0;
            pwfs.wave=ones(floor(resolution));
            pwfs.alpha=pi/2;
            pwfs.c=2;
            pwfs.camera=detector(2*pwfs.c*floor(resolution));
            pwfs.isInitialized=false;
            
            validSlopesMap = logical(padarray(utilities.piston(resolution/2), [resolution/4 resolution/4], 'both'));
            pwfs.validSlopes = [validSlopesMap validSlopesMap];

        end
        
        
        %% Gets and sets
        % The get functions are not necessary here, because all properties
        % are defined as public with respect to GetAccess
        
        %%Get and set alpha
        function alpha = getalpha(pwfs)
            alpha = pwfs.alpha;
        end
        function pwfs = setalpha(pwfs,alpha)
            pwfs.alpha = alpha;
        end
        
        %%Get and set c
        function c = getc(pwfs)
            c = pwfs.c;
        end
        function pwfs = setc(pwfs,c)
            pwfs.c = c;
            'Warning, setting c to any value other than 2 will probably give wrong results'
        end
        
        %%Get and set modulation
        function modulation = getmodulation(pwfs)
            modulation = pwfs.modulation;
        end
        function pwfs = setmodulation(pwfs,modulation)
            pwfs.modulation = modulation;
        end
        
        %%Get and set wave
        function wave = getwave(pwfs)
            wave = pwfs.wave;
        end
        function pwfs = setwave(pwfs,wave)
            if mod(length(wave),2)==0
                pwfs.wave = wave;
            else
                display('WARNING The wave must be a square matrix of even dimension. The new wave setting has been disregarded')
            end
        end
        
        
        %% Relay
        % Method that allows compatibility with the overloaded mtimes 
        % operator, allowing things like source=(source.*tel)*wfs
        function relay(pwfs, src)
            pwfs.wave=src.wave;
            pwfs.pyramidTransform;
            pwfs.camera.frameGrabber=pwfs;
            grab(pwfs.camera,pwfs.lightmap);
        end
        
        
        %% Propagate the wavefront transformed by the pyramid WFS to the detector 
        function pwfs = pyramidTransform(pwfs) 
            n = length(pwfs.wave);
            px_side  = length(pwfs.wave)*2*pwfs.c;
            q = zeros(px_side);
            u = 1+n*(2*pwfs.c-1)/2:n*(2*pwfs.c+1)/2;
            q(u,u) = pwfs.wave;

            [fx,fy] = freqspace(px_side,'meshgrid');
            fx = fx.*floor(px_side/2);
            fy = fy.*floor(px_side/2);
            pym = zeros(px_side);

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
            pym   = pym./sum(abs(pym(:)));

            I4Q = abs(fft2(fft2(q).*fftshift(pym))).^2;

%             figure()
%             h = imagesc(I4Q);
%             axis square
%             colorbar
%             % drawnow

            if pwfs.modulation>0
                [u,v] = ndgrid((0:(px_side-1))./px_side);
                [o,r] = cart2pol(u,v);
                nTheta = round(4*pi*pwfs.modulation);
                I4Q = zeros(px_side,px_side,nTheta);
                fpym = fftshift(pym);
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
            pwfs.lightmap=I4Q;
        end
        
        function INIT(pwfs)
            %This function computes the slopes corresponding to a flat
            %incoming light wave.
            px_side  = length(pwfs.camera.frame);
            half = floor(px_side/2);
            xc = half+1;
            I4Q=pwfs.camera.frame;
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
            

            I = (I1+I2+I3+I4);
            I = sum(I(pwfs.validSlopes(:,1:end/2)))*ones(size(I));
            Sy = (I1-I2+I4-I3)./I;
            Sx = (I1-I4+I2-I3)./I;
            pwfs.referenceSlopes=[Sx,Sy];
            pwfs.isInitialized = true;
            
        end         
        
        function uplus(pwfs)
            %This function updates the detector with whatever light there
            %is in the pyramid, and then procedes to process the output of
            %the detector in order to obtain the slopes
            pwfs.camera.frameGrabber=pwfs;
            grab(pwfs.camera,pwfs.lightmap);
            dataProcessing(pwfs);
        end
        
        function dataProcessing(pwfs)
            px_side  = length(pwfs.camera.frame);
            half = floor(px_side/2);
            xc = half+1;
            I4Q=pwfs.camera.frame;
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
            
            I = (I1+I2+I3+I4); 
            I = sum(I(:))*ones(size(I));
            Sy = (I1-I2+I4-I3)./I;
            Sx = (I1-I4+I2-I3)./I;
            if pwfs.isInitialized == true
                pwfs.slopes=[Sx,Sy]-pwfs.referenceSlopes;
            else
                pwfs.slopes=[Sx,Sy];
            end
        end
        
        function slopesDisplay(pwfs)
            imagesc(pwfs.slopes);
        end
        
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




















