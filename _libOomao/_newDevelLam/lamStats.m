classdef lamStats
    % Zernike statistics static class
    
    methods (Static)
        
        function out = multipleSpectra(nu,atm,zern,tel)
            
            nMode = zern.nMode;
            out = zeros(nMode, length(nu),atm.nLayer);
            parfor kmode=1:nMode
                out(kmode,:,:)  = lamStats.temporalSpectrum(nu,atm,tel,kmode);
            end
            
        end
        
        
        function out = temporalSpectrum(nu,atm,tel,kmode)
            %% SPECTRUM Phase power spectrum density
            %
            % out = phaseStats.spectrum(f,atm) computes the phase power
            % spectrum density from the spatial frequency f and an
            % atmosphere object
            %
            % See also atmosphere
            
            out = zeros(size(nu,2),atm.nLayer);
            for kLayer = 1:atm.nLayer
                nPx = size(atm.layer(kLayer).phase,1);
                pupil = logical(utilities.piston(nPx));
                D = tel.diameterAt(atm.layer(kLayer).altitude);
                
                zern = zernike(kmode+1,'resolution',nPx,'pupil',pupil,'D',D);
                
                atmSlab = slab(atm,kLayer);
                [vx,vy] = pol2cart(atmSlab.layer.windDirection,atmSlab.layer.windSpeed);
                for k=1:numel(nu)
                    if abs(vx)>eps(atmSlab.layer.windSpeed)
                        out(k,kLayer) = out(k,kLayer) + quadgk( @integrandFy , -Inf, Inf);
                    else
                        out(k,kLayer) = out(k,kLayer) + quadgk( @integrandFx , -Inf, Inf);
                    end
                end
                if vx == 0;
                    signvx = 1;
                else
                    signvx = sign(vx);
                end
                if vy == 0;
                    signvy = 1;
                else
                    signvy = sign(vy);
                end
                out(:,kLayer) = out(:,kLayer)*signvx*signvy;
                %%%%
                % --- apply wind direction using Conan95, Eq. 30
            end
            
            function int = integrandFy(fy)
                fx = (nu(k) -fy*vy)/vx;
                %fx = (nu(k))/vx;
                int = zernikeStats.spectrum( hypot(fx,fy) , atan2(fy,fx), atmSlab , zern)/vx;
            end
            
            function int = integrandFx(fx)
                fy = (nu(k) -fx*vx)/vy;
                int = zernikeStats.spectrum( hypot(fx,fy) , atan2(fy,fx), atmSlab , zern)/vy;
            end
        end
    end
end