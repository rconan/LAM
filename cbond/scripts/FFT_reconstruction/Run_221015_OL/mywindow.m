function [w] = mywindow(wname,N,vargin)

    [X,Y] = meshgrid(1:N,1:N);
    
    n0 = (N+1)/2;
    
    R = sqrt((X-n0).^2 + (Y-n0).^2);
    
    w = zeros(N,N);
    [idx] = find(R<=N/2);
    
    if strcmp(wname,'barthann')
        w(idx) = 0.62 - 0.48*abs(R(idx)/N) + 0.38*cos(2*pi*R(idx)/N);
    elseif strcmp(wname,'bartlett')
        w(idx) = 1-2*R(idx)/N;
    elseif strcmp(wname,'bohman')
        rho = R/(N/2);
        w(idx) = (1-rho(idx)).*cos(pi*rho(idx)) + (1/pi)*sin(pi*rho(idx));
    elseif strcmp(wname,'flattop')
        a0 = 1-0.21557895;
        a1 = 0.41663158;
        a2 = 0.277263158;
        a3 = 0.083578947;
        a4 = 0.006947368;
        w(idx) = a0 + a1*cos(2*pi*R(idx)/N) - a2*cos(4*pi*R(idx)/N) ...
            + a3*cos(6*pi*R(idx)/N) - a4*cos(8*pi*R(idx)/N);
    elseif strcmp(wname,'hamming')
        w(idx) = 0.46 + 0.46*cos(2*pi*R(idx)/N);
    elseif strcmp(wname,'hann')
        w(idx) = 0.5 + 0.5*cos(2*pi*(R(idx))/(N-1));
    elseif strcmp(wname,'squareHann')
        w = 0.25*(1-cos(2*pi*(X-1)/(N-1))).*(1-cos(2*pi*(Y-1)/(N-1)));
    elseif strcmp(wname,'nuttall')
        a0 = 1-0.3635819;
        a1 = 0.4891775;
        a2 = 0.1365995;
        a3 = 0.0106411;
        w(idx) = a0 + a1*cos(2*pi*R(idx)/N) - a2*cos(4*pi*R(idx)/N) ...
            + a3*cos(6*pi*R(idx)/N);
    elseif strcmp(wname,'parzen')
        % ???? Is this the right equation?
        w(idx) = 1-(R(idx)/(N/2)).^3;
    elseif strcmp(wname,'tukey')
        r_tukey = vargin*N/2;
        w(idx) = 1/2 * (cos(2*pi/(2*(N/2-r_tukey))*(R(idx)-r_tukey))+1);
        i = find(R<=r_tukey);
        w(i) = 1;
    elseif strcmp(wname,'squareTukey')
        r_tukey = vargin*N/2;
        d = N/2-r_tukey;
        wx = w;
        wy = w;
        wx(abs(X-n0)<=r_tukey) = 1;
        wx(X<=d) = 0.5*(1-cos(pi*(X(X<=d)-1)/(d)));
        wx(X>=N-d) = 0.5*(cos(pi*(X(X>=N-d)-(N-d))/(d))+1);
        wy(abs(Y-n0)<=r_tukey) = 1;
        wy(Y<=d) = 0.5*(1-cos(pi*(Y(Y<=d)-1)/(d)));
        wy(Y>=N-d) = 0.5*(cos(pi*(Y(Y>=N-d)-(N-d))/(d))+1);
        %wy = 1/2 * (cos(pi*(X-1+d)/(d-1))+1);    %(cos(2*pi/(2*(N/2-r_tukey+1/2))*(X-n0-r_tukey))+1);
        %wy = 1/2 * (cos(pi*(Y-1+d)/(d-1))+1);    %(cos(2*pi/(2*(N/2-r_tukey+1/2))*(Y-n0-r_tukey))+1);
        %
        %wy(abs(Y-n0)<=r_tukey) = 1;
        w = zeros(N,N,2);
        w(:,:,1) = wx;
        w(:,:,2) = wy;
    end
    
    
    
    

end

