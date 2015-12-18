% Version 1.3
% 20150728: corrected bug with sign convetion of Zernike tilt for time series data

function PS = ESO_GeneratePS_v1_3(File,Parameters)

Segment = 1;

%% design parameters
SegmentDiameter = 1.41849;
NSegments = 798;
Mask.Spider.LegWidth = 0.5;
Mask.Spider.LegNumber = 6;
Mask.Spider.Theta = 60/180*pi;
Mask.Spider.Theta0 = 30/180*pi;
Mask.Diameter.Out = 38.542;
Mask.Diameter.In = 11.067;

%% load the data
Data  = LoadData(File.ErrorItem);
TimeSeries = isstruct(Data);
if TimeSeries
    Time = Data.PTTtimeDataOptical(:,1);
    PTT = Data.PTTtimeDataOptical(:,2:end);
    h1 = figure;
    plot(Time, PTT(:,(Segment-1)*3+(1:3)));
    legend('Piston','Tip','Tilt');
    drawnow
    title(sprintf('PTT of segment #%.d',Segment))
    R = input('Three options:\n [1] instantaneous time\n [2] time period\n [3] time period with specified duration\n please choose (type 1, 2 or 3) :\n');
    if R == 1
        fprintf(2,'Specify time in the plot using one mouse click...\n');
        t = ginput(1);
        t = t(1);
    elseif R == 2
        fprintf(2,'Specify start and end time in the plot using two mouse clicks...\n');
        t = ginput(2);
        t = t(:,1);
    elseif R == 3
        fprintf(2,'Specify start time in the plot using one mouse click...\n');
        t = ginput(1);
        t = t(1);
        duration = str2double(input(sprintf('Specify the time duration in s ...: '),'s'));
        t(2) = t(1) + duration;
    else
        error('not supported')
    end
else
    error('this is no times series')
end
SegmentVertices = LoadData(File.Segmentation);
Vertices = SegmentVertices(:,1:2);

%% create a grid
N = round(Parameters.Size/Parameters.Resolution);
x = ((1:N) - (1+N)/2) * Parameters.Resolution;
[X,Y] = meshgrid(x);
PixelPositions = [X(:),Y(:)];

%% map pixels to segments
SegmentNumber = FindHexagon(PixelPositions, Vertices, SegmentDiameter);

%% build phase screen
PS.Logical = false(size(X));
PS.Logical(:) = logical(SegmentNumber);
PS.Positions = zeros(sum(PS.Logical(:)),2);
PS.Positions(:,1) = X(PS.Logical(:));
PS.Positions(:,2) = Y(PS.Logical(:));
PS.SegmentNumber = SegmentNumber(logical(SegmentNumber));
    %tt = [X(PS.Logical(:)), Y(PS.Logical(:))];
    
%% compute the spider shadow
Spider = false(numel(X,1),1);
for leg = 1:Mask.Spider.LegNumber
    phi = (leg-1)*Mask.Spider.Theta + Mask.Spider.Theta0;
    P = Mask.Diameter.In/2 * [cos(phi), sin(phi)];
    par = [cos(phi); sin(phi)];
    per = [-sin(phi); cos(phi)];
    in = [X(:),Y(:)] * par - P * par > 0;
    thickness = abs([X(:),Y(:)] * per - P * per) - Mask.Spider.LegWidth/2 < 0;
    Spider = Spider | (in.*thickness);
end
R = sqrt(X.^2+Y.^2);
Obscuration = false(size(X));
Obscuration(logical(Spider)) = true;
Obscuration(R(:) > Mask.Diameter.Out/2 | R(:) < Mask.Diameter.In/2) = true;
PS.Obscuration = Obscuration;
Pupil = PS.Logical .* ~PS.Obscuration;

%% convert the tabulated Zernike coefficients in wavefront
if numel(t) == 1 % single frame
    PS.WFE = zeros(sum(PS.Logical(:)),1);
    dTime = Time-t;
    ind = find(abs(dTime) == min(abs(dTime)));
    ptt = reshape(PTT(ind,:),3,NSegments)';
    for s = 1:NSegments
        LocInd = PS.SegmentNumber == s;
        LocPos = [PS.Positions(LocInd,1),PS.Positions(LocInd,2)];
        LocPos = bsxfun(@minus, LocPos, Vertices(s,:));
        Z = ComputeZernikeModes(LocPos,(1:3),SegmentDiameter);
        if TimeSeries
            Z(:,3) = -Z(:,3);
        end
        PS.WFE(LocInd) = Z * ptt(s,:)';
    end 
    h2 = figure;
    RMS = sqrt(mean(PS.WFE.^2));
    ps = zeros(size(X));
    ps(PS.Logical(:)) = PS.WFE;
    imagesc(ps), axis equal tight
    title(sprintf('RMS WFE: %.0f nm', RMS*1e9))
    if File.Save.Flag% save the phase screen
        if isempty(File.Save.File)
            error('you need to specify an output file')
        end
        switch File.Save.Format
            case 'fits'
                fitswrite(ps,[File.Save.File,'.fits']);
                fitswrite(Pupil,[File.Save.File,'_Pupil','.fits']);
            case {'pdf', 'png'}
                saveas(h2, [File.Save.File,File.Save.Format])
            otherwise
                error('format not supported')
        end
    end
else % time series of frames 
    dt = t(2)-t(1);
    NInc = dt/mean(diff(Time));
    R = input(sprintf('you are generating %.d phase screens, do you want to proceed (y/n): ',ceil(NInc)),'s');
    if ~strcmpi(R,'y')
        error('...')
    end
    ind = Time > t(1) & Time < t(2);
    LocTime = Time(ind);
    fitswrite(Pupil,[File.Save.File,'_Pupil','.fits']);
    h2 = figure;
    wb = waitbar(0,'computing phase screens...');
    for ii = 1:sum(ind)
        PS.WFE = zeros(sum(PS.Logical(:)),1);
        dTime = Time-LocTime(ii);
        locind = find(abs(dTime) == min(abs(dTime)));
        ptt = reshape(PTT(locind,:),3,NSegments)';
        for s = 1:NSegments
            LocInd = PS.SegmentNumber == s;
            LocPos = [PS.Positions(LocInd,1),PS.Positions(LocInd,2)];
            LocPos = bsxfun(@minus, LocPos, Vertices(s,:));
            Z = ComputeZernikeModes(LocPos,(1:3),SegmentDiameter);
            if TimeSeries
                Z(:,3) = -Z(:,3);
            end
            PS.WFE(LocInd) = Z * ptt(s,:)';
        end
        RMS = sqrt(mean(PS.WFE.^2));
        ps = zeros(size(X));
        ps(PS.Logical(:)) = PS.WFE;
        imagesc(ps), axis equal tight
        title(sprintf('RMS WFE: %.0f nm', RMS*1e9))
        if File.Save.Flag % save the phase screen
            if isempty(File.Save.File)
                error('you need to specify an output file')
            end
            switch File.Save.Format
                case 'fits'
                    fitswrite(ps,[File.Save.File,'_',num2str(ii),'.fits']);
                case {'pdf', 'png'}
                    saveas(h2, [File.Save.File,'_',num2str(ii),File.Save.Format])
                otherwise
                    error('format not supported')
            end
        end
        waitbar(ii/sum(ind))
    end
    delete(wb)
end

%%  FONCTIONS  UTILISEES

function Data = LoadData(File)

Pos = max(strfind(File,'.'));
Extension = File(Pos+1:end);
switch Extension
    case 'txt'
        Data = importdata(File);
    case 'mat'
        Data = load(File);
end

%

function SegmentNumber = FindHexagon(PixelPositions, Vertices, SegmentDiameter)

R2_Pixels = sum(PixelPositions.^2,2);
R2_Vertices = sum(Vertices.^2,2);
R2_max = (sqrt(max(R2_Vertices))+(SegmentDiameter/2))^2;
in = find(R2_Pixels<=R2_max);
xy = (PixelPositions(in,:));
Angles =  pi/3*[0:2]';
XY2UV = 4/3*[cos(Angles),sin(Angles)]/SegmentDiameter;
uvs = (XY2UV * Vertices')';
uvsr = round(uvs);
ReconstructedRadius = uvsr*pinv(XY2UV)';
ReconstructedRadius = sqrt(sum(ReconstructedRadius.^2,2));
order = 3;
radius = sqrt(sum(Vertices.^2,2));
m = zeros(numel(radius),order+1);
for u = 0:order
    m(:,u+1) = radius.^u;
end
im = (m'*m)\m';
c = im*ReconstructedRadius;
radius = sqrt(sum(xy.^2,2));
m = zeros(numel(radius),order+1);
for u = 0:order
    m(:,u+1) = radius.^u;
end
xy = xy.*repmat((m*c)./radius,1,size(xy,2));
uv = xy*XY2UV';
uvr = round(uv);
frac = abs(uvr - uv);
[~,rejected] = max(frac,[],2);
loc = zeros(size(uv,1),1);
for k = 1:3
    rej = find(rejected==k);
    range = [mod(k,3)+1,mod(k+1,3)+1];
    [~,loc(rej)] = ismember(uvr(rej,range),uvsr(:,range),'rows');
end
SegmentNumber = zeros(size(PixelPositions,1),1);
SegmentNumber(in) = loc;

%

function Z_array = ComputeZernikeModes(LocPos,Orders,SegmentDiameter)

Positions = bsxfun(@rdivide,LocPos,SegmentDiameter/2);
One = Orders == 1;
% compute the R Positions
R_list = sqrt(sum(Positions.^2,2));
theta_list = atan2(Positions(:,2),Positions(:,1));
MaxOrder = max(Orders);
Z_array = zeros(size(Positions,1),numel(Orders));
[N,M] = FindNM(1:MaxOrder);
M = max(M);
N = max(N);
cosi = zeros(numel(theta_list),M+1);
sinu = zeros(numel(theta_list),M+1);
for mode = 0:M
    cosi(:,mode+1) = cos(mode*theta_list);
    sinu(:,mode+1) = sin(mode*theta_list);
end
Rs(:,1) = ones(size(R_list));
factoriel = ones(N+1,1);
for mode = 1:N
    Rs(:,mode+1) = Rs(:,mode).*R_list;
    factoriel(mode+1) = factoriel(mode)*mode;
end
for mode = 1:numel(Orders)
    Z_array(:,mode) = Zernike_polar(...
        Orders(mode),theta_list,R_list,...
        'cos',cosi,...
        'sin',sinu,...
        'radii',Rs,...
        'fact',factoriel);
end
OPD = Z_array;
OPD(:,~One) = OPD(:,~One)-repmat(mean(OPD(:,~One)),size(OPD,1),1);
Z_array = OPD;
for mode = 1:numel(Orders)-1
    m = OPD(:,mode);
    im = m'/sum(m.^2);
    Coeff = im*OPD(:,mode+1:end);
    OPD(:,mode+1:end) = OPD(:,mode+1:end) - OPD(:,mode)*Coeff;
    Z_array(:,mode+1:end) = Z_array(:,mode+1:end) - Z_array(:,mode)*Coeff;
end
NormalizationFactors = ones(1,numel(Orders));
NormalizationFactors(One) = mean(OPD(:,One));
D = bsxfun(@minus,OPD(:,~One),mean(OPD(:,~One)));
NormalizationFactors(~One) = sqrt(mean(D.^2));
Z_array = bsxfun(@rdivide,Z_array,NormalizationFactors);

%

function Zern = Zernike_polar(mode,theta,r,varargin)

if isempty(varargin)
    if r>1
        Zern=0;
    else
        [n,m] = FindNM(mode);
        R = 0;
        for s=0:(n-m)/2
            R = R + (-1).^s.*prod(1:(n-s)).*r.^(n-2.*s)./(prod(1:s).*prod(1:((n+m)/2-s)).*prod(1:((n-m)/2-s)));
        end
        radial = sqrt(n+1).*R;
        if m==0
            Zern = radial;
        else
            if rem(mode,2)==0
                Zern = cos(m*theta).*radial*sqrt(2);
            else
                Zern = sin(m*theta).*radial*sqrt(2);
            end
        end
    end
else
    for u = 1:2:length(varargin)
        switch varargin{u}
            case 'cos'
                cosi = varargin{u+1};
            case 'sin'
                sinu = varargin{u+1};
            case 'radii'
                radii = varargin{u+1};
            case 'fact'
                fact = varargin{u+1};
        end
    end
    if r>1
        Zern=0;
    else
        [n,m] = FindNM(mode);
        R = 0;
        for s=0:(n-m)/2
            R = R + (-1)^s*fact(1+(n-s))*radii(:,n-2*s+1)/(fact(1+s)*fact(1+((n+m)/2-s))*fact(1+((n-m)/2-s)));
        end
        radial = sqrt(n+1).*R;  
        if m==0
            Zern = radial;
        else
            if rem(mode,2)==0
                Zern = cosi(:,m+1).*radial*sqrt(2);
            else
                Zern = sinu(:,m+1).*radial*sqrt(2);
            end
        end
    end
end

%

function [nf,mf,Sign] = FindNM(mode)

nf = ceil(sqrt(2*mode+0.25)-1.5);
mf = mode - nf.*(nf+1)/2;
odd = logical(mod(nf,2));
% even = find(~mod(nf,2));
mf(odd)  = 2*floor((mf(odd )+1)/2)-1;
mf(~odd) = 2*floor(mf(~odd)/2);
Sign = mod(mode - (nf-1).*nf/2,2);
