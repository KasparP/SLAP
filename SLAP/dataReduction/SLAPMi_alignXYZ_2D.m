function [scandata, P, refIM]= SLAPMi_alignXYZ_2D (scandata, refIM, opts)
%finds the optimal offsets in X,Y,and Z, relative to the reference
%image for a 2D SLAPmi acquisition

if nargin<3
    opts.doMotion = false; %estimate per-frame motion or simply global offset?
end

is3D = any(diff(scandata.frames(1).Z)); %is the scan 3D?
centerXYOffset = [0 0]; %expected offset of scanimage raster scan
maxlag = 300; %maximum sample motion in lxls
refIM.M.pixelsize = diff(refIM.M.coords.X(1:2));

%preprocess the reference image
if ~isfield(refIM, 'data') && isfield(refIM, 'IM')
    IM = refIM.IM;
else
    IM = double(squeeze(refIM.data(:,:,1,:)));
end
if scandata.metadata.aperture
    IM = max(0, IM- prctile(IM(:), 80)); %this makes it easier to align y to yE
end
scandata.refIMcoords = refIM.M.coords;
IM = imtranslate(IM, centerXYOffset);

%get a PSF at the same spacing as the reference image
P = linePSF_full(scandata);

%the Z locations we'll consider
Zrange = 2; %must be an integer
res=  1;
nsteps = floor((Zrange)/res); %number of planes to evaluate in addition to the center plane
Zs = linspace(-nsteps*res, nsteps*res, 2*nsteps+1);
Pz = length(P.coords{3});
IMpadded = zeros(size(IM,1), size(IM,2), size(IM,3)+2*Pz);
IMpadded(:,:,Pz+1:Pz+size(IM,3)) = IM;

%the actual projection data, average over frames  ??Median??
y = mean([scandata.frames.pmtData],2); 
y(isnan(y)) = 0;

%calculate expected projection for each integer Z; we can interpolate
%these for superresolution localization
yE = expectedProjection(P.P, Pz, IMpadded, Zs);
yE = mean(y) .* yE;

%align the grand time-averaged signal
%[bestXY, bestZ] = alignYtoYE_DTW (y, yE, scandata, maxlag);
[bestXY, bestZ] = alignYtoYE (y, yE, scandata, maxlag);

%transform from lxl-space shift to a voltage-space shift
Xtform = [diff(scandata.Vx(1:2)) diff(scandata.Vy(1:2))]; %unit lxl vector in voltage space
Ytform = [-diff(scandata.Vx(end-1:end)) diff(scandata.Vy(end-1:end))];

%shift in pixel space
shiftPX = (bestXY* [Xtform ; Ytform]) ./refIM.M.pixelsize;

%display
IMshifted = imtranslate(IMpadded, -shiftPX);
figure('name', 'shifted image'), imshow3D(IMshifted)
SLAPMi_backproject(scandata,P);
drawnow;
%for single-plane imaging, aligning on each axis tends to identify the center Z
%plane of the PSF with the brightest plane of the sample regardless of offset. 
% This is because correlation doesn't take amplitude into account and
%favors the peakier signal. Once we've corrected in XY we should search in Z again so relative brightness
%information from the different angles can be taken into account
if ~is3D
    res = 0.5; nsteps = floor((Zrange)/res); %number of planes to evaluate in addition to the center plane
    Zs = linspace(-nsteps*res, nsteps*res, 2*nsteps+1);
    yE = expectedProjection(P.P, Pz, IMshifted, Zs);
    yE = mean(y) .* yE;
    Zcorrs = nan(1,size(yE,2));
    for plane = 1:length(Zcorrs)
        Zcorrs(plane) = corr(yE(:,plane), y);
    end
    bestZix = centroid(Zcorrs);
    bestZ = interp1(Zs, bestZix);
end

if opts.doMotion %now fit the location for each timepoint
    keyboard
else
   motion = [shiftPX -bestZ]; keyboard %%%%not sure is this should be +bestZ or -bestZ%%%%
   
   %compare the new expected projection to the data
   Pm = applymotion(P, -motion([2 1 3]), size(IM,3));
   yE = Pm.P*IM(:);
   yE = mean(y)* yE./mean(yE);
   figure('name', 'Mean Actual vs Expected projections after alignment'), plot(yE); hold on, plot(y)
end

end

function [bestXY, bestZ] = alignYtoYE (y, yE, scandata, maxlag)
%crosscorrelate signal with expected.
for line = 4:-1:1
    line_ixs = scandata.line==line;
    yline = y(line_ixs);
    yEline = yE(line_ixs, :);
    CC{line} = normxcorr2_general(yline,yEline); %a handy fairly-fast FEX implementation
    if line==2 || line==3
        F = griddedInterpolant({sqrt(2)*(1-length(yline):length(yline)-1), 1:size(yE,2)}, CC{line});
        CC{line} = F({-maxlag:maxlag, 1:size(yE,2)});
    else
        CC{line} =  CC{line}(length(yline)+ (-maxlag:maxlag),:);
    end
end
%for each pair of Xshift and Yshift, compute the sum of the
%relevant xcorr functions, to produce a 2d objective function.
%coordinate convention for lines (assumes that scan angles are
%separated by precisely 45 degrees... could get this from PSF)
%line1 is '+X'.
%line2 is '-X -Y'.
%line3 is '-X +Y'.
%line4 is '+Y'
[Xshift, Yshift] = meshgrid(-maxlag/2:maxlag/2);
REG = CC{1}(Xshift(:)+maxlag+1,:) + CC{2}(-Yshift(:) - Xshift(:) +maxlag+1,:) + CC{3}(-Xshift(:)+Yshift(:)+maxlag+1,:) + CC{4}(-Yshift(:)+maxlag+1,:);
REG = reshape(REG, [size(Xshift) size(yE,2)]);
[~, maxix] = max(REG(:));
[m, n, bestZ] = ind2sub(size(REG), maxix);
bestXY = [m n] - maxlag/2 - 1;
figure, imshow3D(REG); hold on, scatter(n,m, 'r')
end

function yE = expectedProjection(P,Pz, IMpadded, Zs)
Zrange = ceil(max(Zs));
IMcenter = ceil(size(IMpadded,3)/2);
yE = nan(size(P,1), 2*Zrange+1);
for dZ = -Zrange:Zrange
    imChunk = IMpadded(:,:,IMcenter+dZ + (-(Pz-1)/2:(Pz-1)/2));
    yE(:, dZ+Zrange+1) = P*imChunk(:);
end
yZint = griddedInterpolant({1:size(yE,1), -Zrange:Zrange}, yE);
yE = yZint({1:size(yE,1), Zs});
yE = yE./mean(yE(:));
end

function peak = centroid(I)
    %find center of mass
    I = max(0, I-median(I));
    peak = sum((1:length(I)).*I)./sum(I);
end
