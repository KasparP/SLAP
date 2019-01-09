function [scandata, P, refIM_new] = SLAPMi_alignDTW_2D (scandata, refIM, optsin)
%aligns a scan using dynamic time warping on each axis

opts.dense = false; %whether the reference image is dense
opts.PSFtype = 'delta'; %delta or full 
maxlag = 300;
maxlag_DTW = 30;
nFrames = 50;

if scandata.metadata.do3D
    error('2D align function was called with 3D data!')
end

%preprocess the reference image
if ~isfield(refIM, 'data') %if we got a segmentation
    error('SLAPMi_align no longer works with segmentation files; Pass your refIM, then segment after aligning!');
end

IM = double(squeeze(refIM.data(:,:,1,:)));
refIM.IM = IM;

if ~opts.dense
    IM = max(0, IM- prctile(IM(:), 80)); %this makes it easier to align y to yE
end

if isfield(scandata.metadata, 'SLM') %mask out the SLM in the reference image
    %TMP for backward compatibility
    if ~isfield(scandata.metadata.calib.SLM, 'map')
        calib = [];
        load('E:\SLAPmiData\Calibration\calibration.cal', '-mat');
    end
    %convert SLM voltages to transmission
    Vslm = scandata.metadata.SLM.pattern;
    SLMhigh  = abs(double(Vslm)-scandata.metadata.calib.SLM.lut(:,:,2))<abs(double(Vslm)-scandata.metadata.calib.SLM.lut(:,:,1));
    Tslm = scandata.metadata.calib.SLM.T(:,:,2).* SLMhigh + scandata.metadata.calib.SLM.T(:,:,1).*~SLMhigh;
    
    %interpolate onto image space
    scanfield = refIM.metadata.rois.RoiGroups.imagingRoiGroup.rois{1, 1}.scanfields(1);
    scanfield.meshgrid = @(varargin)(meshgrid(refIM.M.coords.X, refIM.M.coords.Y));
    pixelToRefTransform = scandata.metadata.SLM.pixelToRefTransform;
    Tim = mapSLMToImage(Tslm, pixelToRefTransform, scanfield);
    
    %smooth the mask
    mask = imgaussfilt(Tim, scanfield.pixelToRefTransform(1)*100)';
    refIM.mask = mask;
end

scandata.refIMcoords = refIM.M.coords;

%get a PSF at the same spacing as the reference image
switch opts.PSFtype
    case 'delta'
        P = linePSF_delta(scandata); %a delta function PSF
    case 'full'
        P = linePSF_full(scandata);  %a psf that incorporates the measured point spread function
    otherwise
        error('Invalid PSF type selected')
end
y = [scandata.frames.pmtData];

%remove nans
y(isnan(y)) = 0; % this should be ok because the nans are only ever the missed regions at the edges of each scan
% nans = any(isnan(y),2);
% y(nans,:) = 0;
% P.P(nans,:) = 0;
% nonans = ~any(isnan(y),2);
% scandata.line = scandata.line(nonans);
% scandata.Vx = scandata.Vx(nonans);
% scandata.Vy = scandata.Vy(nonans);
% for f = 1:length(scandata.frames)
%     %scandata.frames(f).Zcenter = scandata.frames(f).Zcenter(nonans);
%     scandata.frames(f).galvoData = scandata.frames(f).galvoData(nonans,:);
%     scandata.frames(f).valid = scandata.frames(f).valid(nonans);
%     scandata.frames(f).pmtData = scandata.frames(f).pmtData(nonans,:);
% end    
% P.P = P.P(nonans,:);
% y = y(nonans,:);

%adjust line brightness of P
disp('AlignDTW: Adjusting P for measured line brightness')
yB = nan(length(scandata.line),1);
for line = 1:4
    line_ixs = scandata.line==line;
    yB(line_ixs) = mean(mean((y(line_ixs,:))));
end
P.P = spdiags(yB,0,length(yB),length(yB))*P.P;

%the Z locations we'll consider
Zrange = 2; %must be an integer
res=  0.5;
nsteps = floor((Zrange)/res); %number of planes to evaluate in addition to the center plane
Zs = linspace(-nsteps*res, nsteps*res, 2*nsteps+1);
Pz = length(P.coords{3});
IMpadded = zeros(size(IM,1), size(IM,2), size(IM,3)+2*Pz);
IMpadded(:,:,Pz+1:Pz+size(IM,3)) = IM.*repmat(mask, 1, 1, size(IM,3));

%the actual projection data, average over frames  ??Median??
disp(['AlignDTW: Using first ' int2str(nFrames) ' frames for reference'])
y = nanmean([scandata.frames(1:nFrames).pmtData],2); y(isnan(y)) = 0;

%calculate expected projection for each integer Z; we can interpolate
%these for superresolution localization
yE = expectedProjection(P.P, Pz, IMpadded, Zs);
yE = mean(y) .* yE;

%compute a DTW distance for each line at each Z plane
CC = cell(1,4);
Ta = cell(4, length(Zs));
D = nan(4, length(Zs));
for line = 1:4
    line_ixs = scandata.line==line;
    yline = filterY(y(line_ixs)); %nonlinear highpass filter
    yEline = filterY(yE(line_ixs, :)); %nonlinear highpass filter
%      yline = y(line_ixs);
%      yEline = yE(line_ixs, :);
    CC{line} = normxcorr2_general(yline,yEline, length(yline)-maxlag);
    for z = 1:length(Zs)
        [~, lagix] = max(CC{line}(:,z));
        lag = lagix-length(yline);
        
        yline2 = [zeros(max(0,lag),1) ; yline ; zeros(max(0,-lag),1)];
        yline2 = smooth(smooth(yline2,3),3); yline2 = yline2./max(yline2); %max because some values can be negative; variance?
        
        yEline2 = [zeros(max(0, -lag),1) ; yEline(:,z) ; zeros(max(0,lag),1)];
        yEline2 = yEline2./max(yEline2); %max because some values can be negative; variance?
        
        [T,D(line,z)] = DTWkp(yline2,yEline2, maxlag_DTW); %Dynamic Time Warping
        Ta{line,z} = T(max(0,-lag) + (1:length(yline))) - max(lag,0);
    end
end
[~,bestZ] = min(sum(D,1));
Ta = Ta(:,bestZ);

%Resample the data and update the projection matrix
pdata = [scandata.frames.pmtData];
for line = 1:4
    line_ixs = scandata.line==line;
    pdata(line_ixs,:) = resample_dtw(pdata(line_ixs,:), Ta{line});
end
for f = 1:length(scandata.frames) 
    scandata.frames(f).pmtData = pdata(:,f); %assign
end

%update the motion field
scandata.motion.dtw = Ta;
scandata.motion.dtwGlobalDone = true;

Zcenter = ceil(size(refIM.IM,3)/2);
refIM_new = interpolate_refIM(refIM, Zcenter+Zs(bestZ), Pz);
IM2 = refIM_new.IM;
IM2 = IM2.*repmat(mask, 1, 1, Pz);

%readjust brightness of P
y2 =[scandata.frames(1:nFrames).pmtData]; y2 = nanmean(y2,2);
notnan = ~isnan(IM2(:));

yE = expectedProjection(P.P, Pz, IMpadded, Zs);
y2_E = P.P(:,notnan)*IM2(notnan);
yB = nan(length(scandata.line),1);
for line = 1:4
     line_ixs = scandata.line==line;
     yB(line_ixs) = nanmean(y2(line_ixs))./ nanmean(y2_E(line_ixs));
end
P.P = spdiags(yB,0,length(yB),length(yB))*P.P;
y2_E = P.P(:,notnan)*IM2(notnan);

%dummy test
figure('Name', ['Alignment for ' scandata.metadata.galvoDataFileName]) 
subplot(2,1,1)
plot(y./nanmean(y)); hold on, plot(yE(:,bestZ)./nanmean(yE(:,bestZ)))
legend({'BEFORE align', 'EXPECTED'}); xlabel('Actual vs expected projection before DTW');
subplot(2,1,2)
plot(y2./nanmean(y2)); hold on, plot(y2_E./mean(y2_E))
legend({'AFTER align', 'EXPECTED'}); xlabel('Actual vs expected projection after DTW')

%validate by backprojecting
% IM_S = S.IM./nanmean(S.IM(:));
% BP = reshape(P.P'*full(y2_E), size(S.IM));
% BP = BP./nanmean(BP(:));
% figure('Name', ['Expected Backprojection for ' scandata.metadata.galvoDataFileName])
% imshow3D(cat(3, nanmean(IM_S,3), nanmean(BP,3)));
% 
% BP = reshape(P.P'*y2, size(S.IM));
% BP = BP./nanmean(BP(:));
% figure('Name', ['Actual Backprojection for ' scandata.metadata.galvoDataFileName])
% imshow3D(cat(3, nanmean(IM_S,3), nanmean(BP,3)));
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

function x2 = filterY(x)
lowpass = imopen(x, strel(ones(100,1)));
for col = 1:size(x,2)
    lowpass(:,col) = smooth(lowpass(:,col),300);
end
x2 = x-lowpass;
end

function peak = centroid(I)
  %find center of mass
    I = max(0, I-median(I));
    peak = sum((1:length(I)).*I)./sum(I);
end

function pmt = resample_dtw (pdata, T)
%similar to interp1 but averages samples between each epoch
Tedges = (T(1:end-1) + T(2:end))/2;
Tedges = [(2*Tedges(1) - Tedges(2)) ; Tedges(:) ; (2*Tedges(end) - Tedges(end-1))];
pdata(isnan(pdata)) = 0;
CSP = cumsum(pdata, 1);
pmt = diff(interp1(1:size(CSP,1), CSP, Tedges),1,1)./repmat(diff(Tedges,1,1), 1, size(CSP,2));
%alternative:
%pdata(line_ixs,:) = interp1(pdata, T);
end