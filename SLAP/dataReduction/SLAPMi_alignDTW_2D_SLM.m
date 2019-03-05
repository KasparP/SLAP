function [scandata, P, refIM, hFig] = SLAPMi_alignDTW_2D_SLM (scandata, refIM, forceAlign, optsin)
%aligns a scan using 2D translation and dynamic time warping while modeling
%the aperture formed by the SLM; this is the richest but most
%computationally expensive alignment model

opts.dense = false; %whether the reference image is dense
opts.PSFtype = 'X'; tooltips.PSFtype = 'delta, full, or X'; %
opts.Zrange = 3; tooltips.Zrange = 'Z range to search for alignment, in pixels. Must be an integer';
opts.Zres=  1; tooltips.Zres = 'Z resolution for alignment, in pixels.can be 0.5 for superresolution';
opts.prior = 0.0005; tooltips.prior = 'strength of spatial prior for center of alignment volume; 0.0005';

maxlag = 300;
maxlag_DTW = 25; %a parameter for how much DTW can warp the measured Ys in time, to match the expected
minSBR = 0.6; %minimum signal-to-background ratio
alignChan = 1;
nIter = 1; %precision of XY alignment; precision = 1/(4^(nIter-1)) pixels

if nargin>3 %UPDATE WITH USER-SPECIFIED OPTIONS
    for field = fieldnames(optsin)'
            opts.(field{1}) = optsin.(field{1});
    end
else
    opts = optionsGUI(opts, tooltips);
end

if scandata.metadata.do3D
    error('2D align function was called with 3D data!')
end

if ~isfield(refIM, 'mask')
    %convert SLM voltages to transmission
    Vslm = scandata.metadata.SLM.pattern;
    SLMhigh  = abs(double(Vslm)-scandata.metadata.calib.SLM.lut(:,:,2))<abs(double(Vslm)-scandata.metadata.calib.SLM.lut(:,:,1));
    Tslm = scandata.metadata.calib.SLM.T(:,:,2).* SLMhigh + scandata.metadata.calib.SLM.T(:,:,1).*~SLMhigh;
    Tslm = max(Tslm, 1e-2); %Sometimes the SLM transmission seems to be underestimated, let's assume it's kind of large
    
    %interpolate onto image space
    scanfield = refIM.metadata.rois.RoiGroups.imagingRoiGroup.rois{1, 1}.scanfields(1);
    scanfield.meshgrid = @(varargin)(meshgrid(refIM.M.coords.X, refIM.M.coords.Y));
    pixelToRefTransform = scandata.metadata.SLM.pixelToRefTransform;
    Tim = mapSLMToImage(Tslm, pixelToRefTransform, scanfield);
    
    %smooth the mask
    mask = imgaussfilt(Tim, 0.005/scanfield.pixelToRefTransform(1))';
    refIM.mask = mask;
end
mask = refIM.mask;
scandata.refIMcoords = refIM.M.coords;

%get a PSF at the same spacing as the reference image
switch opts.PSFtype
    case 'delta'
        P = linePSF_delta(scandata); %a delta function PSF
    case 'full'
        P = linePSF_full(scandata);  %a psf that incorporates the measured point spread function
    case 'X'
        P = linePSF_X(scandata); %a delta-like PSF that compensates for the difference between the line and raster PSFs, without requiring a measured point spread function
    otherwise
        error('Invalid PSF type selected')
end
nChan = size(scandata.frames(1).pmtData,2);
nFrames = length(scandata.frames);
y = [scandata.frames.pmtData];
y = permute(reshape(y, [], nChan, nFrames), [1 3 2]);
y = y(:,:,alignChan); %for multicolor data we align on only one channel for now

%remove nans 
y(isnan(y)) = 0;  %TMP this is temporary, we should make the pipeline deal with missing data

%adjust line brightness of P
disp('AlignDTW: Adjusting P for measured line brightness')
yB = nan(length(scandata.line),1);
for line = 1:4
    line_ixs = scandata.line==line;
    yB(line_ixs) = mean(mean((y(line_ixs,:,alignChan))));
end
P.P = spdiags(yB,0,length(yB),length(yB))*P.P;

%the Z locations we'll consider
nsteps = floor((opts.Zrange)/opts.Zres); %number of planes to evaluate in addition to the center plane
Zs = linspace(-nsteps*opts.Zres, nsteps*opts.Zres, 2*nsteps+1);
Pz = length(P.coords{3});

IMpadded = zeros(size(refIM.IM,1), size(refIM.IM,2), size(refIM.IM,3)+2*Pz);
if opts.dense
    IMpadded(:,:,Pz+1:Pz+size(refIM.IM,3)) = refIM.IM;
else
    IMpadded(:,:,Pz+1:Pz+size(refIM.IM,3)) = max(0, refIM.IM- median(refIM.IM(:)))+ mean(refIM.IM(:))/1000;
end
mask3D = repmat(mask, 1, 1, size(IMpadded,3));

%the actual projection data, averaged over frames
disp('AlignDTW: Using mean of observations for reference')
y = [scandata.frames.pmtData];
y = permute(reshape(y,[], nChan, nFrames), [1 3 2]);
y = yRef(y(:,:,alignChan));

if nargin>2 && ~isempty(forceAlign) %the alignment has been given to us
    deltas = 0;
    [gX,gY,gZ,gD,gT] = alignMultiScale(forceAlign(1),forceAlign(2), deltas);
    dX = gX; dY = gY; dZ=Zs(gZ); bestZ=gZ; bestD = gD; Ta = gT;
else %multiscale alignment; scale of 4 and 1 iteration gives up to 10 pixels of alignment shift
    scale = 4;
    deltas = scale*(-3:3);
    [bestX,bestY,bestZ, bestD,bestT] = alignMultiScale(0,0, deltas);
    for iter = 1:nIter
        scale = scale/4; deltas = scale*(-2:2);
        [gX,gY,gZ,gD,gT] = alignMultiScale(bestX,bestY, deltas);
        if gD<bestD
            bestX = gX; bestY = gY; bestZ = gZ; bestD = gD; bestT = gT;
        end
    end
    dX = bestX; dY = bestY; dZ = Zs(bestZ); Ta = bestT;
end

disp('Alignment offsets:')
disp(['X: ' num2str(dX) ' px'])
disp(['Y: ' num2str(dY) ' px'])
disp(['Z: ' num2str(dZ) ' px'])

%Expected data for original alignment
[yE, ~] = expectedProjection(P.P, Pz, IMpadded, mask3D, Zs,y);

%Resample the data
pdata = [scandata.frames.pmtData];
pdata = permute(reshape(pdata,[], nChan,  nFrames), [1 3 2]);
for line = 1:4
    line_ixs = scandata.line==line;
    for ch = 1:nChan
        pdata(line_ixs,:,ch) = resample_dtw(pdata(line_ixs,:,ch), Ta{line});
    end
end
for ch = 1:nChan
    for f = 1:length(scandata.frames)
        scandata.frames(f).pmtData(:,ch) = pdata(:,f,ch); %assign
    end
end

%resample refIM to recenter and match number of planes with P
ZplanesP = ceil(size(refIM.IM,3)/2) + (-(Pz-1)/2:(Pz-1)/2);
refIM.M.coords.Z = interp1(refIM.M.coords.Z, ZplanesP -dZ);
refIM.IM = imtranslate(refIM.IM, [dX dY -dZ]);
refIM.IM = refIM.IM(:,:,ZplanesP);
if isfield(refIM, 'labels')
    refIM.labels = imtranslate(refIM.labels, [dX dY -dZ]);
    refIM.labels = refIM.labels(:,:,ZplanesP);
end
for ch = 1:size(refIM.data,4)
    refIM.data(:,:,:,ch) = imtranslate(refIM.data(:,:,:,ch), [dX dY -dZ]);
end
refIM.data = refIM.data(:,:,ZplanesP,:);

%update segmentation
if isfield(refIM, 'seg')
    [XX,YY,ZZ] = ndgrid(1:size(refIM.bw,1), 1:size(refIM.bw,2), 1:size(refIM.bw,3));
    XX = XX-dY; YY = YY-dX; ZZ = ZZ+dZ;
    outOfRange = XX<1 | XX>size(refIM.bw,2) | YY<1 | YY>size(refIM.bw,1) | ZZ<1 | ZZ>size(refIM.bw,3);
    refIM.seg(~outOfRange,:) = refIM.seg(sub2ind(size(refIM.bw), XX(~outOfRange), YY(~outOfRange), ZZ(~outOfRange)),:);
    refIM.seg(outOfRange,:) = 0;
    
    select = false(size(refIM.bw)); select(:,:,ZplanesP) = true;
    refIM.seg = refIM.seg(select(:),:);
    
    refIM.bw = imtranslate(refIM.bw, [dX dY -dZ]);
    refIM.bw = refIM.bw(:,:,ZplanesP);
end

%update the motion field
scandata.motion.dtw = Ta;
scandata.motion.dtwGlobalDone = true;
scandata.motion.dXYZ = [dX,dY,dZ];
scandata.motion.error = bestD;

imTMP = imtranslate(IMpadded, [dX dY -dZ]).*mask3D;
IM2 = imTMP(:,:,ceil(size(imTMP,3)/2) + (-(Pz-1)/2:(Pz-1)/2));

%readjust brightness of P for new alignment
y2 = yRef([scandata.frames.pmtData]);
IM2(isnan(IM2(:))) = 0;

[y2_E, b] = expectedProjection(P.P, Pz, IM2, mask3D(:,:,1:size(IM2,3)), 0,y);
coeff = b(1)+b(2).*minSBR;
SBRest = (b(1)+minSBR.*b(2))./b(2);
% notnan = ~isnan(IM2(:));
% y2_E = P.P(:,notnan)*IM2(notnan);
yB = nan(length(scandata.line),1);

for line = 1:4 %adjust P for measured brightness, allowing for slow intensity variations
     line_ixs = find(scandata.line==line);
     w = y(line_ixs);
     v = (y2(line_ixs)./y2_E(line_ixs)); 
     v = max(0.67, min(1.5, v));
     v([1:10 end-10:end]) = nanmean(v); 
     w([1:10 end-10:end]) = 1;
     alpha = csaps(line_ixs,v,1e-6,line_ixs,w);
     %alpha =  nanmean(y2(line_ixs))./ nanmean(y2_E(line_ixs));
     yB(line_ixs) = coeff.*alpha;
     y2_E(line_ixs) = y2_E(line_ixs).*alpha;
end
P.P = spdiags(yB,0,length(yB),length(yB))*P.P;

%display
hFig = figure('Name', ['dX: ' num2str(dX) ' dY:' num2str(dY) ' dZ: ' num2str(dZ) ' px; SBR~=' num2str(SBRest,3) ' ' scandata.metadata.galvoDataFileName]);
subplot(2,1,1)
plot(y./nanmean(y)); hold on, plot(yE(:,bestZ)./nanmean(yE(:,bestZ)))
legend({'BEFORE align', 'EXPECTED'}); xlabel('Actual vs expected projection before alignment');
subplot(2,1,2)
plot(y2./nanmean(y2)); hold on, plot(y2_E./mean(y2_E))
legend({'AFTER align', 'EXPECTED'}); xlabel('Actual vs expected projection after alignment')
drawnow

%cut down P and S
Psz = [length(P.coords{1}), length(P.coords{2}), length(P.coords{3})];
switch opts.PSFtype
    case 'X'
        selectCenter = false(Psz); selectAbove =  false(Psz); selectBelow =  false(Psz);
        selectCenter(:,:,P.coords{3}==0) = true;
        selectAbove(:,:,P.coords{3}==1.5) = true;
        selectBelow(:,:,P.coords{3}==-1.5) = true;

        centerplane = (refIM.data(:,:,P.coords{3}==0,:).*P.opts.weights(1) + P.opts.weights(2).*refIM.data(:,:,P.coords{3}==0.75,:) + P.opts.weights(2).*refIM.data(:,:,P.coords{3}==-0.75,:))./(P.opts.weights(1)+2*P.opts.weights(2));
        refIM.data = cat(3, refIM.data(:,:,P.coords{3}==-1.5,:), centerplane,  refIM.data(:,:,P.coords{3}==1.5,:));
        
        centerplane = (refIM.IM(:,:,P.coords{3}==0).*P.opts.weights(1) + P.opts.weights(2).*refIM.IM(:,:,P.coords{3}==0.75) + P.opts.weights(2).*refIM.IM(:,:,P.coords{3}==-0.75))./(P.opts.weights(1)+2*P.opts.weights(2));
        refIM.IM = cat(3, refIM.IM(:,:,P.coords{3}==-1.5), centerplane, refIM.IM(:,:,P.coords{3}==1.5));
        
        if isfield(refIM, 'seg') %if we are using a segmentation; the usual case
            centerplane = max(refIM.labels(:,:,P.coords{3}<=0.8), [],3);
            refIM.labels = cat(3,refIM.labels(:,:,P.coords{3}==-1.5), centerplane, refIM.labels(:,:,P.coords{3}==1.5));

            selectOneAbove = false(Psz); selectOneBelow = false(Psz);
            selectOneAbove(:,:,P.coords{3}==0.75) = true;
            selectOneBelow(:,:,P.coords{3}==-0.75) = true;
            refIM.seg(selectCenter,:) = (P.opts.weights(1).*refIM.seg(selectCenter,:) + P.opts.weights(2).*(refIM.seg(selectOneAbove,:)+refIM.seg(selectOneBelow,:)))./(P.opts.weights(1)+2*P.opts.weights(2));
            refIM.seg = refIM.seg(selectCenter(:) | selectAbove(:) | selectBelow (:),:);
            refIM.bw = cat(3, refIM.bw(:,:,P.coords{3}==-1.5), any(refIM.bw(:,:,abs(P.coords{3})<0.8),3), refIM.bw(:,:,P.coords{3}==1.5));
        end
        if abs(abs(diff(refIM.M.coords.Z(1:2)))-0.75)>0.1
            warning('Your image was acquired at incorrect Z spacing for the X PSF; it should be acquired at 0.75um resolution. Continuing...')
        end
        refIM.M.coords.Z = [-1.5 0 1.5];
        
        P.P = [P.P(:,selectBelow) P.P(:,selectCenter).*(P.opts.weights(1)+2*P.opts.weights(2)) P.P(:,selectAbove)];
        P.coords{3} = [-1.5 0 1.5];
        
    otherwise
        sumP = sum(sum(reshape(full(sum(P.P,1)),Psz)));
        selectPlanes = sumP>max(sumP)/2;
        select = false(Psz); select(:,:,selectPlanes) = true;
        P.P = P.P(:, select(:));
        if isfield(refIM, 'seg')
            refIM.seg = refIM.seg(select(:),:);
            refIM.bw = refIM.bw(select);
        end
        P.coords{3} = P.coords{3}(selectPlanes);
        refIM.data = refIM.data(:,:,selectPlanes,:);
        refIM.IM = refIM.IM(:,:,selectPlanes);
        refIM.labels = refIM.labels(:,:,selectPlanes);
end


%END MAIN FUNCTION

    function [dX,dY,dZ,bestD, Ta] = alignMultiScale(Xcenter,Ycenter,deltas, Zsamps)
       tic;
       if nargin<4
           Zsamps = Zs;
       end
           
       T = cell(length(deltas),length(deltas),length(Zsamps),4);
       D = nan(length(deltas),length(deltas),length(Zsamps),4);
       for x_ix = 1:length(deltas)
           disp(['Multiscale alignment, ' int2str(x_ix) ' of ' int2str(length(deltas))])
           for y_ix = 1:length(deltas)
               %translate mask
                imTMP = imtranslate(IMpadded, [Xcenter+deltas(x_ix) Ycenter+deltas(y_ix)]);
                [yE, ~] = expectedProjection(P.P,Pz, imTMP, mask3D, Zsamps, y); %superresolution interpolated expected projection
                
%                 yE = mean(y) .* yE; 
                for line = 1:4
                    line_ixs = scandata.line==line;
                    [yline, ~] = filterY(y(line_ixs)); %highpass filter
                    [yEline, ~] = filterY(yE(line_ixs, :)); %highpass filter
                    CC = normxcorr2_general(yline,yEline, length(yline)-maxlag); %get the global offset
                    for z_ix = 1:length(Zsamps)
                        [~, lagix] = max(CC(:,z_ix)); %correct for average shift
                        lag = lagix-length(yline);
                        [Ta,D(x_ix,y_ix,z_ix,line)] = DTWkp(interp1(yline, (1:length(yline))-lag, 'linear', 0),yEline(:,z_ix), maxlag_DTW); %Dynamic Time Warping
                        T{x_ix,y_ix,z_ix,line} = Ta - lag;
                    end
                end
           end
       end
       Dsum = sum(D,4);
       
       %add a weak prior for the center of the alignment volume
       [priorA,priorB,priorC] = ndgrid(deltas, deltas, Zs*4);
       prior = 1 + opts.prior*(priorA.^2 + priorB.^2 + priorC.^2);
       
       [bestD,minix] = min(Dsum(:) .* prior(:));
       [Xix, Yix, Zix] = ind2sub(size(Dsum), minix);
       dX = Xcenter+deltas(Xix); dY = Ycenter+deltas(Yix); dZ = Zix;
       Ta = squeeze(T(Xix,Yix,Zix,:));
       toc;
    end


function [yE, b] = expectedProjection(P,Pz, IMpadded, mask, Zs,y)
IMpadded = IMpadded.*mask;
Zrng = ceil(max(Zs));
IMcenter = ceil(size(IMpadded,3)/2);
yE = nan(size(P,1), 2*Zrng+1);
for DZ = -Zrng:Zrng
    imChunk = IMpadded(:,:,IMcenter+DZ + (-(Pz-1)/2:(Pz-1)/2));
    yE(:, DZ+Zrng+1) = P*imChunk(:);
end
if length(Zs)>1
    yZint = griddedInterpolant({1:size(yE,1), -Zrng:Zrng}, yE);
    yE = yZint({1:size(yE,1), Zs});
end
yBG = P*mask(1:size(P,2))'; yBG = yBG  .*  mean(yE(:))./mean(yBG(:));

%regression
for z_ix = 1:length(Zs)
    C = [yE(:,z_ix) yBG+yE(:,z_ix).*minSBR ones(size(yBG))];
    b = lsqnonneg(C, y);
    yE(:,z_ix) = C*b;
end
end

end


function [x2, lowpass] = filterY(x)
lowpass = imgaussfilt(x,[25 eps], 'padding', 'symmetric');
%lowpass = imopen(x, strel(ones(500,1)));
%lowpass = min(lowpass, imgaussfilt(lowpass,[100 eps], 'padding', 'symmetric'));
x2 = x-lowpass;
end

function pmt = resample_dtw (pdata, T)
%similar to interp1 but averages samples between each epoch
Tedges = (T(1:end-1) + T(2:end))/2;
Tedges = [(2*Tedges(1) - Tedges(2)) ; Tedges(:) ; (2*Tedges(end) - Tedges(end-1))];
nans = isnan(pdata);
pdata(isnan(pdata)) = 0;
CSP = cumsum(pdata, 1);
CSP(nans) = nan;
pmt = diff(interp1(1:size(CSP,1), CSP, Tedges),1,1)./repmat(diff(Tedges,1,1), 1, size(CSP,2));
%alternative:
%pdata(line_ixs,:) = interp1(pdata, T);
end

function yOut = yRef(yIn)
    %get an estimate of the F0 line intensities
    yOut = imgaussfilt(yIn,[eps 150], 'padding', 'symmetric');
    yIn(isnan(yIn)) = nanmin(yOut(:));
    yOut = imgaussfilt(yIn,[eps 150], 'padding', 'symmetric');
    yOut = min(yOut,[],2);
    
%     nCut = round(size(yIn,2)/4);
%     [~, sortorder] = sort(nansum(yIn,1));
%     yOut = nanmean(yIn(:, sortorder(10:10+nCut)), 2);
%     yOut(sum(~isnan(yIn),2)<30) = nan;
%     %replace nans with prediction from surrounding measurements
%     yOut(isnan(yOut)) = max(0, min(yOut)); %, inpaint_nans(yOut));
end