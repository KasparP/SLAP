function refIM = SLAPMi_refIM(optsin)
%produces a reference image from a set of scanimage stacks

%options:
opts.noiseThresh = 30; tooltips.noisethresh = 'threshold for noise in the image';
%opts.estimatePhase:    forces phase estimation even if it has previously been cached
%opts.doMotion:         do motion correction
%opts.alignchan:        channel to use for alignment, 2 by default

%FORMAT FOR A REFIM:
%the scanimage stack should consist of 2 ROIs of the same square area,
%oriented at 90 degrees to each other (a 'fast x' and a 'fast y' image);
%The 'Lock to #frames' checkbox should be enabled
%At the moment, #frames per slice must be 1. This could be fixed easily with a
%bit of added code

basedir = 'E:\SLAPMidata\';
if ~nargin || ~isfield(optsin, 'fns')
    [fns, dr] = uigetfile([basedir filesep '*.tif'], 'multiselect', 'on');
    if ~iscell(fns)
        fns = {fns};
    end
    fns = sort_nat(fns);
else
    fns = optsin.fns;
    dr = optsin.dr;
    optsin =rmfield(optsin, {'fns', 'dr'});
end

opts.blackLevel = 3;                    tooltips.blackLevel = 'percentage of pixels to set to Black in reference image';
opts.alignChan = 1;                    tooltips.alignChan = 'which channel to use for alignment';
%opts.estimatePhase = true;             tooltips.estimatePhase = 'force bidi phase and warp estimation. if FALSE, will used a hashed value';
opts.doMotion = false;                 tooltips.doMotion = 'registed Motion during raster stacks';
opts.verbose = false;                  tooltips.verbose = 'Show figures?';
opts.blocks = false;                    tooltips.blocks = 'Were your stacks acquired in blocks that might need different phase corrections?'; 
if nargin %UPDATE WITH USER-SPECIFIED OPTIONS
    for field = fieldnames(optsin)'
            opts.(field{1}) = optsin.(field{1});
    end
    if ~isfield(optsin, 'doMotion') %force the options GUI if motion correction wasn't specified
        opts = optionsGUI(opts, tooltips);
    end
else
    opts = optionsGUI(opts, tooltips);
end
savefn = [fns{1}(1:end-10) '_REF.tif'];
savedr = dr;

%options not in the GUI
opts.gamma = 0.6;                      tooltips.gamma = 'gamma for saving 8bit images for ilastik; does not affect .mat files';
opts.satMultiplier = 1;                %increase this number if the 8-bit '.tif' reference images are dim. You will have to retrain Ilastik classifiers
opts.estimatePhase = true;

%Sanity check
for fnum = 1:length(fns)
    reader = ScanImageTiffReader([dr filesep fns{fnum}]);
    metadata = parseMetaData(reader.metadata());
    numChannels = length(metadata.hChannels.channelSave);
    alignChan = find(metadata.hChannels.channelSave==opts.alignChan,1);
    if isempty(alignChan)
        if numChannels==1
            alignChan=1;
            warning(['You selected alignment CH: ' int2str(opts.alignChan) ' but only CH:' int2str(metadata.hChannels.channelSave) ' was recorded. Continuing...'])
        else
            error(['The alignment channel ' int2str(opts.alignChan) ' was not recorded']);
        end
    end
    zs =  metadata.hStackManager.zs;
    %Make sure files are consistent:
    if metadata.hRoiManager.mroiEnable
        rois = metadata.rois.RoiGroups.imagingRoiGroup.rois;
        assert(all(rois{1}.scanfields.sizeXY == rois{2}.scanfields.sizeXY));
        assert(all(rois{1}.scanfields.centerXY == rois{2}.scanfields.centerXY));
        assert(rois{1}.scanfields.rotationDegrees==0 && rois{2}.scanfields.rotationDegrees==90);
        assert(all(rois{1}.scanfields.pixelResolutionXY == rois{2}.scanfields.pixelResolutionXY));
        assert(diff(rois{1}.scanfields.pixelResolutionXY) ==0);
        
        scanfield = rois{1}.scanfields;
        xres = scanfield.pixelResolutionXY(1);
        yres = scanfield.pixelResolutionXY(2);
        [xx,yy] = meshgrid(linspace(1/(xres*2),1-1/(xres*2),xres),linspace(1/(yres*2),1-1/(yres*2),yres));
        [xx,yy] = scanimage.mroi.util.xformMesh(xx,yy,scanfield.affine);
        scannerToRefTransform = metadata.hScan2D.scannerToRefTransform;
        [xx,yy] = scanimage.mroi.util.xformMesh(xx,yy,inv(scannerToRefTransform));
        X = xx(1,:); Y = yy(:,1)';
        if fnum==1
            M.coords.X = X;
            M.coords.Y = Y;
            M.coords.Z = zs;
        else
            assert(all(M.coords.X == X & M.coords.Y == Y));
            %assert(all(M.coords.Z == zs));
        end
    else
        keyboard %Not an MROI image
    end
end

%accumulate data
tic;
disp('Loading data...')
Lx = length(M.coords.X);
IMfastY = nan(length(M.coords.X), length(M.coords.X), length(M.coords.Z), length(fns), numChannels);
IMfastX = nan(length(M.coords.X), length(M.coords.X), length(M.coords.Z), length(fns), numChannels);
for fnum = 1:length(fns)
    fileName=[dr filesep fns{fnum}];
    reader = ScanImageTiffReader(fileName);
    vol=double(reader.data());
    for chan = 1:numChannels
        IMfastY(:,:,:, fnum, chan) = vol(:,1:Lx,chan:numChannels:end);
        IMfastX(:,:,:, fnum, chan) = rot90(vol(:, Lx + (1:Lx),chan:numChannels:end));
    end
end
toc;

%correct black level; this influences phase error adjustment
for ch = 1:numChannels
    tmp = IMfastX(:,:,:,:,ch); 
    thresh = prctile(tmp(randsample(numel(tmp), 1e6)) , opts.blackLevel);
    IMfastY(:,:,:,:,ch) = IMfastY(:,:,:,:,ch)-thresh;
    IMfastX(:,:,:,:,ch) = IMfastX(:,:,:,:,ch)-thresh;
end

%adjustment for scan phase errors
roiHash = DataHash({metadata.rois.RoiGroups.imagingRoiGroup.rois,  metadata.hRoiManager.scanFrameRate metadata.hScan2D.linePhase});
if exist([basedir '\GalvoTraces\RasterPhase\' roiHash '.mat'], 'file') && ~opts.estimatePhase
    load([basedir '\GalvoTraces\RasterPhase\' roiHash '.mat']);
else
    disp('Unknown mROI settings, fitting scan phase errors')
    nBlocks = 1;
    blockIDs = ones(1,length(fns))';
    if opts.blocks %get the blocks
        nBlocks = 0;
        done=false(size(fns));
        while any(~done)
            nBlocks = nBlocks+1;
            s = listdlg('PromptString',['Select files in block: ' int2str(nBlocks)],...
                'ListString',fns(~done));
            fdo= find(~done); fdo = fdo(s);
            blockIDs(fdo) = nBlocks;
            done(fdo) = true;
        end
    end
    shiftX = cell(1,nBlocks); shiftY = cell(1,nBlocks);
    for block = 1:nBlocks
        blockSelect = find(blockIDs ==block);
        IMX = IMfastX(:,:,:,blockSelect, alignChan); IMY = IMfastY(:,:,:,blockSelect,alignChan);
        
        B = squeeze(sum(sum(IMX,1),2));
        select = B>prctile(B(:), 100*max(1-(numel(blockSelect)/numel(B)), 0.85));
        IMX = IMX(:,:,select); IMY = IMY(:,:,select); %collapse across all frames here
        
        %find phase
        shiftX{block} = findPhase(IMX, opts.noiseThresh);
        shiftY{block} = findPhase(permute(IMY, [2 1 3]), opts.noiseThresh); %IMfastY{alignchan} = permute(IMfastY, [2 1 3]);
    end
end

% if ~opts.doMotion
%     medFX =  squeeze(nanmedian(IMfastX,4));
%     medFY =  squeeze(nanmedian(IMfastY,4));
%     figure('Name', 'fastX before BiDi correction'), imshow3D(medFX(:,:,:,1));
%     figure('Name', 'fastY before BiDi correction'), imshow3D(medFY(:,:,:,1));
% end

%fix phase for all images
tic;
disp('Correcting bidi phase...')
for fnum = 1:length(fns)
    for ch = 1:numChannels
        IMfastY(:,:,:,fnum,ch) = permute(fixPhase(permute(IMfastY(:,:,:,fnum, ch), [2 1 3]), shiftY{blockIDs(fnum)}), [2 1 3]);
        IMfastX(:,:,:,fnum,ch) = fixPhase(IMfastX(:,:,:,fnum,ch), shiftX{blockIDs(fnum)});
    end
end
toc;

if opts.doMotion
    tic;
    disp('Aligning Motion Y...')
    alignOpts.noiseThresh = opts.noiseThresh;
    IMfastY = alignMotion(IMfastY, alignChan, alignOpts);
    disp('Aligning Motion X...')
    IMfastX = alignMotion(IMfastX, alignChan, alignOpts);
    toc;
else
    %no motion correction
end

%Combine images
indCut = min(floor(length(fns)/4), ceil(length(fns)/20));
if ~indCut
    medFX = squeeze(nanmean(IMfastX,4));
    medFY = squeeze(nanmean(IMfastY,4));
else %trimmed mean
    IMfastX_s = sort(IMfastX,4);
    medFX = squeeze(nanmean(IMfastX_s(:,:,:,1+indCut:end-indCut,:),4));
    clear IMfastX_s
    
    IMfastY_s = sort(IMfastY,4);
    medFY = squeeze(nanmean(IMfastY_s(:,:,:,1+indCut:end-indCut,:),4));
    clear IMfastY_s
end
medFX(isnan(medFX(:))) = nanmedian(medFX(:)); %get rid of nans for easier registration, imtranslate
medFY(isnan(medFY(:))) = nanmedian(medFY(:)); %get rid of nans for easier registration, imtranslate

%realign each individual image to the reference for that direction, by strips
if opts.doMotion
    disp('Aligning fast motion...')
    res = [size(IMfastX,1) size(IMfastX,2)];
    [C,R] = meshgrid(1:res(1),1:res(2));
%     IMfastX2 = nan(size(IMfastX));
%     IMfastY2 = nan(size(IMfastY));
    for Z = 1:size(IMfastX,3)
        disp(['Z=' int2str(Z) ' of ' int2str(size(IMfastX,3))])
        [RSx,CSx] = getTranslationByStrips(squeeze(IMfastX(:,:,Z,:,alignChan)), medFX(:,:,Z,alignChan));
        [RSy, CSy] = getTranslationByStrips(permute(squeeze(IMfastY(:,:,Z,:,alignChan)),[2,1,3]), medFY(:,:,Z,alignChan)');
        %CSx = CSx-mean(CSx,2); CSy = CSy-mean(CSy,2);
        for ch = 1:size(IMfastX,5)
            for imN = 1:size(IMfastX,4)
                CC = C + repmat(CSx(:,imN),1, res(2));
                RR = R + repmat(RSx(:,imN), 1, res(2));
                IMfastX(:,:,Z,imN,ch) = reshape(interp2(IMfastX(:,:,Z,imN,ch), CC(:), RR(:)), res(1), res(2));
                CC = C + repmat(RSy(:,imN)',res(1), 1);
                RR = R + repmat(CSy(:,imN)', res(1), 1);
                IMfastY(:,:,Z,imN,ch) = reshape(interp2(IMfastY(:,:,Z,imN,ch), CC(:), RR(:)), res(1), res(2));
            end
        end
    end
%Combine images again
if ~indCut
    medFX = squeeze(nanmean(IMfastX,4));
    medFY = squeeze(nanmean(IMfastY,4));
else %trimmed mean
    IMfastX_s = sort(IMfastX,4);
    medFX = squeeze(nanmean(IMfastX_s(:,:,:,1+indCut:end-indCut,:),4));
    clear IMfastX_s
    
    IMfastY_s = sort(IMfastY,4);
    medFY = squeeze(nanmean(IMfastY_s(:,:,:,1+indCut:end-indCut,:),4));
    clear IMfastY_s
end
medFX(isnan(medFX(:))) = nanmedian(medFX(:)); %get rid of nans for easier registration, imtranslate
medFY(isnan(medFY(:))) = nanmedian(medFY(:)); %get rid of nans for easier registration, imtranslate
end
    
%register median X and median Y images to each other
prior.strength = 0.5;
prior.width = 30;
prior.sharpness = 50;
for Z = 1:size(medFX,3)
    output = dftregistration_wprior(fft2(imgaussfilt(medFX(:,:,Z,alignChan),2)), fft2(imgaussfilt(medFY(:,:,Z, alignChan),2)), 1, prior);
    rowshift = output(3);
    colshift = output(4);
    for ch = 1:numChannels
        medFX(:,:,Z,ch) = imtranslate(medFX(:,:,Z,ch), -[colshift/2 rowshift/2]);
        medFY(:,:,Z,ch) = imtranslate(medFY(:,:,Z,ch), [colshift/2 rowshift/2]);
    end
end

%fix warp
[warpX, warpY] = findWarp(medFX(:,:,:,alignChan),medFY(:,:,:,alignChan));
for ch = 1:numChannels
    medFX(:,:,:,ch) = fixWarp(medFX(:,:,:,ch),warpX);
    medFY(:,:,:,ch) = permute(fixWarp(permute(medFY(:,:,:,ch), [2 1 3]), warpY), [2 1 3]);
end
if opts.verbose
    figure('Name', 'Difference Image'), imshow3D(medFX(:,:,:,alignChan)-medFY(:,:,:,alignChan),[])
end

%assemble rigid image
rigidImage = nan([size(vol,1), size(vol,2)/2, length(M.coords.Z),numChannels]);
for Z = 1:size(medFX,3)
    for ch = 1:numChannels
        %correct low frequency noise
        correction = zeros(size(medFX,1),1);
        for row = 1:size(medFX,1)
            %valid = ~isnan(medFX(row,:,Z,ch)+concensusY(row,:,Z,ch));
            correction(row) = trimmean(medFX(row,:,Z,ch) - medFY(row,:,Z,ch), 20);
            %COULD DO CORRECTION VIA HISTOGRAM MATCHING;
            %i.e. call hist, find maximum, fit centroid, subtract differences.
        end
        concensusX = medFX(:,:,Z,ch) - repmat(correction, 1, size(medFX,2));
        
        correction = zeros(1, size(medFY,2));
        for col = 1:size(medFY,2)
            %valid = ~isnan(concensusX(:,col)+concensusY(:,col));
            correction(col) = trimmean(medFY(:,col,Z,ch) - medFX(:,col,Z,ch), 20);
            %COULD DO CORRECTION VIA HISTOGRAM MATCHING;
            %i.e. call hist, find maximum, fit centroid, subtract differences.
        end
        concensusY = medFY(:,:,Z,ch) - repmat(correction, size(medFY,1),1);
        
        rigidImage(:,:,Z,ch) = (concensusX+concensusY)/2; %arithmetic mean  - geometric mean: %concensus = sqrt(concensusX.*concensusY);
    end
end

% %make the black level of each plane the same?; i.e. match the XXth
% %percentile
%no longer valid since we take reference images without a mask!
% for Z = 1:size(rigidImage,4)
%     for ch = 1:size(rigidImage,3)
%         plane = rigidImage(:,:,Z,ch);
%         rigidImage(:,:,Z,ch) = plane - prctile(plane(:), 10);
%     end
% end
for ch = 1:numChannels
    rigidImage(:,:,:,ch) = max(0, rigidImage(:,:,:,ch)- prctile(reshape(rigidImage(:,:,:,ch), 1,[]), opts.blackLevel));  %set zero level
end

refIM.data = uint16(rigidImage.*((2^16-1)/max(rigidImage(:))));
refIM.metadata = metadata;
refIM.M = M;
refIM.M.pixelSizeV = diff(refIM.M.coords.X(1:2));
refIM.IM = double(refIM.data(:,:,:,alignChan)); %a black and white double precision image useful for quick operations
refIM.alignChan = alignChan; %which of the channels in refIM is used for alignment
refIM = alignPlanes(refIM);

%SAVE
figure, imshow3D(refIM.data(:,:,:,alignChan)); drawnow;
disp('Saving...')
save([savedr savefn(1:end-4)], 'refIM', '-v6') %mat file
options.overwrite = true;
for ch = 1:numChannels
    %save image
    errorcode = saveastiff(uint8(opts.satMultiplier*(rigidImage(:,:,:,ch).^opts.gamma)), [savedr savefn(1:end-4) '_Ch' int2str(ch) '.tif'], options); %tiff file for easy reading
    if errorcode
        keyboard
    end
end

% if strcmpi(input('Save these shifts for future image registrations? Y|[N] >>', 's'), 'y')
%     disp('Saving...')
%     hashdata = {metadata.rois.RoiGroups.imagingRoiGroup.rois,  metadata.hRoiManager.scanFrameRate metadata.hScan2D.linePhase}; %#ok<NASGU>
%     save( [basedir '\GalvoTraces\RasterPhase\' roiHash '.mat'], 'shiftX', 'shiftY', 'warpX', 'warpY', 'hashdata');
% end

end

function yy = findPhase(IMset, noiseThresh)
edgecut = 15;
blocksize = 5;

shifts = -2:0.5:2;
%shifts = -4:0.5:4;
sorted1 = nth_element(IMset(:), ceil(numel(IMset)*0.5));
sorted2 = nth_element(IMset(:), ceil(numel(IMset)*0.01));
noiseamp = sorted1(ceil(numel(IMset)*0.5)) - sorted2(ceil(numel(IMset)*0.01)) +40; %do this better
%noiseamp = nanmedian(IMset(:))-prctile(IMset(:), 1) + 40; %3*sqrt(estimatenoise(reshape(IMset(:,:,1), 1,[]))); %estimate of digitizer noise; this should be less than the signal of a single photon
noiseamp = noiseThresh;

H2 = floor((size(IMset,1)-1)/2);
cols = edgecut+max(shifts)+1:blocksize:size(IMset,2)-edgecut + min(shifts);
rowshift = nan(length(cols), size(IMset,3));
GOF = nan(size(rowshift));

for frameix = 1:size(IMset,3)
    IM1 = medfilt2(IMset(1:2:2*H2,:,frameix),[3 3]);
    IM1 = max(0, min(IM1, prctile(IM1(:), 99.995))-noiseamp);
    IM2 = medfilt2(IMset(2:2:2*H2,:,frameix),[3 3]);
    IM2 = max(0, min(IM2, prctile(IM2(:), 99.995))-noiseamp);
    IM3 = [IM1(2:end,:); zeros(1,size(IM1,2))];
    try
        CC = normxcorr2_general(IM1,IM2, 0.9*numel(IM1));
    catch ME
       disp(['Frameix ' int2str(frameix) ': image was too dim to perform phase correction and dewarping']);
       continue
    end
    [~,maxix] = max(CC(:)); [~,j] = ind2sub(size(CC), maxix);
    meanShift = j-size(IM1,2);
    
    %show an interlaced image of the corrected mean shift
    
    %IM1int = griddedInterpolant(IM1);
    IM2int = griddedInterpolant(IM2);
    for col_ix = 1:length(cols)
        col = cols(col_ix);
        E1 = nan(1, length(shifts));
        E3 = E1;
        for shiftix = 1:length(shifts)
            IM2col = IM2int({1:size(IM2,1), col+meanShift-shifts(shiftix) + (-blocksize:blocksize)});
            E1(shiftix) = sum(sum(abs(IM1(:,col + (-blocksize:blocksize))-IM2col)));
            E3(shiftix) = sum(sum(abs(IM3(:,col+ (-blocksize:blocksize))-IM2col)));
            %better distance metric?
        end
        if sum(~isnan(E1))>3 && sum(~isnan(E3))>3
            
            [minE, minIX] = min(E1+E3);
            rowshift(col_ix, frameix) = -meanShift + shifts(minIX);
            GOF(col_ix,frameix) = mean(E1+E3)-minE;
            
%             [minE1, minIX1] = min(E1);
%             [minE3, minIX3] = min(E3);
%             rowshift(col_ix, frameix) = -meanShift + (shifts(minIX1) + shifts(minIX3))/2;
%             GOF(col_ix,frameix) = (mean(E1)-minE1) + (mean(E3)-minE3);
        end
    end
end

%WEIGHT EACH ROWSHIFT MEASUREMENT BY THE QUALITY OF FIT
weight = nansum(GOF,2);
rowshift = nansum((rowshift+meanShift).*GOF,2)./(weight + nanmax(weight)/50) -meanShift;
meanShift = -(nansum(rowshift.*weight)./nansum(weight));
rowshift(isnan(rowshift)) = -meanShift;
rowshift = medfilt2(rowshift, [3 1]); %filter out brief inconsistent spikes
cs = csaps([1 cols size(IMset,2)],[-meanShift rowshift' -meanShift],1e-8,[], [nanmax(weight)/10; weight ;nanmax(weight)/10]);
yy = fnval(cs, 1:size(IMset,2));
% figure('name', 'Phase Correction'), scatter(cols, rowshift')
% hold on, plot(yy)
end

function [IMcorrected, yy] = fixPhase(IMset, yy)
%apply correction
IMcorrected = IMset;
grid = 1:size(IMset,2);
for ix = 1:size(IMset,3)
    IMcorrected(2:2:end,:,ix) = interp1(grid', IMset(2:2:end,:,ix)', grid-yy)';  %IM2int({IM2int.GridVectors{1}, grid-yy});
end
end

function [xx, yy] = findWarp(imX,imY)
edgecut = 15;
shifts = -3:0.5:3;
nZ = size(imX,3);

sorted1 = nth_element(imX(:), ceil(numel(imX)*0.5));
sorted2 = nth_element(imX(:), ceil(numel(imX)*0.01));
noiseamp = sorted1(ceil(numel(imX)*0.5)) - sorted2(ceil(numel(imX)*0.01)) +15;
%noiseamp = nanmedian(imX(:))-prctile(imX(:),1)+15; %sqrt(estimatenoise(reshape(IMset(:,:,1), 1,[]))); %estimate of digitizer noise; this should be less than the signal of a single photon

warpX = nan(nZ,size(imY,2));
GOFx = warpX;
warpY = nan(nZ,size(imX,1));
GOFy = warpY;
for Z = 1:nZ
    medX = max(0, imX(:,:,Z)-noiseamp);
    medY = max(0, imY(:,:,Z)-noiseamp);

    Xint = griddedInterpolant(medX(edgecut+1:end-edgecut, :));
    for col = edgecut+max(shifts)+1:size(medY,2)-edgecut + min(shifts)
        E = nan(1, length(shifts));
        Ycol = medY(edgecut+1:end-edgecut,col);
        for shiftix = 1:length(shifts)
            Xcol =  Xint({1:length(Ycol),col+shifts(shiftix)});
            E(shiftix) =  -sum(abs(Ycol-Xcol));%corr(Ycol, Xcol);
        end
        [maxE, maxIX] = max(E);
        warpX(Z,col) = shifts(maxIX);
        GOFx(Z,col) = maxE-mean(E); 
    end
    
    Yint = griddedInterpolant(medY(:,edgecut+1:end-edgecut));
    for row = edgecut+max(shifts)+1:size(medX,1)-edgecut + min(shifts)
        E = nan(1, length(shifts));
        Xrow = medX(row,edgecut+1:end-edgecut)';
        for shiftix = 1:length(shifts)
            Yrow =  Yint({row+shifts(shiftix),1:length(Xrow)})'; %interp1(medX(edgecut+1:end-edgecut,:)', col-shifts(shiftix))';
            E(shiftix) =  -sum(abs(Yrow-Xrow)); %corr(Yrow, Xrow);
        end
        [maxE, maxIX] = max(E);
        warpY(Z,row) = shifts(maxIX);
        GOFy(Z,row) = maxE-mean(E); 
    end
end

warpX = nansum(warpX.*GOFx,1)./nansum(GOFx,1);
weight = nansum(GOFx,1);
warpX(isnan(warpX)) = 0;
weight([1 end]) = nanmax(weight)/100;
csX = csaps(1:length(warpX),warpX',1e-5,[], weight);
xx = fnval(csX, 1:length(warpX));
% figure('name', 'X Warp Correction'), scatter(1:length(warpX), warpX)
% hold on, plot(xx)

warpY = nansum(warpY.*GOFy,1)./nansum(GOFy,1);
weight = nansum(GOFy,1);
warpY(isnan(warpY)) = 0;
weight([1 end]) = nanmax(weight)/100;
csY = csaps(1:length(warpY),warpY',1e-5,[], weight);
yy = fnval(csY, 1:length(warpY));
% figure('name', 'Y Warp Correction'), scatter(1:length(warpY), warpY)
% hold on, plot(yy)
end

function IMout = fixWarp(IM, xx)
IMout = nan(size(IM));
for ix = 1:size(IM,3)
    IMout(:,:,ix) = interp1(IM(:,:,ix)', (1:size(IM,2))+xx)';
end
end

function  SI = parseMetaData(Mstring)
SI = [];
Mstring = strsplit(Mstring, '\n');
try
    for i = 1:length(Mstring)
        eval([Mstring{i} ';']);
    end
catch
    i = i-1;
end
SI.rois = loadjson(strjoin(Mstring(i+1:end), '\n'));
end