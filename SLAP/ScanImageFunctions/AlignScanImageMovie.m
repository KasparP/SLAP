function dataset = AlignScanImageMovie
%default options:
    opts.motionChan = 1;
    opts.activityChan = [1 2];  tooltips.activityChan = 'activity channels, e.g. [1 2] or [1]';
    opts.threshDense = true; tooltips.threshDense = 'Set intensity thresholds for dense imaging?';
    opts = optionsGUI(opts, tooltips);
    
    opts.forceStage = false(1, 100);
if ~nargin
    disp('Starting Stage 1:')
    %get a list of movies
    [fns, dr] = uigetfile('*.tif', 'Select SI Activity movies', 'MultiSelect', 'On');
    if ~iscell(fns)
        fns = {fns};
    end
    drsave = dr;
    fnsave = [fns{1}(1:end-16) '_rasterDataset'];
    
    phaseChan = opts.motionChan; %channel to use for phase correction
    
    %for each movie, perform bidi phase correction
    for fnum = 1:length(fns)
        disp(['Loading and Phase-correcting file: ' int2str(fnum) ' of ' int2str(length(fns))])
        
        reader = ScanImageTiffReader([dr filesep fns{fnum}]);
        IM = permute(double(reader.data()), [2 1 3]);

        if fnum==1
            dataset.fns = fns;
            dataset.metadata = parseMetaData(reader.metadata());
            dataset.nChan= length(dataset.metadata.hChannels.channelSave);
            dataset.IM = nan(size(IM,1), size(IM,2), size(IM,3)/dataset.nChan, dataset.nChan, length(fns), 'single');
        end
        
        phaseOffset = findPhase(IM(:,:,phaseChan:dataset.nChan:end)); %find phase with channel 1
        for ch = 1:dataset.nChan
            dataset.IM(:,:,:,ch, fnum) = fixPhase(IM(:,:,ch:dataset.nChan:end), phaseOffset);
        end
    end
    
    %save data
    dataset.stage = 1;
    dataset.filename = [drsave filesep fnsave];
    
    save(dataset.filename, 'dataset', '-v7.3') 
end

if dataset.stage<2 || opts.forceStage(2)

    
    disp('Starting Stage 2: Alignment')
    alignOpts.window = 30;  %maximum movement;
    alignOpts.nRefs = 3; %number of reference values to produce; we will use the most common one
    alignOpts.nFramesRef = 300; %create a reference from first nFramesRef frames
    alignOpts.doplot = false;
    alignOpts.prior.strength = 0.3;
    alignOpts.prior.width = 30;
    alignOpts.prior.sharpness = 50;
    alignOpts.chunkSize = 30;
    %align
    disp('aligning...')
    if all(opts.activityChan==opts.motionChan)
        [alignedA, ~,~] = align2D_notemplate(reshape(dataset.IM(:,:,:,opts.motionChan,:), size(dataset.IM,1), size(dataset.IM, 2), []), alignOpts); %ALIGN
    else %align other channels
        [~, shifts, ~] = align2D_notemplate(reshape(dataset.IM(:,:,:,opts.motionChan,:), size(dataset.IM,1), size(dataset.IM, 2), []), alignOpts); %ALIGN
        alignedA = zeros(size(dataset.IM,1), size(dataset.IM,2), size(shifts,2),size(dataset.IM,4), 'single');
        for ch = 1:length(opts.activityChan)
            imageSet= reshape(dataset.IM(:,:,:,opts.activityChan(ch),:), size(dataset.IM,1), size(dataset.IM, 2), []);
            imageSet(isnan(imageSet)) = nanmean(imageSet(randi(numel(imageSet),1, min(numel(imageSet), 1e6))));
            for frame = 1:size(imageSet,3)
                alignedA(:,:,frame,ch) = imtranslate(imageSet(:,:,frame), shifts(:,frame)', 'FillValues', prctile(reshape(imageSet(:,:,frame),1,[]), 5));
            end
            clear imageSet
        end
    end
    dataset = rmfield(dataset, 'IM');
    
    %correct slow noise
    disp('correcting slow noise...')
    
    for ch = 1:size(alignedA,4)
        smoothA = mean(alignedA(:,:,1:min(300,end),ch),3);
        threshA = prctile(smoothA(:), 50)+2;
        dimRef  = imerode(imdilate(smoothA<threshA, ones(3)), ones(7));
        for f = 1:size(alignedA,3)  %for each frame
            if mod(f,300)==151 && f<(size(alignedA,3)-150)
                smoothA = mean(alignedA(:,:,f-150:min(f+150,end), ch),3);
                threshA = prctile(smoothA(:), 50)+2;
                dimRef  = imerode(imdilate(smoothA<threshA, ones(3)), ones(7));
            end
            
            %correct row noise from slow digitizer fluctuations
            for row = 1:size(smoothA,1)
                select = dimRef(row,:);
                if sum(select)>30
                    alignedA(row,:,f,ch) = alignedA(row,:,f,ch) - (mean(alignedA(row,select,f,ch) - smoothA(row,select)));
                else
                    disp('too few pixels for correction')
                end
            end
        end
    end
    
    %register by strips
    dostrips = true;
    if dostrips
        inds = round(linspace(1,0.99*size(alignedA,3),256));
        ref = trimmean(alignedA(:,:,:,1),60,3);
        noiselevel = prctile(ref(:), 50)+5;
        tmp = registerByStrips(alignedA(:,:,inds,1),ref,noiselevel);
        ref = trimmean(tmp,60,3);
        alignedA = registerByStrips(alignedA,ref,noiselevel);
    end
    %metric1 = ones(1,size(smoothA,3)); metric2 = ones(1,size(smoothA,3)); 
    metric3 = ones(1,size(smoothA,3));
    dataset.ref = ref;
    refHI = dataset.ref(:)>(noiselevel*1.5);
    %ptile = 100*sum(refHI(:))./numel(refHI);
    for f = 1:size(alignedA,3)
        maskF = imgaussfilt(alignedA(:,:,f,1), 0.6);
        %metric1(f) = -corr(maskF(:), dataset.ref(:));
        %metric2(f) = -corr(maskF(:)>prctile(maskF(:), ptile), refHI);
        metric3(f) = -corr(maskF(:)>(noiselevel*1.5), refHI);
    end
    motionMetric = metric3;
    
    dataset.motionMetric = motionMetric;
    dataset.aligned = reshape(alignedA, size(alignedA,1),size(alignedA,2), [], length(dataset.fns), size(alignedA,4));
    dataset.stage = 2;
    disp(['Saving dataset to: ' fnsave]);
    save([drsave filesep fnsave], 'dataset', '-v7.3') 
    disp('done save')
end
end

function yy = findPhase(IMset)
edgecut = 15;
blocksize = 5;
shifts = -2:0.5:2;
sorted1 = nth_element(IMset(:), ceil(numel(IMset)*0.5));
sorted2 = nth_element(IMset(:), ceil(numel(IMset)*0.01));
noiselevel = sorted1(ceil(numel(IMset)*0.5)) + (sorted1(ceil(numel(IMset)*0.05)) - sorted2(ceil(numel(IMset)*0.01)));
%noiseamp = nanmedian(IMset(:))-prctile(IMset(:), 1) + 40; %3*sqrt(estimatenoise(reshape(IMset(:,:,1), 1,[]))); %estimate of digitizer noise; this should be less than the signal of a single photon

H2 = floor((size(IMset,1)-1)/2);
cols = edgecut+max(shifts)+1:blocksize:size(IMset,2)-edgecut + min(shifts);
rowshift = nan(length(cols), size(IMset,3));
GOF = nan(size(rowshift));

for frameix = 1:size(IMset,3)
    IM1 = medfilt2(IMset(1:2:2*H2,:,frameix),[3 3]);
    IM1 = max(0, IM1-noiselevel);
    IM2 = medfilt2(IMset(2:2:2*H2,:,frameix),[3 3]);
    IM2 = max(0, min(IM2, max(IM2(:))-noiselevel));
    IM3 = [IM1(2:end,:); zeros(1,size(IM1,2))];
    
    CC = normxcorr2_general(IM1,IM2, 0.9*numel(IM1));
    [~,maxix] = max(CC(:)); [~,j] = ind2sub(size(CC), maxix);
    meanShift = j-size(IM1,2);
    
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
            [minE1, minIX1] = min(E1);
            [minE3, minIX3] = min(E3);
            rowshift(col_ix, frameix) = -meanShift + (shifts(minIX1) + shifts(minIX3))/2;
            GOF(col_ix,frameix) = (mean(E1)-minE1) + (mean(E3)-minE3);
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

function imX = registerByStrips(imX,ref, noiseamp)
smoothing = 8.3;
stripheight = 16;
ref = ref-min(ref(:));
%refFilt = ref - imgaussfilt(ref, stripheight/2);
edgecut = stripheight;
medX = max(0, imX(:,:,:,1)-noiseamp);
rowCs = stripheight+1:stripheight/2:size(imX,1)-stripheight;
outsize = [2*stripheight+1 size(ref,2)]+[stripheight+1 size(medX,2)-(2*edgecut)]-1;
for rowIX = length(rowCs):-1:1
    a = ref(rowCs(rowIX)+(-stripheight:stripheight),:);
    fftArot(:,:,rowIX) = fft2(rot90(a,2),outsize(1),outsize(2));
end
lsA = nan([outsize length(rowCs)]); lsA2 = lsA;
for T = 1:size(imX,3)
    T
    warpX = nan(1,size(imX,1));
    warpY = nan(1,size(imX,1));
    GOF = warpX;
    for rowIX = 1:length(rowCs)
        rowC = rowCs(rowIX);
        b = medX(rowC+(-stripheight/2:stripheight/2),edgecut+1:end-edgecut,T);
        if T==1
            [C , lsA(:,:,rowIX), lsA2(:,:,rowIX)] = normxcorr2_general_FFT(b, ref(rowC+(-stripheight:stripheight),:), fftArot(:,:,rowIX), [],[],numel(b));
            C = rot90(C,2);
        else
            C = rot90(normxcorr2_general_FFT(b, ref(rowC+(-stripheight:stripheight),:), fftArot(:,:,rowIX), lsA(:,:,rowIX), lsA2(:,:,rowIX),numel(b)),2);
        end
        %C = freqxcorr(fftArot(:,:,rowIX),b,outsize);
        %C2 = normxcorr2_general(b, ref(rowC+(-stripheight:stripheight),:), 100);
        %C = xcorr2(ref(rowC+(-stripheight:stripheight),:), medX(rowC+(-stripheight/2:stripheight/2),edgecut+1:end-edgecut));
        [~, op] = max(C(:)); %minC = min(min(C(ceil((size(C,1)-stripheight)/2):floor((size(C,1)+stripheight)/2), ceil((size(C,2)-stripheight)/2):floor((size(C,2)+stripheight)/2))));
        [opX,opY] = ind2sub(size(C), op);
        warpX(rowC) = -(opX-ceil(size(C,1)/2));
        warpY(rowC) = -(opY-ceil(size(C,2)/2));
        GOF(rowC) = sum(b(:)); %maxC-minC;
    end
    GOF = GOF-0.9*nanmin(GOF);
    discard = isnan(warpX) | isnan(warpX) | abs(warpX)>=(stripheight/2-1) | abs(warpY)>=(stripheight/2-1);
    warpX(discard) = 0;
    warpY(discard) = 0;
    GOF(discard) = 0;
    weight = GOF;
    weight([1 end]) = nanmax(weight)/100;
    
    csX = csaps(1:length(warpX),warpX',10^(-smoothing),[], weight);
    xx = fnval(csX, 1:length(warpX));
    csY = csaps(1:length(warpY),warpY',10^(-smoothing),[], weight);
    yy = fnval(csY, 1:length(warpY));
    
    D = cat(3, -repmat(yy', 1, size(imX,2)), -repmat(xx', 1, size(imX,2)));
    for ch = 1:size(imX,4)
        imX(:,:,T,ch) = imwarp(imX(:,:,T,ch),D);
    end
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
