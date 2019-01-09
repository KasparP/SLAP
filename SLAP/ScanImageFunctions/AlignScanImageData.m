function dataset = AlignScanImageData(dataset,optsin)
opts.forceStage = false(1, 100);
if ~nargin
    disp('Starting Stage 1:')
    %get a list of movies
    [fns, dr] = uigetfile('*.tif', 'Select SI Activity movies', 'MultiSelect', 'On');
    drsave = dr;
    fnsave = [fns{1}(1:end-16) 'dataset'];
    
    phaseChan = 1; %channel to use for phase correction
    
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
        
        %stimulus timing
        desc = reader.descriptions();
        ix = strfind(desc{1},  'epoch = ');
        eval([desc{1}(ix:end-1) ';']);
        dataset.stimTime(fnum) = datenum(epoch);
        for ix = 1:length(desc)/2
            loc = strfind(desc{2*ix},  'frameTimestamps_sec =');
            dataset.frameTimes(ix) = str2double(desc{2*ix}(loc+22:loc+29));
        end
    end

    %STIMULUS
    stimfiles = dir([dr  filesep '..' filesep '*Timings*']);
    if ~isempty(stimfiles)
        stimList = [];
        for sf = 1:length(stimfiles)
            SF = load([dr filesep '..' filesep stimfiles(sf).name]);
            stimList = [stimList ; SF.time_stamp]; %#ok<AGROW>
        end
    else
        error('no timing fiiles in folder');
    end
    dist = repmat(stimList(:,1), 1, length(dataset.stimTime)) - repmat(dataset.stimTime, size(stimList,1),1);
    distP = dist; distM = dist;
    distP(distP<0) = nan; distM(distM>0) = nan;
    mindistP = nanmin(distP, [], 1);
    mindistM = nanmax(distM,[],1);
    if var(mindistM)<var(mindistP)
        dataset.stimulus.stimDelay = nanmedian(mindistM);
        dataset.stimulus.delayVariance = nanvar(mindistM);
    else
        dataset.stimulus.stimDelay = nanmedian(mindistP);
        dataset.stimulus.delayVariance = nanvar(mindistP);
    end
    dataset.stimulus.stim = nan(1,length(fns));
    for fnum = 1:length(fns)
        [dataset.stimulus.timeError(fnum), minIx] = min(abs(dataset.stimTime(fnum) - stimList(:,1) + dataset.stimulus.stimDelay));
        dataset.stimulus.stim(fnum, 1:8) = stimList(minIx:minIx+7,2);
        dataset.stimulus.stimTime(fnum, 1:8) = stimList(minIx:minIx+7,1);
    end
    
    %save data
    dataset.stage = 1;
    save([drsave filesep fnsave], 'dataset', '-v7.3') 
end

if dataset.stage<2 || opts.forceStage(2)
    motionChan = 2; %use red channel for alignment; little/no activity in this channel?
    activityChan = 1;
    
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
    if activityChan==motionChan
    [alignedA, shifts, errors] = align2D_notemplate(reshape(dataset.IM(:,:,:,motionChan,:), size(dataset.IM,1), size(dataset.IM, 2), []), alignOpts); %ALIGN
    %align other channel
    else
        [~, shifts, ~] = align2D_notemplate(reshape(dataset.IM(:,:,:,motionChan,:), size(dataset.IM,1), size(dataset.IM, 2), []), alignOpts); %ALIGN
        imageSet= reshape(dataset.IM(:,:,:,activityChan,:), size(dataset.IM,1), size(dataset.IM, 2), []);
        imageSet(isnan(imageSet)) = nanmean(imageSet(randi(numel(imageSet),1, min(numel(imageSet), 1e7))));
        alignedA = zeros(size(imageSet), 'single');
        for frame = 1:size(imageSet,3)
            alignedA(:,:,frame,1) = imtranslate(imageSet(:,:,frame), shifts(:,frame)');
        end
        clear imageSet
    end
    dataset = rmfield(dataset, 'IM');
    
    %correct slow noise
    disp('correcting slow noise...')
    %smoothA = reshape(imgaussfilt(reshape(alignedA, size(alignedA,1)*size(alignedA,2), size(alignedA,3)),[eps 50]), size(alignedA));
    
    
    smoothA = mean(alignedA(:,:,1:min(300,end)),3);
%    smoothM = smoothA;
%     if motionChan~=activityChan
%         %smoothM =reshape(smooth2KP(reshape(alignedM, size(alignedA,1)*size(alignedA,2), size(alignedA,3)),50), size(alignedA));
%         smoothM = mean(alignedM(:,:,1:min(300,end),3));
%     end
    threshA = prctile(smoothA(:), 85); %threshM = prctile(smoothM(:), 85); 
    motionMetric = ones(1,size(smoothA,3));
    for f = 1:size(alignedA,3)  %for each frame
        if mod(f,300)==151
            smoothA = mean(alignedA(:,:,f-150:min(f+150,end)),3);
%            smoothM = smoothA;
%             if motionChan~=activityChan
%                 smoothM = mean(alignedM(:,:,f-150:min(f+150,end)),3);
%             end
        end
        dimRef  = imerode(imdilate(smoothA(:,:,f)<threshA, ones(3)), ones(7));
       %correct row noise from slow digitizer fluctuations
        for row = 1:size(smoothA,1)
            select = dimRef(row,:);
            if sum(select)>50
                alignedA(row,:,f) = alignedA(row,:,f) - (mean(alignedA(row,select,f) - smoothA(row,select)));
%                 if activityChan~=motionChan
%                     alignedM(row,:,f) = alignedM(row,:,f) - (mean(alignedM(row,select,f) - smoothM(row,select)));
%                 end
            end
        end
    end
    
    dodemons = true;
    if dodemons
        ref = mean(alignedA,3);
        for f = 1:size(alignedA,3)
            disp(['Nonrigid Align Frame:' int2str(f) ' of ' int2str(size(alignedA,3))])
            [~, alignedA(:,:,f)] = imregdemons(alignedA(:,:,f), ref, 20, 'AccumulatedFieldSmoothing',2, 'PyramidLevels', 1, 'DisplayWaitBar', false);
        end
    end
    
    dataset.ref = ref;
    ref = min(prctile(ref(:),99), max(ref, prctile(ref(:),1)));
    for f = 1:size(alignedA,3)
        motionMetric(f) = -corr(reshape(alignedA(:,:,f),[],1), reshape(ref,[],1));
    end
    
    dataset.motionMetric = motionMetric;
    dataset.aligned = reshape(alignedA, size(alignedA,1),size(alignedA,2), [], length(dataset.fns));
    dataset.stage = 2;
    disp(['Saving dataset to: ' fnsave]);
    save([drsave filesep fnsave], 'dataset', '-v7.3') 
    disp('done save')
end

%censor motion and zero out baseline
keyboard
censor = imtophat(dataset.motionMetric, ones(1,251))>0.033;
censor = censor | medfilt1(double(censor), 9);
sz = size(dataset.aligned);
dataset.aligned = reshape(dataset.aligned, sz(1),sz(2),sz(3)*sz(4)); 
dataset.aligned(:,:,censor) = nan; %nan out motion frames
%ref = nanmean(dataset.aligned,3);
% dataset.aligned= dataset.aligned - prctile(ref(:),0.1);
% ref = max(0, ref- prctile(ref(:),0.1));
dataset.aligned = reshape(dataset.aligned, sz);


%rearrange movies into stimuli
for mov = 1:size(dataset.aligned,4)
    
    
end

respAvg = nanmean(dataset.aligned, 4);

ref = min(respAvg,[], 3); 
ref = max(0, ref-prctile(ref(:),0.1));

refD = imdilate(ref, strel('disk', 1,0));

lambda = 10*prctile(ref(:),20);
respAvg = (respAvg - nanmin(respAvg,[],3))./sqrt(ref+lambda);
figure('name', 'Average Response'); imshow3D(respAvg);

%average response for each stimulus
stimuli = unique(dataset.stimulus.stim);
resp = nan(size(dataset.aligned,1), size(dataset.aligned,2), size(dataset.aligned,3), length(stimuli));
for s = 1:length(stimuli)
    select = dataset.stimulus.stim==stimuli(s);
    resp(:,:,:,s) = nanmean(dataset.aligned(:,:,:,select), 4);
    resp(:,:,:,s) = (resp(:,:,:,s) - min(resp(:,:,:,s),[],3))./sqrt(refD+lambda);
end

%Make HSV plot
rad = linspace(0, 2*pi, length(stimuli)+1);
rad = reshape(rad(1:end-1), 1,1,1,[]);
rad = repmat(rad, size(resp,1), size(resp,2),size(resp,3));
%hue is the preferred direction
H = max(0,min(1, (pi + circ_mean(rad,resp,4))/(2*pi)));

%saturation is the degree of preference
pref = circ_r(rad,resp,pi/4,4);
S = pref.*sum(resp,4);
S = min(1, S./(0.8*max(S(:))));
SD = imdilate(S,strel('disk', 1,0)); %dilated

V = min(1, respAvg./(0.8*max(respAvg(:))));

respRGB =  reshape(hsv2rgb([H(:) SD(:) V(:)]), [size(H) 3]); 
figure, imshow3D(respRGB);


colors = hsv(8);
resp = reshape(resp, size(resp,1)*size(resp,2), size(resp,3), size(resp,4));
abort = false;
while ~abort
    figure('name', 'Select a region of interest'), imshow(sqrt(max(0,dataset.ref)),[])
    bw = roipoly;
    figure,
    for i = 1:8
        hold on,
        plot(squeeze(mean(resp(bw(:), :,i),1)), 'color', colors(i,:))
    end
    abort = strcmpi(input('Type q to end>>', 's'), 'q');
end 
end
function yy = findPhase(IMset)
edgecut = 15;
blocksize = 5;
shifts = -2:0.5:2;
sorted1 = nth_element(IMset(:), ceil(numel(IMset)*0.5));
sorted2 = nth_element(IMset(:), ceil(numel(IMset)*0.01));
noiseamp = sorted1(ceil(numel(IMset)*0.5)) - sorted2(ceil(numel(IMset)*0.01)) + 40;
%noiseamp = nanmedian(IMset(:))-prctile(IMset(:), 1) + 40; %3*sqrt(estimatenoise(reshape(IMset(:,:,1), 1,[]))); %estimate of digitizer noise; this should be less than the signal of a single photon

H2 = floor((size(IMset,1)-1)/2);
cols = edgecut+max(shifts)+1:blocksize:size(IMset,2)-edgecut + min(shifts);
rowshift = nan(length(cols), size(IMset,3));
GOF = nan(size(rowshift));

for frameix = 1:size(IMset,3)
    IM1 = medfilt2(IMset(1:2:2*H2,:,frameix),[3 3]);
    IM1 = max(0, min(IM1, prctile(IM1(:), 99.99))-noiseamp);
    IM2 = medfilt2(IMset(2:2:2*H2,:,frameix),[3 3]);
    IM2 = max(0, min(IM2, prctile(IM2(:), 99.99))-noiseamp);
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
rowshift = nansum(rowshift.*GOF,2)./(weight + nanmax(weight)/50);
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