function dataset = SI_calcTuning(dataset)

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
    end
    
    %save data
    dataset.stage = 1;
    save([drsave filesep fnsave], 'dataset', '-v7.3') 
end

if dataset.stage<2
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
    tic
    if activityChan==motionChan
    [alignedA, shifts, errors] = align2D_notemplate(reshape(dataset.IM(:,:,:,motionChan,:), size(dataset.IM,1), size(dataset.IM, 2), []), alignOpts); %ALIGN
    %align other channel
    else
        disp('Using Motion Channel for Alignment.')
        [~, shifts, ~] = align2D_notemplate(reshape(dataset.IM(:,:,:,motionChan,:), size(dataset.IM,1), size(dataset.IM, 2), []), alignOpts); %ALIGN
        imageSet= reshape(dataset.IM(:,:,:,activityChan,:), size(dataset.IM,1), size(dataset.IM, 2), []);
        imageSet(isnan(imageSet)) = nanmean(imageSet(randi(numel(imageSet),1, min(numel(imageSet), 1e7))));
        alignedA = zeros(size(imageSet), 'single');
        for frame = 1:size(imageSet,3)
            alignedA(:,:,frame,1) = imtranslate(imageSet(:,:,frame), shifts(:,frame)');
        end
        clear imageSet
    end
    toc
    dataset = rmfield(dataset, 'IM');
    
    %correct slow noise
    disp('correcting slow noise...')
    tic
    motionMetric = ones(1,size(alignedA,3));
    
    for f = 1:size(alignedA,3)  %for each frame
        if mod(f,150)==1
            smoothA = mean(alignedA(:,:,max(1,f-150):min(f+150,end)),3);
            threshA = prctile(smoothA(:), 85);
            dimRef  = imerode(imdilate(smoothA<threshA, ones(3)), ones(7));
        end
        
       %correct row noise from slow digitizer fluctuations
        for row = 1:size(smoothA,1)
            select = dimRef(row,:);
            if sum(select)>50
                alignedA(row,:,f) = alignedA(row,:,f) - (mean(alignedA(row,select,f) - smoothA(row,select)));
            end
        end
    end
    toc
    
    ref = mean(alignedA,3);
    
    dodemons = false;
    if dodemons
        tic
        disp('Demons Registration...')    
        for f = 1:size(alignedA,3)
            disp(['Nonrigid Align Frame:' int2str(f) ' of ' int2str(size(alignedA,3))])
            [~, alignedA(:,:,f)] = imregdemons(alignedA(:,:,f), ref, 20, 'AccumulatedFieldSmoothing',2, 'PyramidLevels', 1, 'DisplayWaitBar', false);
        end
        toc
    end
    
    dataset.ref = ref;
    ref = min(prctile(ref(:),99), max(ref, prctile(ref(:),1)));
    for f = 1:size(alignedA,3)
        motionMetric(f) = -corr(reshape(alignedA(:,:,f),[],1), reshape(ref,[],1));
    end
    

    dataset.motionMetric = motionMetric;
    dataset.aligned = reshape(alignedA, size(alignedA,1),size(alignedA,2), [], length(dataset.fns));
    dataset.stage = 2;
    for ch = 1:size(dataset.aligned,4)
        dataset.aligned(:,:,:,ch) = dataset.aligned(:,:,:,ch) - prctile(reshape(dataset.aligned(:,:,1:10:end,ch),1,[]), 25);
    end
    
    disp(['Saving dataset to: ' fnsave]);
    save([drsave filesep fnsave], 'dataset', '-v7.3') 
    disp('done save')
end

%plot a mean-subtracted average movie
mmov = nanmean(dataset.aligned,4);
mmov = mmov- nanmean(mmov,3);
figure, imshow3D(mmov);

%select ROIS
f = figure('Name', 'DRAW ROIS'); ax=axes; im=imagesc(dataset.ref);
while(isgraphics(f))
    h = imrect(ax);
    pos = wait(h);
    delete(h);
    if isempty(pos)
        break
    end
    
    selected = dataset.aligned(round(pos(2)):round(pos(2)+pos(4)), round(pos(1)):round(pos(1)+pos(3)), :,:);
    %plot the average over space
    total = squeeze(nansum(nansum(selected,1),2));
    total = total./nanmean(total,1);
    
    figure, plot(total), hold on, plot(mean(total,2), 'k', 'linewidth', 2)
    times = floor(5*(1:8) * 3.4);
    hold on, scatter(times, ones(size(times)), 'r', 'linewidth', 2)
end

% %censor motion and zero out baseline
% figure('name', 'Select a censor threshold for motion'), plot(imtophat(dataset.motionMetric, ones(1,251)))
% xlabel('frame'), ylabel('Motion Metric: Poorly aligned frames have large values')
% [~,y] = ginput;
% if isempty(y)
%     thresh = 0.033;
% else
%     thresh = y(end);
% end
% 
% censor = imtophat(dataset.motionMetric, ones(1,251))>thresh;
% censor = censor | medfilt1(double(censor), 9);
% sz = size(dataset.aligned);
% dataset.aligned = reshape(dataset.aligned, sz(1),sz(2),sz(3)*sz(4)); 
% dataset.aligned(:,:,censor) = nan; %nan out motion frames
% 
% stdIM = nanvar(dataset.aligned,1, 3)./(nanmean(dataset.aligned, 3)-min(dataset.aligned(:)));



%show a mean-subtracted average movie



keyboard
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