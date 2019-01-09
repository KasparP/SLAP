function dataset = analyze_SLAPMiTuning

%select files
[fns, dr] = uigetfile('*.mat', 'Select your saved RECON files for analysis', 'multiselect', 'on');
if ~iscell(fns)
    fns = {fns};
end
fns = sort_nat(fns);


[fnRef, drRef] = uigetfile([dr filesep '*.mat'], 'select your reference Segmentation _SEG file');

%LOAD STIMULUS
stimfiles = dir([dr  filesep '*Timings*']);
if ~isempty(stimfiles)
    stimList = [];
    for sf = 1:length(stimfiles)
        SF = load([dr filesep filesep stimfiles(sf).name]);
        stimList = [stimList ; SF.time_stamp]; %#ok<AGROW>
    end
else
    error('no timing files in folder');
end

%load each file, identify stimulus time
dataset.stimTime = nan(1,length(fns));
for fnum = 1:length(fns)
    load([dr fns{fnum}]);
    if fnum==1
        dataset.X = nan([size(sys_recon.output.F) length(fns)]);
    end
    valid = sys_recon.input.ref.fusedInto>0;
    dataset.X(valid,:,fnum) = sys_recon.output.F(sys_recon.input.ref.fusedInto(valid), :);
    
    %stimulus timing
    dataset.stimTime(fnum) = sys_recon.input.scandata.metadata.timeNow;
end

%time of stimulus within the recording
prePeriod = [51:1000 3500:size(dataset.X,2)-50];
stimPeriod = 1001:3000;
postPeriod = 3001:(size(dataset.X,2)-50);

F0 = repmat(min(dataset.X(:, prePeriod,:),[],2), 1, size(dataset.X,2),1);
dataset.dFF = max(0, (dataset.X-F0)./(F0+1));
%p = 0.5; %norm for per-movie, per-seed normalization, should be in (0-1]
%dFFnorm = sum(abs(dataset.dFF).^p, 2).^(1/p);
%dataset.dFF = dataset.dFF./repmat(dFFnorm, 1, size(dataset.dFF,2),1); %normalize

%get stimulus IDs
dist = repmat(stimList(:,1), 1, length(dataset.stimTime)) - repmat(dataset.stimTime, size(stimList,1),1);
distP = dist; distM = dist;
distP(distP<0) = nan; distM(distM>0) = nan;
mindistP = nanmin(distP, [], 1);
mindistM = nanmax(distM,[],1);
if nanvar(mindistM)<nanvar(mindistP)
    dataset.stimulus.stimDelay = nanmedian(mindistM);
    dataset.stimulus.delayVariance = nanvar(mindistM);
else
    dataset.stimulus.stimDelay = nanmedian(mindistP);
    dataset.stimulus.delayVariance = nanvar(mindistP);
end
dataset.stimulus.stim = nan(length(fns),1);
dataset.stimulus.timeError = nan(length(fns),1);
dataset.stimulus.stimTime = nan(length(fns),1);
for fnum = 1:length(fns)
    [dataset.stimulus.timeError(fnum), minIx] =  min(abs(dataset.stimTime(fnum) - stimList(:,1) + dataset.stimulus.stimDelay));
    dataset.stimulus.timeError(fnum) = 86400*dataset.stimulus.timeError(fnum); %convert to seconds
    dataset.stimulus.stim(fnum) = stimList(minIx,2);
    dataset.stimulus.stimTime(fnum) = stimList(minIx,1);
end
figure,  scatter(86400*(dataset.stimulus.stimTime-min(dataset.stimulus.stimTime)), dataset.stimulus.stim);
xlabel('Stimulus time (s)');
ylabel('Stimulus ID');
title(['mean absolute stimulus timing error = ' num2str(mean(abs(dataset.stimulus.timeError))) ' seconds']);

%mean response for each stimulus
stims = unique(dataset.stimulus.stim);
nStim = length(stims);
meanResp = nan(size(dataset.dFF,1), size(dataset.dFF,2),nStim);
meanAmp = nan(size(dataset.dFF,1),nStim);
amp(:,:,1) = squeeze(nanmean(dataset.dFF(:,prePeriod, :),2));
amp(:,:,2) = squeeze(nanmean(dataset.dFF(:,stimPeriod, :),2));
amp(:,:,3) = squeeze(nanmean(dataset.dFF(:,postPeriod, :),2));

for stimNum = 1:nStim
    sel = dataset.stimulus.stim==stims(stimNum);
    meanResp(:,:,stimNum) = nanmean(dataset.dFF(:,:, sel),3);
    meanAmp(:,stimNum) = nanmean(amp(:,sel,2)./ amp(:,sel,1),2);
    semAmp(:,stimNum) = nanstd(amp(:,sel,2)./ amp(:,sel,1),0,2)./sqrt(sum(sel));
end

respAvg = nanmean(meanResp,3);
totalAmp = nanmean(meanAmp,2);

%compute an RGB color for each segment
%hue is the preferred direction    
rad = linspace(0, 2*pi, 8+1);
rad = reshape(rad(1:end-1), 1,[]);
rad = repmat(rad, size(meanAmp,1), 1);

H = max(0,min(1, (pi + circ_mean(rad,meanAmp,2))/(2*pi)));

%saturation is the degree of preference
pref = circ_r(rad,meanAmp,pi/4,2);
S = pref.*nansum(meanAmp,2);
S = min(1, S./(0.8*max(S(:))));

%V is the average amplitude
V = min(1, totalAmp./(0.8*max(totalAmp(:))));

RGB = hsv2rgb([H(:) S(:) V(:)]);

rVol = zeros([size(sys_recon.input.ref_image) 3]);
for seed = 1:size(sys_recon.input.S,2)-2
    if ~any(sys_recon.input.S(:,seed))
       continue 
    end
    rVol = rVol + reshape(sys_recon.input.S(:,seed)*RGB(seed,:), size(rVol));
end

thresh = prctile(rVol(rVol>0), 99);
rVol = min(1, rVol./thresh);

%A draggable point to measure tuning
hF = figure;
hIm = imshow( reshape(rVol, size(rVol,1), size(rVol,2)*size(rVol,3), size(rVol,4)));
hAxIm = get(hIm, 'parent'); set(hAxIm, 'pos', [0 0.5 1 0.5]);
for s =1:8
    hAxR(s) = axes('Pos', [0.1, (s-1)*0.5/8, 0.4, 0.48/8]);
end
hAxT = axes('Pos', [0.55, 0.05, 0.45, 0.45]);
%hT = title('Segment');
colors = hsv(nStim);
hPt = impoint(hAxIm); hPt.setColor('y');
done = false;
while ~done
    figure(hF);
    while ~waitforbuttonpress
        if ~ishandle(hF)
            break
        end
    end
    xy = getPosition(hPt);
    x = xy(1); y = xy(2);
    ind = sub2ind([size(rVol,1), size(rVol,2)*size(rVol,3)], round(y), round(x));
    [maxval, bestSeg] = max(refIM.seg(ind,1:end));
    if maxval>0
        set(hF,'name', ['Segment:' int2str(bestSeg)]),
        %subplot(nStim,2,[2 2*nStim])
        cla(hAxT);
        hold(hAxT,'on')
        errorbar(hAxT,meanAmp(bestSeg,:), semAmp(bestSeg,:)); %tuning curve
        set(hAxT, 'xlim', [0.5 nStim+0.5])
        for stimNum = 1:nStim
            %subplot(nStim,2,2*stimNum-1)%plot trials of each individual stim
            cla(hAxR(stimNum))
            hold(hAxR(stimNum), 'on')
            imagesc(hAxR(stimNum), squeeze(dataset.dFF(bestSeg,:, dataset.stimulus.stim==stims(stimNum)))');
            set(hAxR(stimNum),'ytick', [])
            if stimNum~=1
                set(hAxR(stimNum), 'xtick', [])
            end
            ylabel(hAxR(stimNum),int2str(stimNum));
        end
        xlabel(hAxR(1), 'Time (ms)')
        
    end
    done= ~ishandle(hF);
end

keyboard