function  dataset = SLAPMiTuning(dataset,refIM)
% [~,dr] = uigetfile('*.mat','Choose the stimulus timing file');
% %LOAD STIMULUS
% stimfiles = dir([dr  filesep '*Timings*']);
% if ~isempty(stimfiles)
%     stimList = [];
%     for sf = 1:length(stimfiles)
%         SF = load([dr filesep filesep stimfiles(sf).name]);
%         stimList = [stimList ; SF.time_stamp]; %#ok<AGROW>
%     end
% else
%     error('no timing files in folder');
% end


if nargin<1
   dataset = [];
   [fn, dr] = uigetfile('*.mat', 'Select a tuning dataset');
   load([dr fn]);
end
if nargin<2
    refIM = dataset.refIM;
end

%%%%%% normalize
% dataset.X = dataset.X./(max(dataset.X,[],2)+1);


%%%%%
prePeriod = 500:1000;
stimPeriod = 1100:3000;
postPeriod = 3100:(size(dataset.X,2)-50);

%%%%%%%%%
insideSLM = repmat(dataset.SLM > 0.9, 1,1,size(refIM.IM,3)); 
%imTemp =  repmat(insideSLM,1,1,3);
% dataset.S = dataset.S.*imTemp(:);
% ss = sum(dataset.S);
% dataset.X(ss == 0,:) = 0;
% 
insideSLMsegs = full(sum(refIM.seg(insideSLM(:),:),1)>0);
% dataset.S = dataset.S(:,insideSLMsegs);
% dataset.X = dataset.X(insideSLMsegs,:,:);


%% normalize the data
[nSeeds,trialLen,~] = size(dataset.X);

F0tmp = min(dataset.X(:, 50:end-50,:),[],2);
dataset.dFF = max(0, (dataset.X-F0tmp)./(F0tmp+1)); %this absorbs any resting F0 into the existing dFF calculation

%divisive normalization
pnorm = 2;
dataset.dFF = dataset.dFF./(nanmean(dataset.dFF(:, [prePeriod postPeriod], :).^pnorm,2).^(1/pnorm));

%normalize to max
% dataset.dFF = 2*dataset.dFF./(max(dataset.dFF(:, [prePeriod postPeriod],:),[],2)+1); %this normalizes the amplitudes (while suppressing very flat traces)
% dataset.dFF(isnan(dataset.X)) = nan;

% Fs = 1016;
% [LX, HX] = decompose(dataset.dFF,Fs);
% dataset.dFF = LX;

figure('units', 'norm', 'pos', [0.05    0.65    0.9    0.2]),  scatter(86400*(dataset.stimulus.stimTime-min(dataset.stimulus.stimTime)), dataset.stimulus.stim);
xlabel('Stimulus time (s)');
ylabel('Stimulus ID');
title(['mean absolute stimulus timing error = ' num2str(mean(abs(dataset.stimulus.timeError))) ' seconds']);
set(gca,'fontsize',20); drawnow;
%mean response for each stimulus
stims = unique(dataset.stimulus.stim);
nStim = length(stims);

% rad = linspace(0,pi,nStim/2+1);

rad = linspace(0,2*pi,nStim/2+1);
rad = rad([1:4 1:4]);

timeconstant = 40; %smoothing time constant for the timecourse plot
indices = 1:nStim;

weights = nan(nSeeds,nStim);
weightsN = nan(nSeeds,nStim);
weightsSTD = nan(nSeeds,nStim);
%weights = zeros(nSeeds,nStim);
% dataset.X(isnan(dataset.X))= 0;

for stimNum = 1:nStim
    sel = dataset.stimulus.stim==stims(stimNum);
    X = nanmean(dataset.dFF(:,:,sel),3);
    %weights(:,stimNum) = nanmean(X(:,stimPeriod),2)-nanmean(X(:,prePeriod),2);
    weights(:,stimNum) = nanmean(X(:,stimPeriod).^pnorm,2).^(1/pnorm)  - nanmean(X(:,[prePeriod postPeriod]).^pnorm,2).^(1/pnorm);
    weightsN(:,stimNum) = sum(~isnan(sum(dataset.dFF(:,stimPeriod,sel),2)),3);
    %Y = squeeze(nanmean(dataset.dFF(:,stimPeriod,sel),2)-nanmean(dataset.dFF(:,prePeriod,sel),2));
    Y = squeeze(nanmean(dataset.dFF(:,stimPeriod,sel),2)-nanmean(dataset.dFF(:,[prePeriod postPeriod], sel).^pnorm,2).^(1/pnorm));
    weightsSTD(:,stimNum) = nanstd(Y,0,2)/sqrt(sum(sel));
end

[H,S,V] = circular_mean(repmat(rad,nSeeds,1),weights(:,indices));
S = S.*(min(weightsN,[],2)./max(min(weightsN,[],2))); %weigh the saturation by the number of observations
% V = min((V-min(V))./(prctile(V,95)-min(V)),1);

H(~insideSLMsegs) = 0.5;
S(~insideSLMsegs) = 1;
V(~insideSLMsegs) = 0;

RGB = hsv2rgb([H(:) S(:) V(:)]);
disp('Calculating volume tuning image...')
rVol = calc_rVol(RGB,refIM);


colors = hsv(nStim);
ysz = [size(refIM.IM,1), size(refIM.IM,2)];
Z = ceil(size(rVol,3)/2); nZ = size(rVol,3);
hF = figure('color', 'w');
hIm = imshow(squeeze(rVol(:,:,Z,:)));
hAxIm = get(hIm, 'parent'); set(hAxIm, 'pos', [0 0 1 1]);
hold(hAxIm,'on');

hF2 = figure('color', 'w');
hAxRaster = axes('Pos', [0.05, 0.05, 0.9, 0.45], 'tickdir', 'out', 'linewidth', 1.5);
hAxTimecourse = axes('Pos', [0.05, 0.55, 0.9, 0.2], 'tickdir', 'out', 'linewidth', 1.5);
hAxTuning = axes('Pos', [0.05, 0.8, 0.9, 0.2], 'tickdir', 'out', 'linewidth', 1.5);
hPt = [];
ptPos = [0 0 0];

hFB = figure('color', 'k'); hAxB = axes('parent', hFB);

set (hF, 'WindowScrollWheelFcn', @scrollFunc);
set (hF, 'WindowButtonDownFcn', @bdf);
set (hF, 'WindowKeyPressFcn', @kpf);
set ([hF hF2], 'CloseRequestFcn', @exit);



    function scrollFunc(obj, evnt)
        UPDN = evnt.VerticalScrollCount;
        Z = max(1,min(nZ, Z+UPDN));
        redrawImage;
    end

    function kpf(src, event)
        
    end

    function bdf(src, event)
        cp = get(hAxIm,'CurrentPoint');
        cp = cp(1, 1:2);
        %If the point is within the axes
        if all(cp>=1 & cp<=ysz)
            ptPos = round([cp(1) cp(2) Z]);
            ind = sub2ind(size(refIM.IM),ptPos(2), ptPos(1),ptPos(3));
            [maxval, bestSeg] = max(refIM.seg(ind,:));
            if ~isempty(bestSeg)
                
                %find all shaft ROIs within XXX microns, and
                %calculate average tuning of those
                %show the centroids of the shaft ROIs
                
                %keyboard
                
                set(hF,'name', ['Segment:' int2str(bestSeg)]),
                set(hF2,'name', ['Segment:' int2str(bestSeg)]),
                
                %tuning plot
                cla(hAxTuning); cla(hAxTimecourse);
                hold(hAxTuning,'on'); hold(hAxTimecourse,'on');
                errorbar(hAxTuning,weights(bestSeg,:), weightsSTD(bestSeg,:), 'linewidth', 2, 'color', 'k', 'linestyle', ':'); %tuning curve
                for stimind= indices
                   scatter(hAxTuning,stimind, weights(bestSeg,stimind), 'sizedata', 120, 'marker', 'o', 'markerFaceColor', colors(stimind,:), 'markerEdgeColor', 'k', 'linewidth', 2);
                end
                set(hAxTuning, 'xlim', [0.5 nStim+0.5])
                set(hAxTimecourse, 'xlim', [0 size(dataset.dFF,2)-50]) 
                
                raster = cell(nStim,1);
                for stimind = indices
                    raster{stimind,1} = squeeze(dataset.dFF(bestSeg,:, dataset.stimulus.stim==stims(stimind)))';
                    if stimind<nStim
                        raster{stimind,2} = nan(1, size(dataset.dFF,2));
                    end
                    
                    %plot the timecourse
                    stimR =imgaussfilt(squeeze(dataset.dFF(bestSeg,:, dataset.stimulus.stim==stims(stimind))),[timeconstant eps]);
                    stimMean = nanmean(stimR,2);
                    stimSTD = nanstd(stimR,0,2)./sqrt(sum(~isnan(stimR),2));

                    hLineTC(stimind) = plot(hAxTimecourse,stimMean, 'color', colors(stimind,:), 'linewidth', 2);
                    hPatchTC(stimind) = patch(hAxTimecourse, [1:size(dataset.dFF,2) fliplr(1:size(dataset.dFF,2))]', [stimMean+stimSTD ; flipud(max(0, stimMean-stimSTD))], colors(stimind,:),'edgecolor', 'none', 'FaceAlpha', 0.3);
                end
                
                cla(hAxRaster); 
                imagesc(hAxRaster, cell2mat(reshape(raster',[],1)));
                hold(hAxRaster, 'on');
                nlines= 0;
                for ix = 1:nStim
                    yt(ix) = nlines + floor((size(raster{ix,1},1)+1)/2);
                    if ix<nStim
                        nlines = nlines+size(raster{ix,1},1)+1;
                        plot(hAxRaster, [1 size(raster{ix,1},2)], [nlines nlines], 'k', 'linewidth', 2.5)
                    end
                end
                plot(hAxRaster, [1016 1016 nan 3*1016 3*1016], [get(hAxRaster, 'ylim') nan get(hAxRaster, 'ylim')], ':r', 'linewidth', 2)
                set(hAxRaster, 'ytick', yt, 'yticklabel', int2str((1:nStim)'), 'Clim', [0 3])
                xlabel(hAxRaster, 'Time (samples)');
                xlabel(hAxTimecourse, 'Time (samples)')
                xlabel(hAxTuning, 'Stimulus #')
                
                
                if ~isvalid(hAxB)
                    hFB = figure('color', 'k'); hAxB = axes('parent', hFB);
                else
                    cla(hAxB)
                end
                base = 1;
                for stim = 8:-1:1
                    for trace = 1:8
                        plot(hAxB, base+raster{stim,1}(trace,:)/3, 'color', colors(stim,:), 'linewidth', 1.5)
                        hold(hAxB, 'on')
                        base = base+1;
                    end
                    
                    base = base+2;
                end
                yl = get(hAxB, 'ylim');
                plot(hAxB, [1016 1016 nan 3*1016 3*1016], [yl nan yl], ':', 'linewidth', 2, 'color', 'w')
                set(hAxB, 'linewidth', 2, 'color', 'k', 'ycolor', 'none', 'xcolor', 'w', 'box', 'off', 'tickdir', 'out', 'fontsize', 12)
                xlabel('time (frames)')

            else
                keyboard
            end
            redrawImage;
        end
    end

    function exit(src, evnt)
        if strcmp(questdlg('Exit SLAPMiTuning?', 'SLAPMiTuning', 'Yes', 'No', 'No'), 'Yes')
            delete([hF hF2]);
        end
    end

    function redrawImage
        set(hIm, 'CData', squeeze(rVol(:,:,Z,:)));
        delete(hPt);
        if ptPos(3)==Z
            hPt = scatter(hAxIm, ptPos(1),ptPos(2), 'sizedata', 50, 'linewidth', 2, 'markeredgecolor', 'r');
        else
            hPt = scatter(hAxIm, ptPos(1),ptPos(2), 'sizedata', 10, 'linewidth', 2, 'markeredgecolor', 'r');
        end
    end

end

function rVol = calc_rVol(RGB, refIM)
RGB(any(isnan(RGB),2),:) = 0;

blockseg = refIM.seg;
rVol = blockseg*RGB(1:end-2,:);
rVol = rVol./repmat(sum(blockseg,2), 1,3);
rVol(~any(blockseg,2),:) = 0;

rVol = reshape(full(rVol), [size(refIM.IM) 3]);
centerplanes = ceil(size(rVol,3)/2)+(-1:1);
selected = max(max(rVol(:,:,centerplanes,:),[],4),[],3);
thresh = prctile(selected(selected>0), 80);
rVol = rVol./max(max(rVol,[],4),thresh);

%weight by refIM intensity?
doweight = true;
if doweight
    IM3Dgamma = sqrt(refIM.IM);
    IM3Dgamma = min(1, IM3Dgamma./prctile(IM3Dgamma(:), 99.9));
    rVol = rVol.*IM3Dgamma;
end
end


function [H,S,V] = circular_mean(rad,w)
A = sum(exp(1i*rad).*w,2);
H = 1/2 + angle(A)/(2*pi);

Sthresh = 0.6;
S = min(1, (abs(A)./sum(abs(w),2))./Sthresh);
S(isnan(A)) = nan;

V = sqrt(sum(w.^2,2));
V = min(1, V./prctile(V,99));
V(isnan(A)) = nan;
end

function [LX, HX] = decompose(X,Fs)
disp('Filtering the data ...')
if nargin<2; Fs = 1016; end
Ts = 1/Fs;
T = size(X,2);
% t = (0:T-1)*Ts;
Y = fft(X,[],2);
f = Fs*linspace(0,1/2,T/2+1);
fc = 10;
if mod(T,2) == 0
   f = [f,fliplr(f(2:end-1))];
else
   f = [f,fliplr(f(2:end))];
end

YLX = Y; YLX(:,f>fc,:)=0;
YHX = Y; YHX(:,f<=fc,:)=0;

LX = real(ifft(YLX,[],2));
HX = real(ifft(YHX,[],2));

end


