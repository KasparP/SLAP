function  SLAPMiTuning_Podgorski(dataset,refIM)
if nargin<1
   dataset = [];
   [fn, dr] = uigetfile('*.mat', 'Select a tuning dataset');
   disp('loading dataset...');
   load([dr fn]);
   disp('dataset loaded.');
end
if nargin<2
    refIM = dataset.refIM;
end

dataset = QCSLAPMiDataset (dataset);


%%%%%% normalize
% dataset.X = dataset.X./(max(dataset.X,[],2)+1);


%%%%%
prePeriod = 500:1000;
bleachPeriod1 = 100:550;
bleachPeriod2 = 551:1000;
stimPeriod = 1100:3000;
postPeriod = []; % 3100:(size(dataset.X,2)-50);
expFitPeriodPre = [50:1000]; expFitPeriodPost = [3550:3950];
timeconstant = 40; %smoothing time constant for the timecourse plot

%%%%%%%%%
insideSLM = repmat(dataset.SLM > 0.9, 1,1,size(refIM.IM,3)); 
insideSLMsegs = full(sum(refIM.seg(insideSLM(:),:),1)>0);

%% normalize the data
[nSeeds,trialLen,ntrial] = size(dataset.X);

stims = unique(dataset.stimulus.stim);
nStim = length(stims);
rad = linspace(0,2*pi,nStim/2+1);
rad = rad([1:4 1:4]);
pnorm = 1;
indices = 1:nStim;

weights = nan(nSeeds,nStim);
weightsN = nan(nSeeds,nStim);
weightsSTD = nan(nSeeds,nStim);

for stimNum = 1:nStim
    sel = dataset.stimulus.stim==stims(stimNum);
    X = nanmean(dataset.dPhotons(:,:,sel),3);
    weights(:,stimNum) = nanmean(X(:,stimPeriod).^pnorm,2).^(1/pnorm)  - nanmean(X(:,[prePeriod postPeriod]).^pnorm,2).^(1/pnorm);  
    weightsN(:,stimNum) = sum(~isnan(sum(dataset.dPhotons(:,stimPeriod,sel),2)),3); % number of valid recordings
    Y = squeeze(nanmean(dataset.dPhotons(:,stimPeriod,sel),2)-nanmean(dataset.dPhotons(:,[prePeriod postPeriod], sel).^pnorm,2).^(1/pnorm));
    weightsSTD(:,stimNum) = nanstd(Y,0,2)/sqrt(sum(sel));
    F0(:,stimNum) = nanmean(dataset.F0_photons(:,sel),2);  
end
XBleach  = nanmean(dataset.dPhotons,3);
Bleach = nanmean(XBleach(:,bleachPeriod1),2)-nanmean(XBleach(:,bleachPeriod2),2);
Bleach = max(Bleach,0);
weights = weights(:,indices)+Bleach;

[H,S,V] = circular_mean(repmat(rad,nSeeds,1),weights);

V = V.*(min(weightsN,[],2)./max(min(weightsN,[],2))); %weigh the value by the number of valid observations

V = max(0, V-prctile(V(V>0), 25)); %

Sthresh = max(prctile(S(~isnan(S)), 60), 0.35);
S = S./Sthresh; S(S>1) = 1;
Vthresh = max(V(:));%prctile(V(~isnan(V)),99.9);
V = V./Vthresh; V(V>1) = 1;
H(isnan(H)) = 0; S(isnan(S)) = 0; V(isnan(V)) = 0;

H(~insideSLMsegs) = 0.5;
S(~insideSLMsegs) = 1;
V(~insideSLMsegs) = 0;

RGB = hsv2rgb([H(:) S(:) V(:)]);

%color balance
sel = any(RGB>0,2);
gamma = 0.25;
for c = 1:3
    factor =(1/mean(RGB(sel,c))).^gamma;
    RGB(:,c) = RGB(:,c).*factor;
end

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
segs = [];

hFB = figure('color', 'k'); hAxB = axes('parent', hFB);

set (hF, 'WindowScrollWheelFcn', @scrollFunc);
set (hF, 'WindowButtonDownFcn', @bdf);
set (hF, 'WindowKeyPressFcn', @kpf);


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
            if strcmpi(src.SelectionType, 'extend')
                ptPos = [ptPos; round([cp(1) cp(2) Z])];
                ind = sub2ind(size(refIM.IM),ptPos(end,2), ptPos(end,1),ptPos(end,3));
                [maxval, bestSeg] = max(refIM.seg(ind,:));
                segs = [segs; bestSeg];
            else
                ptPos = round([cp(1) cp(2) Z]);
                ind = sub2ind(size(refIM.IM),ptPos(2), ptPos(1),ptPos(3));
                [maxval, bestSeg] = max(refIM.seg(ind,:));
                segs = bestSeg;
            end
            if ~isempty(segs)
                set(hF,'name', ['Segment:' int2str(bestSeg)]),
                set(hF2,'name', ['Segment:' int2str(bestSeg)]),
                
                %tuning plot
                F0seg = sum(F0(segs,:),1);
                cla(hAxTuning); cla(hAxTimecourse);
                hold(hAxTuning,'on'); hold(hAxTimecourse,'on');
                errorbar(hAxTuning,sum(weights(segs,:),1)./F0seg, sqrt(sum(weightsSTD(segs,:).^2,1))./F0seg, 'linewidth', 2, 'color', 'k', 'linestyle', ':'); %tuning curve - errors are being approximated; should fix
                for stimind= indices
                   scatter(hAxTuning,stimind, sum(weights(segs,stimind),1)./F0seg(stimind), 'sizedata', 120, 'marker', 'o', 'markerFaceColor', colors(stimind,:), 'markerEdgeColor', 'k', 'linewidth', 2);
                end
                set(hAxTuning, 'xlim', [0.5 nStim+0.5])
                set(hAxTimecourse, 'xlim', [0 size(dataset.dPhotons,2)-50]) 
                    
                raster = cell(nStim,1);
                for stimind = indices
                    raster{stimind,1} = squeeze(sum(dataset.dPhotons(segs,:, dataset.stimulus.stim==stims(stimind)),1))';
                    if stimind<nStim
                        raster{stimind,2} = nan(1, size(dataset.dPhotons,2));
                    end
                    
                    %plot the timecourse
                    stimR =imgaussfilt(squeeze(sum(dataset.dPhotons(segs,:, dataset.stimulus.stim==stims(stimind)),1)),[timeconstant eps]);
                    stimMean = nanmean(stimR,2);
                    stimMean = stimMean - mean(stimMean(prePeriod));
                    if any(~isnan(stimMean))
                        atb = fit_exp([expFitPeriodPre expFitPeriodPost],[stimMean(expFitPeriodPre)' min(0, stimMean(expFitPeriodPost))']);
                        stimMean = stimMean - (atb(3)+atb(1)*exp(-(1:size(stimMean,1))'/atb(2)));
                    end
                    %TCtmp(stimind) = mean(stimMean(stimPeriod));
                    stimSTD = nanstd(stimR,0,2)./sqrt(sum(~isnan(stimR),2));
                    
                    F0stim = squeeze(sum(nanmean(dataset.F0_photons(segs,dataset.stimulus.stim==stims(stimind)),2),1));
                    stimMean = stimMean./F0stim;
                    stimSTD = stimSTD./F0stim;
                    
                    hLineTC(stimind) = plot(hAxTimecourse,stimMean, 'color', colors(stimind,:), 'linewidth', 2);
                    hPatchTC(stimind) = patch(hAxTimecourse, [1:size(dataset.dPhotons,2) fliplr(1:size(dataset.dPhotons,2))]', [stimMean+stimSTD ; flipud(stimMean-stimSTD)], colors(stimind,:),'edgecolor', 'none', 'FaceAlpha', 0.3);
                end
                %figure, plot(TCtmp)
                
                cla(hAxRaster);
                Raster = cell2mat(reshape(raster',[],1));
                RasterScale = max(prctile(Raster(:),99),1);
                
                imagesc(hAxRaster, Raster);
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
                set(hAxRaster, 'ytick', yt, 'yticklabel', int2str((1:nStim)'), 'Clim', [0 RasterScale])
                xlabel(hAxRaster, 'Time (samples)');
                xlabel(hAxTimecourse, 'Time (samples)')
                xlabel(hAxTuning, 'Stimulus #')
                
                
                if ~isvalid(hAxB)
                    hFB = figure('color', 'w'); hAxB = axes('parent', hFB);
                else
                    cla(hAxB)
                end
                base = 1;
                for stim = 8:-1:1
                    for trace = 1:size(raster{stim,1},1)
                        plot(hAxB, base+raster{stim,1}(trace,:)/RasterScale, 'color', colors(stim,:), 'linewidth', 1.5)
                        hold(hAxB, 'on')
                        base = base+1;
                    end
                    
                    base = base+2;
                end
                yl = get(hAxB, 'ylim');
                plot(hAxB, [1016 1016 nan 3*1016 3*1016], [yl nan yl], ':', 'linewidth', 2, 'color', 'w')
                set(hAxB, 'linewidth', 2, 'color', 'w', 'ycolor', 'none', 'xcolor', 'w', 'box', 'off', 'tickdir', 'out', 'fontsize', 12)
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
        hPt = [];
        for ptix = 1:size(ptPos,1)
        if ptPos(ptix,3)==Z
            hPt(end+1) = scatter(hAxIm, ptPos(ptix,1),ptPos(ptix,2), 'sizedata', 50, 'linewidth', 2, 'markeredgecolor', 'r');
        else
            hPt(end+1) = scatter(hAxIm, ptPos(ptix,1),ptPos(ptix,2), 'sizedata', 10, 'linewidth', 2, 'markeredgecolor', 'r');
        end
        end
    end

end

function [rVol, rVol2D] = calc_rVol(RGB, refIM)
RGB(any(isnan(RGB),2),:) = 0;

blockseg = refIM.seg;
select = any(blockseg,1);
blockseg(:,select) = blockseg(:,select)./sum(blockseg(:,select),1);
rVol = blockseg*RGB(1:end-2,:);
rVol(~any(blockseg,2),:) = 0;
rVol = reshape(full(rVol), [size(refIM.IM) 3]);


rVol2D = squeeze(sum(rVol, 3));
thresh = prctile(rVol2D(rVol2D>0), 99.9);
rVol2D = rVol2D./max(max(rVol2D(3), thresh));
figure, imshow(rVol2D);


centerplanes = ceil(size(rVol,3)/2)+(-1:1);
selected = max(max(rVol(:,:,centerplanes,:),[],4),[],3);
thresh = prctile(selected(selected>0), 99);
rVol = rVol./max(max(rVol,[],4),thresh);
end


function [H,S,V] = circular_mean(rad,w)
% w = w - min(w,[],2);
A = sum(exp(1i*rad).*w,2);
H = 1/2 + angle(A)/(2*pi);
S = abs(A)./sum(abs(w),2);
% w = w - nanmean(w,2);
V = sqrt(sum(w.^2,2));
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


function atb_est = fit_exp(x,y)
F = @(atb,xdata) atb(1)*exp(-xdata/atb(2)) + atb(3);
atb0 = [max(y(1)-y(end),1),500,y(end)];
atb_est = lsqcurvefit(F,atb0,x,y,[0 0 -100]',[100,10000,0]', optimset('Display', 'off'));
end
