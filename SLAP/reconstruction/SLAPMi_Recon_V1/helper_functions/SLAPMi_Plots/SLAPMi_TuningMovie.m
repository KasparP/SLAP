function SLAPMi_TuningMovie(dataset,optsin)
if nargin<1
   dataset = [];
   [fn, dr] = uigetfile('*.mat', 'Select a tuning dataset');
   disp('loading dataset...');
   load([dr fn]);
   disp('dataset loaded.');
end

hF = [];
opts.skipFrames = 5;  tooltips.skipFrames = 'Draw every Nth frame';
opts.appendname = ''; tooltips.skipFrames = 'Append this to the movie filename. If empty, solver params will be appended';
opts.MovieType = 'dPhotons'; %dff or dPhotons
opts.verbose = 100;

if nargin>1 %UPDATE WITH USER-SPECIFIED OPTIONS
    if ischar(optsin)
        opts.appendname = optsin;
    else
        for field = fieldnames(optsin)'
            opts.(field{1}) = optsin.(field{1});
        end
    end
else
      opts = optionsGUI(opts, tooltips);
end

prePeriod = 500:1000;
expFitPeriodPre = 50:1000;
expFitPeriodPost = 3550:3950;
%%%%%%%%%


%% normalize the data
[nSeeds,trialLen,ntrial] = size(dataset.X);
pnorm = 1;
%clipping to first 8 repetitions
select = 1:64;
dataset.dFF = dataset.dPhotons(:,:,select);
dataset.X = dataset.X(:,:,select);
dataset.filenames = dataset.filenames(select);
dataset.stimulus.stim = dataset.stimulus.stim(select);
dataset.stimulus.timeError = dataset.stimulus.timeError(select);
dataset.stimulus.stimTime = dataset.stimulus.stimTime(select);
dataset.Z = dataset.Z(select);
dataset.stimTime = dataset.stimTime(select);
dataset.spikes = dataset.spikes(select);

stims = unique(dataset.stimulus.stim);
nStim = length(stims);

rad = linspace(0,2*pi,nStim/2+1);
rad = rad([1:4 1:4]);

timeconstant = 80; %smoothing time constant for the timecourse plot
X1 = squeeze(isnan(dataset.X(:,1,:)));
for stimix = 1:nStim
    Xn(:,stimix) = sum(X1(:,dataset.stimulus.stim==stims(stimix)),2);
end
OutsideSLMSegs = any(Xn > length(dataset.stimulus.stim)/nStim/2,2);
% stimR = imgaussfilt3(squeeze(dataset.dFF),[eps timeconstant eps]);
stimR = smoothdata(dataset.dFF,2,'gaussian',timeconstant);
stimR(isnan(dataset.dFF)) = nan;
 
stimMean = nan(size(stimR,1), size(stimR,2), nStim);
for stimix = 1:nStim
    stimMean(:,:,stimix) = nanmean(stimR(:,:,dataset.stimulus.stim==stims(stimix)),3);
    weightsN(:,stimix) = sum(~isnan(stimR(:,1,dataset.stimulus.stim==stims(stimix))),3); 
%     atb = fit_exp(expFitPeriod,squeeze(stimMean(:,expFitPeriod,stimix)));     
%     stimMean(:,:,stimix) = stimMean(:,:,stimix) - (atb(3)+atb(1)*exp(-(1:size(stimMean,2))/atb(2)))';
    for segind = 1:size(stimR,1)
        if isnan(stimMean(segind,1,stimix))
            continue
        else
        stimMean(segind,:,stimix) = stimMean(segind,:,stimix) - mean(stimMean(segind,prePeriod,stimix),2);
        atb = fit_exp([expFitPeriodPre expFitPeriodPost],...
            [squeeze(stimMean(segind,expFitPeriodPre,stimix)),...
            min(0,squeeze(stimMean(segind,expFitPeriodPost,stimix))) ]);     
        stimMean(segind,:,stimix) = stimMean(segind,:,stimix) - (atb(3)+atb(1)*exp(-(1:size(stimMean,2))/atb(2)));
        end
    end
end


drawframes = 1:opts.skipFrames:trialLen;
sz = size(dataset.refIM.IM);
dataset.Zcenter = 6;
drawPixels = ((dataset.Zcenter-2)*sz(1)*sz(2)+1):((dataset.Zcenter+1)*sz(1)*sz(2));
valid3D = reshape(full(dataset.refIM.seg(drawPixels,:)),sz(1), sz(2),3,[]);
S2D = reshape(squeeze(sum(valid3D,3)),sz(1)*sz(2),[]);
valid3D = any(valid3D(:,:,:,~OutsideSLMSegs),4);
valid2D = any(valid3D,3);

ref2D = dataset.refIM.IM(:,:,dataset.Zcenter+(-1:1)); ref2D(~valid3D) = 0; 
ref2D = sum(ref2D,3);
Svalid = S2D(valid2D(:),~OutsideSLMSegs);
% Svalid = Svalid./(sum(Svalid,2));
Svalid = Svalid./(sum(Svalid));
Svalid(isnan(Svalid)) = 0;
%%
ix = 1;
H = nan(nSeeds,length(drawframes));
S = H; V = H;
for frame = drawframes
    [H(:,ix),S(:,ix),V(:,ix)] = circular_mean(repmat(rad,nSeeds,1),squeeze(stimMean(:,frame,:)));
    S(:,ix) = S(:,ix).*(min(weightsN,[],2)./max(min(weightsN,[],2))); %weigh the saturation by the number of observations
    % Could do the same with V
    ix = ix + 1;
end
Sthresh = 0.6;
S = S./Sthresh; S(S>1) = 1;
Vthresh = prctile(V(:),99.9);
factor = sqrt(Vthresh);
% Vthresh = max(V);
V = V./Vthresh; V(V>1) = 1;
H(OutsideSLMSegs,:) = 0.5;
S(OutsideSLMSegs,:) = 1;
V(OutsideSLMSegs,:) = 0;

ix = 1;
rVol = nan(size(Svalid,1),3*length(drawframes));
for frame = drawframes
    RGB = hsv2rgb([H(:,ix) S(:,ix) V(:,ix)]);
    rVol(:,(ix-1)*3+1:3*ix) = calc_rVol(RGB);
    ix = ix+1;
end
% rVolFactor = max(0.8,max(rVol,[],2));
rVolFactor = prctile(max(rVol(:,50*3:end),[],2),95);
rVol = rVol./max(rVolFactor,max(rVol(:,50*3:end),[],2));
%%
try
mkdir(dataset.dr, 'TuningMovies');
catch
    disp('Could not make TuningMovies directory !!')
end

VideoName = [dataset.dr filesep 'TuningMovies' filesep dataset.filenames{1}(1:end-16) '_TuningMovie'];
try
    W = VideoWriter(VideoName);
catch
    warning(['!! Could not write video in default directory, writing to PWD instead: ' pwd ' !!'])
    W = VideoWriter([dataset.filenames{1}(1:end-16) '_TuningMovie' opts.appendname]);
end
W.FrameRate = 30;
open(W);

f = figure('pos', [ 1 41 1280 1280]); hF(end+1) = f;
ax = axes(f, 'units', 'normalized');
h_im = imshow(0.5*ones([sz(1) sz(2) 3]), 'parent', ax, 'Border','tight');
colormap(ax, hsv2rgb([ones(100,1), linspace(0,1,100)', 0.8*ones(100,1)]))
h_cb = colorbar(ax);
axis image
set(f, 'resize', 'off')

set(h_cb, 'position', [0.96 0.06 0.01 0.87], 'axislocation', 'in', 'tickdir', 'out', 'ticklength', 0, 'ticks', [0.01 0.99], 'ticklabels', num2str([0 ; Sthresh], 3), 'ycolor', 'w', 'fontsize', 14)
h_cb.Title.String = 'ORI Index';
h_cb.Title.Color = 'w'; h_cb.Title.FontSize = 14;
h_cb.Label.Rotation = 0;
h_txt = text(ax,1,1, 'Frame X. t= XXX ms', 'units', 'characters', 'verticalalignment', 'baseline', 'horizontalalignment', 'left', 'color', 'w');
set(h_txt, 'position', [1 1])
h_txt.FontSize = 14;
linelength = (10/(sz(1)/5))/2; %10 microns. The whole FOV is 1280/5 = 256um
hs = annotation('line', [.96-linelength 0.96+linelength], [0.02, 0.02]);
hs.LineWidth = 5;
hs.Color = 'w';
hScalebarLabel = text(ax,0.96,0.03, '10 {\mu}m', 'verticalalignment', 'baseline', 'units', 'normalized', 'horizontalalignment', 'center', 'color', 'w');
hScalebarLabel.FontSize = 14;
ix = 1;
for frame = drawframes
    h_txt.String = ['  Frame #' int2str(frame) '       ' num2str((frame-1)/1.016, '%.1f') ' ms'];
    if mod(frame, opts.verbose)<opts.skipFrames
        disp(['Writing frame: ' int2str(frame)])
    end
    %sys.control_params.frame = frame;
    imFrame = min(1,repmat(sqrt(reshape(sum(ref2D,3),[],1)), 1, 3)./factor);
    imFrame(valid2D, :) = rVol(:,(ix-1)*3+1:3*ix);
    imFrame(~valid2D,:) = 0;
    F = reshape(imFrame, [sz(1),sz(2), 3]);
    set(h_im, 'Cdata', F); drawnow;
    writeVideo(W, getframe(ax));
    ix = ix+1;
end
disp('Done writing movie.');
disp(['Video Location:' VideoName]);
close(W)
%%
    function rVol = calc_rVol(RGB)
    RGB(any(isnan(RGB),2),:) = 0;
    rVol = Svalid*RGB(~OutsideSLMSegs,:);
%     rVol = rVol.*ref2D(valid2D);
    end
end

function [H,S,V] = circular_mean(rad,w)
% w = w - min(w,[],2);
A = sum(exp(1i*rad).*w,2);
H = 1/2 + angle(A)/(2*pi);
S = abs(A)./sum(abs(w),2);
% w = w - nanmean(w,2);
V = sqrt(sum(w.^2,2));
end

function atb_est = fit_exp(x,y)
F = @(atb,xdata) atb(1)*exp(-xdata/atb(2)) + atb(3);
atb0 = [max(y(1)-y(end),1),500,y(end)];
atb_est = lsqcurvefit(F,atb0,x,y,[0 0 -100]',[100,10000,0]',optimset('display','off'));

% I = any(isnan(y),2);
% y_temp = y(~I,:);
% 
% [n,T] = size(y_temp);
% 
% F = @(atb,xdata) atb(:,1).*exp(-xdata./atb(:,2)) + atb(:,3);
% atb0 = [max(y(:,1)-y(:,end),1),500*ones(1,size(y,1))',y(:,end)];
% atb_est = lsqcurvefit(F,atb0,x,y,zeros(n,3),repmat([100;10000;100],n,1),optimset('display','off'));
% 
% coeffs = nan(size(y,1),3);
% coeffs(~I,:) = atb_est;

end
