function [hF,mRaw,mDFF] = SLAPMi_Plot(sys, optsin)
if isfield(sys.output,'SLM_mask')
sys.opts.SLM_threshold = prctile(sys.output.SLM_mask(:),5);
end
if isfield(sys.output,'SLM_mask')
    InsideSLM = sys.output.SLM_mask > sys.opts.SLM_threshold;
    nPlanes = size(sys.input.ref_image,3);
    SLM_extended = repmat(InsideSLM,1,1,nPlanes);
    S_temp = (sys.input.S./sum(sys.input.S,1)).*~SLM_extended(:);
    sys.output.OutsideSLMSeeds = sum(S_temp,1)>1e-2;
end
    
    sys.output.F(sys.output.OutsideSLMSeeds,:) = 0;
    sys.output.spikes(sys.output.OutsideSLMSeeds,:) = 0;
    disp('Removing segments outside SLM');
    
hF = [];
opts.skipFrames = 5;  tooltips.skipFrames = 'Draw every Nth frame';
opts.appendname = ''; tooltips.appendname = 'Append this to the movie filename. If empty, solver params will be appended';
opts.t0 = 1016; tooltips.t0 = 'optional; Frame # corresponding to stimulus onset';
opts.t_end =  3200;   tooltips.t_end = 'optional; Last frame of visual motion stimulus';

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

opts.MovieType = 'min(dFF,Z)'; %opts.MovieType = {'''min(dFF,Z)'''; '''dPhotons'''; '''dFF'''; '''Zscore'''};

if ~nargin ||   isempty(sys)
    [fns, dr] = uigetfile('*.mat', 'Select you saved RECON files for plotting', 'multiselect', 'on');
    if ~iscell(fns)
        fns = {fns};
    end
    sys_recon = 0;
    for fnum = 1:length(fns)
        load([dr filesep fns{fnum}]);
        hF = SLAPMi_Plot(sys_recon, opts);
        try
            delete(hF);
        end
    end
    return
end
if isempty(opts.appendname)
    opts.appendname = AppendName(sys);
end
if isfield(sys.opts,'StimOnset')    
    StimOnset = sys.opts.StimOnset;
% elseif isfield(opts,'StimOnset')
%     StimOnset = opts.StimOnset;
else
    StimOnset = [];
end


%% Summary of the Reconstructions

hF(end+1) = figure;
subplot(2,2,[1 2]);

plot(sys.output.F(:, 1:end-20)');

title('Delta F/F for each segment');
subplot(2,2,3);
plot(sys.control_params.objective);
title('Objective function'); xlabel('Iteration')
ylabel('Objective value')
recon_error = (sys.input.y-sys.output.PS*(sys.output.F) - sys.output.additive_baseline - sys.output.stim_artifact)./sys.input.y;
subplot(2,2,4); imagesc(recon_error,[-3 3]);
title('Reconstruction Error')
drawnow;


%drawframes = [1:opts.skipFrames:399 400:699 700:opts.skipFrames:1990]; %for 2-spot uncaging
drawframes = [1:opts.skipFrames:opts.t0 opts.t0+1:opts.t_end opts.t_end:opts.skipFrames:(size(sys.output.F,2)-10)];
%drawframes = 1:opts.skipFrames:(size(sys.output.F,2)-10);

%% Calculating Reconstructed Movies
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch opts.MovieType
    case 'dFF'
        [mRaw,mDFF,valid2D] = calc_Movies(sys,drawframes,[],StimOnset);
        dffMax = round(2*max(1, min(5, prctile(mDFF(:), 99.5))))/2;
    case 'dPhotons'
        [mRaw,mDFF,valid2D] = calc_dPhotons(sys,drawframes,[],StimOnset);
        dffMax = max(5,prctile(mDFF(:), 99.9));
    case 'Zscore'
        [mRaw,mDFF,valid2D] = calc_Zscore(sys,drawframes,[],StimOnset);
        dffMax = round(2*max(2.5,prctile(mDFF(:), 99.9)))/2;
    case 'min(dFF,Z)'
        [mRaw1,mDFF1,valid2D1] = calc_Movies(sys,drawframes,[],StimOnset);
        dffMax = round(2*max(1, min(5, prctile(mDFF1(:), 99.5))))/2  %1.5
        
        [mRaw2,mDFF2,valid2D2] = calc_Zscore(sys,drawframes,[],StimOnset);
        ZMax = round(2*max(2.5,prctile(mDFF2(:), 99.9)))/2 %3
        
        mRaw = min(mRaw1, mRaw2.*(dffMax./ZMax));
        mDFF = min(mDFF1, mDFF2.*(dffMax./ZMax));
        valid2D = valid2D1 & valid2D2;
    otherwise
        keyboard
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gamma correction to better see dim features
factor = sqrt(prctile(mRaw(:), 99.9));
hsvV = min(1, sqrt(mRaw(:))./factor);
hsvS = min(1,mDFF(:)./dffMax);
hsvH = ones(size(mDFF(:)));
rgb_image = hsv2rgb([hsvH hsvS hsvV]);
rgb_image = reshape(rgb_image, [size(mRaw) 3]);
try
mkdir(sys.opts.dr, 'ReconMovies');
mkdir([sys.opts.dr filesep 'ReconMovies'], opts.MovieType);
catch
    disp('Error: Could not make ReconMovies directory! Attempting to Continue')
end
VideoName = [sys.opts.dr filesep 'ReconMovies' filesep opts.MovieType sys.opts.filename(1:end-11) 'dffMovie' opts.appendname];
opts.videoMode = 'MPEG-4';
try
    W = VideoWriter(VideoName, opts.videoMode);
catch
    warning(['!! Could not write video in default directory, writing to PWD instead: ' pwd ' !!'])
    W = VideoWriter([sys.opts.filename(1:end-16) '_' opts.MovieType '_' opts.appendname], opts.videoMode);
end
W.FrameRate = 30;
open(W);

f = figure('pos', [ 1 41 1280 1280]); hF(end+1) = f;
ax = axes(f, 'units', 'normalized');
h_im = imshow(0.5*ones([sys.opts.dim 3]), 'parent', ax, 'Border','tight', 'initialmag', 300);
colormap(ax, hsv2rgb([ones(100,1), linspace(0,1,100)', 0.8*ones(100,1)]))
h_cb = colorbar(ax);
axis image
set(f, 'resize', 'off')
set(h_cb, 'position', [0.96 0.06 0.01 0.87], 'axislocation', 'in', 'tickdir', 'out', 'ticklength', 0, 'ticks', [0.01 0.99], 'ticklabels', num2str([0 ; dffMax], 3), 'ycolor', 'w', 'fontsize', 14)
switch opts.MovieType
    case 'dPhotons'
        h_cb.Title.String = '{\Delta}Photons  ';
    case 'Zscore'
        h_cb.Title.String = 'Z score  ';
    case 'dFF'
        h_cb.Title.String = '{\Delta}F/F_0 ';
    case 'min(dFF,Z)'
        h_cb.Title.String = '{\Delta}F/F_0 ';
    otherwise
        keyboard
end
h_cb.Title.Color = 'w'; h_cb.Title.FontSize = 12.5;
h_cb.Label.Rotation = 0;
h_txt = text(ax,1,1, 'Frame X. t= XXX ms', 'units', 'characters', 'verticalalignment', 'baseline', 'horizontalalignment', 'left', 'color', 'w');
set(h_txt, 'position', [1 1])
h_txt.FontSize = 14;
linelength = (10/(sys.opts.dim(1)/5))/2; %10 microns. The whole FOV is 1280/5 = 256um
hs = annotation('line', [.93-linelength 0.93+linelength], [0.01, 0.01]);
hs.LineWidth = 5;
hs.Color = 'w';
hScalebarLabel = text(ax,0.93,0.02, '10 {\mu}m', 'verticalalignment', 'baseline', 'units', 'normalized', 'horizontalalignment', 'center', 'color', 'w');
hScalebarLabel.FontSize = 12;
if opts.t_end
    htxtStimOn= text(ax,1,1,'      Blank Screen', 'units', 'characters', 'verticalalignment', 'baseline', 'horizontalalignment', 'left', 'color', [0.5 0.5 0.5]);
    set(htxtStimOn, 'position', [70 1])
    htxtStimOn.FontSize = 14;
end
for frame = 1:length(drawframes)
    h_txt.String = ['  Frame #' int2str(drawframes(frame)) '       t= ' num2str((drawframes(frame)-1-opts.t0)/1.016, '%.1f') ' ms'];
    if any(frame==round(linspace(1,length(drawframes), 40)))
        disp(['Writing frame: ' int2str(drawframes(frame))])
    end
    %sys.control_params.frame = frame;
    imFrame = min(1,repmat(sqrt(reshape(sum(sys.input.ref_image,3),[],1)), 1, 3)./factor);
    imFrame(valid2D, :) = squeeze(rgb_image(:,frame,:));
    imFrame(~valid2D,3) = imFrame(~valid2D,3)* 0.6 + 0.05;
    imFrame(~valid2D,1:2) = imFrame(~valid2D,1:2) * 0.5;
    %imFrame= imFrame.*repmat(sqrt(sys.SLMmask(:)), 1, 3);  %black out the background
    F = reshape(imFrame, [sys.opts.dim 3]);
    set(h_im, 'Cdata', F);
    
    %plot stimulus time
    if opts.t_end
        if drawframes(frame)>opts.t0
            phase = round(360*mod(2*floor((drawframes(frame)-opts.t0)/16.6)./floor((opts.t_end-opts.t0)/16.6),1));
            if drawframes(frame)<opts.t_end
                htxtStimOn.String = ['Grating Phase: ' int2str(phase) '�'];
                set(htxtStimOn, 'color', [0.8 0.8 1])
            else
                htxtStimOn.String = ['      Blank Screen'];
                set(htxtStimOn, 'color', [0.5 0.5 0.5])
            end
        end
    end
    written = false;
    while ~written
    try
        writeVideo(W, getframe(ax));
        written = true;
    catch
        disp('Error writing frame. Retrying...'); pause(0.1);
    end
    end
end
disp('Done writing movie.');
disp(['Video Location:' VideoName]);
close(W)

end
