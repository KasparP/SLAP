function [hF,mRaw,mDFF] = SLAPMi_ZPlot(sys, optsin)
if isfield(sys.output,'SLM_mask') && ~isfield(sys.output,'OutsideSLMSeeds')
    InsideSLM = sys.output.SLM_mask < 0.8;
    nPlanes = size(sys.input.ref_image,3);
    SLM_extended = repmat(InsideSLM,1,1,nPlanes);
    S_temp = sys.input.S.*SLM_extended(:);
    sys.output.OutsideSLMSeeds = sum(S_temp) > 0;
    sys.output.F(sys.output.OutsideSLMSeeds,:) = 0;
    sys.output.spikes(sys.output.OutsideSLMSeeds,:) = 0;
end
    
% sys.output.F(end,:) = 0;
sys.opts.dr = ['/Volumes' sys.opts.dr(18:end)];
hF = [];
opts.skipFrames = 2;  tooltips.skipFrames = 'Draw every Nth frame';
opts.appendname = ''; tooltips.skipFrames = 'Append this to the movie filename. If empty, solver params will be appended';
opts.inVivo = false;  tooltips.inVivo = 'Abbas'' flag for in vivo data';

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
if ~isfield(sys.opts,'StimOnset')
    StimOnset = [];
else
    StimOnset = sys.opts.StimOnset;
end


%% Summary of the Reconstructions
% hF(end+1) = figure;

% if isfield(opts,'inVivo') && opts.inVivo
%     A = sys.output.F./max(sys.output.F,[],2);
%     normalizedSpikes = filter(sys.opts.theta,1,A,[],2);
%     F = (A+cumsum([0;max(A(1:end-1,:),[],2)-min(A(1:end-1,:),[],2)]));
%     plot(F'); hold on; line([1000 1000],[0 max(F(:))]); line([3000 3000],[0 max(F(:))]); 
% end
% hF(end+1) = figure;
% subplot(2,2,[1 2]);
% 
% plot(sys.output.F');

% title('Total Seed Intensities');
% subplot(2,2,3);
% plot(sys.control_params.objective);
% title('Objective function'); xlabel('Iteration')
% ylabel('Objective value')
% recon_error = (sys.input.y-sys.output.PS*(sys.output.F) - sys.output.additive_baseline - sys.output.stim_artifact)./sys.input.y;
% subplot(2,2,4); imagesc(recon_error,[-3 3]);
% title('Reconstruction Error')
% drawnow;


%drawframes = [1:opts.skipFrames:399 400:699 700:opts.skipFrames:1990]; %for 2-spot uncaging
drawframes = 1:opts.skipFrames:(size(sys.output.F,2)-10);
% StimOnset = ceil(StimOnset/opts.skipFrames);

%% Calculating Reconstructed Movies
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [mRaw,mDFF,valid2D] = calc_ZMovies(sys,drawframes,[],StimOnset);
[mRaw,mDFF,valid2D] = calc_ZMovies_Bootstrap(sys,drawframes,[],StimOnset);
% dffMax = max(1, min(5, prctile(mDFF(:), 99.5)));
dffMax =  min(5,prctile(mDFF(:), 99.5));
% dff_threshold = 0.01;
% dffMax = min(5, prctile(mDFF(mDFF>dff_threshold), 99));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gamma correction to better see dim features
factor = sqrt(prctile(mRaw(:), 99.9));
hsvV = min(1, sqrt(mRaw(:))./factor);
hsvS = min(1,mDFF(:)./dffMax);
hsvH = ones(size(mDFF(:)));
% hsvH = max(0,min(1,ZF(:)/1.96));
rgb_image = hsv2rgb([hsvH hsvS hsvV]);
rgb_image = reshape(rgb_image, [size(mRaw) 3]);

% mkdir(sys.opts.dr, 'ReconMovies');
VideoName = [sys.opts.dr filesep 'ReconMovies' filesep sys.opts.filename(1:end-11) 'ZMovie' opts.appendname];
try
    W = VideoWriter(VideoName);
catch
    warning(['!! Could not write video in default directory, writing to PWD instead: ' pwd ' !!'])
    W = VideoWriter([sys.opts.filename(1:end-16) '_ZMovie' opts.appendname]);
end
W.FrameRate = 30;
open(W);

f = figure('pos', [ 1 41 1280 1280]);
ax = axes(f, 'units', 'normalized');
h_im = imshow(0.5*ones([sys.opts.dim 3]), 'parent', ax, 'Border','tight');
colormap(ax, hsv2rgb([ones(100,1), linspace(0,1,100)', 0.8*ones(100,1)]))
h_cb = colorbar(ax);
axis image
set(f, 'resize', 'off')
set(h_cb, 'position', [0.96 0.06 0.01 0.87], 'axislocation', 'in', 'tickdir', 'out', 'ticklength', 0, 'ticks', [0.01 0.99], 'ticklabels', num2str([0 ; dffMax], 3), 'ycolor', 'w', 'fontsize', 14)
h_cb.Title.String = '{\Delta}F/F_0 '; h_cb.Title.Color = 'w'; h_cb.Title.FontSize = 14;
h_cb.Label.Rotation = 0;
h_txt = text(ax,1,1, 'Frame X. t= XXX ms', 'units', 'characters', 'verticalalignment', 'baseline', 'horizontalalignment', 'left', 'color', 'w');
set(h_txt, 'position', [1 1])
h_txt.FontSize = 14;
linelength = (10/(1280/5))/2; %10 microns. The whole FOV is 1280/5 = 256um
hs = annotation('line', [.96-linelength 0.96+linelength], [0.02, 0.02]);
hs.LineWidth = 5;
hs.Color = 'w';
hScalebarLabel = text(ax,0.96,0.03, '10 {\mu}m', 'verticalalignment', 'baseline', 'units', 'normalized', 'horizontalalignment', 'center', 'color', 'w');
hScalebarLabel.FontSize = 14;
for frame = 1:length(drawframes)
    h_txt.String = ['  Frame #' int2str(drawframes(frame)) '       ' num2str((drawframes(frame)-1)/1.016, '%.1f') ' ms'];
    if mod(drawframes(frame), sys.solver_params.verbose)<opts.skipFrames
        disp(['Writing frame: ' int2str(drawframes(frame))])
    end
    %sys.control_params.frame = frame;
    imFrame = min(1,repmat(sqrt(reshape(sum(sys.input.ref_image,3),[],1)), 1, 3)./factor);
    imFrame(valid2D, :) = squeeze(rgb_image(:,frame,:));
%     imFrame(~valid2D,3) = imFrame(~valid2D,3) + 0.15;
    imFrame(~valid2D,3) = imFrame(~valid2D,3)* 0.6 + 0.05;
    imFrame(~valid2D,1:2) = imFrame(~valid2D,1:2) * 0.5;
    
    %imFrame= imFrame.*repmat(sqrt(sys.SLMmask(:)), 1, 3);  %black out the background
    F = reshape(imFrame, [sys.opts.dim 3]);
    set(h_im, 'Cdata', F);
    writeVideo(W, getframe(ax));
end
disp('Done writing movie.');
disp(['Video Location:' VideoName]);
close(W)

end

