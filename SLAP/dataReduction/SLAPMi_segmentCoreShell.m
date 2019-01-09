function S = SLAPMi_segmentCoreShell(refIM,optsin)
%Segment a SLAPMi reference image that has already been classified
%    (or provide your own pixel labels for segmentation via input opts)

opts.distThreshUM = 1.4;
opts.distThreshSpineUM = (0.8/1.4)*opts.distThreshUM;
opts.maxseeds = Inf;
opts.visualize = false;
if nargin>1 && ~isempty(optsin)
    for field = fieldnames(optsin)'
        opts.(field{1}) = optsin.(field{1});
    end
end

%READ DATA
if (~nargin || isempty(refIM))
    if ~isfield(opts, 'fn')
        [fn, dr] = uigetfile('*.mat', 'Select refIM_ILSTK data');
        opts.fn = [dr fn];
    end
    [dr, fn] = fileparts(opts.fn); dr = [dr filesep];
    load([dr fn]);
end
if isfield(opts, 'labels') && ~isempty(opts.labels)
    refIM.labels = opts.labels;
end

if max(refIM.labels(:))<4
    %this is not a core-shell segmentation!
    keyboard
else
    opts.XYscale = refIM.metadata.objectiveResolution * median(diff(refIM.M.coords.X)); %pixel size in microns
    refIM.IM(isnan(refIM.IM)) = 0;
    
    S = segment3D(refIM,opts);
    %S = segmentCoreShell3D(refIM,XYscale, opts);
end
%visualize S
if opts.visualize
    if isfield(opts, 'savefile')
        visualize_S(S, opts.savefile);
    else
        visualize_S(S);
    end
end