function S = SLAPMi_segmentOnce(refIM,optsin)
%Segment a SLAPMi reference image that has already been classified
%    (or provide your own pixel labels for segmentation via input opts)

opts.distThreshUM = 1.5;
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

if max(refIM.labels(:))>3
    %segment by labels; the segmentation was already performed, e.g. in Takashi's data
    keyboard
else
    
    XYscale = refIM.metadata.objectiveResolution * median(diff(refIM.M.coords.X)); %pixel size in microns
    opts.dist_thresh = round(opts.distThreshUM/XYscale);
    refIM.IM(isnan(refIM.IM)) = 0;
    
    SS = strel('diamond', 1);
    SSbig = strel('diamond', 5);
    bwSoma = refIM.labels==3;
    bwSoma = imclose(imopen(bwSoma, SSbig), SSbig);
    bwSoma = bwareaopen(bwSoma, 400);
    bwSoma = imfill(bwSoma, 4, 'holes');
    refIM.labels(bwSoma) = 3; %reclassify anything set to soma
    refIM.labels(refIM.labels==3 & ~bwSoma) = 2; %reclassify rejected somata as dendrite
    bw = refIM.labels==2; %2 is the label for dendrites;  1=outside 3=soma
    bw = imdilate(bw, SS);
    bw = bwareaopen(bw, 20, 4);%remove regions of size less than X pixels
    bw = imerode(bw, SS);
    bw = bwareaopen(bw, 5, 4);%remove regions of size less than X pixels
    refIM.labels(refIM.labels==2 & ~bw) = 1; %reclassify rejected dendrite as outside
    
    seg3D = segmentSkel3D(bw,refIM.IM,XYscale, opts);
    nS = size(seg3D.seg,2);
    S = refIM;
    S.bw = seg3D.bw;
    
    %add in somata
    [L, nSoma] = bwlabeln(bwSoma);
    if nSoma
        %i,pixels %j,seeds, %v, values
        i = cell(nSoma,1); j = i; v = i;
        for Lix = 1:nSoma
            i{Lix} = find(L==Lix);
            j{Lix} = ones(size(i{Lix}))*(nS+Lix);
            v{Lix} = refIM.IM(L==Lix);
            if any(isnan(v{Lix}))
                error('The segmentation matrix contains nans');
            end
        end
        
        [I,J,V] = find(seg3D.seg);
        S.seg = sparse([I ; cell2mat(i)], [J ; cell2mat(j)], [V ; cell2mat(v)], size(seg3D.seg,1), nS+nSoma);
    else
        S.seg = seg3D.seg;
    end
    clear seg3D;
    
    S.segLabel = [2*ones(1, nS) 3*ones(1,nSoma)];   
    
    %remove '0' segments
    select = any(S.seg,1);
    S.seg = S.seg(:,select);
    
    %ensure that the segments add up to the reference image within the support
    valid = full(any(S.seg,2));
    S.seg(valid,:) = S.seg(valid,:).*repmat(S.IM(valid)./(sum(S.seg(valid,:),2)), 1, size(S.seg,2));
end
%visualize S
if opts.visualize
    if isfield(opts, 'savefile')
        visualize_S(S, opts.savefile);
    else
        visualize_S(S);
    end
end