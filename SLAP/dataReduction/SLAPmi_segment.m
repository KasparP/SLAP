function S = SLAPMi_segment(refIM,optsin)
%Segment a SLAPMi reference image that has already been classified
%    (or provide your own pixel labels for segmentation via input opts)

opts = [];
if nargin>1
    opts = optsin;
end

%READ DATA
if (~nargin || isempty(refIM))
    if ~isfield(optsin, 'fn')
        [fn, dr] = uigetfile('*.mat', 'Select refIM_ILSTK data');
        optsin.fn = [dr fn];
    end
    [dr, fn] = fileparts(optsin.fn); dr = [dr filesep];
    load([dr fn]);
end
if nargin>1 && isfield(optsin, 'labels') && ~isempty(optsin.labels)
    refIM.labels = optsin.labels;
end

if max(refIM.labels(:))>3
    %segment by labels; the segmentation was already performed, e.g. in Takashi's data
    keyboard
else
    
    XYscale = refIM.metadata.objectiveResolution * median(diff(refIM.M.coords.X)); %pixel size in microns
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
    
    %OPTIONS for segmentation
    if ~isfield(opts, 'maxseeds')
        opts.maxseeds = 1000; %maximum number of seeds per plane
    end
    opts.mask = refIM.mask;
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
end
%visualize S
if nargin && isfield(opts, 'visualize') && opts.visualize
    if isfield(opts, 'savefile')
        visualize_S(S, opts.savefile);
    else
        visualize_S(S);
    end
end