function S = fuseSeeds(S, optsin)
%reduce the number of seeds in the segmentation, to make sure no seeds are
%much dimmer than the average, and fuse highly overlapping seeds

disp('Fusing dim seeds...')
opts.fuseThresh = 0.5/1000; %threshold, as fraction of mean intensity, below which to call a seed dim; should be in range (0.1-1)/1000
if nargin>1 %UPDATE WITH USER-SPECIFIED OPTIONS
    for field = fieldnames(optsin)'
         opts.(field{1}) = optsin.(field{1});
    end
else

npx = size(S.IM,1)*size(S.IM,2);
nZ = size(S.seg,1)/npx;

mask3 = repmat(S.mask, 1,1,nZ);
S.fusedInto = 1:size(S.seg,2);

%sort seeds by intensity
segMasked = S.seg .* mask3(:);
intensities = full(sum(segMasked,1));
S.fusedInto(intensities==0) = 0;
[sorted, sortorder] = sort(intensities, 'ascend');

%2D projection of segmentation for measuring overlap
seg2D =S.seg(1:npx,:);
for Z = 2:nZ
    seg2D = seg2D+ S.seg((Z-1)*npx + (1:npx),:);
end

%merge seeds that overlap significantly (in 2D!)
nfused = 0;
for S1 = sortorder %start with the dimmest seeds
    if S.fusedInto(S1)==0
        continue
    end
    iS1 = sum(seg2D(:,S1));
    %find the seed that captures the most intensity of this one
    overlap = nan(1,size(seg2D,2));
    for S2 = 1:size(seg2D,2)
        if S2==S1
            continue
        end
        iBoth = sum(min(seg2D(:,S1), seg2D(:,S2)));
        overlap(S2) = iBoth./iS1;
    end
    [maxval,maxix] = nanmax(overlap);
    if maxval>0.5
        nfused = nfused+1;
        S.seg(:,maxix) = S.seg(:,maxix)+S.seg(:, S1);
        S.seg(:, S1) = 0;
        seg2D(:,maxix) = seg2D(:,maxix)+seg2D(:, S1);
        seg2D(:,S1) = 0;
        S.fusedInto(S.fusedInto==S1) = maxix;
%     elseif maxval >0.1
%         figure, imshow(cat(3, reshape(full(seg2D(:,maxix)), size(S.IM,1), size(S.IM,2)), reshape(full(seg2D(:,S1)), size(S.IM,1), size(S.IM,2)), reshape(full(seg2D(:,S1)), size(S.IM,1), size(S.IM,2)))/1000);
    end
end
disp(['Fused ' int2str(nfused) ' overlapping seeds']);

%re-sort seeds by intensity
segMasked = S.seg .* mask3(:);
intensities = full(sum(segMasked,1));

%discard seeds with intensity less than THRESH times the average
%numValid = sum(sum(S.seg,1)>0);
intensities = intensities./sum(intensities); %normalize to mean
[sorted, sortorder] = sort(intensities, 'ascend');

%don't fuse the seeds on the border of the SLM, these should be kept small
%to better deal with model mismatch and motion.
maskmask = S.mask(:)>0.8;
onBorder = full(any(seg2D(maskmask,:),1) & any(seg2D(~maskmask,:),1));
nfused = 0; ndiscard = 0;
for int_ix = 1:find(sorted<opts.fuseThresh,1,'last')
    fuse_ix = sortorder(int_ix);
    if intensities(fuse_ix)>0 && intensities(fuse_ix)<opts.fuseThresh && ~onBorder(fuse_ix)
        %find the seed with highest overlap with the current seed
        support = any(seg2D((seg2D(:,fuse_ix)>0),:), 1);
        support(fuse_ix) = false;
        support(intensities>2*opts.fuseThresh) = false;
        support(onBorder) = false;
        support = find(support);
        if ~isempty(support)
            nfused = nfused+1;
            [~,minix] = min(intensities(support)); %merge to the dimmest neighbour
            minix = support(minix);
            
%           figure, imshow(cat(3, reshape(full(seg2D(:,minix)), size(S.IM,1), size(S.IM,2)), reshape(full(seg2D(:,fuse_ix)), size(S.IM,1), size(S.IM,2)), reshape(full(seg2D(:,fuse_ix)), size(S.IM,1), size(S.IM,2)))/1000); 
            S.seg(:,minix) = S.seg(:,minix)+S.seg(:, fuse_ix);
            S.seg(:, fuse_ix) = 0;
            seg2D(:,minix) = seg2D(:,minix)+seg2D(:, fuse_ix);
            seg2D(:,fuse_ix) = 0;
            
            S.fusedInto(S.fusedInto==fuse_ix) = minix;
            
        elseif intensities(sortorder(int_ix))<opts.fuseThresh/2
            ndiscard = ndiscard+1;
            S.seg(:, fuse_ix) = 0;
            seg2D(:, fuse_ix) = 0;
            S.fusedInto(S.fusedInto==fuse_ix) = 0;
        end
    end
end

disp(['Fused ' int2str(nfused) ' dim seeds']);
disp(['Discarded ' int2str(ndiscard) ' dim seeds']);
end