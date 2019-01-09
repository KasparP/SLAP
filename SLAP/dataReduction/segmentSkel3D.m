function S = segmentSkel3D (bw, image, XYscale, optsin)
%Segment a 3D image by skeletonization

%skeletonize each plane
%starting with the most prominent skeleton points (i.e. brightest in a
%tophat filtered image):
%remove nearby points in xy
%if the planes directly above or below are in the mask,
%place points there as part of this seed,
%and remove enarby points in xy
opts.mask = ones(size(image));
opts.scramble = false;
opts.maxseeds = Inf; %maximum number of elements
opts.sharpness = 1.1; %how sharp the edges between seeds should be; a nonnegative number. Higher is sharper
opts.dist_thresh = round(1.8/XYscale);  %minimum distance between seeds, in pixels: typical is 1.8 um
if nargin>3 %UPDATE WITH USER-SPECIFIED OPTIONS
    for field = fieldnames(optsin)'
        opts.(field{1}) = optsin.(field{1});
    end
end

h = fspecial('gaussian', ceil(2/XYscale), 0.5/XYscale);

i1 = cell(size(bw,3),1);
i2 = cell(size(bw,3),1);
i3 = cell(size(bw,3),1);
B =  cell(size(bw,3),1);

image = double(image);
image_masked = image.*opts.mask;
for z  = 1:size(bw,3)
    I_gauss = image_masked(:,:,z);
    %I_gauss = imfilter(image_masked(:,:,z), h);%filtered image will be used to find centers of dendrites
    I_gauss = imtophat(I_gauss,strel('disk',ceil(1/XYscale))); %radius of 1 micron
    
    CC = bwconncomp(bw(:,:,z));
    bwZ = false(size(I_gauss));
    for Lix = 1:length(CC.PixelIdxList) %for each object in bw
        sel = false(size(I_gauss));
        sel(CC.PixelIdxList{Lix}) = true;
        imvals = I_gauss(sel);
        threshN = max(20, length(imvals)*0.5);
        CC1 = bwconncomp(sel); CC2 = CC1;
        while CC1.NumObjects==1 && sum(sel(:))>threshN
            CC2 = CC1;
            [~,minix] = min(imvals);
            sel(CC.PixelIdxList{Lix}(minix)) = false;
            imvals(minix) = inf;
            CC1 = bwconncomp(sel);
        end
        bwZ(CC2.PixelIdxList{1}) = true;
    end
    
    %skeletonize
    skel = bwmorph(bw(:,:,z),'skel',Inf);
    
    %break circles in the skeleton
    holes = imfill(skel,4,'holes') & (~skel);
    while(any(holes(:)))
        [L,num] = bwlabel(holes); SE = strel('diamond',1);
        for Lix = 1:num %for each hole
            %get pixels on perimeter
            [is,js] = find(imdilate(L==Lix, SE) & skel);
            %remove the dimmest one
            vs = I_gauss(sub2ind(size(L),is,js));
            [~, minix] = min(vs);
            skel(is(minix),js(minix)) = false;
        end
        holes = imfill(skel,4,'holes') & (~skel);
    end

    ends = bwmorph(skel, 'endpoints'); %TODO: prioritize endpoints for segmentation?
    [i1{z},i2{z}] = find(skel);
    i3{z} = z*ones(size(i2{z}));
    B{z} = I_gauss(sub2ind(size(bw(:,:,z)),i1{z},i2{z})); %brightness of skeleton points; we'll keep the brightest ones
end

disp('Segmenting image...')
i1 = cell2mat(i1); i2 = cell2mat(i2); i3 = cell2mat(i3);
B = cell2mat(B);
[B, sortorder] = sort(B, 'descend');
B0 = B(ceil(length(B)/3));
if opts.scramble %temporary code for simulating alternate segmentations
    perm = randperm(length(sortorder));
    B  = B(perm);
    sortorder = sortorder(perm);
end
i1 = i1(sortorder);
i2 = i2(sortorder);
i3 = i3(sortorder);

repeat=true;
while repeat
added = cell(size(i3));
exclude = false(length(B),1);
for ix = 1:length(B) %for every point, in order of prominence
    if ~exclude(ix)
        dist_thresh = opts.dist_thresh .* max(1, min(2, sqrt(B0/B(ix))));
        
        N1 = i3==i3(ix); %N1: neighbours in this plane in the original skeletonization
        
        %flood-fill to exclude nearby points along the skeleton within this plane
        N2 = sqrt((i2(ix)-i2(N1)).^2 + (i1(ix)-i1(N1)).^2)<dist_thresh; %of the candidate points to remove, only remove the ones within the distance threshold
        i1N1 = i1(N1); i2N1 = i2(N1);
        fillDist = squareform(pdist([i2N1(N2) i1N1(N2)]));
        filled = (i2N1(N2)==i2(ix) & i1N1(N2)==i1(ix))';
        for iter = 1:ceil(dist_thresh)
            filled = filled | any(fillDist(filled,:)<1.5,1); %flood fill
        end
        N2(N2) = filled;
        N1(N1) = N2; 
        
        for z2 = i3(ix)-1:-1:max(1,i3(ix)-2) %extend the segment vertically within the mask
            if bw(i1(ix),i2(ix),z2)
                N3 = i3==z2;
                filled = sqrt((i2(ix)-i2(N3)).^2 + (i1(ix)-i1(N3)).^2)<dist_thresh;
                N3(N3) = filled;
                N1 = N1 | N3;
                %N1 = N1 | i3==z2;
                added{ix} = [added{ix} z2];
            else
                break
            end
        end
        for z2 = i3(ix)+1:min(size(bw,3), i3(ix)+2) %extend the segment vertically within the mask
            if bw(i1(ix),i2(ix),z2)
                N3 = i3==z2;
                filled = sqrt((i2(ix)-i2(N3)).^2 + (i1(ix)-i1(N3)).^2)<dist_thresh;
                N3(N3) = filled;
                N1 = N1 | N3;
                %N1 = N1 | i3==z2;
                added{ix} = [added{ix} z2];
            else
                break
            end
        end

        N1(1:ix) = false; %don't remove points we've previously assigned
%         N2 = sqrt((i2(ix)-i2(N1)).^2 + (i1(ix)-i1(N1)).^2)<dist_thresh; %of the candidate points to remove, only remove the ones within the distance threshold
%         N1(N1) = N2;
        
        exclude = exclude | N1;
    end
end
numPointsInMask =  sum(opts.mask(sub2ind(size(opts.mask), i1(~exclude),i2(~exclude)))>0.5);
if numPointsInMask>opts.maxseeds
    repeat = true;
    opts.dist_thresh = max(opts.dist_thresh+1, round(opts.dist_thresh*(opts.maxseeds/numPointsInMask)));
    disp(['Too many seeds were selected, increasing distance to ' num2str(opts.dist_thresh*XYscale, 4) ' um'])
else
    repeat = false;
    disp(['Number of seeds within mask: ' int2str(numPointsInMask)])
end
end
i1 = i1(~exclude); i2 = i2(~exclude); i3 = i3(~exclude); added = added(~exclude);

%we are going to calculate influence of seeds over nearby pixels with a
%'flood fill' operation
nh_size = 3*opts.dist_thresh; %maximum radius of seed influence, in pixels. This should be >2*dist_thresh

S.bw = bw;
I = cell(4*length(i1),1); J = I; V = I;
cellIX = 1;
for Z = 1:size(bw,3)
    i_ixs = find(i3==Z  |  cellfun(@(x)any(x==Z), added));
    bw_padded = false([size(bw,1), size(bw,2)] + 2*nh_size);
    bw_padded(nh_size+1:end-nh_size, nh_size+1:end-nh_size) = bw(:,:,Z);
    
    nhood = zeros(2*nh_size+1,2*nh_size+1,length(i_ixs));
    influence = zeros(2*nh_size+1,2*nh_size+1,length(i_ixs));
    influence(nh_size+1,nh_size+1,:) = 1;
    for seed = 1:length(i_ixs)
        nhood(:,:,seed) = bw_padded((0:2*nh_size)+i1(i_ixs(seed)), (0:2*nh_size)+i2(i_ixs(seed)));
    end
    
    SE = [0.5 1 0.5; 1 1 1 ; 0.5 1 0.5];
    %Flood fill:
    for i = 1:(3*opts.dist_thresh+1)
        influence = influence + (convn(influence, SE, 'same') .* nhood);
    end
    influence = influence./sum(sum(influence,1),2); %so seeds at edges don't get swamped
    
    %accumulate indices for sparse matrix, #voxels x #seeds
    for seed = 1:length(i_ixs)
        [row, col, val] = find(influence(:,:,seed));
        seg_ixs = sub2ind(size(bw), row+i1(i_ixs(seed))-(nh_size+1),col+i2(i_ixs(seed))-(nh_size+1),Z*ones(size(col)));
        %accumulators for sparse matrix:
        I{cellIX} = seg_ixs; J{cellIX} = i_ixs(seed)*ones(size(seg_ixs)); V{cellIX} = val; cellIX = cellIX+1;
    end
end

%error checking: make sure every pixel is influenced by a seed
no_seed = bw;
no_seed(unique(cell2mat(I))) = false; 
num = 0;
if any(no_seed(:))
    no_seed = bwareaopen(no_seed, 4, 4);%discard isolated pixels
    [L,num] = bwlabeln(no_seed);
    for Lix = 1:num
        I{cellIX} = find(L(:)==Lix); J{cellIX} = (length(i1)+Lix)*ones(size(I{cellIX})); V{cellIX} = ones(size(I{cellIX}));
        cellIX = cellIX+1;
    end
end
%I:pixel J:seed V:value
i = cell2mat(I); j = cell2mat(J); v = cell2mat(V);
%apply boundary sharpening
v = v.^opts.sharpness;
sumIM = sum(sparse(i,j,v, numel(bw), length(i1)+num),2);
select = v>(0.05.*sumIM(i));
i = i(select); j = j(select); v=v(select);
v = v./sumIM(i);
S.seg = sparse(i,j,v.*image(i), numel(bw), length(i1)+num);
segMasked = sparse(i,j,v.*image_masked(i), numel(bw), length(i1)+num);

%discard seeds with intensity less than THRESH times the average
disp('Fusing dim seeds...')
thresh = 0.2*max(1, size(S.seg,2)/1000);
intensities = full(sum(segMasked,1)); intensities = intensities./mean(intensities);
[sorted, sortorder] = sort(intensities, 'ascend');

%don't fuse the seeds on the border of the SLM, these should be kept small
%to better deal with model mismatch and motion.
maskmask = opts.mask(:)>0.8;
onBorder = any(S.seg(maskmask,:),1) & any(S.seg(~maskmask,:),1);

for int_ix = 1:find(sorted<thresh,1,'last')
    if intensities(sortorder(int_ix))<thresh && ~onBorder(sortorder(int_ix))
        %find the seed with highest overlap with the current seed
        support = any(S.seg((S.seg(:,sortorder(int_ix))>0),:), 1);
        support(sortorder(int_ix)) = false;
        support(intensities>2*thresh) = false;
        support(onBorder) = false;
        support = find(support);
        if ~isempty(support)
            [~,minix] = min(intensities(support)); %merge to the dimmest neighbour
            minix = support(minix);
            
            S.seg(:,minix) = S.seg(:,minix)+S.seg(:, sortorder(int_ix));
            S.seg(:, sortorder(int_ix)) = 0;
        elseif intensities(sortorder(int_ix))<thresh/2
            S.seg(:, sortorder(int_ix)) = 0;
        end
    end
end

disp('Done segmentation.');
end