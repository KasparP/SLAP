function refIM = segmentCoreShell3D (refIM, XYscale, optsin)
%Segment a 3D image by skeletonization

%skeletonize each plane
%starting with the most prominent skeleton points (i.e. brightest in a
%tophat filtered image):
%remove nearby points in xy
%if the planes directly above or below are in the mask,
%place points there as part of this seed,
%and remove enarby points in xy
image = refIM.IM;
image(isnan(image)) = 0;

opts.mask = ones(size(image));
opts.maxseeds = Inf; %maximum number of elements
opts.sharpness = 1.1; %how sharp the edges between seeds should be; a nonnegative number. Higher is sharper
opts.distThreshUM = 1.5;  %minimum distance between shaft seeds, in pixels: typical is 1.8 um
opts.distThreshSpineUM = 0.8; %minimum lateral distance between a spine and another spine in the plane above or below it, which doesn't also appear in the same plane
if nargin>3 %UPDATE WITH USER-SPECIFIED OPTIONS
    for field = fieldnames(optsin)'
        opts.(field{1}) = optsin.(field{1});
    end
end

opts.distThresh = round(opts.distThreshUM/XYscale);
opts.distThreshSpine = round(opts.distThreshSpineUM/XYscale);

%sometimes Ilastik can make gross errors at edges, which we don't care about
%anyways
refIM.labels([1:8 end-7:end],:,:) = 0;
refIM.labels(:,[1:8 end-7:end],:) = 0;

bw = refIM.labels>1;
SS = strel('diamond', 1);

bw = imerode(bw, SS);
bw = imdilate(bw, SS);
bw = bwareaopen(bw, 16, 4);%remove regions of size less than X pixels

core = refIM.labels == 2;
%shell = refIM.labels== 3;
shaft = refIM.labels==4;
bw = bw|core|shaft;
refIM.bw = bw;

%we will assign seeds to spine heads (capital indexes) and shafts
%(lowercase) separately:
i1 = cell(size(bw,3),1); I1 = cell(size(bw,3),1);
i2 = cell(size(bw,3),1); I2 = cell(size(bw,3),1);
i3 = cell(size(bw,3),1); I3 = cell(size(bw,3),1);
b =  cell(size(bw,3),1); B =  cell(size(bw,3),1);

image = double(image);
for z  = 1:size(bw,3)
    I_gauss = image(:,:,z);
    %I_gauss = imfilter(image(:,:,z), h);%filtered image will be used to find centers of dendrites
    I_gauss = imtophat(I_gauss,strel('disk',ceil(1/XYscale))); %background subtract with radius of 1 micron
    
    %add the spine heads
    CoreZ = core(:,:,z);
    skelCoreZ = false(size(core(:,:,z)));
    SEveto = strel('disk', 4);
    CC = bwconncomp(core(:,:,z));
    [~,sortorder] = sort(cellfun(@(x)(numel(x)), CC.PixelIdxList), 'descend');
    for Lix = sortorder %for each spine head
        if any(CoreZ(CC.PixelIdxList{Lix}))
            ints = I_gauss(CC.PixelIdxList{Lix});
            [~,maxix] = max(ints);
            skelCoreZ(CC.PixelIdxList{Lix}(maxix)) = true;
            
            veto = false(size(CoreZ));
            veto(CC.PixelIdxList{Lix}) = true;
            veto = imdilate(veto, SEveto);
            CoreZ(veto) = false;
        end
    end
    
    %skelShaftZ =  bwmorph(shaft(:,:,z),'skel',Inf);
    skelShaftZ =  bwmorph(bw(:,:,z),'skel',Inf);
    
    [I1{z},I2{z}] = find(skelCoreZ);
    I3{z} = z*ones(size(I2{z}));
    B{z} = I_gauss(sub2ind(size(bw(:,:,z)),I1{z},I2{z})); %brightness of skeleton points; we'll keep the brightest ones
    
    [i1{z},i2{z}] = find(skelShaftZ);
    i3{z} = z*ones(size(i2{z}));
    b{z} = I_gauss(sub2ind(size(bw(:,:,z)),i1{z},i2{z})); %brightness of skeleton points; we'll keep the brightest ones
end

disp('Segmenting image...')
I1 = cell2mat(I1); I2 = cell2mat(I2); I3 = cell2mat(I3);
B = cell2mat(B);
[B, sortorderB] = sort(B, 'descend');

i1 = cell2mat(i1); i2 = cell2mat(i2); i3 = cell2mat(i3);
b = cell2mat(b);
[b, sortorderb] = sort(b, 'descend');
b0 = b(ceil(length(b)/3));

i1 = i1(sortorderb); I1 = I1(sortorderB);
i2 = i2(sortorderb); I2 = I2(sortorderB);
i3 = i3(sortorderb); I3 = I3(sortorderB);

exclude = false(length(b),1);
EXCLUDE = false(length(B),1); %Capitals represent the core skeletonization, lowercase is the rest
ADDED = cell(size(I3));
ADDEDTF = false(length(I3), size(refIM.IM,3));
for ix = 1:length(B)
    if ~EXCLUDE(ix)
        %exclude nearby points in the shaft skeletonization
        n1 = i3==I3(ix) & sqrt((I2(ix)-i2).^2 + (I1(ix)-i1).^2)<opts.distThreshSpine;
        N1 = I3==I3(ix) & sqrt((I2(ix)-I2).^2 + (I1(ix)-I1).^2)<opts.distThreshSpine;
        N1(1:ix) = false;
        exclude = exclude | n1;
        EXCLUDE = EXCLUDE | N1;
        
        for z2 = I3(ix)-1:-1:max(1,I3(ix)-2) %extend the segment vertically within the mask
            n4 = I3==z2; n4(ix:end) = false;
            n4 = n4 | ADDEDTF(:,z2); %potential neighbours that veto the extension
            mindist = min(sqrt((I2(ix)-I2(n4)).^2 + (I1(ix)-I1(n4)).^2)); %minimal distance to an already-added point
            if (isempty(mindist) || mindist>opts.distThreshSpine) && bw(I1(ix),I2(ix),z2)
                ADDED{ix} = [ADDED{ix} z2];
                ADDEDTF(ix,z2) = true;
                
                n2 = i3==z2;
                filled = sqrt((I2(ix)-i2(n2)).^2 + (I1(ix)-i1(n2)).^2)<opts.distThreshSpine;
                n2(n2) = filled;
                n1 = n1 | n2;
                
                N2 = I3==z2;
                filled = sqrt((I2(ix)-I2(N2)).^2 + (I1(ix)-I1(N2)).^2)<opts.distThreshSpine;
                N2(N2) = filled;
                N2(1:ix) = false;
                EXCLUDE = EXCLUDE | N2;
            else
                break
            end
        end
        
        for z2 = I3(ix)+1:min(size(bw,3), I3(ix)+2) %extend the segment vertically within the mask
            n4 = I3==z2; n4(ix:end) = false;
            n4 = n4 | ADDEDTF(:,z2); %potential neighbours that veto the extension
            mindist = min(sqrt((I2(ix)-I2(n4)).^2 + (I1(ix)-I1(n4)).^2)); %minimal distance to an already-added point
            if (isempty(mindist) || mindist>opts.distThreshSpine) && bw(I1(ix),I2(ix),z2)
                ADDED{ix} = [ADDED{ix} z2];
                ADDEDTF(ix,z2) = true;
                
                n2 = i3==z2;
                filled = sqrt((I2(ix)-i2(n2)).^2 + (I1(ix)-i1(n2)).^2)<opts.distThreshSpine;
                n2(n2) = filled;
                n1 = n1 | n2;
                
                N2 = I3==z2;
                filled = sqrt((I2(ix)-I2(N2)).^2 + (I1(ix)-I1(N2)).^2)<opts.distThreshSpine;
                N2(N2) = filled;
                N2(1:ix) = false;
                EXCLUDE = EXCLUDE | N2;
            else
                break
            end
        end
        n1(1:ix) = false; %don't remove points we've previously assigned
        exclude = exclude | n1;
    end
end

% for z = 1:max(I3)
%     select = (I3==z | cellfun(@(x)(any(x==z)), ADDED)) & ~EXCLUDE;
%     D = squareform(pdist([I1(select) I2(select)]));
%     D(logical(eye(size(D)))) = inf;
%     if any(D<opts.distThreshSpine)
%         keyboard
%     end
% end


added = cell(size(i3));
addedTF = false(length(i3), size(refIM.IM,3));
for ix = 1:length(b) %for every point, in order of prominence
    if ~exclude(ix)
        dist_thresh = opts.distThresh;
        %dist_thresh = opts.distThresh .* max(1, min(2, sqrt(b0/b(ix))));
        
        n1 = i3==i3(ix); %n1: neighbours in this plane in the original skeletonization
        %exclude nearby points within this plane
        n2 = sqrt((i2(ix)-i2(n1)).^2 + (i1(ix)-i1(n1)).^2)<dist_thresh;
        n1(n1) = n2; 
        
        for z2 = i3(ix)-1:-1:max(1,i3(ix)-2) %extend the segment vertically within the mask
            n4 = i3==z2; n4(ix:end) = false;
            n4 = n4 | addedTF(:,z2);
            N4 = I3==z2 | ADDEDTF(:,z2);
            mindist = min([inf ; sqrt((i2(ix)-i2(n4)).^2 + (i1(ix)-i1(n4)).^2)]);
            MINDIST = min([inf ; sqrt((i2(ix)-I2(N4)).^2 + (i1(ix)-I1(N4)).^2)]);
            if (MINDIST>opts.distThreshSpine && mindist>opts.distThresh) && bw(i1(ix),i2(ix),z2)
                n3 = i3==z2;
                filled = sqrt((i2(ix)-i2(n3)).^2 + (i1(ix)-i1(n3)).^2)<dist_thresh;
                n3(n3) = filled;
                n1 = n1 | n3; %add these to the excluded list
                added{ix} = [added{ix} z2];
                addedTF(ix,z2) = true;
            else
                break
            end
        end
        for z2 = i3(ix)+1:min(size(bw,3), i3(ix)+2) %extend the segment vertically within the mask
            n4 = i3==z2; n4(ix:end) = false;
            n4 = n4 | addedTF(:,z2);
            N4 = I3==z2 | ADDEDTF(:,z2);
            mindist = min([inf ; sqrt((i2(ix)-i2(n4)).^2 + (i1(ix)-i1(n4)).^2)]);
            MINDIST = min([inf ; sqrt((i2(ix)-I2(N4)).^2 + (i1(ix)-I1(N4)).^2)]);
            if (MINDIST>opts.distThreshSpine && mindist>opts.distThresh) && bw(i1(ix),i2(ix),z2)
                n3 = i3==z2;
                filled = sqrt((i2(ix)-i2(n3)).^2 + (i1(ix)-i1(n3)).^2)<dist_thresh;
                n3(n3) = filled;
                n1 = n1 | n3; %add these to the excluded list
                added{ix} = [added{ix} z2];
                addedTF(ix,z2) = true;
            else
                break
            end
        end

        n1(1:ix) = false; %don't remove points we've previously assigned
        exclude = exclude | n1;
    end
end

%merge the two sets of points
I1 = I1(~EXCLUDE); I2 = I2(~EXCLUDE); I3 = I3(~EXCLUDE); ADDED = ADDED(~EXCLUDE);
i1 = i1(~exclude); i2 = i2(~exclude); i3 = i3(~exclude); added = added(~exclude);

%seglabel = [ones(length(I1),1) ; 2*ones(length(i1,1))];
i1 = [I1;i1]; i2 = [I2;i2]; i3 = [I3;i3]; added = [ADDED;added];

% for z = 1:max(i3)
%     select = i3==z | cellfun(@(x)(any(x==z)), added);
%     D = squareform(pdist([i1(select) i2(select)]));
%     D(logical(eye(size(D)))) = inf;
%     if any(D(:)<opts.distThreshSpine)
%         keyboard
%     end
% end

%we are going to calculate influence of seeds over nearby pixels with a
%'flood fill' operation
nh_size = 3*opts.distThresh; %maximum radius of seed influence, in pixels. This should be >2*dist_thresh
I = cell(4*length(i1),1); J = I; V = I;
cellIX = 1;
for Z = 1:size(bw,3)
    i_ixs = find(i3==Z  |  cellfun(@(x)any(x==Z), added));
    bw_padded = false([size(bw,1), size(bw,2)] + 2*nh_size);
    bw_padded(nh_size+1:end-nh_size, nh_size+1:end-nh_size) = bw(:,:,Z);
    IMpadded = zeros([size(bw,1), size(bw,2)] + 2*nh_size);
    IMpadded(nh_size+1:end-nh_size, nh_size+1:end-nh_size) = refIM.IM(:,:,Z);
    
    nhood = false(2*nh_size+1,2*nh_size+1,length(i_ixs));
    nhoodIM = zeros(2*nh_size+1,2*nh_size+1,length(i_ixs));
    for seed = 1:length(i_ixs)
        nhood(:,:,seed) = bw_padded((0:2*nh_size)+i1(i_ixs(seed)), (0:2*nh_size)+i2(i_ixs(seed)));
        nhoodIM(:,:,seed) = IMpadded((0:2*nh_size)+i1(i_ixs(seed)), (0:2*nh_size)+i2(i_ixs(seed)));
    end
    nhoodIM = nhoodIM./repmat(max(max(nhoodIM,[],1),[],2), size(nhood,1),size(nhood,2),1);
    influence = zeros(2*nh_size+1,2*nh_size+1,length(i_ixs));
    influence(nh_size+1,nh_size+1,:) = nhoodIM(nh_size+1,nh_size+1,:);
    
    SE = ones(3);
    %Flood fill:
    for i = 1:(3*opts.distThresh+1)
        influence = (influence + convn(influence, SE, 'same')).*nhoodIM.* nhood;
    end
    influence = influence./sum(sum(influence,1),2); %so seeds at edges don't get swamped
    influence(isnan(influence)) = 0;
    
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
refIM.seg = sparse(i,j,v.*image(i), numel(bw), length(i1)+num);
segMasked = sparse(i,j,v.*image(i), numel(bw), length(i1)+num);

%discard seeds with intensity less than THRESH times the average
disp('Fusing dim seeds...')
thresh = 0.2*max(1, size(refIM.seg,2)/1000);
intensities = full(sum(segMasked,1)); intensities = intensities./mean(intensities);
[sorted, sortorderb] = sort(intensities, 'ascend');

%don't fuse the seeds on the border of the SLM, these should be kept small
%to better deal with model mismatch and motion.
maskmask = opts.mask(:)>0.8;
onBorder = any(refIM.seg(maskmask,:),1) & any(refIM.seg(~maskmask,:),1);
nPixels = sum(refIM.seg>(refIM.IM(:)/2), 1); npxThr =50; 


%fuse seeds if they are below the intensity threshold, if they have fewer
%than XX pixels, or if they have too many neighbours
for int_ix = 1:find(sorted<thresh,1,'last')
    overlap = sum(refIM.seg((refIM.seg(:,sortorderb(int_ix))>0),:)>0, 1);
    if (intensities(sortorderb(int_ix))<thresh || nPixels(sortorderb(int_ix))<npxThr || sum(overlap>0)>5) && ~onBorder(sortorderb(int_ix))
        %find the seed with highest overlap with the current seed
        support = overlap>0;
        support(sortorderb(int_ix)) = false;
        support(intensities>4*thresh) = false;
        support(onBorder) = false;
        if any(support)
            support(overlap<max(overlap(support))) = false;
            support = find(support);
            [~,minix] = min(intensities(support)); %merge to the dimmest neighbour
            minix = support(minix);
            
            refIM.seg(:,minix) = refIM.seg(:,minix)+refIM.seg(:, sortorderb(int_ix));
            refIM.seg(:, sortorderb(int_ix)) = 0;
            
            nPixels(minix) = sum(refIM.seg(:,minix)>(refIM.IM(:)/2), 1); %update nPixels
            intensities(minix) = intensities(minix) + intensities(sortorderb(int_ix));
        elseif intensities(sortorderb(int_ix))<thresh/2
            refIM.seg(:, sortorderb(int_ix)) = 0;
        end
    end
end

%remove '0' segments
select = any(refIM.seg,1);
refIM.seg = refIM.seg(:,select);

%ensure that the segments add up to the reference image within the support
valid = full(any(refIM.seg,2));
refIM.seg(valid,:) = refIM.seg(valid,:).*(refIM.IM(valid)./(sum(refIM.seg(valid,:),2))); %, 1, size(refIM.seg,2));

disp('Done segmentation.');
end