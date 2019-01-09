function refIM = segment3D (refIM, optsin)
%Segment a 3D image by skeletonization

%skeletonize each plane
%starting with the most prominent skeleton points (i.e. brightest in a
%tophat filtered image):
%remove nearby points in xy
%if the planes directly above or below are in the mask,
%place points there as part of this seed,
%and remove enarby points in xy
image = double(refIM.IM);
image(isnan(image)) = 0;

opts.XYscale = 0.2;
opts.mask = ones(size(image));
opts.maxseeds = Inf; %maximum number of elements
opts.sharpness = 1.1; %how sharp the edges between seeds should be; a nonnegative number. Higher is sharper
opts.distThreshUM = 1.2;  %minimum distance between shaft seeds, in pixels: typical is 1.8 um
if nargin>1 %UPDATE WITH USER-SPECIFIED OPTIONS
    for field = fieldnames(optsin)'
        opts.(field{1}) = optsin.(field{1});
    end
end

distThresh = round(opts.distThreshUM/opts.XYscale);

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
bwThin = bw;

i1 = cell(size(bw,3),1); 
i2 = cell(size(bw,3),1); 
i3 = cell(size(bw,3),1); 
b =  cell(size(bw,3),1);

for z  = 1:size(bw,3)
    I_gauss = image(:,:,z);
    I_gauss = imtophat(I_gauss,strel('disk',ceil(1/opts.XYscale))); %background subtract with radius of 1 micron
    
    %thin bw to cut apart merged spines
    I_thin = I_gauss./imgaussfilt(I_gauss, 80); %compensate for dimmer/brighter regions
    thresh = 4; %prctile(I_thin(bw(:,:,z)), 90)/3;
    bwZ = I_thin>thresh & bw(:,:,z);
    bwZ = bwareafilt(bwZ,[6 inf]);
    bwThin(:,:,z) = bwZ;
    
    %Adjust BW to prevent very large areas from being segmented 'in'
    bwThresh = bw(:,:,z) & bwareafilt(bw(:,:,z) & (I_thin>0.8), [10 Inf]); 
    holes = ~bwThresh & imfill(bwThresh, 'holes');
    smallholes = holes & ~bwareaopen(holes, 10);
    bw(:,:,z) = bw(:,:,z) & imdilate(bwThresh | smallholes, strel('diamond', 1));

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
    
    %skelCoreZ =  bwmorph(coreZ,'skel',Inf);
    skelShaftZ =  bwmorph(bwZ,'skel',Inf);
    
    [i1{z},i2{z}] = find(skelShaftZ | skelCoreZ);
    i3{z} = z*ones(size(i2{z}));
    b{z} = I_gauss(sub2ind(size(bw(:,:,z)),i1{z},i2{z})); %brightness of skeleton points; we'll keep the brightest ones
end
refIM.bw = bw;

i1 = cell2mat(i1); i2 = cell2mat(i2); i3 = cell2mat(i3); 
pointType = refIM.labels(sub2ind(size(refIM.labels), i1,i2,i3));
b = cell2mat(b);
b(pointType==4 | pointType==2) = 2*b(pointType==4 | pointType==2); %favor putting points on certain pixels first?
[b, sortorderb] = sort(b, 'descend');
b0 = b(ceil(length(b)/3));

i1 = i1(sortorderb);
i2 = i2(sortorderb); 
i3 = i3(sortorderb); 
pointType = pointType(sortorderb);

disp('Selecting seed points...')
exclude = false(length(b),1);
added = cell(size(i3));
addedTF = false(length(i3), size(refIM.IM,3));
for ix = 1:length(b) %for every point, in order of prominence
    if ~exclude(ix)
        %exclude nearby points within this plane
        distThreshOther = distThresh; distThreshShaft = distThresh;
        if pointType(ix)==4 %if shaft, exclude shaft points over a wider distance
            distThreshShaft = 1.5*distThresh;
        elseif pointType==2
            distThreshShaft = distThresh/1.5; distThreshOther = distThresh/1.5;
        end
        
        n1 = false(size(i1));
        select = pointType==4 & i3==i3(ix);
        n1(select) =  sqrt((i2(ix)-i2(select)).^2 + (i1(ix)-i1(select)).^2)<distThreshShaft;
        select = pointType~=4 & i3==i3(ix);
        n1(select) =  sqrt((i2(ix)-i2(select)).^2 + (i1(ix)-i1(select)).^2)<distThreshOther;

        for z2 = i3(ix)-1:-1:max(1,i3(ix)-3) %extend the segment vertically within the mask
            n4 = i3==z2; n4(ix:end) = false; 
            n4 = (n4 | addedTF(:,z2)) & ~exclude;
            
            select = n4 & pointType==4;
            d1 =  any(sqrt((i2(ix)-i2(select)).^2 + (i1(ix)-i1(select)).^2)<distThreshShaft);
            select = n4 & pointType~=4;
            d2 =  any(sqrt((i2(ix)-i2(select)).^2 + (i1(ix)-i1(select)).^2)<distThreshOther);
            if  ~(d1 || d2) && bw(i1(ix),i2(ix),z2)
                veto = false(size(n1));
                select = i3==z2 & pointType == 4;
                veto(select) = sqrt((i2(ix)-i2(select)).^2 + (i1(ix)-i1(select)).^2)<distThreshShaft;
                select = i3==z2 & pointType ~= 4;
                veto(select) = sqrt((i2(ix)-i2(select)).^2 + (i1(ix)-i1(select)).^2)<distThreshOther;
                n1 = n1 | veto; %add vetoes to the excluded list
                
                added{ix} = [added{ix} z2];
                addedTF(ix,z2) = true;
            else
                break
            end
        end
        
        for z2 = i3(ix)+1:min(size(bw,3), i3(ix)+3) %extend the segment vertically within the mask
            n4 = i3==z2; n4(ix:end) = false; 
            n4 = (n4 | addedTF(:,z2)) & ~exclude;
            
            select = n4 & pointType==4;
            d1 =  any(sqrt((i2(ix)-i2(select)).^2 + (i1(ix)-i1(select)).^2)<distThreshShaft);
            select = n4 & pointType~=4;
            d2 =  any(sqrt((i2(ix)-i2(select)).^2 + (i1(ix)-i1(select)).^2)<distThreshOther);
            if  ~(d1 || d2) && bw(i1(ix),i2(ix),z2)
                veto = false(size(n1));
                select = i3==z2 & pointType == 4;
                veto(select) = sqrt((i2(ix)-i2(select)).^2 + (i1(ix)-i1(select)).^2)<distThreshShaft;
                select = i3==z2 & pointType ~= 4;
                veto(select) = sqrt((i2(ix)-i2(select)).^2 + (i1(ix)-i1(select)).^2)<distThreshOther;
                n1 = n1 | veto; %add vetoes to the excluded list
                
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

i1 = i1(~exclude); i2 = i2(~exclude); i3 = i3(~exclude); added = added(~exclude);
pointType = pointType(~exclude);

%calculate influence of seeds over nearby pixels with a
%'flood fill' operation
disp('Filling Volume...')
nh_size = 3*distThresh; %maximum radius of seed influence, in pixels. This should be >2*distThresh
I = cell(4*length(i1),1); J = I; V = I;
cellIX = 1;
for zz = 1:size(bw,3)
    i_ixs = find(i3==zz  |  cellfun(@(x)any(x==zz), added));
    bw_padded = padarray(bw(:,:,zz), [nh_size nh_size]);
    shaft_padded = padarray(shaft(:,:,zz), [nh_size nh_size]);
    nonshaft_padded = bw_padded & ~shaft_padded;
    IMpadded = padarray(refIM.IM(:,:,zz), [nh_size nh_size]);
    
    nhood = false(2*nh_size+1,2*nh_size+1,length(i_ixs));
    nhoodbw = false(2*nh_size+1,2*nh_size+1,length(i_ixs));
    nhoodIM = zeros(2*nh_size+1,2*nh_size+1,length(i_ixs));
    for seed = 1:length(i_ixs)
        if shaft(i1(i_ixs(seed)), i2(i_ixs(seed)), i3(i_ixs(seed))) %i3, not zz here
            nhood(:,:,seed) = shaft_padded((0:2*nh_size)+i1(i_ixs(seed)), (0:2*nh_size)+i2(i_ixs(seed)));
        else
            nhood(:,:,seed) = nonshaft_padded((0:2*nh_size)+i1(i_ixs(seed)), (0:2*nh_size)+i2(i_ixs(seed)));
        end
        nhoodbw(:,:,seed) = bw_padded((0:2*nh_size)+i1(i_ixs(seed)), (0:2*nh_size)+i2(i_ixs(seed)));
        nhoodIM(:,:,seed) = IMpadded((0:2*nh_size)+i1(i_ixs(seed)), (0:2*nh_size)+i2(i_ixs(seed)));
    end
    nhoodIM = nhoodIM./repmat(max(max(nhoodIM,[],1),[],2), size(nhood,1),size(nhood,2),1);
    influence = zeros(2*nh_size+1,2*nh_size+1,length(i_ixs));
    influence(nh_size+1,nh_size+1,:) = nhoodIM(nh_size+1,nh_size+1,:);
    
    SE = ones(3)/9;
    %Flood fill:
    for i = 1:2*distThresh
        influence = influence + convn(influence.*nhoodIM, SE, 'same').*nhood;
    end
    for i = 1:distThresh
        influence = influence + convn(influence.*nhoodIM, SE, 'same').*nhoodbw;
    end
    influence = influence./sum(sum(influence,1),2); %so seeds at edges don't get swamped
    influence(:,:,pointType(i_ixs)==4) = 1.5*influence(:,:,pointType(i_ixs)==4); %the shaft points are bigger; compensate
    influence(isnan(influence)) = 0;
    
    %accumulate indices for sparse matrix, #voxels x #seeds
    for seed = 1:length(i_ixs)
        [row, col, val] = find(influence(:,:,seed));
        seg_ixs = sub2ind(size(bw), row+i1(i_ixs(seed))-(nh_size+1),col+i2(i_ixs(seed))-(nh_size+1),zz*ones(size(col)));
        %accumulators for sparse matrix:
        I{cellIX} = seg_ixs; J{cellIX} = i_ixs(seed)*ones(size(seg_ixs)); V{cellIX} = val; cellIX = cellIX+1;
    end
end

%error checking: make sure every pixel is influenced by a seed
disp('Filling unlabeled areas...')
no_seed = bw;
no_seed(unique(cell2mat(I))) = false; 
num = 0;
if any(no_seed(:))
    no_seed = bwareaopen(no_seed, 4, 4);%discard isolated pixels
    [L,num] = bwlabeln(no_seed);
    i1B = nan(num,1); i2B = nan(num,1); i3B = nan(num,1); pTB = nan(num,1);
    for Lix = 1:num
        L_inds = find(L(:)==Lix);
        [~,maxind] = max(refIM.IM(L_inds));
        [i1B(Lix), i2B(Lix), i3B(Lix)] = ind2sub(size(refIM.IM), L_inds(maxind));
        pTB(Lix) = refIM.labels(maxind);
        
        I{cellIX} = L_inds; J{cellIX} = (length(i1)+Lix)*ones(size(I{cellIX})); V{cellIX} = ones(size(I{cellIX}));
        cellIX = cellIX+1;
    end
    i1 = [i1 ; i1B]; i2 = [i2 ; i2B]; i3 = [i3 ; i3B]; pointType = [pointType ; pTB]; 
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

%FUSION CRITERIA
disp('Fusing dim seeds...')
thresh = 0.2*max(1, size(refIM.seg,2)/1000);
intensities = full(sum(refIM.seg,1)); intensities = intensities./mean(intensities);
[sorted, sortorderb] = sort(intensities, 'ascend');
%nPixels = sum(refIM.seg>(refIM.IM(:)/2), 1); npxThr =50; 
for int_ix = find(sorted>0,1,'first'):find(sorted<thresh,1,'last')
    %overlap = sum(refIM.seg((refIM.seg(:,sortorderb(int_ix))>0),:)>0, 1); %how much each other seed overlaps with this one
    if (intensities(sortorderb(int_ix))<thresh) % || (nPixels(sortorderb(int_ix))<npxThr)  %&& ~onBorder(sortorderb(int_ix))
        %find the seed with highest overlap with the current seed
        %support = overlap>0;
        support = any(refIM.seg((refIM.seg(:,sortorderb(int_ix))>0),:), 1);
        support(sortorderb(int_ix)) = false;
        support(intensities>5*thresh) = false;
        if any(support)
            support = find(support);
            overlap = sum(min(repmat(refIM.seg(:,sortorderb(int_ix)), 1, length(support)), refIM.seg(:,support)),1);
            [~,maxix] = max(overlap);
            select= support(maxix);
            
%             [~,minix] = min(intensities(support)); %merge to the dimmest neighbour
%             minix = support(minix);
            
            refIM.seg(:,select) = refIM.seg(:,select)+refIM.seg(:, sortorderb(int_ix));
            refIM.seg(:, sortorderb(int_ix)) = 0;
            
            intensities(select) = intensities(select) + intensities(sortorderb(int_ix));
        elseif intensities(sortorderb(int_ix))<thresh/2
            refIM.seg(:, sortorderb(int_ix)) = 0;
        end
    end
end

%remove '0' segments
select = any(refIM.seg,1);
refIM.seg = refIM.seg(:,select);
i1 = i1(select);
i2 = i2(select);
i3 = i3(select);
pointType = pointType(select);

%ensure that the segments add up to the reference image within the support
valid = full(any(refIM.seg,2));
refIM.seg(valid,:) = refIM.seg(valid,:).*(refIM.IM(valid)./(sum(refIM.seg(valid,:),2))); %, 1, size(refIM.seg,2));
refIM.pts = [i1 i2 i3];
refIM.ptType = pointType;
refIM.segopts = opts;

disp('Done segmentation.');
end