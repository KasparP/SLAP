function refIM = SLAPMi_ROI_seg (ROIs, refIM)
% merge segments in a reference image according to provided ROIs
XYscale = 0.2;

%ensure that ROIs are nonoverlapping
if any(any (all(ROIs,3)))
    error('ROIs overlap; please correct');
end

if isfield(refIM,'seg')
    for rix = 1:size(ROIs,3)
        ROI3 = repmat(ROIs(:,:,rix), 1, 1, size(refIM.bw,3));
        inside = any(refIM.seg(ROI3(:),:),1); %inside = ~any(S.seg(~ROI3(:),:),1);
        
        refIM.seg = [sum(refIM.seg(:, inside),2) refIM.seg(:,~inside)];
    end
    refIM = rmfield(refIM, {'pts', 'ptType'});
else
    image = double(refIM.IM);
    image(isnan(image)) = 0;
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
    
    for z  = 1:size(bw,3)
        I_gauss = image(:,:,z);
        I_gauss = imtophat(I_gauss,strel('disk',ceil(1/XYscale))); %background subtract with radius of 1 micron
        
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
    end
    refIM.bw = bw;
    
    for rix = 1:size(ROIs,3)
        ROI3 = repmat(ROIs(:,:,rix), 1, 1, size(refIM.bw,3));
        refIM.seg(:,rix) = ROI3(:).*refIM.bw(:).*refIM.IM(:);
    end
    refIM.seg(:,end+1) = (refIM.IM(:) - sum(refIM.seg,2)).*refIM.bw(:);
end

end