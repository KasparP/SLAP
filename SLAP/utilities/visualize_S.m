function S_block = visualize_S(S, filepath)
nS = size(S.seg,2);

RGB = rand(3,nS); RGB = RGB./repmat(sum(RGB,1), 3,1);
S_RGB = sqrt([S.seg*RGB(1,:)' S.seg*RGB(2,:)' S.seg*RGB(3,:)']);
S_RGB = 1.5* S_RGB./max(S_RGB(:));
S_RGB = reshape(full(S_RGB), [size(S.bw,1) size(S.bw,2) size(S.bw,3) 3]);
figure, imshow3D(S_RGB);

IM = full(sum(S.seg,2));
for i = 1:size(S.seg,2)
    valid = S.seg(:,i)>0;
    S.seg(valid,i) = min(1,S.seg(valid,i)./IM(valid));
end

S_block = [S.seg*RGB(1,:)' S.seg*RGB(2,:)' S.seg*RGB(3,:)'];
S_block = reshape(full(S_block), [size(S.bw,1) size(S.bw,2) size(S.bw,3) 3]);

if ~isfield(S,'mask')
    S.mask = ones(size(S.bw,1), size(S.bw,2));
end
M = repmat(0.25*S.mask,1,1,size(S_block,3),size(S_block,4));
M(S_block>0.01) = 0;
%figure, imshow3D(S_block./prctile(S_block(:), 99.9));
%EXAMPLE FIGURE
%figure, imshow(max(squeeze(S_RGB(:,:,7,:)),0.2*mask) .*max(mask,.5))
figure, imshow3D((S_block./prctile(S_block(:), 99.9)+M).*repmat(max(S.mask,.4),1,1,size(S_block,3),size(S_block,4)));
if nargin>1
    if isempty(filepath)
        [dr, fn] = uiputfile('*.tif', 'Select a filename to save the segmentation');
        filepath = [dr filesep fn(1:end-4)];
    end
    options.overwrite = true; options.color = true;
    errorcode1 = saveastiff(permute(uint8(255*S_RGB), [1 2 4 3]), [filepath '_RGB.tif'], options); %tiff file for easy reading
    errorcode2 = saveastiff(permute(uint8(255*S_block), [1 2 4 3]), [filepath '_block.tif'], options); %tiff file for easy reading
    if errorcode1 || errorcode2
        keyboard
    end
end
end