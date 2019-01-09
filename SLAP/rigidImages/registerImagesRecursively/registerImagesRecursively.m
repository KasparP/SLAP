function [images, data] = registerImagesRecursively(images,data,forcemerge)
%%
% inputs
%   images: cell array of images to register; all images need to be of the
%           same size and data type
%   data: internal variable for recursion
%   forcemerge: if true, the merge criterion is ignored and all images merged
%
%
% outputs:
%   images: registered image
%   data: internal variable

%% prepare data
if nargin < 2 || isempty(data)
    data = struct('fftimages',{{}},'weights',{{}},'corrtable',[],'paddedimages',{{}},'merges',{{}});
    
    numImages = length(images);
    data.corrtable = nan(numImages,numImages,3);
    data.iterations = 0;
    
    for idx = 1:numImages
        if ~isfloat(images{idx})
            images{idx} = single(images{idx});
        end
    
        IM = images{idx};
        sz = ceil(size(IM)./4).*8;
        
        data.paddedimages{idx} = zeros(sz,'like',IM);
        data.paddedimages{idx}(sz(1)/4+1:sz(1)/4+size(IM,1), sz(2)/4+1:sz(2)/4+size(IM,2)) = IM;
        
        data.weights{idx} = zeros(sz);
        data.weights{idx}(sz(1)/4+1:sz(1)/4+size(IM,1), sz(2)/4+1:sz(2)/4+size(IM,2)) = 1;
        
        data.fftimages{idx} = fft2(data.paddedimages{idx});
        
        images{idx} = data.paddedimages{idx};
        
        data.V{idx} = zeros(sz); %The variance associated with each pixel % initialize to zero for math to be correct
        %data.C{idx} = Inf;  %The merge criterion
    end
end

if nargin < 3 || isempty(forcemerge)
    forcemerge = false;
end

%% recursion criterium
err = data.corrtable(:,:,1);
if all(isinf(err(~logical(eye(size(err))))))
    return % stop recursion
else
    data.iterations = data.iterations + 1;
end

%% image analysis
[data] = buildCorrelationTable(data);

% find two images with minimum error in correlation table
[~,errminidx] = min(reshape(data.corrtable(:,:,1),[],1));
[i,j] = ind2sub(size(data.corrtable),errminidx);

% set the image with more total pixels as reference
if sum(data.weights{i}(:)) >= sum(data.weights{j}(:))
    keepme = i;
    deleteme = j;
else
    keepme = j;
    deleteme = i;
end

if keepme==deleteme
    keyboard
end

IM1 = data.paddedimages{keepme};
IM2 = data.paddedimages{deleteme};
W1 = data.weights{keepme};
W2 = data.weights{deleteme};
V1 = data.V{keepme};
V2 = data.V{deleteme};

row_shift = data.corrtable(keepme,deleteme,2);
col_shift = data.corrtable(keepme,deleteme,3);

% apply motion correction
IM2 = shiftim(IM2,row_shift,col_shift);
W2 = shiftim(W2,row_shift,col_shift);
V2 = shiftim(V2,row_shift,col_shift); % is it okay to linearly interpolate the Variance?

IM = IM1+IM2;
W = W1 + W2;

M1 = IM1./W1; M1(W1<1e-6) = 0; %mean images
M2 = IM2./W2; M2(W2<1e-6) = 0; %mean images
Mnew = (IM1+IM2)./W; Mnew(W<1e-6) = 0; %new mean image
 
Vnew = (max(0,(W1-1)).*V1 + max(0,(W2-1)).*V2 + W1.*(M1-Mnew).^2 + W2.*(M2-Mnew).^2)./max(1,(W-1));
Vnew(W<=1) = 0;

Wvalid_1 = W1>=2 & W2>=1;
Wvalid_2 = W2>=2 & W1>=1;

Cnew_1 = sum(sqrt(Vnew(Wvalid_1)./(W(Wvalid_1)-1)));
C1 = sum(sqrt(V1(Wvalid_1)./(W1(Wvalid_1)-1)));

Cnew_2 = sum(sqrt(Vnew(Wvalid_2)./(W(Wvalid_2)-1)));
C2 = sum(sqrt(V2(Wvalid_2)./(W2(Wvalid_2)-1)));

if ~forcemerge && (Cnew_1>C1 || Cnew_2>C2)
    figure, subplot(1,2,1)
    Vdiff1 = (Vnew-V1)./W1; Vdiff1(~Wvalid_1) = 0;
    imshow(Vdiff1,[]); title('Variance difference- 1')
    subplot(1,2,2)
    Vdiff2 = (Vnew-V2)./W2; Vdiff2(~Wvalid_2) = 0;
    imshow(Vdiff2,[]); title('Variance differnce- 2')

    figure()
    subplot(2,2,1)
    imagesc(images{deleteme});
    axis image
    title('Reject Image');
    subplot(2,2,2)
    imagesc(data.weights{deleteme});
    axis image
    title('Reject Weights');
        
    subplot(2,2,3)
    imagesc(images{keepme});
    axis image
    title('Keep Image');
    subplot(2,2,4)
    imagesc(data.weights{keepme});
    axis image
    title('Keep Weights');
    drawnow
    
    [row_shift col_shift]
    [Cnew_1 C1 Cnew_2 C2]
%     
%     keyboard
    
    %reject the merge, delete the group with fewer images (worse criterion?)
%     if Cj<Ci
%         deleteme = i; keepme = j;
%     else
%         deleteme = j; keepme = i;
%     end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this keeps both images in the pool, but flags them as unmergeable
    data.corrtable(keepme,deleteme,1) = Inf;
    data.corrtable(deleteme,keepme,1) = Inf;
% alternative: effectively remove deletme from the pool
%    data.corrtable(deleteme,:,1) = Inf;
%    data.corrtable(:,deleteme,1) = Inf;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    keepme = [];
    deleteme = [];
end

if ~isempty(keepme) && ~isempty(deleteme)
    data.merges{end+1} = struct('keepme',keepme,'deleteme',deleteme,'row_shift',row_shift,'col_shift',col_shift); %,'Cnew',Cnew
end

if keepme==deleteme
    keyboard
end
if ~isempty(keepme)
    %accept the merge
    IM_normalized = IM./W;
    IM_normalized(W<1e-6) = 0;
    images{keepme} = IM_normalized;
    data.paddedimages{keepme} = IM;
    data.fftimages{keepme} = fft2(IM_normalized);
    data.weights{keepme} = W;
    data.corrtable(keepme,:,:) = NaN; % reset row in correlation table
    data.corrtable(:,keepme,:) = NaN; % reset column in correlation table
    data.V{keepme} = Vnew;
    %data.C{keepme} = Cnew;
end

if ~isempty(deleteme)
    images(deleteme) = [];
    data.paddedimages(deleteme) = [];
    data.fftimages(deleteme) = [];
    data.weights(deleteme) = [];
    data.corrtable(deleteme,:,:) = []; % delete row in correlation table
    data.corrtable(:,deleteme,:) = []; % delete column in correlation table
    data.V(deleteme) = [];
    %data.C(deleteme) = [];
end

% recurse
[images, data] = registerImagesRecursively(images,data,forcemerge);
end

function [data] = buildCorrelationTable(data)
numImages = numel(data.paddedimages);
K = ceil(30/numImages); %number of target images to compare each candidate to; algorithm complexity it O(K^2)
for idx = 1:numImages
    for jdx = idx+1:min(numImages, idx+K)
        if isnan(data.corrtable(idx,jdx,1))
            output = dftregistration(data.fftimages{idx},data.fftimages{jdx},2);
            error = output(1);
            
            if isnan(error)
                error('Error making correlation table: this is probbaly due to images on the GPU (?)- images sometimes come up as NaNs!')
            end

%            diffphase = output(2);
            net_row_shift = output(3);
            net_col_shift = output(4);
            
            % corrtable(:,:,1) is error between images
            % corrtable(:,:,2) is net_row_shift between images
            % corrtable(:,:,3) is net_col_shift between images
            data.corrtable(idx,jdx,:) = shiftdim([error net_row_shift net_col_shift],-1);
            data.corrtable(jdx,idx,:) = shiftdim([error -net_row_shift -net_col_shift],-1);
        end
    end
end
end

function im = shiftim(im,row_shift,col_shift)
[xx,yy] = meshgrid(1:size(im,2),1:size(im,1));
xx = xx-col_shift;
yy = yy-row_shift;
im(:,:) = interp2(im(:,:),xx,yy,'linear',0);
end
