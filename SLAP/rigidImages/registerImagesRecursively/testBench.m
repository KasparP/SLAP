numImages = 50;
maxMotion = 10; % in pixel
sizeImages = [200,200];
noiseAmplitude = 50;
numNoiseImages = 10;  % images that are added to the pool, which contain of noise only

%% prepare source images
close all
S = load('durer','X'); % source image
sourceIm = single(S.X);

images = {};
true_motions = {};
for idx = 1:numImages
    imcenter = [249 379];
    rect = floor([imcenter-sizeImages./2 sizeImages]);
    
    motion =  2*maxMotion*(rand(1,2)-0.5); % generate random motion vector
    rect(1:2) = rect(1:2) + motion;
    images{idx} = interp2(sourceIm, linspace(rect(2),rect(2)+rect(4)-1, rect(4)), linspace(rect(1),rect(1)+rect(3)-1, rect(3))','linear');
    images{idx} = images{idx} + rand(size(images{idx}))*2*noiseAmplitude-noiseAmplitude; % add noise to image;
    true_motions{idx} = motion;
end

for idx = 1:numNoiseImages
    images{end+1} = std(sourceIm(:)) .* randn(sizeImages) + mean(sourceIm(:));
    true_motions{end+1} = [0 0];
end

% shuffle random images into deck
idxs = randperm(length(images));
images = images(idxs);
true_motions = true_motions(idxs);

%mean subtract images before zero-padding
for idx = 1:length(images)
    images{idx} = images{idx} - mean(images{idx}(:));
end

%% perform registration
disp('Registration started.')

start = tic();
[registeredImage,data] = registerImagesRecursively(images);
stop = toc(start);

% tic();
% [registeredImage2,data2] = registerImagesLoop(images);
% toc();

fprintf('Image registration took %fs\n',stop);
fprintf('Registration finished in %d iterations\n',data.iterations);

%% verification
%%% pad images and create weights mask
paddedimages = {};
weights = {};
motiontable = nan(length(images),length(images),2);
for idx = 1:length(images) 
    sourceIm = images{idx};
    sz = ceil(size(sourceIm)./4).*8;
    paddedimages{idx} = zeros(sz,'like',sourceIm);
    paddedimages{idx}(sz(1)/4+1:sz(1)/4+size(sourceIm,1), sz(2)/4+1:sz(2)/4+size(sourceIm,2)) = sourceIm;
    
    weights{idx} = zeros(sz);
    weights{idx}(sz(1)/4+1:sz(1)/4+size(sourceIm,1), sz(2)/4+1:sz(2)/4+size(sourceIm,2)) = 1;
    
    for jdx=idx+1:length(images)
        relative_motion = shiftdim(true_motions{jdx} - true_motions{idx},-1);
        motiontable(idx,jdx,:) = relative_motion;
        motiontable(jdx,idx,:) = -relative_motion;
    end
end

%%% shift padded images and weights mask back
motion_eval_true = cell(1,length(images));
motion_eval_registered = cell(1,length(images));
for mergeidx = 1:length(data.merges)
    merge = data.merges{mergeidx};
    keepme = merge.keepme;
    deleteme = merge.deleteme;
    
    if ~isempty(keepme)
        IM1 = paddedimages{keepme};
        IM2 = paddedimages{deleteme};
        
        col_shift = motiontable(keepme,deleteme,2);
        row_shift = motiontable(keepme,deleteme,1);
        
        sz = size(IM1);
        [xx,yy] = meshgrid(1:sz(2),1:sz(1));
        xx = xx - col_shift;
        yy = yy - row_shift;
        paddedimages{deleteme} = interp2(paddedimages{deleteme},xx,yy,'linear',0);
        weights{deleteme} = interp2(weights{deleteme},xx,yy,'linear',0);
        
        paddedimages{keepme} = paddedimages{keepme} + paddedimages{deleteme};
        weights{keepme} = weights{keepme} + weights{deleteme};
        
        if isempty(motion_eval_true(keepme))
            motion_eval_true{keepme} = [];
            motion_eval_registered{keepme} = [];
        end
        motion_eval_true{keepme}(end+1,:) = [row_shift col_shift];
        motion_eval_registered{keepme}(end+1,:) = [merge.row_shift merge.col_shift];
    end
    
    if ~isempty(deleteme)
        paddedimages(deleteme) = [];
        weights(deleteme) = [];
        motiontable(:,deleteme,:) = [];
        motiontable(deleteme,:,:) = [];
        motion_eval_true(deleteme) = [];
        motion_eval_registered(deleteme) = [];
    end
end

eval_image = {};
eval_weights = {};
for idx = 1:length(paddedimages)
    eval_image{idx} = paddedimages{idx};
    eval_weights{idx} = weights{idx};
    eval_image{idx} = eval_image{idx} ./ eval_weights{idx};
    eval_image{idx}(eval_weights{idx} < 1e-6) = 0;
end

%% plot error

for idx = 1:length(registeredImage)
figure
subplot(2,3,1);
imagesc(registeredImage{idx});
title(sprintf('Registered Image Group %d',idx))
axis image; axis off;
subplot(2,3,2);
imagesc(eval_image{idx});
axis image; axis off;
title(sprintf('Evaluation Image Group %d',idx))
subplot(2,3,3);
imagesc(registeredImage{idx}-eval_image{idx});
axis image; axis off;
title('Error')

subplot(2,3,4);
imagesc(data.weights{idx});
title('Registration samples per pixel')
axis image; axis off;
subplot(2,3,5);
imagesc(eval_weights{idx});
axis image; axis off;
title('Evaluation samples per pixel')
subplot(2,3,6);
imagesc(data.weights{idx}-eval_weights{idx});
axis image; axis off;
title('Error')
end


%% plot true motion vs registered motion
for idx = 1:length(motion_eval_true)
figure
subplot(2,1,1);
plot([motion_eval_true{idx}(:,1),motion_eval_registered{idx}(:,1)])
legend('true row shift','registered row shift');
title(sprintf('Motion correction error for each merge Group %d',idx));

subplot(2,1,2);
plot([motion_eval_true{idx}(:,2), motion_eval_registered{idx}(:,2)])
legend('true col shift','registered col shift')
end


