function [imageSet, shifts, errors, D] = align2D_notemplate(imageSet, optsin)
%input: a 3-dimensional matrix of 2D images to align to each other

%disp('Registering motion without a template...')

%DEFAULT OPTIONS
opts.window = 30;  %maximum movement;
opts.nRefs = 3; %number of reference values to produce; we will use the most common one
opts.nFramesRef = 300; %create a reference from nFramesRef frames
opts.doplot = false;
opts.prior.strength = 0.5;
opts.prior.width = 30;
opts.prior.sharpness = 50;
if nargin>1 %UPDATE WITH USER-SPECIFIED OPTIONS
    for field = fieldnames(optsin)'
        opts.(field{1}) = optsin.(field{1});
    end
end
if ~isfield(opts, 'K')
    opts.K = opts.nRefs; %how far apart to consider values to correlate against each other
end


nFrames = size(imageSet,3);
nFramesRef = min(nFrames, opts.nFramesRef);
selectFrames = 1:floor(nFrames/nFramesRef):nFrames; 

%build correlation table
FF = imageSet(:,:, selectFrames);  %The merged images will be computed within this matrix
thresh = nanmedian(reshape(FF(:,:,1:ceil(nFrames/101):end),1,[]))+ max(1,nanstd(FF(FF<prctile(reshape(FF(:,:,1:ceil(nFrames/101):end),1,[]),90))));
FF = max(0, FF-thresh); %sparsify to get rid of low-amplitude plaiding
FF = FF-repmat(nanmean(nanmean(FF,1),2), size(FF,1), size(FF,2), 1);
FF(isnan(FF)) = 0;

FFTs = cell(nFramesRef, 1);
corrtable = nan(nFramesRef,nFramesRef);
shift = nan(2,nFramesRef,nFramesRef);
N = ones(nFramesRef,1);

updatecorrtable(1:nFramesRef); %populate the correlation table

while length(N) > opts.nRefs
    %find max element of corrtable
    [i,j] = find(corrtable==min(corrtable(:)),1, 'last'); %ensures j>i
    %apply shift and average together, considering N
    combine(i,j);
    updatecorrtable(max(1, i-opts.K):min(length(N), i+opts.K));
end

%select the best reference image
[maxN, ref_ix] = max(N);

refFFT = FFTs{ref_ix};
shifts = nan(2,nFrames);
errors = nan(nFrames,1);
%align each frame to the reference
for frame = 1:nFrames
    f = max(0, imageSet(:,:,frame)-thresh);
    f = f-nanmean(f(:));
    f(isnan(f)) = 0;
    
    %dft registration
    fFFT = fft2(f);
    output = dftregistration_wprior(refFFT,fFFT,1, opts.prior);
    shifts(:,frame) = [output(4) output(3)];
    errors(frame) = output(1);
end
%apply shifts
imageSet(isnan(imageSet)) = 0;
shifts = shifts - repmat(round(mean(shifts,2)), 1, nFrames);
for frame = 1:nFrames
    imageSet(:,:,frame) = imtranslate(imageSet(:,:,frame), shifts(:,frame)', 'fillvalues', median(reshape(imageSet(:,:,frame),1,[])));
end

% if isfield(opts, 'chunkSize') && opts.chunkSize>0
%     ref = trimmean(imageSet,20,3);
%     nRows = size(ref,1);
%     chunkStarts = linspace(1, nRows-opts.chunkSize, ceil(nRows/opts.chunkSize));
%     shifts = nan(2, size(imageSet,3), length(chunkStarts));
%     for c_ix = 1:length(chunkStarts)
%         %get a chunk of the ref
%         refFFT = fft2(ref(chunkStart:chunkStart+opts.chunkSize,:));
%         for frame = 1:size(imageSet,3)
%             fFFT = fft2(imageSet(chunkStart:chunkStart+opts.chunkSize,:,frame));
%             output = dftregistration_wprior(refFFT,fFFT,1, opts.prior);
%             %get a chunk of the frame
%             shifts(:, frame, c_ix) =  [output(4) output(3)];
%         end
%     end
% end

% if isfield(opts, 'demons') && opts.demons
%     ref = trimmean(imageSet,20,3);
%     D = cell(1,size(imageSet,3));
%     for frame = 1:size(imageSet,3)
%        [imageSet(:,:,frame)  D{frame}] = imregdemons(imageSet(:,:,frame),ref); 
%     end
% end

%diagnostic dummy check
if opts.doplot
    figure('name', 'After motion alignment')
    imshow3D(imageSet);
end

%%%%%%END OF MAIN FUNCTION

    function combine(i,j)
        W = N(i)/(N(i)+N(j)); 
        IMi = imtranslate(FF(:,:,i), -ceil([shift(2,i,j) shift(1,i,j)]*(1-W)), 'FillValues', nan);
        IMj = imtranslate(FF(:,:,j), floor([shift(2,i,j) shift(1,i,j)]*W), 'FillValues', nan);
        FF(:,:,i) = nansum(cat(3, N(i).*IMi, N(j)*IMj), 3)./(N(i).*~isnan(IMi)+N(j).*~isnan(IMj)); %average together measurements by weight
        FF(isnan(FF)) = 0;
        
        N(i) = N(i) + N(j);
        shift(:,i,:) = nan;
        shift(:,:,i) = nan;
        corrtable(i,:,:) = nan;
        corrtable(:,i,:) = nan;
        
        %delete j
        FFTs(j,:) = [];
        N(j) = [];
        FF(:,:,j) = [];
        shift(:,:,j) = [];
        shift(:,j,:) = [];
        corrtable(j,:,:) = [];
        corrtable(:,j,:) = [];
    end


    function updatecorrtable (doframes)
        %requiredoverlap = sum(line_ixs) - opts.window;
        for f1 = doframes %update FFTS
            FFTs{f1} = fft2(FF(:,:,f1));
        end
        for f1 = doframes %update correlations
            for f2 = max(1,f1-opts.K) : min(size(corrtable,1),f1+opts.K)
                if isnan(corrtable(f1,f2)) && f1~=f2
                    output = dftregistration(FFTs{f1},FFTs{f2},1);
                    corrtable(f1,f2) = output(1);
                    corrtable(f2,f1) = corrtable(f1,f2);
                    shift(:, f1,f2) = [output(3) output(4)];
                    shift(:,f2,f1) = -[output(3) output(4)];
                end
            end
        end
    end
end