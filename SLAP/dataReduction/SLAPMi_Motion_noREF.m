function [scandata, hF] = SLAPMi_Motion_noREF (scandata, optsin)
%align scandata without a reference image to get rid of motion

if isfield(scandata.frames(1), 'motion') && isfield(scandata.frames(1).motion, 'noREF')
   disp('Fast motion was previously registered, skipping SLAPMi_Motion_noREF')
   return
end

disp('Registering fast motion without a template...')

%DEFAULT OPTIONS
opts.alignChan = 1;
opts.window = 30;  %maximum movement;
opts.K = 6; %how far apart to consider values to correlate against each other
opts.nRefs = floor(opts.K-1); %number of reference values to produce; we will use the most common one
opts.nFramesRef = 500; %create a reference from first nFramesRef frames
opts.doplot = true;
if nargin>1 %UPDATE WITH USER-SPECIFIED OPTIONS
    for field = fieldnames(optsin)'
         opts.(field{1}) = optsin.(field{1});
    end
end

%prior for alignment
prior.strength = 0.5;
prior.width = opts.window;
prior.sharpness = 20;

nchan = size(scandata.frames(1).pmtData,2);
if opts.doplot
   P = [scandata.frames.pmtData];
   P = reshape(P, [],length(scandata.frames));
   Ppre = P;
   hF = figure('name', 'Fast Motion Alignment');
   ax1 = subplot(1,3,1);
   imagesc(P);
   xlabel('Before Motion Alignment')
end

nframes = min(opts.nFramesRef, length(scandata.frames));
FF = permute(reshape([scandata.frames(1:nframes).pmtData], [], nchan, nframes), [1 3 2]); %has size (#lxls x #frames x #channels)
FF(isnan(FF)) = 0;

%build correlation table
FFTs = cell(nframes, 4);
corrtable = nan(nframes,nframes, 4);
shift = nan(nframes,nframes, 4);
N = ones(nframes,1);

updatecorrtable(1:nframes); %populate the correlation table

while length(N) > opts.nRefs
    %find max element of corrtable
    C2 = sum(corrtable,3); %use Fisher's Z transform
    [i,j] = find(C2==min(C2(:)),1, 'last'); %ensures j>i
    
    %apply shift and average together, considering N
    combine(i,j);
    updatecorrtable(max(1, i-opts.K):min(length(N), i+opts.K));
end

%select the best reference image
[maxN, ref_ix] = max(N);
ref = FF(:,ref_ix,:);

%align each frame to the reference
for L = 4:-1:1
    line_ixs = scandata.line==L;
    refFFT =  fft2(ref(line_ixs,:));
    r = (1:sum(line_ixs))';
    shifts = zeros(1,length(scandata.frames));
    errors = shifts;
    for frame = 1:length(scandata.frames) %determine shifts
       f = scandata.frames(frame).pmtData(line_ixs,:);
       f(isnan(f)) = 0;
       f_FFT = fft2(f);
       output = dftregistration_wprior(refFFT,f_FFT,2, prior);
       shifts(frame) = output(3);
       errors(frame) = min(output(1), Inf); %nan errors are infinite
    end
    shifts = smooth(shifts, 35/length(shifts), 'loess'); %polynomial smoothing to avoid cutting off peaks
    shifts = shifts-mean(shifts);
    for frame = 1:length(scandata.frames) %apply shifts
        for ch = 1:nchan
            f = scandata.frames(frame).pmtData(line_ixs,ch);
            scandata.frames(frame).pmtData(line_ixs,ch) = qinterp1(r, f, r - shifts(frame));
        end
    end
    scandata.motion.noREF.shift(L,:) = shifts;
    scandata.motion.noREF.error(L,:) = errors;
end

%correct motion-associated brightness changes using PCA?
pcaCorrect = false;
if pcaCorrect %THIS IS CURRENTLY BROKEN!! Need an unbiased estimator for z given Poiss(z+x), x
    [facs,~] = pca(scandata.motion.noREF.shift);
    P = [scandata.frames.pmtData];
    P = reshape(P, [],length(scandata.frames));
    P2 = P;
    P2(any(isnan(P),2),:) = 0; P2 = P2-mean(P2,2); %mean subtracted
    %regress intensity against top two components
    b = nan(2,size(P2,1));
    for lxl = size(P2,1):-1:1
        b(:,lxl) = regress(P2(lxl,:)', facs(:,1:2));
    end
    P = P - (facs(:,1:2)*b)'; P(P<0) = 0;
    P = reshape(P, [],nchan, length(scandata.frames));
    for frame = 1:length(scandata.frames) %apply correction
        scandata.frames(frame).pmtData = P(:,:, frame);
    end
end

%diagnostic dummy check
if opts.doplot
   P = [scandata.frames.pmtData];
   P = reshape(P, [],length(scandata.frames));
   figure(hF)
   ax2 = subplot(1,3,2);
   imagesc(P); xlabel('After Motion Alignment')
   ax3 = subplot(1,3,3);
   D = P-Ppre; D(isnan(D)) = 0;
   imagesc(D); xlabel('Difference')
   set([ax1 ax2], 'tickdir', 'out')
   linkaxes([ax1 ax2 ax3]);
end

%%%%%%END OF MAIN FUNCTION

    function combine(i,j)
        for line = 1:4
            line_ixs = scandata.line==line;
            r = (1:sum(line_ixs))';
            for chan = 1:nchan
                Fi = qinterp1(r, FF(line_ixs,i,chan), r + shift(i,j,line)/2);
                Fj = qinterp1(r, FF(line_ixs,j,chan), r - shift(i,j,line)/2);
                FF(line_ixs,i,chan) = nansum([N(i)*Fi  N(j)*Fj], 2)./(N(i).*~isnan(Fi)+N(j).*~isnan(Fj)); %average together measurements by weight
            end
        end
        FF(isnan(FF)) = 0;
        
        N(i) = N(i) + N(j);
        shift(i,:,:) = nan;
        shift(:,i,:) = nan;
        corrtable(i,:,:) = nan;
        corrtable(:,i,:) = nan;
        
        %delete j
        FFTs(j,:) = [];
        N(j) = [];
        FF(:,j,:) = [];
        shift(j,:,:) = [];
        shift(:,j,:) = [];
        corrtable(j,:,:) = [];
        corrtable(:,j,:) = [];
    end


    function updatecorrtable (doframes)
        for line = 1:4
            line_ixs = scandata.line==line;
            %requiredoverlap = sum(line_ixs) - opts.window;
            F = FF(line_ixs,:,:);
            for f1 = doframes %update FFTS
                FFTs{f1, line} = fft2(squeeze(F(:,f1,:)));
            end
            for f1 = doframes %update correlations
                for f2 = max(1,f1-opts.K) : min(size(corrtable,1),f1+opts.K)
                    if isnan(corrtable(f1,f2, line)) && f1~=f2
                        output = dftregistration(FFTs{f1,line},FFTs{f2,line},2);
                        %xc = normxcorr2_general(F(:,f1),F(:,f2), requiredoverlap);
                        %xc = xcorr(F(:,f1),F(:,f2), opts.window);
                        corrtable(f1,f2, line) = output(1); 
                        corrtable(f2,f1, line) = corrtable(f1,f2, line);
                        shift(f1,f2,line) = output(3); %change this line for 3D!
                        shift(f2,f1,line) = -output(3);
                        %shift(f1,f2,line) = s - size(F,1);
                    end
                end
            end
        end
    end
end