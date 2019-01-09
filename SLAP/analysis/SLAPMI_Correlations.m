function SLAPMI_Correlations(dataset)
tau = 10;
theta = [1 -exp(-1/tau)];

%temporary, because there was a bug in how datasets were formed that had 0s where nans should be:
%dataset.dPhotons(isnan(dataset.X)) = nan;

prePeriod = 500:1000;
stimPeriod = 1100:3000;
postPeriod = 3050:3999;

%select ROIs for analysis
%minimum response
meanAmp = nanmean(nanmean(dataset.dPhotons(:,stimPeriod, :),2),3) - nanmean(nanmean(dataset.dPhotons(:,prePeriod, :),2),3);
ampThresh = 0.3;
stims = unique(dataset.stimulus.stim);
nStim = length(stims);
for stimix = nStim:-1:1
    Xn(:,stimix) = sum(isnan(dataset.X(:,1,dataset.stimulus.stim==stims(stimix))),3);
end
Xn = max(Xn,[],2);
planeSize = size(dataset.refIM.IM,1)*size(dataset.refIM.IM,2);
CenterPlaneIndices = false(size(dataset.refIM.seg,1),1);
CenterPlaneIndices((dataset.Zcenter-2)*planeSize+(1:3*planeSize)) = true;
sumCenterPlanes = sum(dataset.refIM.seg(CenterPlaneIndices,:));
sumOuterPlanes = sum(dataset.refIM.seg(~CenterPlaneIndices,:));
segSel = meanAmp(1:end-2)>ampThresh  & Xn(1:end-2)<8;
% & (sumCenterPlanes > prctile(sumCenterPlanes,95))' & (sumCenterPlanes>5*sumOuterPlanes)';
nSeedSel = sum(segSel);
% %visualize selected segments
% figure('name', [int2str(nSeedSel) 'selected segments']),
% imshow3D(reshape(full(sum(dataset.refIM.seg(:, segSel),2)), size(dataset.refIM.IM)))

%% Calculate seed distances
[distances,BinnedDistances,unitDistance] = find_distances(dataset,segSel,nSeedSel);
%% Apply selection
[dataset,LPF,HPF] = trimDataset(dataset,segSel,theta);
%% Calculate correlations


Period = stimPeriod;
Len = length(Period);
analysisMode = 7;
nPC2Remove = 10;
Period2Remove = postPeriod;

switch analysisMode
    % 1: signal correlations
    case 1

c_sig = nan(nSeedSel,nSeedSel,nStim);
X = nan(nSeedSel,Len,nStim);  
for stimix = 1:nStim
    Y = nanmean(dataset.X(:,Period2Remove,dataset.stimulus.stim == stimix),3);
    Coeffs = pca(Y'); ProjMtx = eye(size(X,1))-Coeffs(:,1:nPC2Remove)*Coeffs(:,1:nPC2Remove)';
    X(:,:,stimix) = ProjMtx*(nanmean(dataset.X(:,Period,dataset.stimulus.stim == stimix),3));
%    X(:,:,stimix) = normr(X(:,:,stimix) - nanmean(X(:,:,stimix),2));
    X(:,:,stimix) = normr(X(:,:,stimix) - nanmean(Y,2));
    c_sig(:,:,stimix) = X(:,:,stimix)*X(:,:,stimix)';
end

    case 2
        % 2: noise correlations
    case 3
        X = nan(nSeedSel,Len,nStim);  
% 
%         c_sig = nan(nSeedSel,nSeedSel,nStim);
    % 3: Event Calling
    stdSpikes = std(dataset.spikes(:,Period2Remove,:),[],2);
    dataset.events = dataset.spikes > 5 * stdSpikes;
    Z = (dataset.spikes-nanmean(dataset.spikes(:,Period2Remove,:),2)).*dataset.events;
    Z = filter(1,theta,Z,[],2);
    for stimix = 1:nStim
        X(:,:,stimix) = nanmean(Z(:,Period,dataset.stimulus.stim == stimix),3);
    end
    [tPeak,c_sig] = findPeakCorr(X);
    case 4
    % 4: smooth data
    winlen = 0:10:100;
    nWin = length(winlen);
    c_sig = nan(nSeedSel,nSeedSel,nStim,nWin);
    for w = 1:nWin
        if w == 1
            X = dataset.spikes;
        else
            X = smoothdata(dataset.spikes(:,Period,:),2,'gaussian',winlen(w));
        end
        for stimix = 1:nStim
            X(:,:,stimix) = nanmean(X(:,:,dataset.stimulus.stim == stimix),3);
            X(:,:,stimix) = normr(X(:,:,stimix) - nanmean(X(:,:,stimix),2));
            c_sig(:,:,stimix,w) = X(:,:,stimix)*X(:,:,stimix)';
        end
    end
    
    case 5
    X = nan(nSeedSel,Len,nStim);  
% 
%   c_sig = nan(nSeedSel,nSeedSel,nStim);
    % 3: Event Calling
    stdSpikes = std(dataset.spikes(:,Period2Remove,:),[],2);
    dataset.events = dataset.spikes-mean(dataset.spikes,2) > 5 * stdSpikes;
    Z = (dataset.spikes-nanmean(dataset.spikes(:,Period2Remove,:),2)).*dataset.events;
    Z = filter(1,theta,Z,[],2);

    [tPeak,c_sig] = findPeakCorr(Z(:,Period,:));
%     for stimix = 1:nStim
%         X(:,:,stimix) = nanmean(Z(:,Period,dataset.stimulus.stim == stimix),3);
%     end

    
    case 7
        % frequency
        nFreqBand = 50;
        fc = 5;
        X = dataset.NormalizedX(:,Period,:);
        [P0,centerFreqs,XLP,XHP] = decompose(X,nFreqBand,fc);
        P0 = P0-nanmean(P0,2);
        %%
        P = P0;
%         P = nan(nSeedSel,nFreqBand,nStim);
%         for stimix = 1:nStim
%         P(:,:,stimix) = nanmean(P0(:,:,dataset.stimulus.stim == stimix),3);
%         end
        
        [tPeak,c_sig] = findPeakCorr(P,0);
        [~,c_sig_HP] = findPeakCorr(dataset.XHP,0);
end

%% plot as a function of distance
for i = 1:size(c_sig_HP,4)
plot_CorrVsDist(c_sig_HP(:,:,:,i),BinnedDistances,unitDistance);
end


%%
figure;
A = dataset.Xden; A = A-nanmean(A); A = A - nanmean(A,2);
% A = dataset.NormalizedX- dataset.Xden;
% A = dataset.XHP;
for i  = 1:7
    for j = 1:8
        B = A(:,:,(i-1)*8+j);
        clim = [0,prctile(B(:),90)];
        subplot(7,8,(i-1)*8+j); imagesc(B,clim);
        
    end
end
% suptitle('Normalized Traces - Deoised Normalized Traces: Saturation 90th percentile, 7 trials of each stimulus')
suptitle('Highpass filtered signal: Saturation 90th percentile, 7 trials of each stimulus')

fvtool(HPF,'Fs',1016);
%%
keyboard
%%
figure;
trialNum = 1;
A = dataset.NormalizedX(:,:,trialNum);
B = dataset.Xden(:,:,trialNum);
subplot(2,2,1); imagesc(A);title('X')
set(gca,'fontsize',20)
subplot(2,2,2); 
plot(flipud(normc(nansum(A,2))),1:size(A,1));

hold on;
plot(flipud(normc(nansum(B,2))),1:size(B,1));
title('sum(X,2) and sum(Xden,2) overlayed')
set(gca,'fontsize',20)
axis tight
subplot(2,2,3);
plot(normr(nansum(A)));
hold on;
plot(normr(nansum(B)));
title('sum(X) and sum(Xden) overlayed')
set(gca,'fontsize',20)
axis tight
subplot(2,2,4); imagesc(B);
 title('deniosed X')
set(gca,'fontsize',20)
end


function [d,Binned,unitDistance] = find_distances(dataset,segSel,nSeedSel)
d = nan(nSeedSel);
dims = size(dataset.refIM.IM);
segs = dataset.refIM.seg(:,segSel);
segs = reshape(full(segs),[dims,nSeedSel]);
S = squeeze(sum(segs(:,:,dataset.Zcenter-1:dataset.Zcenter+1,:),3));
S = reshape(S,[dims(1)*dims(2),nSeedSel]);

%% find segment centers
centers = nan(nSeedSel,2);

for seed = 1:nSeedSel
    [cdummyX, cdummyY] = ind2sub(dims(1:2),find(S(:,seed)>0));
    centers(seed,:) = [mean(cdummyX),mean(cdummyY)];
end
centers = 0.2*centers; %in um


disp('Calculating segment distances ...')
for i = 1:nSeedSel
    for j = i+1:nSeedSel
            d(i,j) = norm(centers(i,:)-centers(j,:));
            d(j,i) = d(i,j);
    end
end
disp('done')
if min(d(:))< 0.001
    disp('There are fused into each other segments!');
    keyboard;
end
unitDistance = 5; %um
nBins = 10;
Binned = false(nSeedSel,nSeedSel,nBins);
for bin = 1:nBins
    if bin ~= nBins
        Binned(:,:,bin) = d >= (bin-1)*unitDitance & d < bin*unitDitance;
    else
        Binned(:,:,bin) = d >= (bin-1)*unitDitance;
    end
end

end
function [dataset,LPF,HPF] = trimDataset(dataset,segSel,theta)
%apply selection
dataset.dPhotons = dataset.dPhotons(segSel,:,:);
dataset.dPhotonSpikes = filter(theta,1,dataset.dPhotons,[],2);
dataset.X = dataset.X(segSel,:,:);
dataset.spikes = dataset.spikes(segSel,:, :);

%clipping to first 8 repetitions
triSel = true(1,64); triSel(1:3) = false;
dataset.stimulus.stim = dataset.stimulus.stim(triSel);
dataset.stimulus.timeError = dataset.stimulus.timeError(triSel);
dataset.stimulus.stimTime = dataset.stimulus.stimTime(triSel);
dataset.filenames = dataset.filenames(triSel);


dataset.dPhotons = dataset.dPhotons(:,:,triSel);
dataset.dPhotonSpikes = dataset.dPhotonSpikes(:,:,triSel);
dataset.X = dataset.X(:,:,triSel);
dataset.Z = dataset.Z(triSel);
dataset.stimTime = dataset.stimTime(triSel);
dataset.spikes = dataset.spikes(:,:, triSel);

dataset.NormalizedX = nan(size(dataset.X));
for i= 1:size(dataset.X,3)
    dataset.NormalizedX(:,:,i) = normr(dataset.X(:,:,i));
end

N = 50;
Fs = 1016;
Fp  = 5;
Ap  = 0.01;
Ast = 40;
Rp  = (10^(Ap/20) - 1)/(10^(Ap/20) + 1);
Rst = 10^(-Ast/20);
LPF = firceqrip(N,Fp/(Fs/2),[Rp Rst],'passedge');
HPF = firceqrip(N,Fp/(Fs/2),[Rp Rst],'high');
disp(['Filtering the data at ' num2str(Fp) ' Hz ...'])

dataset.XHP = filter(HPF,1,dataset.NormalizedX,[],2);

dataset.Xden = filter(LPF,1,dataset.NormalizedX,[],2);
for j =1:size(dataset.Xden,3)
dataset.Xden(:,:,j) = fliplr(filter(LPF,1,fliplr(dataset.Xden(:,:,j)),[],2));
end

disp('done')
% dataset.Xden = nan(size(dataset.X));
% disp('wavelet denoising the normalized traces')
% for i = 1:size(dataset.X,1)
%     for j = 1:size(dataset.X,3)
%         dataset.Xden(i,:,j) = cmddenoise(dataset.NormalizedX(i,:,j),'db3',1);
%     end
% end
% disp('done')
end
function meanC = plot_CorrVsDist(C,BinnedDistances,unitDistance)
C(C == 0) = nan;
nPlots = size(C,3);
nDist = size(BinnedDistances,3);
legendStr = {};
figure; hold on;
meanC = nan(nPlots,nDist);
for p = 1:nPlots
    for d = 1:nDist
        cdummy = C(:,:,p); cdummy = cdummy(BinnedDistances(:,:,d)); 
        meanC(p,d) = nanmean(cdummy);
    end
    plot((1:nDist)*unitDistance,meanC(p,:),'linewidth',3);
    legendStr{end+1} = ['Stimulus # ' num2str(p)];
end

xlabel('segment distance (\mum)')
ylabel('mean correlation');
legend(legendStr);
set(gca,'fontsize',20);

end
function [tPeak,cPeak] = findPeakCorr(X,tmax)
if nargin<2
    tmax = 100;
end
disp('Finding peak correlations ...')
tPeak = zeros(size(X,1),size(X,1),size(X,3));
cPeak = zeros(size(tPeak));
% X = X - nanmean(X,2);
% disp('Subtracting spatial mean from X')
% X = X - nanmean(X);
disp('Subtracting temporal mean from X')
X = X - nanmean(X,2);

Ct = nan(size(cPeak)); 
for t = -tmax:tmax
    for ix = 1:size(X,3)
        X1 = normr(X(:,abs(t)+1:end,ix));
        Y1 = normr(X(:,1:end-abs(t),ix));
%         X1 = (X(:,abs(t)+1:end,ix));
%         Y1 = (X(:,1:end-abs(t),ix));
%         X1 = normr(X(:,abs(t)+1:end,ix)-nanmean(X(:,abs(t)+1:end,ix),2));
%         Y1 = normr(X(:,1:end-abs(t),ix)-nanmean(X(:,1:end-abs(t),ix),2));
        if t>=0
        Ct(:,:,ix) = X1*Y1';
        else
        Ct(:,:,ix) = Y1*X1';
        end
    end
    Indices = abs(cPeak) < abs(Ct);
    tPeak(Indices) = t;
    cPeak(Indices) = Ct(Indices);
    
    c(t+101) = normr(X(3,abs(t)+1:end,1))*normr(X(3,1:end-abs(t),1))';
%   c(t+101) =  Ct(3,4,1);
    d(t+101) = Ct(4,3,1);
end
% figure; plot(-100:100,c); hold on; plot(-100:100,d);
disp('done');

end
function [P,centerFreqs,XLP,XHP] = decompose(X,nFreqBand,fc)
disp('Filtering the data ...')
Fs = 1016;
fmax = 100;
f1 = 1.3333;
winlen = floor(Fs/f1*2);
Y = fft(X,winlen,2);
f = Fs*linspace(0,1/2,winlen/2);
f = [f,fliplr(f)];
Bands = [0, f1/2, linspace(3*f1/2,fmax,nFreqBand-1)];
P = nan(size(Y,1),nFreqBand,size(Y,3));
    for FreqBand = 1:nFreqBand
        YBand = Y(:,f >= Bands(FreqBand) & f < Bands(FreqBand+1),:);
        P(:,FreqBand,:) = squeeze(mean(abs(2*YBand/winlen).^2,2));
    end
    centerFreqs = Bands(1:end-1)+diff(Bands/2);

Y = fft(X,[],2);
f = Fs*linspace(0,1/2,floor(size(X,2)/2));
f = [f,fliplr(f)];
YLP = Y; YLP(:,f > fc,:) = 0;
YHP = Y; YHP(:,f<= fc | f > fmax,:) = 0;
XLP = real(ifft(YLP,[],2));
XHP = real(ifft(YHP,[],2));
    
end
function transform_X(X)
Y = X;
for i = 1:size(X,3)
    Y(:,:,i) = normr(X(:,:,i));
end

end