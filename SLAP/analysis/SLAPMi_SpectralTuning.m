function  SLAPMi_SpectralTuning(dataset)
tau = 10;
theta = [1 -exp(-1/tau)];

%temporary, because there was a bug in how datasets were formed that had 0s where nans should be:
%dataset.dPhotons(isnan(dataset.X)) = nan;

prePeriod = 500:1000;
% prePeriod = 3050:3999;
stimPeriod = 1100:3000;

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
segSel = meanAmp>ampThresh  & Xn<3; 
segSel = segSel(1:end-2);
nSeedSel = sum(segSel);
%visualize selected segments
figure('name', [int2str(nSeedSel) 'selected segments']),
imshow3D(reshape(full(sum(dataset.refIM.seg(:, segSel),2)), size(dataset.refIM.IM)))

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


rad = linspace(0,2*pi,nStim/2+1);
rad = rad([1:4 1:4]);


Fs = 1016; 
nFreqBand = 20;
stimR = dataset.dPhotons;
[ORI, centerFreqs] = getSpectralTuning(stimR);


% random shuffling ?
doShuffle = 10;
if doShuffle
    ORI_shuff = nan([size(ORI) doShuffle]);
    for rep = 1:doShuffle
        disp(['Shuffle Repeat: ' int2str(rep)])
        ORI_shuff(:,:,rep) = getSpectralTuning(stimR(:,:,randperm(size(stimR,3))));
    end 
end
if doShuffle
    ORI_shuff_spikes = nan([size(ORI) doShuffle]);
    for rep = 1:doShuffle
        disp(['Shuffle Spikes Repeat: ' int2str(rep)])
        dPhotonSpikes = randpermND(dataset.dPhotonSpikes,stimPeriod);
        stimRshuff = filter(1,theta,dPhotonSpikes,[],2);
        %stimRshuff = randpermND(stimR,stimPeriod);
        ORI_shuff_spikes(:,:,rep) = getSpectralTuning(stimRshuff);
    end 
end


figure, plot(centerFreqs, nanmean(ORI,1)), hold on, plot(centerFreqs,nanmean(nanmean(ORI_shuff,1),3))
xlabel('frequency (Hz)');
ylabel('mean ORI index');
legend({'Measured', 'Shuffled Trials'})

figure, plot(centerFreqs, nanmean(ORI,1)), hold on, plot(centerFreqs,nanmean(nanmean(ORI_shuff_spikes,1),3))
xlabel('frequency (Hz)');
ylabel('mean ORI index');
legend({'Measured', 'Shuffled Spikes within Trial'})

keyboard


    function [S, centerFreqs] = getSpectralTuning(stimR)
        [PStim,centerFreqs] = decompose(stimR(:,stimPeriod,:),Fs,nFreqBand);
        %PPre = decompose(stimR(:,prePeriod,:),Fs,nFreqBand);
        stimMeanFreq = nan(size(stimR,1), nStim,nFreqBand);
        for stimix = 1:nStim
%             stimMeanFreq(:,stimix,:) = nanmean(PStim(:,dataset.stimulus.stim==stims(stimix),:),2)...
%                 - nanmean(PPre(:,dataset.stimulus.stim==stims(stimix),:),2);
             stimMeanFreq(:,stimix,:) = nanmean(PStim(:,dataset.stimulus.stim==stims(stimix),:),2);
        end
        
        stimMeanFreq = stimMeanFreq - min(stimMeanFreq,[],2);  %ABBAS SEE THIS LINE
        H = nan(nSeedSel,nFreqBand);
        S = H; V = H;
        for freqix = 1:nFreqBand
            [H(:,freqix),S(:,freqix),V(:,freqix)] = circular_mean(repmat(rad,nSeedSel,1),stimMeanFreq(:,:,freqix));
        end
    end
end


function X = randpermND(X,period)
    for i = 1:size(X,1)
        for j = 1:size(X,3)
            X(i,period,j) = X(i,period(randperm(length(period))),j);
        end
    end
end

function [H,S,V] = circular_mean(rad,w)
% w = w - min(w,[],2);
A = sum(exp(1i*rad).*w,2);
H = 1/2 + angle(A)/(2*pi);
S = abs(A)./sum(abs(w),2);
% w = w - nanmean(w,2);
V = sqrt(sum(w.^2,2));
end








% function [LX, HX] = decompose(X,Fs)
% disp('Filtering the data ...')
% if nargin<2; Fs = 1016; end
% Ts = 1/Fs;
% T = size(X,2);
% % t = (0:T-1)*Ts;
% Y = fft(X,[],2);
% f = Fs*linspace(0,1/2,T/2+1);
% fc = 10;
% if mod(T,2) == 0
%    f = [f,fliplr(f(2:end-1))];
% else
%    f = [f,fliplr(f(2:end))];
% end
% 
% YLX = Y; YLX(:,f>fc,:)=0;
% YHX = Y; YHX(:,f<=fc,:)=0;
% 
% LX = real(ifft(YLX,[],2));
% HX = real(ifft(YHX,[],2));
% 
% end

function [P,centerFreqs] = decompose(X,Fs,nFreqBand)
disp('Filtering the data ...')

Ts = 1/Fs;
T = size(X,2);
fc = 10;
fmax = 50;
% f1 = 3.41/4;%data since 10/1
f1 = 1.3333;

% t = (0:T-1)*Ts;
winlen = floor(Fs/f1*2);
Y = fft(X,winlen,2);
f = Fs*linspace(0,1/2,winlen/2);

f = [f,fliplr(f)];

Bands = [0, f1/2, linspace(3*f1/2,fmax,nFreqBand-1)];
%XB = nan([size(X),nFreqBand]);
% varargout = cell(nFreqBand,1);
P = nan(size(Y,1),size(Y,3),nFreqBand);
    for FreqBand = 1:nFreqBand
%         Y0 = Y;
%         Y0(:,f < Bands(FreqBand) | f > Bands(FreqBand+1),:) = 0;
        YBand = Y(:,f >= Bands(FreqBand) & f < Bands(FreqBand+1),:);
%         XB(:,:,:,FreqBand) = real(ifft(Y0,[],2));
        P(:,:,FreqBand) = squeeze(mean(abs(YBand).^2,2));
%           P(:,:,FreqBand) = squeeze(sqrt(sum(abs(YBand/winlen).^2,2)));
    %     Yband = Y;
    %     Yband(:,~(f> Bands(FreqBand) & f <= Bands(FreqBand+1)),:) = 0;
    %     varargout{FreqBand} = real(ifft(Yband,[],2));
    end
    centerFreqs = Bands(1:end-1)+diff(Bands/2);
end