function [d,cc,peak_times] = SLAPMi_Corr_vs_Dist(TuningDataset,FusedInfoDataset)
dr = pwd;
[fnSeg, drSeg] = uigetfile([dr filesep '*.mat'], 'select your reference Segmentation _SEG file');
load([drSeg filesep fnSeg])

S0 = refIM.seg;
dims = size(refIM.IM);
S = 0;
for plane = (dims(3)+1)/2-1:(dims(3)+1)/2+1
    S = S + S0((plane-1)*dims(1)*dims(2)+(1:dims(1)*dims(2)),:);
end
uniqueSegs = unique(FusedInfoDataset.fusedInto(:)); uniqueSegs(1) = [];

Spikes = TuningDataset.spikes(1:end-2,:,:);
GaussianWidth = 10;
% find segment distances;
nSegs = size(S,2);
nStim = size(Spikes,3);
xspacing = .2; % in um
yspacing = .2;
maxlag = 100;

%filtering
filtering = 'lpf';

switch filtering
    case 'none'
X = imgaussfilt3(Spikes,[eps,GaussianWidth,eps],'Padding','symmetric'); %nans shouldn't be a problem ?
X(:,1:2*GaussianWidth,:) = [];
X(:,end-2*GaussianWidth+1:end,:)=[];

    case 'lpf'
%% Filter the data Now
% LPF
fs = 1016; %Hz
f1 = 2; %Hz
order = 2;
framelen = 100*f1+1; 
X = sgolayfilt(Spikes,order,framelen,[],2);
X = permute(reshape(X,nSegs,nStim,[]),[1 3 2]);
X(:,1:5*GaussianWidth,:) = [];
X(:,end-5*GaussianWidth+1:end,:)=[];

%%%%%%%%
    case 'hpf'
% HPF
[HD HR] = wfilters('db45','h');

    otherwise 
        error('Filter type not specified')
end


%%
disp('Calculating segment centers ...')
centers = nan(nSegs,2);
for i = 1:length(uniqueSegs)
    seg = uniqueSegs(i);
    [cdummyX, cdummyY] = ind2sub(dims(1:2),find(S(:,seg)>0));
    centers(seg,:) = [mean(cdummyX),mean(cdummyY)];
end
disp('done');
centers = centers.*[xspacing, yspacing];

d = nan(nSegs,nSegs); 
cc = nan(nSegs,nSegs,nStim); 
peak_times = nan(nSegs,nSegs,nStim);

disp('Calculating segment distances ...')
for i = 1:length(uniqueSegs)
    seg1 = uniqueSegs(i);
    for j = i+1:length(uniqueSegs)
        seg2 = uniqueSegs(j);
        if seg2 > seg1
            d(seg1,seg2) = norm(centers(seg1,:)-centers(seg2,:));
            d(seg2,seg1) = d(seg1,seg2);
        end
    end
end
disp('done')



X = X - nanmean(X,2);
V = sqrt(size(X,2))*nanstd(X,[],2);
X0 = X./V;
% X0 = nan(size(X));
% for stim = 1:nStim
%     X0(:,:,stim) = normr(X(:,:,stim));
% end
%%
for stim = 1:nStim
    disp(['Calculating peak correlations and timing for stimulus ' num2str(stim) '...']);
    indices = find(FusedInfoDataset.fusedInto(stim,:)>0);
    indices = unique(indices);
    cc(indices,indices,stim) = X0(indices,:,stim)*X0(indices,:,stim)';
    cc(indices,indices,stim) = cc(indices,indices,stim) - diag(nan*diag(cc(indices,indices,stim)));
    peak_times(indices,indices,stim) = 0;
    for lag = -maxlag:GaussianWidth:maxlag
        if lag<0
            cc_temp = (X0(:,1-lag:end,stim))*(X0(:,1:end+lag,stim))';
        elseif lag>0
            cc_temp = (X0(:,1:end-lag,stim))*(X0(:,1+lag:end,stim))';
        end
%         cc_temp = cc_temp - diag(nan*diag(cc_temp)); not necessary
        cc_temp(isnan(cc(:,:,stim))) = nan;
        OverWriteIndices = cc(:,:,stim)< cc_temp;
        cc(:,:,stim) = max(cc_temp,cc(:,:,stim));
        
        peak_times_temp = peak_times(:,:,stim); 
        peak_times_temp(OverWriteIndices) = lag;
        peak_times(:,:,stim) = peak_times_temp;
    end
end
%%
avg_peak_times = nanmean(peak_times,3);
avg_cc = nanmean(cc,3);
% binSize = 1; % in um
nBinsD = 100;
I = zeros(2*maxlag/GaussianWidth+1,nBinsD);
dBinned = linspace(0,max(d(:)),nBinsD+1);
for bin = 1:nBinsD
    Indices = d > dBinned(bin) & d <= dBinned(bin+1) & ~isnan(avg_peak_times);
    I(floor(avg_peak_times(Indices)/GaussianWidth)+maxlag/GaussianWidth+1,bin) = avg_cc(Indices);
end

figure;
imagesc(I');
set(gca, 'XTick', 1:21, 'XTickLabel', -maxlag:GaussianWidth:maxlag) % 10 ticks 
set(gca, 'YTick', 1:5:nBinsD, 'YTickLabel', dBinned(1:5:nBinsD)) % 20 ticks
xlabel('delay (ms)')
ylabel('distance (um)')
set(gca,'fontsize',20);




end
