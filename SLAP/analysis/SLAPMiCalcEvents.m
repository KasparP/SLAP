function dataset = SLAPMiCalcEvents(dataset,prePeriod,stimPeriod,insideSLMsegs)

%% calc F0
X = dataset.X;
X(~insideSLMsegs,:,:) = nan;

% H = zeros(size(X,1),size(X,3));
% H(isnan(squeeze(sum(X,2)))) = nan;
% for i=1:size(X,1)
%     for j=1:size(X,3)
%         if ~isnan(H(i,j))
%         H(i,j) = ttest(X(i,1000,j),0,'Alpha',0.01);
%         end
%     end
% end

stdX = squeeze(std(X,[],2));
goodSeeds = sum(isnan(stdX),2) == 0;
X(~goodSeeds,:,:) = nan;
% X = X(~goodSeeds,:,:);
% stdX =stdX(goodSeeds,:);

Y = reshape(permute(X,[1 3 2]),size(X,1)*size(X,3),[]);
Z = isnan(sum(Y,2));
YY = Y(~Z,:);

sn = estimate_noise(YY(:,stimPeriod));
prc = prctile(YY(:,stimPeriod),5,2);

r = reshape(sn./prc,[],size(X,3));


totalF = squeeze(nansum(X(:,stimPeriod,:)));

Xs = filter(ones(20,1)/20,1,X,[],2);
fstd = squeeze(nanstd(Xs,[],1)./nanmean(Xs,1));

% baseline = estimate_base(YY(:,stimPeriod),noise);
% sn2m = sn./mean(YY(:,stimPeriod),2);


keyboard
%% find bad trials
badTrials = nan(size(Z));
badTrials(~Z) = sn2m;
badTrials = reshape(badTrials,size(X,1),size(X,3));

keyboard;


keyboard
%%
YS = zeros(size(Y));
tc = 15;
for i = 1:size(YY,1)
    i
    YS(i,:) = NND (YY(i,:)', tc);
end
YSmoothed = nan(size(Y));
YSmoothed(~Z,:) = YS;
Spikes = permute(reshape(YSmoothed,size(X,1),size(X,3),size(X,2)),[1,3,2]);
% sys = genFadeSys(YY,stimPeriod,tc);
% sys_fade = FADE(sys);
keyboard;
%%




F0tmp = repmat(min(dataset.X(:, prePeriod,:),[],2), 1, size(dataset.X,2),1);
dataset.dFF = max(0, (dataset.X-F0tmp)./(F0tmp+1)); %this absorbs any resting F0 into the existing dFF calculation
dataset.dFF = 2*dataset.dFF./(max(dataset.dFF,[],2)+1); %this normalizes the amplitudes (while suppressing very flat traces)
keyboard

end

function sys = genFadeSys(YY,stimPeriod,tc)
sys.y = YY(:,stimPeriod);
sys.Order = 1;
sys.theta = [1 -exp(-1/tc)];
sys.num_iters = 50;
sys.min_iters = 30;
end

function sigma = estimate_noise(y)
T = size(y,2);
range_ff = [0.25; 0.5];
[pyy,ff]=pwelch(y',round(T/8),[],1000,1);
ind = ff > range_ff(1) & ff < range_ff(2);
sigma = sqrt(exp(mean(log(pyy(ind,:)/2))))';
end
