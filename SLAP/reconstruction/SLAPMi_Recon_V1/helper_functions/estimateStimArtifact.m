function [stimArtifact, stimStart] = estimateStimArtifact(scandata,stimLength)
%how long was the artifact (roughly)?
t = 2*stimLength; %number of frames

%rough timing
y = [scandata.frames.pmtData];
ysum = nanmedian(y,1);

nans = isnan(y);
y(nans) = 0;
A = imtophat(ysum, ones(1,t*2));
nA = length(A);
select = medfilt2(A>prctile(A, 100*(1-(t-1)/nA))*0.8, [1 ceil(t/4)*2+1]);

%subtract the trend from y
stimStart = find(select,1,'first')-1; stimEnd = find(select,1,'last')+1;
pre = mean(y(:, stimStart-79:stimStart-1),2);
post= mean(y(:, stimEnd+1:stimEnd+79),2);
nulls = interp1([0 1], [pre post]', linspace(0,1, stimEnd-stimStart+80+1))';
obs = y(:, stimStart-40:stimEnd+40);
DD = obs-nulls;
nanrows = all(obs==0,2);

%median filter the mean subtracted to get a 'support region'
support = medfilt2(medfilt2(medfilt2(DD, [100,1], 'symmetric')>0.6, [100 1], 'symmetric'),[3 3], 'symmetric');
support = imdilate(support, ones(200,1));
support(nanrows,:) = false;
%average the rate within the support region for each frame
SA = zeros(size(support));
for col = 1:size(support,2)
    dcol = DD(support(:,col), col);
    dcol(dcol>prctile(dcol, 90)) = prctile(dcol, 90);
    dcol(dcol<prctile(dcol, 10)) = prctile(dcol, 10);
    dcol = imgaussfilt(dcol, [200 1], 'padding', 'symmetric');
    SA(support(:,col),col) = dcol;
end

figure('name', 'Estimated Stimulus Artifact'), imagesc(SA);
stimArtifact = zeros(size(y));
stimArtifact(:, stimStart-40:stimEnd+40) = SA*1.01;
