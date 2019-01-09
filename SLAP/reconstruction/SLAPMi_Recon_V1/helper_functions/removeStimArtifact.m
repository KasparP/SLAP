function scandata = removeStimArtifact(scandata,stimLength)

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

y(:,select) = max(0, y(:,select) - repmat(A(select), size(y,1),1));

mfilt = medfilt2(y, [1 3]);
replace = (select ~= [false select(1:end-1)])  | (select ~= [select(2:end) false]);
y(:,replace) = mfilt(:,replace);
y(nans) = nan;

scandata.stimOnset = find(select,1,'first');

for f = 1:length(scandata.frames)
    scandata.frames(f).pmtData = y(:,f);
end