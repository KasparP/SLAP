function [mRaw,mDFF,valid2D,S2D,Svalid] = calc_Zscore(sys,drawframes,drawPixels,StimOnset)
fprintf('Calculating raw movie... ');
if  ~exist('drawframes') || isempty(drawframes)
    drawframes = 1:sys.opts.T;
end

%% Reconstructed Movie
sz = size(sys.input.ref_image);
SOutsideMask = sys.output.OutsideSLMSeeds;
valid3D = reshape(full(sum(sys.input.S(:,1:end-2),2)>eps), sz)  & repmat(sys.output.SLM_mask, [1 1 sz(3)])>0.8;
ref2D = sys.input.ref_image; ref2D(~valid3D) = 0; 
ref2D = sum(ref2D,3);
if ~exist('drawPixels') || isempty(drawPixels)
    drawPixels = true(size(ref2D));
end

S2D = 0*sys.input.S(1:(sz(1)*sz(2)),:);
for plane = 1:sz(3)
    S2D = S2D +  sys.input.S((plane-1)*(sz(1)*sz(2)) + (1:(sz(1)*sz(2))),:).*reshape(valid3D(:,:,plane), [],1);
end
valid2D = any(valid3D,3) & drawPixels; %pixels within the mask
Svalid = S2D(valid2D(:),~SOutsideMask);

% normalization

%%

if ~exist('StimOnset') || isempty(StimOnset)
    baselineX =  prctile(sys.output.F(:, 20:end-10),5,2);
else
    baselineX =  mean(sys.output.F(:, StimOnset-40:StimOnset-1),2); %for 2-spot uncaging
end

DFFX = (sys.output.F-baselineX)./(baselineX+1);
DFFX = max(DFFX,0);

valid = sys.input.ref.fusedInto>0;
PS = full(sum(sys.output.PS)');
dPhotons = zeros(size(sys.output.F)); Zscore = dPhotons;
dPhotons(valid,:) = PS(valid).*DFFX(sys.input.ref.fusedInto(valid), :);
Zscore(valid,:) = dPhotons(valid,:)./sqrt(3/8+PS(valid).*(baselineX(sys.input.ref.fusedInto(valid))+1));

normalizeSvalid = (Svalid./sum(Svalid,2));
normalizeSvalid(isnan(normalizeSvalid)) = 0;
mDFF = normalizeSvalid*Zscore(~SOutsideMask,:);
% mRaw = repmat(ref2D(valid2D),1,sys.opts.T)*sys.opts.baselineMult + Svalid*(sys.output.F(~SOutsideMask,:));
mRaw = repmat(ref2D(valid2D),1,sys.opts.T)*sys.opts.baselineMult + Svalid*(baselineX(~SOutsideMask,:));
mRaw = mRaw(:,drawframes);
mDFF = max(0, mDFF);
mDFF = mDFF(:,drawframes);
    
fprintf('done\n')
    
end