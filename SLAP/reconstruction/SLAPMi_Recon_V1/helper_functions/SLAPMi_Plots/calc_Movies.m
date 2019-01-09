function [mRaw,mDFF,valid2D,S2D,Svalid] = calc_Movies(sys,drawframes,drawPixels,StimOnset)
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

mDFF = Svalid*sys.output.F(~SOutsideMask,:);
if ~exist('StimOnset') || isempty(StimOnset)
    baselineX =  prctile(sys.output.F(:, 20:end-10),5,2);
else
    baselineX =  mean(sys.output.F(:, StimOnset-40:StimOnset-1),2); %for 2-spot uncaging
%     baselineX =  sys.output.F(:, StimOnset); %for 2-spot uncaging
end
baselineX = repmat(baselineX,1,size(sys.output.F,2));
baseline = Svalid*(baselineX(~SOutsideMask,:));
    
if all(sys.output.additive_baseline(:)<=(sys.opts.dark_noise+eps))
    %as a baseline we are using the average of the measured baseline and the expected baseline from the reference image
    mDFF = (mDFF-baseline)./((repmat(sum(Svalid(:,1:end-2),2),1,sys.opts.T) + baseline)/2);
else
    %under the additive baseline model the additive baseline roughly corresponds to the sum of the reference image, hence:
    mDFF = (mDFF-baseline)./(repmat(sum(Svalid(:,1:end-2),2),1,sys.opts.T) + baseline);
end

mRaw = repmat(ref2D(valid2D),1,sys.opts.T)*sys.opts.baselineMult + Svalid*(sys.output.F(~SOutsideMask,:));

mRaw = mRaw(:,drawframes);
mDFF = max(0, mDFF);
mDFF = mDFF(:,drawframes);

% X = sys.output.F(~SOutsideMask,drawframes);
% 
% if sys.opts.perSeedBaselineFlag
%     mRaw = Svalid*((1+X).*sys.F0(~SOutsideMask,drawframes));
% else
%     % Calculate Raw Movie
%     mRaw = Svalid*X + ref2D(valid2D)*sys.opts.baselineMult;   
% end

%     if exist('StimOnset')
%         baselineX = X(:,StimOnset);
%         F0 = mRaw(:,StimOnset);
%         F0 = repmat(F0,1,size(mRaw,2));
%     %   F0 = Svalid*((1+baselineX).*sys.F0(~SOutsideMask,drawframes));
%         mDFF = max(0,mRaw./F0-1);
% 
%     else
%         F0 = min(mRaw,[],2);
%         mDFF = max(0,mRaw./F0-1);
%     end
    
fprintf('done\n')
    
end