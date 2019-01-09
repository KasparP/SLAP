function [mRaw,mDFF,valid2D,S2D,Svalid] = calc_ZMovies(sys,drawframes,drawPixels,StimOnset)
fprintf('Calculating raw movie... ');
if  ~exist('drawframes') || isempty(drawframes)
    drawframes = 1:sys.opts.T;
end

%% Reconstructed Movie
sz = size(sys.input.ref_image);
if sys.opts.offFocusSpike
    SOutsideMask = full(any((repmat(sys.output.SLM_mask(:)<0.8, [sz(3) 1])).*sys.input.S(:,1:end-2), 1));
else
    SOutsideMask = full(any((repmat(sys.output.SLM_mask(:)<0.8, [sz(3) 1])).*sys.input.S(:,1:end-1), 1));
end
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
%%%%
% HessF = Hessian(sys.output.PS,sys.input.y,sys.output.F,sys.control_params.rate);
% sys.output.F = HessF;
%%%%

mDFF = Svalid*sys.output.F(~SOutsideMask,:);

if ~exist('StimOnset') || isempty(StimOnset)
    baselineX =  min(sys.output.F(:, 20:end-10),[],2);
else
    baselineX =  mean(sys.output.F(:, StimOnset-40:StimOnset-1),2); %for 2-spot uncaging
%     baselineX =  sys.output.F(:, StimOnset); %for 2-spot uncaging
end
baselineX = repmat(baselineX,1,size(sys.output.F,2));
baseline = Svalid*(1+baselineX(~SOutsideMask,:));
mDFF = (mDFF-baseline)./(1+baseline);

mRaw = Svalid*(1+sys.output.F(~SOutsideMask,:));

mRaw = mRaw(:,drawframes);
mDFF = max(0, mDFF);
mDFF = mDFF(:,drawframes);


fprintf('done\n')
    
end
function [HessF,Sym] = Hessian(PS,y,F,rate)
A = mean(y,2)./mean(rate,2).^2;
S = diag(A);
Q = PS'*S*PS;
[V,D] = eig(Q);
% D = diag(D)+0.01;
% d = diag(1./D/norm(1./D));
D = sqrt(diag(D));
d = diag(D/norm(D));
Sym = V*d*V';
Sym = (Sym+Sym')/2;
HessF = Sym*F;
HessF(HessF<0) = 0;
% HessF = sqrt(PS'*(y./(rate.^2).*(PS*F.^2)));
% HessF = PS'*(sqrt(y)./rate.*(PS*F));
end



