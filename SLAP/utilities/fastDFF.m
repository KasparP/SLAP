function [F, dFF] = fastDFF (rawF, tau,binTo)
%does bleaching correction and deltaF/F calculation, 

lxls = size(rawF,1);
window1 = min(201, size(rawF,2));
%window2 = 21;
if nargin <3
binTo = 600;
end
if nargin<2 || isempty(tau)
    tau = 20;
end
kernel = exp(-(0:6*tau)/tau);
binfactor = floor(lxls/binTo);
binned = ceil(lxls/binfactor);

dFF = nan(binned, size(rawF,2));
F0 = nan(binned, size(rawF,2));
F = nan(binned, size(rawF,2));
for i = 1:binned
    Z = nanmean(rawF((i-1)*binfactor+1:min(end, i*binfactor),:));
    Z(isnan(Z)) = 0;
    baseline = smoothfast2(Z, window1);
    baseline(1:window1) = max(baseline(1:window1), baseline(window1));
    baseline= cummin(baseline);
    F0(i,:) = max(baseline, 1); %5e-1);
%     baseline = cummin(medfilt1([inf(1,floor(window1/3)) Z], window1));
%     baseline = baseline(floor(window1/3)+1:end);
%     F0(i,:) = smoothfast2(baseline, window2);
    V = (Z-F0(i,:)) ./ F0(i,:);
    V = conv(NN_KP(V',tau), kernel);  
    %V = filter(1,[1 -kernel(2)],NN_KP(V',tau), [],2);
    dFF(i,:) = V(1:size(rawF,2));
    F(i,:) = Z;
end



end

function Y = smoothfast2 (X, m) %#ok<DEFNU> %this algorithm is slightly faster than smoothfast1 but might have numeric issues and does suboptimal things at the data edges
mm = 2*m+1;
Y = [repmat(X(1),m,1) ; X(:) ; repmat(X(end),m,1)] ;
Y = [0 ; cumsum(Y)] ;
Y = (Y(mm+1:end)-Y(1:end-mm)) / mm ;
end