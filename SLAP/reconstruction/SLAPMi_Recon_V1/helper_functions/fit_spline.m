function [F0, DFF]= fit_spline(F,Fitmethod,StimOnset)
if ~exist('Fitmethod') || isempty(Fitmethod)
    Fitmethod = 'pchip';
end

if ~exist('StimOnset') || isempty(StimOnset)
    StimOnset = 1;
end

p = size(F,1);
F0 = zeros(size(F));
for i = 1:p
    F0(i,:) = fitSplineRow(F(i,:),Fitmethod);
end
%% Two-sided
% F00 = F0;
% F01 = F0;
% for i = 1:p
%     F00(i,:) = fitSplineRow(F(i,:),type);
%     F01(i,:) = fliplr(fitSplineRow(fliplr(F(i,:)),type));
% end
% F0 = max(F00,F01);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if StimOnset >1
F0 = max(F0,F0(:,StimOnset));
end
DFF =  1 - F0./F;
DFF(isnan(DFF)) = 0;
DFF = max(0,DFF);
end

function est_baseline = fitSplineRow(F,type)

    if sum(F) < 0.001
        est_baseline = zeros(size(F));
    else

    T = length(F);
    t = 1:T;
    k0 = convhull(t,F);
    k1 = filter([1 -1],1,k0)>0;
    bp = k0(k1);
    est_baseline = interp1(bp,cummin(F(bp)),t,type);

    end  
end