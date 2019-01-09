function [baseline, y_sub] = SubBase( y,h )
 N = size(y,1);
 y_sub =  zeros(size(y));
 y0 = zeros(size(y));
 for n = 1:N
     if nargin < 2
        [y_sub(n,:), y0(n,:)] = sub_baseline_row(y(n,:));
     else
        [y_sub(n,:), y0(n,:)]= sub_baseline_row(y(n,:),h);
     end
 end
baseline = y0-y_sub;
end

function [y_sub,y0] = sub_baseline_row(y,h)
    if nargin < 2
        y0 = y;
    else 
        y0 = filter(h,1,y);
        y0(1:length(h)) = y(1:length(h));
    end
    y_sub = zeros(size(y0));
    T = length(y0);
    t = 1;
    while t<T
       a_temp = [inf(1,t-1),cummin(y0(t:end))];
       [~,I] = min(a_temp);
       y_sub(t:I) = y0(t:I)-a_temp(t:I);
       t = I+1;
    end
end


