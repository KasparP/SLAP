function ysmooth = smoothCircle(y, width)
%accepts a column vector

% ysmooth = zeros(size(y,1),2*width+1);
% for i = 1:2*width+1
%     ysmooth(:,i) = circshift(y, [i-width-1 0]);
% end
% ysmooth = mean(ysmooth,2);

ysmooth = y;
for i = 1:width
    ysmooth = ysmooth + circshift(y, [i 0]) + circshift(y, [-i 0]);
end
ysmooth = ysmooth/(2*width+1);



% [X,Y] = ndgrid(1:length(y), -width:width);
% ixs = mod(X+Y-1, length(y))+1;
% ysmooth = mean(y(ixs),2);