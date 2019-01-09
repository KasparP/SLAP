function Plot_Uncaging_LoogerLab_DFF(sys)
StimOnset = sys.opts.StimOnset;
tmin = StimOnset; %400
tmax = 550;

sz = size(sys.input.ref_image);
valid3D = reshape(full(sum(sys.input.S(:,1:end-1),2)>eps), sz)  & repmat(sys.output.SLM_mask, [1 1 sz(3)])>0.8;
ref2D = sys.input.ref_image; ref2D(~valid3D) = 0; ref2D = sum(ref2D,3);

S2D = 0*sys.input.S(1:(sz(1)*sz(2)),:);
for plane = 1:sz(3)
    S2D = S2D +  sys.input.S((plane-1)*(sz(1)*sz(2)) + (1:(sz(1)*sz(2))),:).*reshape(valid3D(:,:,plane), [],1);
end
valid2D = any(valid3D,3); %pixels within the mask
Svalid = S2D(valid2D(:),:);

fprintf('Calculating raw movie ... ');
drawframes = 1:size(sys.output.F,2);
 
mDFF = Svalid*sys.output.F;
baselineX = nan(size(sys.output.F));
for i = 1:size(sys.output.F,1)
    baselineX(i,:) =  mean(sys.output.F(i, 400), 2); %for 2-spot uncaging
end
baseline = single(Svalid*baselineX);

mDFF = (mDFF-baseline)./(repmat(full(sum(Svalid(:,1:end),2)),1,sys.opts.T) + baseline);
mDFF(mDFF<0.2) = 0;
mDFF = max(0, mDFF);
mDFF = mDFF(:,drawframes);
clear baseline
fprintf('done. \n');

%%
maxDFF = max(mDFF(:, tmin:tmax),[],2);
minDFF = min(mDFF(:, tmin:tmax),[],2);
MM = max(maxDFF);
mm = min(minDFF);
% minDFF = 0*mDFF(:,tmin);
% maxDFF_time = inf(size(mDFF,1),1);
% for i = 1:sum(valid2D(:))
%     T= find(mDFF(i,tmin:tmax)>= minDFF(i) + 0.67*(maxDFF(i)-minDFF(i)),1,'first');
%     if ~isempty(T)
%        maxDFF_time(i) = T;
%     end
% end

%%
%Prepare data
% maxDFF_time = maxDFF_time-min(maxDFF_time);
% tmin = 35; tmax = 100;
% % cmin = 0.4; 
% cmin = 0.6;
% cmax = 1;
% H = cmin + (cmax-cmin)*max(0,min(1, (maxDFF_time-tmin)./(tmax-tmin))); %10ms->0.4 50ms->1
% H(H>=cmin) = 1+cmin-H(H>=cmin);
H = ones(size(maxDFF));
S = min(1, maxDFF/1);
% S(S>0) = 1;
V = Svalid*(sys.output.F(:,1)+1);
ff = max(V)/4; 
V = min(1, (V./ff));

rgb_image = hsv2rgb([H, S, V]);

rgb_im = repmat(full(((S2D(:,end)*(sys.output.F(end,1)+1))./ff)),1,3);
rgb_im(valid2D,:) = rgb_image;
rgb_im = reshape(rgb_im, sys.opts.dim(1), sys.opts.dim(2), 3);
rgb_im = rgb_im.*repmat((sys.output.SLM_mask), 1, 1, 3);
%-------%
% for 5-14-FOV1-10 ms
uncagingx = round([791.2340  506.0005]);
uncagingy = round([508.4363  801.8780]);
% for 5-14-FOV2-10 ms
% uncagingx = round([529  780]);
% uncagingy = round([813  479]);
inds = sub2ind(sys.opts.dim, uncagingy, uncagingx);
Spts = 0*sys.input.S(inds,:);
for plane = 1:sz(3)
    Spts = Spts +  sys.input.S((plane-1)*(sz(1)*sz(2))+inds,:);
end
ptDFF = Spts*sys.output.F;
B = [ptDFF(1,410)  ; ptDFF(2,400) ];
ptDFF = max(0,(ptDFF-B)./(sum(Spts,2)+B));

%%% figure 1
% figure, plot((-399:1000)*(1/1.016), ptDFF(:, 1:1400)', 'linewidth', 2)
% legend(strcat({'Uncaging Location '}, int2str([1:length(uncagingx)]')));

%%%%%%%%%%%%%%%
% [dffMax,tMax] = max(ptDFF,[],2); 
% Stim2Rise1 = ptDFF(1,401:tMax(1));
% Stim2Rise2 = ptDFF(2,401:tMax(2));
% df1 = abs((Stim2Rise1-dffMax(1)/2));
% df2 = abs((Stim2Rise2-dffMax(2)/2));
% [m1,halfRiseTime1] = min(df1);
% [m2,halfRiseTime2] = min(df2);
%%%%%%%%
% time2peakDelay = halfRiseTime1-halfRiseTime2
% indices = 376:500;
%%%% Zoomed
%%% figure 2
% figure, plot((indices-400)*(1/1.016), ptDFF(:, indices)', 'linewidth', 2);
% hold on; 
% line([halfRiseTime1, halfRiseTime1],[1 0]);
% line([halfRiseTime2, halfRiseTime2],[1 0]);
% legend(strcat({'Uncaging Location '}, int2str([1:length(uncagingx)]')));

%%
normalizedDFF = (ptDFF-min(ptDFF,[],2))./(max(ptDFF,[],2)-min(ptDFF,[],2));
% %%% figure 3
% figure, plot((indices-400)*(1/1.016), normalizedDFF(:, indices)', 'linewidth', 2)
% legend(strcat({'Uncaging Location '}, int2str([1:length(uncagingx)]')));


%-------%
%kymograph

%axis of kymograph
% for 5-14-FOV1-10 ms
% axisX =   [581.2963 441.9052];
% axisY =  [578.3 1016];
% % for 5-14-FOV2-10 ms
% % axisX =   [620 430];
% % axisY =  [744 862];
% axisDir = [diff(axisX) diff(axisY)];

%points for kymograph
% for 5-14-FOV1-10 ms
x = [579.4538  572.8897  547.9462  530.8795  521.6897  511.1872  495.4333  490.1821  480.9923  470.4897  458.6744  444.2333];
y = 1e3*[0.5906    0.6221    0.6773    0.7153    0.7521    0.7941    0.8295    0.8978    0.9333    0.9490    0.9845    1.0120];
% for 5-14-FOV2-10 ms
% x = [449 465 478 491 505 522 535 546 557 568 581 595];
% y = [853 847 840 837 828 815 808 802 792 784 769 761];

xs = [];
ys = [];
dists = sqrt(diff(x).^2 + diff(y).^2);
for ix = 1:length(x)-1
    tmpx = round(linspace(x(ix),x(ix+1), ceil(dists(ix)/3)));
    tmpy = round(linspace(y(ix),y(ix+1), ceil(dists(ix)/3)));
    xs = [xs tmpx(2:end)];
    ys = [ys tmpy(2:end)];
end
% % for 5-14-FOV1-10 ms
xs = xs(25:125);
ys = ys(25:125);
% for 5-14-FOV2-10 ms
% do nothing

%get the traces for the points:
% inds = sub2ind(sys.opts.dim, round(ys), round(xs));
% Spts = 0*sys.input.S(inds,:);
% for plane = 1:sz(3)
%     Spts = Spts +  sys.input.S((plane-1)*(sz(1)*sz(2))+inds,:);
% end
% ptDFF = Spts*sys.output.F;
% B = mean(ptDFF(:,400),2);
% ptDFF = max(0,(ptDFF-B)./(sum(Spts,2)+B));

%project onto the position axis
% axisPos = nan(1, length(xs));
% for p = 1:length(xs)
%    axisPos(p) = dot([xs(p) ys(p)], axisDir)./norm(axisDir);
% end
% uncagingCenter = 620;
% axisPos = (axisPos-uncagingCenter)./5;
%%% figure 4
%plot kymograph
% figure('name', 'Kymograph', 'Color', 'w')
% imagesc([-400 1500]*(1/1.016), [min(axisPos) max(axisPos)], ptDFF(:,1:1901)) 
% ylabel('Position (um)')
% xlabel('Time (ms)')
% 
% hold on,
% plot([5 5], [min(axisPos) max(axisPos)], 'r:', 'linewidth', 2)
% set(gca, 'tickdir', 'out', 'linewidth', 1.5, 'box', 'off')
% 
% hcb = colorbar;
% ylabel(hcb, '{\Delta}F/F')

%traces
% for 5-14-FOV1-10 ms
% selected = [55 59 64 67 71 78 82 89 99];
% for 5-14-FOV2-10 ms
% selected = 5:5:50;

% cmp = flipud(parula(length(selected)+2));
% distances = axisPos(selected); distances = distances-min(distances);
% DFF = ptDFF(selected, 300:1400);
%%% figure 5
% figure, 
% for ix = 1:size(DFF,1)
%     plot(-100:1000, DFF(ix,:), 'linewidth', 2, 'color', cmp(2+ix,:))
%     hold on
% end
% xlim([-100 1000]);
% xlabel('time (ms)')
% ylabel('{\Delta}F/F_0','Interpreter', 'TeX')
% hl = legend(strcat(num2str(distances', 3), {' um'}));
% title(hl, 'Distance')
%%
rgb_ref = repmat(squeeze(sum(sys.input.ref_image,3)),1,1,3);
rgb_ref(repmat(sys.output.SLM_mask>0.1,1,1,3)) = 0;
rgb_ref(repmat(sys.output.SLM_mask<0.1,1,1,3)) = 0.00002*rgb_ref(repmat(sys.output.SLM_mask<0.1,1,1,3));
rgb_ref(rgb_im>0.0001) = 0;

%% Raster Image
%%% figure 6
f1 = figure('color','w'); ax = axes; hIM = imshow(rgb_im+rgb_ref, 'parent', ax);
% f1 = figure('color','w'); ax = axes; hIM = imshow(rgb_im, 'parent', ax);
cmaprange = linspace(mm,MM,100);
% cmap1 = (hsv(length(cmaprange)));
cmap1 =  flipud(myColorMap(length(cmaprange)));
ax1 = axes('pos', [2000 2000 1 1]);
% colormap(ax1,cmap1(100*cmin:100*cmax,:));
colormap(ax1,cmap1);
h_cb1 = colorbar(ax1, 'axislocation', 'in','fontsize',20);
set(h_cb1, 'pos', [0.83    0.10    0.01   0.83], 'tickdir', 'out', 'ticklength', 0, 'ticks', linspace(cmaprange(1),cmaprange(end),3), 'ticklabels', int2str(linspace(cmaprange(1),cmaprange(end),3)'), 'fontsize', 12)
h_cb1.Label.String = 'Saturation = Max {\Delta}F/F_0';
h_cb1.Label.FontSize = 14;
% cmap2 = (hsv2rgb([0.4*ones(length(cmaprange),1), linspace(0,1,length(cmaprange))', 0.95*ones(length(cmaprange),1)]));
% ax2 = axes('pos', [2000 2000 1 1]);
% h_cb2 = colorbar(ax2, 'axislocation', 'in','fontsize',20);
% colormap(ax2, cmap2);
% set(h_cb2, 'pos', [0.9    0.10    0.01   0.83], 'tickdir', 'out', 'ticklength', 0, 'ticks', [0 1], 'ticklabels', int2str([0 1]'), 'fontsize', 12)
% h_cb2.Label.String = 'Saturation = Max {\Delta}F/F_0';
% h_cb2.Label.FontSize = 14;
%kymograph axis
hold(ax, 'on')
hKL = plot(ax, [xs(1) xs(end)]-40, [ys(1) ys(end)]-25, 'y', 'linewidth', 2);
%uncaging locations
% uncagingx = round([791.2340  506.0005]);
% uncagingy = round([508.4363  801.8780]);
axes(ax)
hA2 = arrow([uncagingx'+22 uncagingy'], [uncagingx'+15 uncagingy'], 'color', 'y', 'length', 15, 'TipAngle', 30);

%%
keyboard
end

function c = myColorMap(n)
c = ones(n,3);
c(:,2) = linspace(0,1,n);
c(:,3) = linspace(0,1,n);
end