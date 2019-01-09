function Plot_Uncaging(sys)
StimOnset = sys.opts.StimOnset;
tmin = StimOnset; %400
tmax = StimOnset + 150;

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

% 
mDFF = Svalid*sys.output.F;
baselineX = nan(size(sys.output.F));
for i = 1:size(sys.output.F,1)
    baselineX(i,:) =  mean(sys.output.F(i, 400), 2); %for 2-spot uncaging
end
baseline = single(Svalid*baselineX);

mDFF = (mDFF-baseline)./(repmat(full(sum(Svalid(:,1:end),2)),1,sys.opts.T) + baseline);
mDFF = max(0, mDFF);
mDFF = mDFF(:,drawframes);
clear baseline
fprintf('done. \n');

%%
maxDFF = max(mDFF(:, tmin:tmax),[],2);
minDFF = 0*mDFF(:,tmin);
maxDFF_time = inf(size(mDFF,1),1);
for i = 1:sum(valid2D(:))
    T= find(mDFF(i,tmin:tmax)>= minDFF(i) + 0.25*(maxDFF(i)-minDFF(i)),1,'first');
    if ~isempty(T)
       maxDFF_time(i) = T;
    end
end

%%
%Prepare data
maxDFF_time = maxDFF_time-min(maxDFF_time);
tmin = 10; tmax = 22;
cmin = 0.8; cmax = 0.79999; %0.6666; %0.5;
H = cmin + (cmax-cmin)*max(0,min(1, (maxDFF_time-tmin)./(tmax-tmin))); %10ms->0.4 50ms->1
S = min(1, maxDFF/1);
V = Svalid*(sys.output.F(:,1)+1);
ff = max(V)/4; 
V = min(1, (V./ff));

rgb_image = hsv2rgb([H, S, V]);

rgb_im = repmat(full(((S2D(:,end)*(sys.output.F(end,1)+1))./ff)),1,3);
rgb_im(valid2D,:) = rgb_image;
rgb_im = reshape(rgb_im, sys.opts.dim(1), sys.opts.dim(2), 3);
rgb_im = rgb_im.*repmat((sys.output.SLM_mask), 1, 1, 3);
%-------%
% uncagingx = round([293  223]);
% uncagingy = round([214  284]);
% uncagingx = round([791.2340  506.0005]);
% uncagingy = round([508.4363  801.8780]);
% inds = sub2ind(sys.opts.dim, uncagingy, uncagingx);
% Spts = 0*sys.input.S(inds,:);
% for plane = 1:sz(3)
%     Spts = Spts +  sys.input.S((plane-1)*(sz(1)*sz(2))+inds,:);
% end
% ptDFF = Spts*sys.output.F;
% B = [ptDFF(1,StimOnset+10)  ; ptDFF(2,StimOnset) ];
% ptDFF = max(0,(ptDFF-B)./(sum(Spts,2)+B));
% figure, plot((-299:500)*(1/1.016), ptDFF(:, StimOnset+(-299:500))', 'linewidth', 2)
% legend(strcat({'Uncaging Location '}, int2str([1:length(uncagingx)]')));

%%%%%%%%%%%%%%%
% [dffMax,tMax] = max(ptDFF,[],2); 
% Stim2Rise1 = ptDFF(1,StimOnset+1:tMax(1));
% Stim2Rise2 = ptDFF(2,StimOnset+1:tMax(2));
% df1 = abs((Stim2Rise1-dffMax(1)/2));
% df2 = abs((Stim2Rise2-dffMax(2)/2));
% [m1,halfRiseTime1] = min(df1);
% [m2,halfRiseTime2] = min(df2);
% %%%%%%%%
% halfRiseTime1-halfRiseTime2
% indices = StimOnset+(-24+100);
% %%%% Zoomed
% figure, plot((indices-StimOnset)*(1/1.016), ptDFF(:, indices)', 'linewidth', 2);
% hold on; 
% line([halfRiseTime1, halfRiseTime1],[1 0]);
% line([halfRiseTime2, halfRiseTime2],[1 0]);
% legend(strcat({'Uncaging Location '}, int2str([1:length(uncagingx)]')));

%%
% normalizedDFF = (ptDFF-min(ptDFF,[],2))./(max(ptDFF,[],2)-min(ptDFF,[],2));
% 
% figure, plot((indices-StimOnset)*(1/1.016), normalizedDFF(:, indices)', 'linewidth', 2)
% legend(strcat({'Uncaging Location '}, int2str([1:length(uncagingx)]')));


%-------%
%kymograph

%axis of kymograph
% axisX =   [581.2963 441.9052];
% axisY =  [578.3 1016];
% axisDir = [diff(axisX) diff(axisY)];

%points for kymograph
%Looger lab
% x = [579.4538  572.8897  547.9462  530.8795  521.6897  511.1872  495.4333  490.1821  480.9923  470.4897  458.6744  444.2333];
% y = 1e3*[0.5906    0.6221    0.6773    0.7153    0.7521    0.7941    0.8295    0.8978    0.9333    0.9490    0.9845    1.0120];
% 11-21 FOV8
% x = [617.5432  607.6617  598.5814  594.0412  583.0914  579.6195   574.2782   568.9368  561.1918  552.1115  546.5031  533.9509];
% y = [948.4773  936.4593  921.5034  893.4613  878.5055  874.4995   865.6862   859.2766  851.7987  840.8489  835.5075  814.1421];

% xs = [];
% ys = [];
% dists = sqrt(diff(x).^2 + diff(y).^2);
% for ix = 1:length(x)-1
%     tmpx = round(linspace(x(ix),x(ix+1), ceil(dists(ix)/3)));
%     tmpy = round(linspace(y(ix),y(ix+1), ceil(dists(ix)/3)));
%     xs = [xs tmpx(2:end)];
%     ys = [ys tmpy(2:end)];
% end
% 
% xs = xs(25:125);
% ys = ys(25:125);
% 
% %get the traces for the points:
% inds = sub2ind(sys.opts.dim, round(ys), round(xs));
% Spts = 0*sys.input.S(inds,:);
% for plane = 1:sz(3)
%     Spts = Spts +  sys.input.S((plane-1)*(sz(1)*sz(2))+inds,:);
% end
% ptDFF = Spts*sys.output.F;
% B = mean(ptDFF(:,StimOnset),2);
% ptDFF = max(0,(ptDFF-B)./(sum(Spts,2)+B));
% 
% %project onto the position axis
% axisPos = nan(1, length(xs));
% for p = 1:length(xs)
%    axisPos(p) = dot([xs(p) ys(p)], axisDir)./norm(axisDir);
% end
% uncagingCenter = 620;
% axisPos = (axisPos-uncagingCenter)./5;

% %plot kymograph
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
% selected = [55 59 64 67 71 78 82 89 99];
% % selected = [38 39 40 41 47 48];
% cmp = flipud(parula(length(selected)+2));
% distances = axisPos(selected); distances = distances-min(distances);
% DFF = ptDFF(selected, StimOnset-100:StimOnset+500);
% figure, 
% for ix = 1:size(DFF,1)
%     plot(-100:500, DFF(ix,:), 'linewidth', 2, 'color', cmp(2+ix,:))
%     hold on
% end
% xlim([-100 1000]);
% xlabel('time (ms)')
% ylabel('{\Delta}F/F_0','Interpreter', 'TeX')
% hl = legend(strcat(num2str(distances', 3), {' um'}));
% title(hl, 'Distance')
%%
rgb_ref = repmat(squeeze(sum(sys.input.ref_image,3)),1,1,3);
rgb_ref(repmat(~(sys.output.SLM_mask<1),1,1,3)) = 0;
rgb_ref = min(1,0.00002*rgb_ref) +0.1;
rgb_ref(repmat(any(rgb_im>0,3),1,1,3)) = 0;

rgb_ref(:,:,1) = 0.5*rgb_ref(:,:,1);

%% Raster Image
hF = figure('color','w');

hAxIm = axes; 
hIm = imshow(rgb_im+rgb_ref, 'parent', hAxIm);
% f1 = figure('color','w'); ax = axes; hIM = imshow(rgb_im, 'parent', ax);
cmap1 = hsv(100); cmap1 = flipud(cmap1(67:end,:));
ax1 = axes('pos', [2000 2000 1 1]);
% colormap(ax1,cmap1(100*cmin:100*cmax,:));
colormap(ax1,cmap1);
h_cb1 = colorbar(ax1, 'axislocation', 'in','fontsize',20);
set(h_cb1, 'pos', [0.83    0.10    0.01   0.83], 'tickdir', 'out', 'ticklength', 0, 'ticks', linspace(0, 1, 11), 'ticklabels', int2str(linspace(tmin,tmax,11)'), 'fontsize', 12)
h_cb1.Label.String = 'Hue = Delay (ms)';
h_cb1.Label.FontSize = 14;
cmap2 = (hsv2rgb([0.8*ones(100,1), linspace(0,1,100)', 0.95*ones(100,1)]));
ax2 = axes('pos', [2000 2000 1 1]);
h_cb2 = colorbar(ax2, 'axislocation', 'in','fontsize',20);
colormap(ax2, cmap2);
set(h_cb2, 'pos', [0.9    0.10    0.01   0.83], 'tickdir', 'out', 'ticklength', 0, 'ticks', [0 1], 'ticklabels', int2str([0 1]'), 'fontsize', 12)
h_cb2.Label.String = 'Saturation = Max {\Delta}F/F_0';
h_cb2.Label.FontSize = 14;
%kymograph axis
% hold(hAxIm, 'on')
% hKL = plot(hAxIm, [xs(1) xs(end)]-40, [ys(1) ys(end)]-25, 'y', 'linewidth', 2);
%uncaging locations
% uncagingx = round([791.2340  506.0005]); Looger Lab
% uncagingy = round([508.4363  801.8780]);

uncagingx = round([822  524]); % New dataset 11-21-FOV8
uncagingy = round([532  853]);

axes(hAxIm)
hA2 = arrow([uncagingx'+12 uncagingy'], [uncagingx'+5 uncagingy'], 'color', 'r', 'length', 15, 'TipAngle', 30);

% hF2 = figure;
% hAxActivity= axes('Pos', [0.1 0.1 0.8 0.8], 'tickdir', 'out', 'linewidth', 1.5);
% set (hF, 'WindowButtonDownFcn', @bdf);
% ysz = size(ref2D);
% Segs = 0;
% for plane = 1:sys.opts.Psize
% Segs = Segs+ sys.input.S((1:numel(ref2D)) + (plane-1)*numel(ref2D),:);
% end
%     function bdf(src, event)
%         cp = get(hAxIm,'CurrentPoint');
%         cp = cp(1, 1:2);
%         %If the point is within the axes
%         if all(cp>=1 & cp<=ysz)
%             ptPos = round([cp(1) cp(2)]);
%             ind = sub2ind(ysz,ptPos(2), ptPos(1));
%             bestSeg = find(Segs(ind,:)>0);
%             if ~isempty(bestSeg)
% 
%                 set(hF,'name', ['Segment:' int2str(bestSeg) ' Position: ' num2str(cp)]),
%                 set(hF2,'name', ['Segment:' int2str(bestSeg) ' Position: ' num2str(cp)]),
%                 
%                 %Activity plot
%                 cla(hAxActivity);
%                 plot(hAxActivity, sys.output.F(bestSeg,:)');
%             end
%         end
%     end


end