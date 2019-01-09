function Red_Green_plot(sys)
%%
% close all;
figure;
imagesc(sys.output.SLM_mask); axis square; title('SLM mask'); drawnow;
G0 = double(max(sys.input.GreenIm-200,0));
R0 = double(sys.input.RedIm);
%% Reference Image
figure; imshow(800*uint16(sqrt((cat(3,R0, G0, 0*G0)))) + 0*uint16(sqrt((cat(3,R0.*sys.output.SLM_mask, G0.*sys.output.SLM_mask, 0*G0)))) );
drawnow;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SLM = sys.output.SLM_mask > 0.8;
Segmented = sum(reshape(full(sum(sys.input.S(:,1:end-1),2)),size(sys.input.ref_image)),3)>0;
figure; imshow(800*uint16(sqrt((cat(3,R0.*Segmented, G0.*Segmented, 0*G0)))) + 0*uint16(sqrt((cat(3,R0.*sys.output.SLM_mask, G0.*sys.output.SLM_mask, 0*G0)))) );


R = R0.*Segmented.*SLM;
G = G0.*Segmented.*SLM;
MostRed = R >2.1*G;
RG = G./R;
Green = medfilt2(RG > prctile(RG(:),95));
Yellow = imdilate(G > prctile(G(:),99),strel('disk',10));
Red = medfilt2(MostRed & ~Yellow);
%%

drawPixels = Red | Green;
StimOnset = 500;
preStimTime = 50; % 50 ms
PostStimTime = 1000; %1 second
drawFrames = StimOnset-preStimTime:StimOnset+PostStimTime-1;

%% Reconstructed Movie
% Calculate Raw Movie
[mRaw,mDFF,valid2D] = calc_Movies(sys,drawFrames,drawPixels);
% [mRaw,mDFF,valid2D] = calc_Movies(sys,drawFrames,drawPixels,preStimTime);
%%
valid2DG = valid2D & Green;
valid2DR = valid2D & Red;

figure;
subplot(1,2,1);
imshow(100*uint16(sqrt((cat(3,R, G, 0*G))))+ 600*uint16(sqrt((cat(3,0*R, G.*valid2DG, 0*G)))) );
title('Selected Green Pixels');
set(gca,'fontsize',16);
subplot(1,2,2);
imshow(100*uint16(sqrt((cat(3,R, G, 0*G))))+ 600*uint16(sqrt((cat(3,R.*valid2DR, 0*G, 0*G)))) );
title('Selected Red Pixels')
set(gcf,'units','normalized','outerposition',[0 0 1 1],'defaulttextinterpreter','latex')
set(gca,'fontsize',16);
drawnow
%%
nPixels = 1000;
dff_threhsold = 0.01;
dffMax = min(5, prctile(mDFF(mDFF>dff_threhsold), 99));
valid = valid2D ;
Rvalid = Red(valid);
Gvalid = Green(valid);
%%
% Option 1
RGvalid = (Rvalid)./(Rvalid+Gvalid);
[s, sortorder] = sort(RGvalid);
s = [s(1:nPixels);s(end-nPixels+1:end)];
sortorder = [sortorder(1:nPixels);sortorder(end-nPixels+1:end)];
sorted_mDFF = mDFF(sortorder,:);
sorted_mRaw = sqrt(mRaw(sortorder,:));

PlotmeanDFF(sorted_mDFF,nPixels);
PlotEvokedDFF(mDFF,mRaw,sorted_mDFF,nPixels,preStimTime,valid2D,Red,R0,G0)
%% Optional: sort according to brightness
[s0,order0] = sort(sum(sorted_mRaw,2),'descend');
s = s(order0);
sorted_mDFF = sorted_mDFF(order0,:);
sorted_mRaw = sorted_mRaw(order0,:);

% Gamma correction
RGB_hue = [s 1-s 0*s];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
factor = sqrt(prctile(mRaw(:), 99));
hue = rgb2hsv(RGB_hue);
H = repmat(hue(:,1), 1, size(sorted_mDFF,2));
S = min(1, sorted_mDFF./dffMax);
V = min(1,sorted_mRaw./factor);
V = sqrt(V);

RGB = hsv2rgb(cat(3, H, S, V));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
f = figure;
ax1 = axes(f);
imshow(RGB, 'parent', ax1);
set(ax1,'ytick',[]); set(gca,'fontsize',20);
set(ax1,  'pos', [0.12    0.1    0.65    0.8]);
axis(ax1, 'normal')

hold on

cmap2 = hsv(300); cmap2 = cmap2(1:100,:); %[linspace(1,0, 100)' linspace(0,1,100)' linspace(0,0,100)'];
ax2 = axes('pos', [2000 2000 1 1]);
h_cb2 = colorbar(ax2, 'axislocation', 'in','fontsize',12);
colormap(ax2, cmap2);
set(h_cb2, 'pos', [0.80    0.1    0.01    0.8], 'tickdir', 'out', 'ticklength', 0, 'ticks', [0 1], 'ticklabels', {'', ''})
h_cb2.Label.String = 'Hue = Ground Truth Pixel Color Ch1/(Ch1+Ch2)';
h_cb2.Label.Position = [1 0.5000 0];

cmap3 = (hsv2rgb([0.18*ones(100,1), linspace(0,1,100)', ones(100,1)]));
ax3 = axes('pos', [2000 2000 1 1]);
h_cb3 = colorbar(ax3, 'axislocation', 'in','fontsize',12);
colormap(ax3, cmap3);
set(h_cb3, 'pos', [0.85    0.1    0.01    0.8], 'tickdir', 'out', 'ticklength', 0, 'ticks', [0 1], 'ticklabels', num2str([0 prctile(sorted_mDFF(:), 99.9)]', 2))
h_cb3.Label.String = 'Saturation = {\Delta}F/F';
h_cb3.Label.Position = [1 0.5000 0];

cmap4 = (hsv2rgb([ones(100,1), zeros(100,1), linspace(0,1,100)']));
ax4 = axes('pos', [2000 2000 1 1]);
h_cb4 = colorbar(ax4, 'axislocation', 'in','fontsize',12);
colormap(ax4, cmap4);
set(h_cb4, 'pos', [0.9    0.1    0.01    0.8], 'tickdir', 'out', 'ticklength', 0, 'ticks', [0 1], 'ticklabels', int2str([0 1]'))
h_cb4.Label.String = 'Value = \surd Pixel intensity (a.u.)';
h_cb4.Label.Position = [1 0.5000 0];

drawnow
set(gca,'fontsize',16);
set(gcf,'units','normalized','outerposition',[0 0 1 1],'defaulttextinterpreter','latex')
%%



end

function PlotmeanDFF(sorted_mDFF,nPixels)
%%
fsize = 16;
figure;
ax1 = subplot(2,1,1); plot(mean(sorted_mDFF(1:nPixels,:)),'g');
hold on;
plot(mean(sorted_mDFF(end-nPixels+1:end,:)),'r');
xlabel('Time (ms)')
ylabel('\Delta F / F');
title('Mean \Delta F / F for Green Pixels')
set(gca,'fontsize',fsize)
axis tight
ax2 = subplot(2,1,2); plot(mean(sorted_mDFF(end-nPixels+1:end,:)),'r');
xlabel('Time (ms)')
ylabel('\Delta F / F');
title('Mean \Delta F / F for Red Pixels')
set(gca,'fontsize',fsize)
axis tight
linkaxes([ax1 ax2],'x');
drawnow
set(gcf,'units','normalized','outerposition',[0 0 1 1],'defaulttextinterpreter','latex')
set(gca,'fontsize',16);
keyboard
end
function PlotEvokedDFF(mDFF,mRaw,sorted_mDFF,nPixels,preStimTime,valid2D,Red,R0,G0)
%%
fsize = 16;
figure;
subplot(1,3,1);
eval_length = 100; % ms after the stimulus
Red_post = mean(sorted_mDFF(end-nPixels+1:end, preStimTime+(1:eval_length)),2);
Green_post= mean(sorted_mDFF(1:nPixels,preStimTime+(1:eval_length)),2);
boxplot([Green_post, Red_post],'Labels',{'Green Pixels ','Red'},'Whisker',2);
set(gca,'fontsize',fsize);axis square
set(gca,'fontsize',16);
drawnow
%%
[mdffMax, mDFFTime] = max(mean(sorted_mDFF(1:nPixels,:)));
subplot(1,3,2) 
EvokedDFF = nan(size(valid2D));
EvokedDFF(valid2D) = mDFF(:,mDFFTime);
imshow(800*uint16(sqrt((cat(3,R0.*EvokedDFF, G0.*EvokedDFF, 0*G0)))));
set(gca,'fontsize',16);
drawnow
%%
mRawRed = nan(size(valid2D));
mRawRed(valid2D) = mean(mRaw,2);
subplot(1,3,3);
imagesc(mRawRed.*Red); axis square
set(gcf,'units','normalized','outerposition',[0 0 1 1],'defaulttextinterpreter','latex')
set(gca,'fontsize',16);
end
