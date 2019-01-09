function Red_Green_Raster(sys)
%%
G0 = sys.input.GreenIm;
R0 = sys.input.RedIm;
%% Reference Image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SLM = sys.output.SLM_mask > 0.8;
Segmented = sum(reshape(full(sum(sys.input.S(:,1:end-1),2)),size(sys.input.ref_image)),3)>0;
R = R0.*Segmented.*SLM;
G = G0.*Segmented.*SLM;

%%

drawPixels = R>prctile(R(:),99.9) | G>prctile(G(:),99.9);
StimOnset = 500;
preStimTime = 50; % 50 ms
PostStimTime = 1000; %1 second
drawFrames = StimOnset-preStimTime:StimOnset+PostStimTime-1;

%% Reconstructed Movie
% Calculate Raw Movie
[mRaw,mDFF,valid2D] = calc_Movies(sys,drawFrames,drawPixels);
%%
nPixels = 1000;
% dff_threhsold = 0.01;
% dffMax = min(5, prctile(mDFF(mDFF>dff_threhsold), 99));
dffMax = 1;
valid = valid2D ;
Rvalid = R(valid);
Gvalid = G(valid);
%%
% Option 1
RGvalid = (Rvalid)./(Rvalid+Gvalid);
[s, sortorder] = sort(RGvalid);
% s = [s(1:nPixels);s(end-nPixels+1:end)];
% sortorder = [sortorder(1:nPixels);sortorder(end-nPixels+1:end)];
sorted_mDFF = mDFF(sortorder,:);
sorted_mRaw = sqrt(mRaw(sortorder,:));
%% Optional: sort according to brightness
% [s0,order0] = sort(sum(sorted_mRaw,2),'descend');
% s = s(order0);
% sorted_mDFF = sorted_mDFF(order0,:);
% sorted_mRaw = sorted_mRaw(order0,:);

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


