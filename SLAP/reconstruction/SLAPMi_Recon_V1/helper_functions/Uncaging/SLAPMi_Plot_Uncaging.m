function SLAPMi_Plot_Uncaging()
dr = '/Volumes/LaCie/InVitro_Recons/';
cd(dr);
filenames = dir('**/*RECON*.mat');
nFiles = length(filenames);
% LED1 = [286, 214]; LED2 = [214,284];
UncagingLoc1Traces = [];
UncagingLoc2Traces = [];
for filenum = 1:nFiles
load([filenames(filenum).folder filesep filenames(filenum).name]);
filenames(filenum).FOV = str2double(filenames(filenum).folder(end));

StimOnset = sys_recon.opts.StimOnset;
filenames(filenum).StimOnset = StimOnset;

if filenum == 1 || filenames(filenum).FOV - filenames(filenum-1).FOV ~= 0
tmin = StimOnset; %400
tmax = StimOnset + 150;

sz = size(sys_recon.input.ref_image);
valid3D = reshape(full(sum(sys_recon.input.S(:,1:end-1),2)>eps), sz)  & repmat(sys_recon.output.SLM_mask, [1 1 sz(3)])>0.8;
ref2D = sys_recon.input.ref_image; ref2D(~valid3D) = 0; ref2D = sum(ref2D,3);

S2D = 0*sys_recon.input.S(1:(sz(1)*sz(2)),:);
for plane = 1:sz(3)
    S2D = S2D +  sys_recon.input.S((plane-1)*(sz(1)*sz(2)) + (1:(sz(1)*sz(2))),:).*reshape(valid3D(:,:,plane), [],1);
end
valid2D = any(valid3D,3); %pixels within the mask
Svalid = S2D(valid2D(:),:);

fprintf('Calculating raw movie ... ');
drawframes = 1:size(sys_recon.output.F,2);

% 
mDFF = Svalid*sys_recon.output.F;
baselineX = nan(size(sys_recon.output.F));
for i = 1:size(sys_recon.output.F,1)
    baselineX(i,:) =  mean(sys_recon.output.F(i, 400), 2); %for 2-spot uncaging
end
baseline = single(Svalid*baselineX);

mDFF = (mDFF-baseline)./(repmat(full(sum(Svalid(:,1:end),2)),1,sys_recon.opts.T) + baseline);
mDFF = max(0, mDFF);
mDFF = mDFF(:,drawframes);
clear baseline
fprintf('done. \n');

%%
maxDFF = max(mDFF(:, tmin:tmax),[],2);
minDFF = 0*mDFF(:,tmin);
maxDFF_time = inf(size(mDFF,1),1);
for i = 1:sum(valid2D(:))
    T= find(mDFF(i,tmin:tmax)>= minDFF(i) + 0.67*(maxDFF(i)-minDFF(i)),1,'first');
    if ~isempty(T)
       maxDFF_time(i) = T;
    end
end

ysz = size(ref2D);
%%
%Prepare data
maxDFF_time = maxDFF_time-min(maxDFF_time);
tmin = 35; tmax = 100;
cmin = 0.4; cmax = 1;
H = cmin + (cmax-cmin)*max(0,min(1, (maxDFF_time-tmin)./(tmax-tmin))); %10ms->0.4 50ms->1
S = min(1, maxDFF/1);
V = Svalid*(sys_recon.output.F(:,1)+1);
ff = max(V)/4; 
V = min(1, (V./ff));

rgb_image = hsv2rgb([H, S, V]);

rgb_im = repmat(full(((S2D(:,end)*(sys_recon.output.F(end,1)+1))./ff)),1,3);
rgb_im(valid2D,:) = rgb_image;
rgb_im = reshape(rgb_im, sys_recon.opts.dim(1), sys_recon.opts.dim(2), 3);
rgb_im = rgb_im.*repmat((sys_recon.output.SLM_mask), 1, 1, 3);

%%
rgb_ref = repmat(squeeze(sum(sys_recon.input.ref_image,3)),1,1,3);
rgb_ref(repmat(sys_recon.output.SLM_mask>0.1,1,1,3)) = 0;
rgb_ref(repmat(sys_recon.output.SLM_mask<0.1,1,1,3)) = 0.00002*rgb_ref(repmat(sys_recon.output.SLM_mask<0.1,1,1,3));
rgb_ref(rgb_im>0.0001) = 0;

%% Raster Image
if filenum == 1
hF = figure('color','w');
else
figure(hF);
end
hAxIm = axes; 
hIm = imshow(rgb_im+rgb_ref, 'parent', hAxIm);
% f1 = figure('color','w'); ax = axes; hIM = imshow(rgb_im, 'parent', ax);
cmap1 = hsv(100);
ax1 = axes('pos', [2000 2000 1 1]);
% colormap(ax1,cmap1(100*cmin:100*cmax,:));
colormap(ax1,cmap1(35:95,:));
h_cb1 = colorbar(ax1, 'axislocation', 'in','fontsize',20);
set(h_cb1, 'pos', [0.83    0.10    0.01   0.83], 'tickdir', 'out', 'ticklength', 0, 'ticks', linspace(0, 1, 11), 'ticklabels', int2str(linspace(tmin,tmax,11)'), 'fontsize', 12)
h_cb1.Label.String = 'Hue = Delay (ms)';
h_cb1.Label.FontSize = 14;
cmap2 = (hsv2rgb([0.4*ones(100,1), linspace(0,1,100)', 0.95*ones(100,1)]));
ax2 = axes('pos', [2000 2000 1 1]);
h_cb2 = colorbar(ax2, 'axislocation', 'in','fontsize',20);
colormap(ax2, cmap2);
set(h_cb2, 'pos', [0.9    0.10    0.01   0.83], 'tickdir', 'out', 'ticklength', 0, 'ticks', [0 1], 'ticklabels', int2str([0 1]'), 'fontsize', 12)
h_cb2.Label.String = 'Saturation = Max {\Delta}F/F_0';
h_cb2.Label.FontSize = 14;

%%
    Segs = 0;
    for plane = 1:sys_recon.opts.Psize
        Segs = Segs+ sys_recon.input.S((1:numel(ref2D)) + (plane-1)*numel(ref2D),:);
    end
    
    disp('Click on the uncgaing centers!')
    [uncagingx,uncagingy] = ginput(2);
    ptPos1 = round([uncagingx(1) uncagingy(1)]);
    ptPos2 = round([uncagingx(2) uncagingy(2)]);
    ind1 = sub2ind(ysz,ptPos1(2), ptPos1(1));
    bestSeg1 = find(Segs(ind1,:)>0);
    ind2 = sub2ind(ysz,ptPos2(2), ptPos2(1));
    bestSeg2 = find(Segs(ind2,:)>0);
    filenames(filenum).uncagingx = uncagingx;
    filenames(filenum).uncagingy = uncagingy;
    filenames(filenum).bestSeg1 = bestSeg1;
    filenames(filenum).bestSeg2 = bestSeg2;
    
end


UncagingLoc1Traces = [UncagingLoc1Traces; Segs(ind1,bestSeg1)*sys_recon.output.F(bestSeg1,:)/sum(Segs(ind1,bestSeg1))];
UncagingLoc2Traces = [UncagingLoc2Traces; Segs(ind2,bestSeg2)*sys_recon.output.F(bestSeg2,:)/sum(Segs(ind2,bestSeg2))];

end

save([dr filesep 'UncagingTraces'],'UncagingLoc1Traces','UncagingLoc2Traces','filenames');
keyboard;
end