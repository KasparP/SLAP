function hFig = scandataPlot2D(scandata, optsin)
%default options
opts.doPCA  = true;         tooltips.doPCA = 'show a plot with NMF';
opts.tau = 20;              tooltips.tau = 'Time constant for DFF filter, in frames';
opts.dffMax = 2;
opts.plotChanIx = 1;

if nargin>1  && ~isempty(optsin) %UPDATE WITH USER-SPECIFIED OPTIONS
    for field = fieldnames(optsin)'
         opts.(field{1}) = optsin.(field{1});
    end
else
   opts = optionsGUI(opts, tooltips);
end

F = scandata.frames;
Nchans = size(F(1).pmtData,2);
P = reshape([F.pmtData], size(F(1).pmtData,1), size(F(1).pmtData,2), []);
if isfield(scandata, 'dt')
    frametimes = scandata.dt*(0:size(P,3)-1);
else
    frametimes = (0:size(P,3)-1);
end
hFig = [];

hFig(end+1) = figure('name', scandata.metadata.galvoDataFileName, 'pos', [-1049 1  1050 1604]); %[1 1201 1600 1124]);  %#ok<AGROW>
[F, dFF] = fastDFF(P(:,opts.plotChanIx:Nchans:end), opts.tau);
subplot(1,2,1), imagesc(F); xlabel('time (frames)'), set(gca, 'ytick', []);
subplot(1,2,2), imagesc(dFF); xlabel('time (frames)'), set(gca, 'ytick', [], 'clim', [0 opts.dffMax]); colorbar;

hFig(end+1) = figure('name', scandata.metadata.galvoDataFileName, 'pos', [9 49 784 1068]);  %#ok<AGROW>
if isfield(scandata.opts, 'diodeCh') && scandata.opts.diodeCh
        %show a denoised average trace
        denoisethresh = 1; %remove Fourier components X log units of power above average background
        s = squeeze(nansum(P,1));
        S = fft(s);
        Spow = log(abs(S));
        Spow = Spow-smooth(Spow,50);
        peaks = Spow>denoisethresh; peaks([1:ceil(0.02*end)+1 floor(0.98*end):end]) = false;
        peaks = find(peaks);
        S(peaks) = 0; S(peaks+1) = 0; S(peaks-1) = 0;
        s2 = ifft(S);
        s2 = s2-smooth(s2, 20*opts.tau) + mean(s2);
        
        subplot(2,1,1); 
        plot(1000*frametimes, squeeze(nansum(P,1)), 'linestyle', ':'); 
        hold on
        plot(1000*frametimes, s2); 
        
        legend({'original', 'denoised'});
        xlabel('time (ms)'), ylabel('Brightness (photons/frame)'); ylim([-Inf; Inf]);
else
    subplot(2,1,1); plot(1000*frametimes, squeeze(nansum(P,1))); xlabel('time (ms)'), ylabel('Brightness (photons/frame)'); ylim([0; Inf]);
end
subplot(2,1,2); plot(squeeze(nanmean(P(:,:,1:min(100,end)),3))); xlabel('lxl'), ylabel('brightness'); ylim([0; Inf]);

if opts.doPCA
    hFig = [hFig SLAPMi_PCA(scandata)]; %#ok<AGROW>
end

end