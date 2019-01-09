function [scandata, hFig] = SLAPMi_reduce(optsin)
import slapmi

%default options
smoothWindow = 33; %smoothing window for galvo positions; set higher if you have more measurement noise
EF_thresh = 2e-2; %absolute pointing error threshold for E and F mirrors, 2e-2 is equivalent to ~300nm in sample space
opts.align = false;  %align to remove motion, good for in vivo data
opts.ignoreSaved = false; %if data has been reduced before, reduce it again anyways
opts.amp_setup = 'MPPC';
opts.doRegression = false;
opts.doNND = false; %nonnegative deconvolution of PMT traces
opts.doPlot = true; %show plots
opts.doPCA  = true; %show a plot with NMF
opts.path = []; % path to file
if nargin %UPDATE WITH USER-SPECIFIED OPTIONS
    for field = fieldnames(optsin)'
         opts.(field{1}) = optsin.(field{1});
    end
end

hFig = [];
basedir = 'E:\SLAPmiData\';

switch opts.amp_setup
    case 'MPPC'
        invert = false; %Whether the amplifier was inverted
        G = 2.3; %time shift parameter; how many samples before the peak slope does the photon response start?
        
        %W: parameters for reducing laser pulses to intensities
        W.tau = 5; W.gamma = exp(-1/W.tau);
        W.photonAmp = 0.0575;
        W.template = []; %[0.0793 ; 0.2799 ; 0.5619 ; 0.6290 ; 0.4122 ; 0.1721 ; 0.0675]; %template of photon response; set below in amp_setup; if empty, it will be estimated
        W.windowSamps = 2:25; %which samples (after onset of pulse) to average as the pixel response
        W.windowSize = max(W.windowSamps);
        W.baselineSamps = 5;
        W.doNND = opts.doNND;
        W.doRegression = opts.doRegression;
    otherwise
        error('Amplifier Setup Unspecified');
end

%SELECT INPUT FILES
if isempty(opts.path)
    [fns, dr] = uigetfile([basedir filesep '*.gdat'], 'multiselect', 'on');
    if isempty(fns)
        return
    end
    if ~iscell(fns)
        fns = {fns};
    end
else
    [dr,name,ext] = fileparts(opts.path);
    fns{1} = [name ext];
end

%PROCESS FILES
for filenum = 1:length(fns)
    scandata = [];
    
    %IO
    name = fns{filenum};
    [pathstr,name,ext] = fileparts([dr filesep name]);
    if ~strcmpi(ext, '.gdat')
        name = [name ext]; %#ok<AGROW>
    end
    if exist([pathstr filesep name '_REDUCED.mat'], 'file') && ~opts.ignoreSaved
        disp([name ': This file has been previously reduced; Loading Results'])
        load([pathstr filesep name '_REDUCED.mat']);
    else
        
        disp('Reading data...')
        tic
        channels= 2 %#ok<NOPRT>
        [metaData, pmtData, galvoData, saveData, Zdata] = SLAPMi.openDataFiles([pathstr filesep name], channels);
        saveData.metaData = metaData; %fold the metadata written at acquisition time into the other saved data
        
        %tmp to reduce disk usage
        if isfield(saveData.calib.SLM, 'map')
            saveData.calib.SLM = rmfield(saveData.calib.SLM, 'map');
        end
        
        if invert %Amplifier was inverting
            pmtData.data = -pmtData.data;
        end
        
        calib = saveData.calib;
        Nsamps = size(pmtData.data,1);
        Nchans = size(pmtData.data,2);
        toc
        
        disp('Finding pixel times...')
        tic
        
        %Identify Laser Frequency
         %fActual = 4.996295164679796e+06;
         samplesPerCycle = 50; % 50.037075825167360;
%to estimate frequency /w fft:
%         dS = Nsamps;
%         fs = 1/pmtData.dt;  %sampling frequency
%         n2 = 2^nextpow2(dS);
%         df = fs / n2; %frequency resolution of fft
%         F = fft(pmtData.data(1:min(dS, end), 1), n2);
%         fSearch = 4.999e6; %expected laser rate
%         fC = round(fSearch/df); fSearch = fC*df;
%         delta = round(100000/df); %difference in Hz to search around the center frequency
%         [~,maxind] = max(abs(F(fC+1 + (-delta:delta)))); %+1 because the first bin of F is 0
%         fActual = fSearch + df*(maxind-(delta+1));
%         samplesPerCycle = fs/fActual;      
        
        %identify phase in first few and last few samples
        chunkSize = 2e3;
        laserframes = (G+1:samplesPerCycle:Nsamps-samplesPerCycle)';
        B = linspace(0, samplesPerCycle, ceil(samplesPerCycle)+1);
        dB = B(2)-B(1);
        B = repmat(B, chunkSize,1);
        A1 =  repmat(laserframes(1:chunkSize), 1, size(B,2));
        A2 =  repmat(laserframes(end-chunkSize+1:end), 1, size(B,2));
        
        IPLNT = griddedInterpolant(pmtData.data(:,1));
        C1 = reshape(IPLNT(B(:)+A1(:)), size(B));
        Cavg1 = nth_element(C1, round(0.98*size(C1,1)));
        Cavg1 = Cavg1(round(0.98*size(C1,1)), :);
        C2 = reshape(IPLNT(B(:)+A2(:)), size(B));
        Cavg2 = nth_element(C2, round(0.98*size(C2,1)));
        Cavg2 = Cavg2(round(0.98*size(C2,1)), :);
        
        %get the centroid of the peak slope
        P1 = laserframes(1)  + centroid(diff(Cavg1))*dB; %timing of first pulse
        P2 = laserframes(end) + centroid(diff(Cavg2))*dB; %timing of last pulse
        N = length(laserframes); %number of pulses
        
        %the laser clock wobbles, so we check for phase precession
        precs = -ceil(1.2*N/1e5):ceil(1.2*N/1e5);
        slope = nan(1, length(precs));
        
        for prec_ix = 1:length(precs)
            laserframes = linspace(P1, P2, N+precs(prec_ix));
            laserframes = laserframes(1:floor(end/1e4):end);
            slope(prec_ix)= sum(IPLNT(laserframes+1)-IPLNT(laserframes)); %
        end
        [maxval, bestix] = max(slope);
        bestP = precs(bestix);
        clear IPLNT;
        
        %the timing of every laser pulse
        laserframes = linspace(P1-G, P2-G, N+bestP); 
        if laserframes(1)<W.baselineSamps+1
            laserframes = laserframes(2:end);
        end
        N = length(laserframes);
        fActual = (N-1) / (pmtData.dt * diff(laserframes([1 end])));
        toc
%         %Uncomment below to visually inspect if the laser pulses are being detected right
%           window1 = 200e6; width = 1e5; select = laserframes>=window1 & laserframes<=window1+width;
%           figure, plot(window1:window1+width, pmtData.data(window1:window1+width,1)), hold on, scatter(laserframes(select), mean(pmtData.data(window1:window1+width,1))*ones(1, sum(select)), 'r');


        disp('Reducing pixeldata...'); tic;
        
        
        
        chunkSize = 1e6; %#laser pulses per chunk (not pmt samples)
        nChunks = ceil(length(laserframes)/chunkSize);
        pixeldata = nan(length(laserframes),Nchans);
        BG2 = nth_element(pmtData.data(floor(laserframes(1:min(end,chunkSize))-1),:), ceil(0.01*min(length(laserframes),chunkSize))); BG2 = BG2(ceil(0.01*min(length(laserframes),chunkSize))) + 1e-3; %Ch2 background
        
        %pixeldata = cell(nChunks,1);
        for chunk = 1:nChunks
            chunkStart = (chunk-1)*chunkSize +1;
            chunkEnd = min(N, chunk*chunkSize);
            
            dataStart = max(1, floor(laserframes(chunkStart))-W.baselineSamps);
            dataEnd = min(length(pmtData.data), ceil(laserframes(chunkEnd)+20));
            dataChunk = pmtData.data(dataStart:dataEnd,:);
            chunkFrames = laserframes(chunkStart:chunkEnd) - dataStart + 1; %the laserframes indexed into the chunk
            
            pixeldata(chunkStart:chunkEnd,1) = getPixelData(dataChunk - BG2, chunkFrames, W);
            %pixeldata{chunk} = getPixelData(D_deconv, chunkFrames, W, template);
%             if Nchans>1
%                 keyboard
%                 BG1 = min(medfilt2(pmtData.data(1:chunkSize,2), [5,1]));
%                 D_deconv = fastDeconv(pmtData.data(:,2)-BG1, gamma);
%                 pixeldata(:,2) = getPixelData(D_deconv, laserframes, W);
%             end
        end
        clear D_deconv dataChunk;
        toc
        tic
        %generate Galvo positions at the pixel data locations
        galvodata =nan(length(laserframes), 6);
        dtPMT = pmtData.dt;
        G_INT = griddedInterpolant(galvoData.t, 1:length(galvoData.t));
        galvoframes = G_INT(laserframes*dtPMT + calib.galvos.AIdelay);
        for Gn = 1:6
            tmp = galvoData.data(:,Gn);
            tmp = smoothfast1(smoothfast1(tmp, smoothWindow), smoothWindow-5);
            tGI = griddedInterpolant(tmp);
            galvodata(:,Gn) = calib.galvos.AI2AO{Gn}(tGI(galvoframes));
        end
        clear tGI tmp galvoData G_INT pmtData;
        toc
        
        disp('Generating Frame Data...'); tic;
        %subdivide data into frames
        framePeriod = saveData.framePeriod;
        sampsPerFrame = framePeriod/dtPMT; %The galvos are retriggerred according to the sample rate!
        nFrames = floor((Nsamps-max(calib.galvos.AIdelay)/dtPMT)/sampsPerFrame);
        pxPerFrame = 4*round(saveData.framePeriod*fActual/4);
        
        %Desired Galvo Signal
        apCalib = calib; %copy of the calibration with the correct aperture information
        if saveData.aperture
            apCalib.galvos.linelength = saveData.PSF.aperture.linelength * saveData.aperture/0.4;
            for line = 1:4
                lname = (['line' int2str(line)]);
                apCalib.galvos.offset.(lname).X = saveData.PSF.aperture.(lname).offsetX;
                apCalib.galvos.offset.(lname).Y = saveData.PSF.aperture.(lname).offsetY;
            end
            galvoDesired = galvo_trace_SLAPmi(pxPerFrame,saveData.tiling,apCalib, true);
        else
            apCalib.galvos.linelength = saveData.PSF.aperture.linelength;
            galvoDesired = galvo_trace_SLAPmi(pxPerFrame,saveData.tiling,apCalib, true);
        end
        
        %figure out the cycle frequency of the galvo signals
        GD = galvodata(~isnan(galvodata(:,1)),1);
        fs = 1;
        freq = 1/(saveData.framePeriod*fActual);
        [pxx, f] = periodogram(GD,[],freq*(1-5e-5:2e-6:1+5e-5), fs);
        [~,maxind] = max(pxx);
        cycleframes = fs/f(maxind);
        
        %the desired X,Y
        line = saveData.lineIDs(ceil((1:pxPerFrame)*length(saveData.lineIDs)/pxPerFrame));
        angle = [-pi/8 -3*pi/8 -pi/8 -3*pi/8];
        desT = nan(length(line)/4, 4); desTd = desT;
        desTedges = nan(1+ length(line)/4, 4);
        VxDes = nan(pxPerFrame,1); VyDes = nan(pxPerFrame,1);
        for lineID = 1:4
            if lineID<=2 %lines 1,2
                desX = galvoDesired(line==lineID, 1);
                desY = galvoDesired(line==lineID, 2);
                
                [desT(:,lineID), desTd(:,lineID)] = V2T(desX, desY, apCalib.galvos.offset.(['line' int2str(lineID)]).X,apCalib.galvos.offset.(['line' int2str(lineID)]).Y,angle(lineID));
                VxDes(line==lineID) = desX; VyDes(line==lineID) = desY;
            elseif lineID<=4 %lines 3,4
                desX = galvoDesired(line==lineID, 3);
                desY = galvoDesired(line==lineID, 4);
                
                [desT(:,lineID), desTd(:,lineID)] = V2T(desX, desY, apCalib.galvos.offset.(['line' int2str(lineID)]).X,apCalib.galvos.offset.(['line' int2str(lineID)]).Y,angle(lineID));
                VxDes(line==lineID) = desX; VyDes(line==lineID) = desY;
            end
            
            dT = diff(desT(1:2,lineID));
            desTedges(:,lineID) =  [desT(1,lineID)-dT/2 ;  desT(:,lineID)+dT/2];
        end
        
        %process Zdata
        Zcenters = zeros(nFrames,1); %default
        nVolumes = nFrames; %default
        if ~isempty(Zdata)
            %tmp
            if ~isfield(calib.piezo, 'AIdelay')
                calib.piezo.AIdelay = 1e-3;
            end
            
            pnAI = length(Zdata.h)-1; %piezo samples per Volume
            saveData.nVolumes = nFrames/saveData.nPlanes; %set this correctly, the actual acquired frames may differ from planned
            nVolumes = ceil(saveData.nVolumes);
            
            Z = smoothfast1(smoothfast1(smoothfast1(Zdata.Zdata,9),9),9);
            Zframes = nan(pnAI,nVolumes); Z = Z(1:min(end, numel(Zframes)));
            Zframes(1:length(Z)) = calib.piezo.AI2AO(Z);
            
            %center Z for each plane
            sampsPerPlane = (pnAI+1)/saveData.nPlanes;
            Zcenters = interp1(Zframes, sampsPerPlane*(0:saveData.nPlanes-1) + (sampsPerPlane)/2 + Zdata.AIrate*calib.piezo.AIdelay, 'linear', 'extrap');
            
            %center Z for each line
            Zline = cell(1,4);
            ltimes = [1 3 2 4];
            for l = 1:4
                Zline{l} = interp1(Zframes, sampsPerPlane*(0:saveData.nPlanes-1) + (ltimes(l)-0.5)*(sampsPerPlane/4) + Zdata.AIrate*calib.piezo.AIdelay, 'linear', 'extrap');
            end
            
            %Z position of every lxl
            pxPerVol = saveData.nPlanes*pxPerFrame;
            lsZsamps = linspace(0,1, pnAI+1);
            lsPsamps = linspace(0,1, pxPerVol+1);  %(Zframes(end,:) + Zframes(1, [2:end end]))/2 ;
            Zlxls = interp1(lsZsamps, [Zframes ;  Zframes(1,[2:end end])], lsPsamps(1:end-1) + calib.piezo.AIdelay/(pnAI/Zdata.AIrate), 'linear','extrap');
        end
        
        %at this point, the galvo data and pmt data are co-referenced
        %The scope timing is calibrated such that the frame should start at sample 1
        frames = struct; frames(saveData.nPlanes, nVolumes).Zcenter = []; %extize
        frameStart = 0;
        for frame = 1:nFrames
            if frameStart+pxPerFrame>size(pixeldata,1)
                disp('There was a mismatch in the frame size, this is likely due to the requirement that pxPerFrame is a multiple of 4')
                break
            end
            
            %get the pmtData for this frame
            F_pmtData = pixeldata(frameStart+1:frameStart+pxPerFrame,:);
            
            frames(frame).galvoData = galvodata((frameStart+1:frameStart+pxPerFrame),:);
            
            %frame data is valid only if the E and F galvos are very close to their
            %target positions
            E_low = abs(frames(frame).galvoData(:,5)-apCalib.galvos.offset.E(1))<EF_thresh;
            E_high = abs(frames(frame).galvoData(:,5)-apCalib.galvos.offset.E(2))<EF_thresh;
            F_low = abs(frames(frame).galvoData(:,6)-apCalib.galvos.offset.F(1))<EF_thresh;
            F_high = abs(frames(frame).galvoData(:,6)-apCalib.galvos.offset.F(2))<EF_thresh;
            
            EOM_switch = line~=[line(end) ; line(1:end-1)];
            if frame==1
                EOM_switch(1:round(1.75e-6*fActual)) = true; %the delay to turn on the pockels cell on the first frame
            end
            EOM_switch = conv(double(EOM_switch), [zeros(1, ceil(3e-6*fActual)-1) ones(1, ceil(3e-6*fActual))], 'same') >0;
            
            %which frames are valid for reconstruction
            frames(frame).valid = (line==1 & E_low) | (line==2 & E_high) | (line==3 & F_low) | (line==4 & F_high);
            frames(frame).valid = frames(frame).valid & ~EOM_switch;
            frames(frame).pmtData = nan(pxPerFrame, Nchans); %extize
            
            %interpolate the data onto the predetermined sample space
            for lineID = 1:4
                linSelect = line==lineID & frames(frame).valid;
                if lineID<=2
                    Vx = frames(frame).galvoData(linSelect, 1);
                    Vy = frames(frame).galvoData(linSelect, 2);
                elseif lineID<=4
                    Vx = frames(frame).galvoData(linSelect, 3);
                    Vy = frames(frame).galvoData(linSelect, 4);
                end
                [T,~] = V2T(Vx, Vy, apCalib.galvos.offset.(['line' int2str(lineID)]).X,apCalib.galvos.offset.(['line' int2str(lineID)]).Y,angle(lineID));
                
                %assemble pmtData
                [sortedT,sortorder] = sort(T);
                sortedP = F_pmtData(linSelect,1); sortedP = sortedP(sortorder);
                CSP = cumsum(sortedP);
                CSP_ix = interp1qr(sortedT,(1:length(sortedT))', desTedges(:,lineID));
                FP = diff(qinterp1(1:length(CSP), CSP, CSP_ix))./diff(CSP_ix);
                frames(frame).pmtData(line==lineID,1) = FP;
                
                %frames(frame).pmtData(line==lineID,1) = interp1(T, F_pmtData(linSelect,1), desT(:,lineID));
                
                if Nchans>1
                    frames(frame).pmtData(:,2) = interp1(T, F_pmtData(linSelect,2), desT(:,lineID));
                end
            end
            
            frames(frame).Zcenter = Zcenters(frame); %Zposition
            if ~isempty(Zdata)
                %frames(frame).Zline = [Zline{1}(frame) Zline{2}(frame) Zline{3}(frame) Zline{4}(frame)];
                frames(frame).Z = Zlxls(frameStart+1:frameStart+pxPerFrame)';
            end
            %move to the next frame
            %frameStart = find(laserframes(frameStart+pxPerFrame+ [-3:3]) <= frame*(sampsPerFrame+3), 1, 'last') + frameStart + pxPerFrame - 4; %-4 compensates for the [-3:3] range of find, to the left
            frameStart = round(frame*cycleframes);
        end
        toc; 
        
        %create a scandata structure; this will be passed to linePSF and image reconstruction methods
        scandata.Vx = VxDes;
        scandata.Vy = VyDes;
        if ~isfield(calib.galvos, 'pixelsizeVperUM'); calib.galvos.pixelsizeVperUM = 1/45.3; end
        scandata.pixelSizeUM = abs(dT)./calib.galvos.pixelsizeVperUM;
        scandata.line = line;
        scandata.frames = frames;
        scandata.metadata = saveData;
        scandata.opts = opts;
        scandata.dt = cycleframes/fActual;
        scandata.laserFreq = fActual;
        if opts.align %align the data being displayed/saved; good to do for in vivo data
            scandata = SLAPMi_Motion_noREF(scandata);
        end
    end
    
    if opts.align && ~scandata.opts.align %align the data being displayed/saved; good to do for in vivo data
        scandata = SLAPMi_Motion_noREF(scandata);
    end
    
    %save the reduced data
    if ~exist([pathstr filesep name '_REDUCED.mat'], 'file') || opts.ignoreSaved
        disp('Saving...')
        save([pathstr filesep name '_REDUCED.mat'], 'scandata', '-v6') %v6 gives faster save but has max file size
    end
    
    %plot
    if opts.doPlot
        disp('Plotting...')
        if scandata.metadata.nPlanes>1
            hFig = scandataPlot3D(scandata);
        else
            hFig = scandataPlot2D(scandata);
        end
        drawnow;
    end
       
    disp('Done!')
end
end

function [D, W] = getPixelData(data, laserframes, W)

if W.doNND %needs to be implemented
    keyboard 
    tic
    data = NN_KP(data, W.tau);
    toc
    
    [AA, BB] = meshgrid(laserframes, 0:W.windowSize-1);
    CGI = griddedInterpolant(data);
    CC = CGI(AA(:)+BB(:));
    CC = reshape(CC,size(AA));
elseif W.doRegression
    keyboard %needs to be reimplemented
    
    %subtract decay
    A = mean(repmat(exp((-W.baselineSamps:-1)/W.tau)', 1, size(AA,2)) .* CC(1:W.baselineSamps,:),1);
    CC = CC(W.baselineSamps+1:end,:) - A*exp(-(0:W.windowSize-1)/W.tau); %decay corrected window samples
    
    if isempty(W.template) %we can feed in a known template for the amplifier, or estimate it, below
        sortedAmps = sort(sum(CC,1));
        nn = round(max(0.5*size(CC,2), min(0.985*size(CC,2), find(sortedAmps<0.06,1,'last'))));
        if isempty(nn); nn = 10; end
        CCmedian = nth_element(CC(:, 1:min(end, 5e6))', nn);
        W.template = CCmedian(nn,:)';
        W.template = W.template / sqrt(sum(W.template.^2)); %Normalize the template
    end
    Trep = repmat(W.template,1,size(CC,2));
    D = sum(Trep.*CC,1)./sum(W.template.^2); %regression
else
    [AA, BB] = meshgrid(laserframes, -W.baselineSamps:W.windowSize-1);
    CGI = griddedInterpolant(data);
    CC = CGI(AA(:)+BB(:));
    CC = reshape(CC,size(AA));
    
    %decay subtraction
    A = mean(repmat(exp((-W.baselineSamps:-1)/W.tau)', 1, size(AA,2)) .* CC(1:W.baselineSamps,:),1);
    F = A*(sum(exp(-(W.windowSamps-1)/W.tau))); %correction factor
    D = sum(CC(W.baselineSamps+W.windowSamps, :), 1)-F;
%     A = mean(CC(1:W.baselineSamps,:),1);
%     F = A*length(W.windowSamps);
%     D = sum(CC(W.baselineSamps+W.windowSamps, :), 1)-F;
end

%Toss negative data
D = D ./ W.photonAmp; %scale to # photons
D = D-median(D(D<0.6));  %subtract offset
D(D<0.5) = 0; %Zero out small values
end

function D = fastDeconv(Dalign, gamma)
%deconvolve an exponential kernel from data
D =conv2(Dalign,[0; 1; -gamma],'same'); %mex'd implementation
end

function peak = centroid(I)
    %find peak in a cyclic signal
    I = I-median(I);
    [~, peak] = max(I); %rough peak
    srnd = max(0, I(mod((peak-2:peak+2)-1, length(I))+1));
    ix = sum((1:5).*srnd)./sum(srnd) - 3;
    peak = mod(peak+ix, length(I));
end

function [T,Td] = V2T(Vx,Vy,Xoff,Yoff,angle)
T = (Vx - Xoff)*cos(angle) - (Vy-Yoff)*sin(angle);
Td = (Vx - Xoff)*sin(angle) + (Vy-Yoff)*cos(angle);
end

function [Vx,Vy] = T2V(T,Td,Xoff,Yoff,angle) %#ok<DEFNU>
Vx = T*cos(angle) + Td*sin(angle) + Xoff;
Vy = -T*sin(angle) + Td*cos(angle) + Yoff;
end

function c = smoothfast1(y, width)
n = length(y); 
c = filter(ones(width,1)/width,1,y);
cbegin = cumsum(y(1:width-2));
cbegin = cbegin(1:2:end)./(1:2:(width-2))';
cend = cumsum(y(n:-1:n-width+3));
cend = cend(end:-2:1)./(width-2:-2:1)';
c = [cbegin;c(width:end);cend];
end

function Y = smoothfast2 (X, m) %#ok<DEFNU> %this algorithm is slightly faster than smoothfast1 but might have numeric issues and does suboptimal things at the data edges
mm = 2*m+1;
Y = [repmat(X(1),m,1) ; X(:) ; repmat(X(end),m,1)] ;
Y = [0 ; cumsum(Y)] ;
Y = (Y(mm+1:end)-Y(1:end-mm)) / mm ;
end
