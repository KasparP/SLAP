function [scandata, hFig] = SLAPMi_reduce(optsin)
import slapmi

%default options
opts.align = false;         tooltips.align = 'align to remove motion, good for in vivo data';
opts.channels = 1;          tooltips.channels = 'which channels to process; if empty do all channels';
opts.ignoreSaved = false;   tooltips.ignoreSaved = 'if data has been reduced before, reduce it again anyways';
opts.doPlot = true;         tooltips.doPlot = 'show plots';
opts.doPCA  = true;         tooltips.doPCA = 'show a plot with NMF';
opts.dffMax = 2;            tooltips.dffMax = 'Max of colorbar for DFF plot';
opts.tau  = 20;             tooltips.tau = 'time constant for DFF plot, in frames';
opts.diodeCh = 0;           tooltips.diodeCh = 'Set this to the channel that is connected to the diode to normalize by excitation power';

if nargin %UPDATE WITH USER-SPECIFIED OPTIONS
    for field = fieldnames(optsin)'
         opts.(field{1}) = optsin.(field{1});
    end
else
   if evalin('base','exist(''hSLAPMi'',''var'')')
      prefPath =  evalin('base', '[hSLAPMi.dataDir filesep hSLAPMi.user filesep ''SLAPMi_reduce_prefs.mat'']');
      if exist(prefPath, 'file')
          load(prefPath);
      end
      opts = optionsGUI(opts, tooltips);
      save(prefPath, 'opts');
   else
      opts = optionsGUI(opts, tooltips);
   end
end

%options not in the GUI
forcesave = false; %This flag will be set if saved data is updated
smoothWindow = 33; %smoothing window for galvo positions; set higher if you have more measurement noise
EF_thresh = 2e-2; %absolute pointing error threshold for E and F mirrors, 2e-2 is equivalent to ~300nm in sample space
dtPMT = 1/(250e6);     tooltips.dtPMT = 'digitizer sample rate';
opts.doNND = false;         tooltips.doNND = 'nonnegative deconvolution of PMT traces';
opts.doRegression = false;  tooltips.doRegression = 'reject dark photons using regression';
opts.amp_setup = 'MPPC';    tooltips.amp_setup = 'supported: ''MPPC''';
hFig = [];
switch opts.amp_setup
    case 'MPPC'
        invert = false; %Whether the amplifier was inverted
        G = 4; %INTEGER; time shift parameter; how many samples before the peak does the photon response start?
        
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
if ~isfield(opts, 'path') || isempty(opts.path)
    basedir = 'E:\SLAPmiData\';
    if evalin('base','exist(''hSLAPMi'',''var'')')
        basedir = [basedir filesep evalin('base', 'hSLAPMi.user')];
    end
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
        tic
        disp('Reading data...')
        try
            [metaData, pmtDataFileName, galvoData, saveData, Zdata] = SLAPMi.openDataFiles([pathstr filesep name]);
        catch ME
            if length(fns)>1
                disp(['Error occurred while loading datafile: ' fns{filenum} '. Continuing...'])
                continue
            else
                rethrow(ME)
            end
        end
        saveData.metaData = metaData; %fold the metadata written at acquisition time into the other saved data
        calib = saveData.calib;
        toc
        
        tic;
        disp('Reducing pixeldata...'); 
        fid = fopen(pmtDataFileName);

        nSamps = metaData.pmtSamplesWritten;
        samplesPerCycle = 50; % samplerate/laserrate
        fActual = 1/(samplesPerCycle*dtPMT);
        channels = opts.channels;
        if isempty(channels) || any(channels>3)
           channels = [1 2]; 
        end
        readChannels = channels;
        if opts.diodeCh
            readChannels = [channels opts.diodeCh];
        end
        chunkSize = 1e6; %~#laser pulses per chunk (not pmt samples)
        nChans = length(channels);
        nOverlap = samplesPerCycle;
        nChunks = ceil((nSamps-nOverlap)/(chunkSize*samplesPerCycle));
        
        framePtr = 0;
        dataStart = 0;
        BG = nan(1,nChans);
        signalChannelIX = find(readChannels~=opts.diodeCh,1);
        
        for chunk = 1:nChunks
            if chunk==1
                dataChunk = [nan(nOverlap,length(readChannels)) ; getDataChunk(fid, metaData, chunkSize*samplesPerCycle, readChannels)]; %get a chunk of data
                
                C = reshape(dataChunk(1:samplesPerCycle*floor(min(chunkSize, 1e4)),signalChannelIX), samplesPerCycle, [])';
                Cavg = nth_element(C, round(0.98*size(C,1))); Cavg = Cavg(round(0.98*size(C,1)), :);
                P = round(centroid(Cavg)); %phase
                laserframes = (mod(P-G-2, samplesPerCycle)+2):samplesPerCycle:nSamps;
                pixeldata = nan(length(laserframes),nChans);
                
                for ch_ix = 1:nChans %estimate background
                    sorted = nth_element(dataChunk(laserframes(1:min(end,chunkSize))-1,ch_ix), ceil(0.01*min(length(laserframes),chunkSize))); 
                    BG(ch_ix) = sorted(ceil(0.01*min(length(laserframes),chunkSize))) + 3e-3;
                    dataChunk(1:nOverlap,ch_ix) = BG(ch_ix);
                end
                
                if opts.diodeCh %if one of the channels is being used for the laser intensity measurement
                    diodeDelay = -6; %the diode produces a signal XXX samples (4ns/sample) before the MPPC
                    diodeData = nan(length(laserframes),1);
                end
            else
               dataChunk = [dataChunk(end-nOverlap+1:end,:) ; getDataChunk(fid, metaData, chunkSize*samplesPerCycle, readChannels)]; 
            end
            
            chunkFrames = laserframes(framePtr+1:min(end, framePtr+chunkSize)) - dataStart + nOverlap; %the laserframes indexed into the chunk
            if length(chunkFrames)<2
                continue
            end
            while chunkFrames(end)>size(dataChunk,1)-W.windowSize
                chunkFrames(end) = [];
            end
            
            for ch_ix = 1:length(channels)
                pixeldata(framePtr+(1:length(chunkFrames)),ch_ix) = getPixelData(dataChunk(:,ch_ix) - BG(ch_ix), chunkFrames, W); %reduce data
            end
            if opts.diodeCh
                diodeData(framePtr+(1:length(chunkFrames))) = dataChunk(chunkFrames+diodeDelay,end) + dataChunk(chunkFrames+diodeDelay-1,end) + dataChunk(chunkFrames+diodeDelay+1,end) - 3*0.005; %- bias
            end
            dataStart = dataStart+chunkSize*samplesPerCycle; %update indexes
            framePtr = framePtr+length(chunkFrames); %update indexes
        end
        fclose(fid);
        clear dataChunk;
        toc
        
        tic
        disp('Reducing galvodata...'); 
        %generate Galvo positions at the pixel data locations
        galvodata =nan(length(laserframes), 6);
        G_INT = griddedInterpolant(galvoData.t, 1:length(galvoData.t));
        galvoframes = G_INT(laserframes*dtPMT + calib.galvos.AIdelay);
        for Gn = 1:6
            tmp = galvoData.data(:,Gn);
            tmp = smoothfast1(smoothfast1(tmp, smoothWindow), smoothWindow-5);
            tGI = griddedInterpolant(tmp);
            galvodata(:,Gn) = calib.galvos.AI2AO{Gn}(tGI(galvoframes));
        end
        clear tGI tmp galvoData G_INT;
        toc
        
        tic;
        disp('Generating Frame Data...'); 
        %subdivide data into frames
        framePeriod = saveData.framePeriod;
        sampsPerFrame = framePeriod/dtPMT; %The galvos are retriggerred according to the sample rate!
        nFrames = floor((nSamps-max(calib.galvos.AIdelay)/dtPMT)/sampsPerFrame);
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
            if opts.diodeCh
                F_diodeData = diodeData(frameStart+1:frameStart+pxPerFrame,:);
            end
            
            frames(frame).galvoData = galvodata((frameStart+1:frameStart+pxPerFrame),:);
            
            %frame data is valid only if the E and F galvos are very close to their
            %target positions
            E_low = abs(frames(frame).galvoData(:,5)-apCalib.galvos.offset.E(1))<EF_thresh;
            E_high = abs(frames(frame).galvoData(:,5)-apCalib.galvos.offset.E(2))<EF_thresh;
            F_low = abs(frames(frame).galvoData(:,6)-apCalib.galvos.offset.F(1))<EF_thresh;
            F_high = abs(frames(frame).galvoData(:,6)-apCalib.galvos.offset.F(2))<EF_thresh;
            if sum(F_low)<length(line)/8 || sum(F_high)<length(line)/8  ||sum(E_low)<length(line)/8 || sum(E_high)<length(line)/8 
                if ~any(diff(frames(frame).galvoData(1:min(end,100),5)))
                    error('The FPGA seems to have stalled. Data is corrupt.')
                else
                    error('The 1D galvos failed to track!');
                end
            end
            EOM_switch = line~=[line(end) ; line(1:end-1)];
            if frame==1
                EOM_switch(1:round(1.75e-6*fActual)) = true; %the delay to turn on the pockels cell on the first frame
            end
            EOM_switch = conv(double(EOM_switch), [zeros(1, ceil(3e-6*fActual)-1) ones(1, ceil(3e-6*fActual))], 'same') >0;
            
            %which frames are valid for reconstruction
            frames(frame).valid = (line==1 & E_low) | (line==2 & E_high) | (line==3 & F_low) | (line==4 & F_high);
            frames(frame).valid = frames(frame).valid & ~EOM_switch;
            frames(frame).pmtData = nan(pxPerFrame, nChans); %extize
            
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
                sortedP = F_pmtData(linSelect,:); sortedP = sortedP(sortorder,:);
                CSP = cumsum(sortedP,1);
                CSP_ix = interp1qr(sortedT,(1:length(sortedT))', desTedges(:,lineID));
                FP = diff(qinterp1(1:size(CSP,1), CSP, CSP_ix))./diff(CSP_ix);
                frames(frame).pmtData(line==lineID,:) = FP;
                if opts.diodeCh
                    sortedD = F_diodeData(linSelect,:); sortedD = sortedD(sortorder,:).^2; %Squared here to reflect two-photon excitation, before interpolation!
                    CSD = cumsum(sortedD,1);
                    FD = diff(qinterp1(1:size(CSD,1), CSD, CSP_ix))./diff(CSP_ix);
                    frames(frame).diodeData(line==lineID,:) = FD;
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
        
        %Normalize PMT data to diode data
        if opts.diodeCh
            P = [frames.pmtData]; %3D not implemented!
            D = [frames.diodeData];
            D = D./repmat(mean(D,2),1,size(D,2)); %normalize to local transmission
            P = P./D;
            %remove average correlation of diode to PMT; this is minor and can be removed
            DD = nanmean(D,1)'-smooth(nanmean(D,1));
            PP = nanmean(P,1)'-smooth(nanmean(P,1));
            b = regress(PP,DD);
            P = P - repmat(b*DD', size(P,1),1);
            for frame = 1:numel(frames)
                frames(frame).pmtData = max(0, P(:,frame));
            end
        end
      
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
    
    %add stimulus
    stimfiles = dir([dr filesep '*Timings*']);
    if ~isempty(stimfiles)
        %add stimulus information
        stimList = [];
        StimDelay = -3.2e-5;
        for sf = 1:length(stimfiles)
            SF = load([dr filesep stimfiles(sf).name]);
            stimList = [stimList ; SF.time_stamp]; %#ok<AGROW>
            if isfield(SF, 'StimDelay')
                StimDelay = SF.StimDelay;
            end
        end
        [scandata.stimulus.timeError, minIx] = min(abs(scandata.metadata.timeNow - stimList(:,1) + StimDelay));
        scandata.stimulus.stim = stimList(minIx,2);
        disp(['Stimulus Information Added. Timing error:' num2str(scandata.stimulus.timeError*86400,2) ' seconds'])
        forcesave = true; %force save
    end
    
    %save the reduced data
    if ~exist([pathstr filesep name '_REDUCED.mat'], 'file') || opts.ignoreSaved || forcesave
        disp('Saving...')
        bytedata = whos('scandata');
        if bytedata.bytes>2e9
            save([pathstr filesep name '_REDUCED.mat'], 'scandata', '-v7.3') %v6 gives faster save but has max file size
        else  
            save([pathstr filesep name '_REDUCED.mat'], 'scandata', '-v6') %v6 gives faster save but has max file size
        end
    end
    
    %plot
    if opts.doPlot
        disp('Plotting...')
        if length(fns)>20 %batch mode; don't leave the windows open
            try %#ok<TRYNC>
                close(hFig) %close the figures from the previous file
            end
        end
        hFig = scandataPlot(scandata, opts);
        if ~scandata.metadata.do3D
            disp('Saving Plots...')
            for figIX = 1:length(hFig)
                try %#ok<TRYNC>
                    saveas(hFig(figIX), [pathstr filesep name '_PLOT' int2str(figIX) '.fig']);
                end
            end
            drawnow;
        end
    end
    disp('Done!')
end
end

function data = getDataChunk(fid, metaData, nsamples, channels)
%process pmt data file
% each element is a u64 packed with 4 i16 values:
% [ch0 sample 1, ch0 sample 2, ch1 sample 1, ch1 sample 2]
dat = fread(fid, nsamples/2, '*uint64');
dat = typecast(dat, 'uint32');
if all(sort(channels) == [1 2]) || (length(channels)==1 && channels>2)
        data = double([typecast(dat(1:2:end), 'int16') typecast(dat(2:2:end), 'int16')]) * (metaData.inputRange/2) / 2^15;
elseif channels == 2
        data = double(typecast(dat(2:2:end), 'int16')) * (metaData.inputRange/2) / 2^15;
elseif channels == 1
        data = double(typecast(dat(1:2:end), 'int16')) * (metaData.inputRange/2) / 2^15;
else
    error('Unsupported channels argument')
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
    CC = data(AA(:)+BB(:));
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
z = median(D(D<0.6)); 
if isnan(z)
    warning('Your sample seems to have become extremely bright');
    z=-1.5;
end
D = D-z;  %subtract offset
D(D<0.5) = 0; %Zero out small values
end

function peak = centroid(I)
    %find peak in a cyclic signal
    I = I-median(I);
    [~, peak] = max(I); %rough peak
    srnd = max(0, I(mod((peak-2:peak+2)-1, length(I))+1));
    ix = sum((1:5).*(srnd+eps))./sum(srnd+eps) - 3;
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
