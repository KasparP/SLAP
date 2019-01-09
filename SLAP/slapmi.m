classdef slapmi < handle
    %% USER PROPS
    properties (SetObservable)
        % general settings
        logData = true;
        retriggerable = false;
        stopGalvosWhileLogging = true;
        logAIs = false;  %a logical array of channels to log from the BeamDAQ on each acquisition. The first element refers to AI1, second to AI2, etc. AI0 is always reserved for logging the piezo.
        fileNameStem;
        fileNameApp;
        pmtDataFileName;
        galvoDataFileName;
        ZDataFileName;cc
        aiDataFileName;
        metaDataFileName;
        metaFileName;
        calib;  %calibration file
        PSFdataFN; %path to a PSF measurement obtained with measurePSF
        PSF;
        refIMFN = [];
        
        framesToCollect = 1000;
        samplesPerLaserClk = 10;
        samplesExpected;        %  #of samples per frame expected from FPGA, = number of laser clocks*samples per clock /2
        inputRange = 5;
        user = 'Kaspar';
        dataDir = 'E:\SLAPmidata\';
        
        useExternalClock = true;
        syncDaqsToFpga = true;
        
        % acquisition plan
        scanLineLength = 8;         % [degrees] length of scan line
        scanCenterPoint = [0 0];    % center point of region to scan
        tiling = 1;                 % tiling factor; how many lines are scanned across FOV
        aperture = 0.55;            %the length of the scan to perform; this is the radius of the the enclosing circle around the center of the SLM
        
        power = 10;
        do3D = false;              %do we use the fast Piezo for 3D imaging?
        analysisFrames = 300;       %number of frames to store in buffer for analysis
        timing = 'hardware';
        nLines = 4;             %number of projection angles to use, generally 4;
        res = 1222;
        framePeriod = nan;           %line time in ms
        volumePeriod = 0.05;         %volume time in ms
        Zmin = 0;
        Zmax = 0;
        dZ = [];
        ZBiDi = false;
        nPlanes = [];
        nVolumes = [];
        acqTime = [];
        binFactor = 1;          %Number of laser pulses to bin into a single pixel for display
        tunePiezo = true;
        
        % hardware settings
        rioDevice = 'RIO1';
        scannerDaqName = 'DEV7';   % scanner daq has 4 galvos and 3 pockels cells
        beamDaqName = 'DEV5';      % power control pockels cell
        piezoDaqName = 'DEV6';     % z piezo
        scannerDaqChannels = [0 4 8 12 16 20 24];
        startTrigTerm = [];
        acqStartPulseTerm = 0;
        laserClockTerm = 0;
        frameClockTerm = 0;
        laserClkDelay = 1.9e-7;
        SLMctrNumber = 1; 
        SLMoutputTerminal = 'PFI14';
        
        % ao waveforms
        powerWaveformRaw;
        powerWaveform;
        galvoWaveform;
        beamWaveform;
        piezoWaveform;
        lineIDs; %which line is being imaged at each point in the galvoWaveform
        
        piezoDesiredX;
        
        %data
        status = 'Idle';
        streamCorrupt = false;
        dispCorrupt = false;
        latestFrameTag;
        latestFrame;
        
        %hardware properties; these could be loaded as Datafiles
        piezo = struct('AOrate', 5e5, 'AIrate', 1e4, 'h', [], 'starting', false, 'verbose', true); 
        pockels = struct('deadtime', 1.5e-6);
        laser = struct('repRate', 5E6, 'lineLength', 200);
        galvos = struct('calFramePeriod', 0, 'calPos', 0, 'AOrate', 1e6);
        SLM = struct('ROIs',[]);
        auxAI = struct('nSamples',0, 'AIrate', 5e4);
        SLMCounter=struct('COrate',1e2,'EveryNFrame',4,'lowTicks',2,'highTicks',2);

        
        sampleRate = 500000; %the minimum sample rate we need to synchronize acquisition to
    end
    
    properties (Hidden)
        hPiezoFig;
        hPiezoAx;
        hPiezoLines;
    end
    
    %% INTERNAL PROPS
    properties (SetObservable, Hidden)
        hSI;
        hGui;
        hFpga;
        hSLM;
        
        hTaskPowerBeam;
        hTaskPowerBeamOnDmd;
        hTaskGalvoBeams;
        hTaskGalvo567OnDmd;
        hTaskPiezo;
        hTaskPiezoFeedback;
        hTaskPiezoLogger;
        hTaskAILogger;
        hTaskFlipMirror;
        hTaskShutter;
        hTaskTriggerOutOnStart;
        hTaskTriggerOutOnDmd;
        hTaskSLMTriggerOut;
        hRouteRegistry;
        
        
        laserClockFullTermInt;
        laserClockFullTermExt;
        frameClockFullTermExt;
        frameClockFullTermInt;
        volClockFullTermInt;
        trgdVolClockFullTermInt;
        refClockTermInt = 'PXI_Trig7';
        
        
        %timers
        hTimerPollData;
        hTimerPiezo;
        
        %fpga session info
        shareSiFpga = false;
        siFpgaRunning = false;
        fpgaRunning = false;
        bitfilePath;
        fpgaSessionId;
        pmtFifoNum;
        galvoFifoNum;
        
        metaHeaderString;
        
        wasRunning = false;
        
        %datalogger params
        % This controls the rolling memory buffer.
        % 1000 pages * 250000 elements/page * 8 bytes/element = 2GB
        % The datalogger will attempt to read one page at a time from the
        % fifo. Increasing the page size will increase how much data is
        % attempted to be read each operation. Increasing the number of
        % pages will increase the size of the buffer without changing how
        % much data is attempted to be read each operation.
        memoryPages = 5000;             %number of blocks to allocate for rotating data buffer
        memoryPageSize = 250000;        %size of each block in number of elements. Each element is 8 bytes
    end
    
    %% LIFECYCLE
    methods
        function obj = slapmi(hSI,simulate)
            % load DAQmx
            [~] = dabs.ni.daqmx.System();
            
            if nargin > 0
                obj.hSI = hSI;
            end
            
            if nargin < 2
                simulate = false;
            end
            
            %Update the minimum samplerate
            obj.sampleRate = min([obj.piezo.AOrate, obj.galvos.AOrate]);
            
            % initialize fpga
            if ~simulate
                obj.bitfilePath = fullfile(fileparts(mfilename('fullpath')), 'FPGA\Bitfiles\KasparScan5170.lvbitx');
                obj.hFpga = dabs.ni.rio.NiFPGA(obj.bitfilePath,false);
                obj.shareSiFpga = most.idioms.isValidObj(obj.hSI) && obj.hSI.fpgaMap.isKey(obj.rioDevice);
                obj.siFpgaRunning = obj.shareSiFpga;
                
                if ~obj.shareSiFpga
                    fprintf('Initializing FPGA... ');
                    dabs.ni.oscope.clearSession();
                    err = dabs.ni.oscope.startSession(obj.rioDevice,obj.bitfilePath);
                    assert(err == 0, 'Error when attempting to connect to NI 517x device. Code = %d', err);
                    
                    err = dabs.ni.oscope.configureSampleClock(obj.useExternalClock,250000000);
                    assert(err == 0, 'Error configuring external reference clock: %d . Is the laser on?', err);
                    
                    obj.hFpga.openSession(obj.rioDevice);
                    obj.hFpga.reset();
                    obj.hFpga.run();
                    obj.fpgaRunning = true;
                    obj.configOscope();
                    disp('Done!');
                end
            end
            
            %Park the laser
            obj.hTaskGalvo567OnDmd = most.util.safeCreateTask('SLAPMIGalvos567OnDmd');
            obj.hTaskGalvo567OnDmd.createAOVoltageChan(obj.scannerDaqName, obj.scannerDaqChannels(end-2:end));
            obj.hTaskGalvo567OnDmd.writeAnalogData([-2 -2 0]); %park galvos
            delete(obj.hTaskGalvo567OnDmd);
            obj.hTaskPowerBeamOnDmd = most.util.safeCreateTask('SLAPMIPowerBeamOnDemand');
            obj.hTaskPowerBeamOnDmd.createAOVoltageChan(obj.beamDaqName, 0);
            obj.hTaskPowerBeamOnDmd.writeAnalogData(0); %turn off beam
            delete(obj.hTaskPowerBeamOnDmd);
            
            %check system ram
            [~, b] = memory();
            gb = b.PhysicalMemory.Total / 2^30;
            if gb < 24
                most.idioms.warn('This system only has %.0fGB of RAM. Buffer size is set to 2GB.',gb);
                obj.memoryPages = 1000;
            end
            
            % initialize acquisition mex
            fprintf('Initializing datalogger... ');
            obj.fileNameStem = [datestr(now,'yyyy-mm-dd') '_'];
            obj.fileNameApp = 0;
            pgs = AcquisitionMex('init',obj);
            disp('Done!');
            pgSz = obj.memoryPageSize*.000000008;
            mem = double(pgs)*pgSz;
            if pgs ~= obj.memoryPages
                most.idioms.warn('Failed to allocate requested memory. Only %d pages (%f GB) allocated.', pgs, mem);
            else
%                 fprintf('Allocated requested memory (%d pages, %.1f GB).\n', pgs, mem);
            end
            
            obj.hTimerPollData = timer();
            obj.hTimerPollData.TimerFcn = @obj.pollFpgaShortData;
            obj.hTimerPollData.ExecutionMode = 'fixedSpacing';
            obj.hTimerPollData.Period =  0.03;
            obj.hTimerPollData.Name = 'SLAPMi data poll timer';
            
            obj.hTimerPiezo = timer();
            obj.hTimerPiezo.StartDelay = 0.5;
            obj.hTimerPiezo.TimerFcn = @(varargin)obj.hTaskPiezoFeedback.start();
            obj.hTimerPiezo.Name = 'SLAPMi piezo feedback timer';
            
            try
                obj.calib = SLAPMi.loadCalibrationData([obj.dataDir filesep 'Calibration\calibration.cal']);
            catch
                most.idioms.warn('Failed to load calibration data.');
            end
            
            %load PSF data
            try
                [FN, PSFdataDR]= uigetfile([obj.dataDir filesep 'PSF\*.psf']);
                obj.PSFdataFN = [PSFdataDR FN];
                S = load([PSFdataDR filesep FN], '-mat');
                obj.PSF = rmfield(S.data, 'testframe');
            catch
                warning('No PSF data loaded. You will have to manually select PSFdata during data reduction!')
            end
        end
        
        function delete(obj)
            obj.endAcq();
            
            obj.fpgaRunning = false;
            delete(obj.hFpga);
            obj.restoreSiFpga();
            if ~obj.shareSiFpga
                dabs.ni.oscope.clearSession();
            end
            
            AcquisitionMex('exit');
            clear('AcquisitionMex');
            most.idioms.safeDeleteObj(obj.hFpga);
            most.idioms.safeDeleteObj(obj.hTimerPollData);
            most.idioms.safeDeleteObj(obj.hTimerPiezo);
        end
        
        function exit(obj)
            if most.idioms.isValidObj(obj.hGui)
                obj.hGui.exit(true);
            else
                delete(obj);
                evalin('base','clear hSLAPMi');
            end
        end
    end
    
    methods
        function takeFpga(obj)
            if obj.siFpgaRunning
                msg = true;
                fprintf('Acquiring FPGA resource from ScanImage... ');
                
                hSiFpga = obj.hSI.fpgaMap(obj.rioDevice).hFpga;
                try
                    hSiFpga.closeSession();
                catch
                end
                obj.siFpgaRunning = false;
            else
                msg = false;
            end
            
            if ~obj.fpgaRunning
                if ~msg
                    msg = true;
                    fprintf('Acquiring FPGA resource from ScanImage... ');
                end
                
                dabs.ni.oscope.clearSession();
                err = dabs.ni.oscope.startSession(obj.rioDevice,obj.bitfilePath);
                assert(err == 0, 'Error when attempting to connect to NI 517x device. Code = %d', err);

                err = dabs.ni.oscope.configureSampleClock(obj.useExternalClock,250000000);
                assert(err == 0, 'Error configuring external reference clock: %d . Is the laser on?', err);

                obj.hFpga.openSession(obj.rioDevice);
                obj.hFpga.reset();
                obj.hFpga.run();
                obj.configOscope();
                obj.fpgaRunning = true;
            end
            
            if msg
                disp('Done!');
            end
        end
        
        function restoreSiFpga(obj)
            if obj.shareSiFpga
                if obj.fpgaRunning
                    msg = true;
                    fprintf('Returning FPGA resource to ScanImage... ');
                    
                    obj.hFpga.closeSession();
                    obj.fpgaRunning = false;
                else
                    msg = false;
                end
                
                if ~obj.siFpgaRunning
                    if ~msg
                        msg = true;
                        fprintf('Returning FPGA resource to ScanImage... ');
                    end
                    
                    dabs.ni.oscope.clearSession();
                    err = dabs.ni.oscope.startSession(obj.rioDevice,obj.hSI.fpgaMap(obj.rioDevice).bitfilePath);
                    assert(err == 0, 'Error when attempting to connect to NI 517x device. Code = %d', err);
                    hSiFpga = obj.hSI.fpgaMap(obj.rioDevice).hFpga;
                    hSiFpga.openSession(obj.rioDevice);
%                     hSiFpga.reset();
%                     hSiFpga.run();
                    obj.siFpgaRunning = true;
                    
                    % ensure oscope settings are restored to what SI wants
                    obj.hSI.hScan2D.channelsInputRanges = obj.hSI.hScan2D.channelsInputRanges;
                end
                
                if msg
                    disp('Done!');
                end
            end
        end
        
        function piezoCallbackSetH(obj,task,evt)
            disp('Setting H...'); tic;
           
            D = [evt.data];
            NAI = round(obj.piezo.AIrate*obj.volumePeriod);
            NAO = length(obj.piezoWaveform);
            
            obj.hTaskPiezoFeedback.stop();
            
            discardN = ceil(0.5/obj.volumePeriod); %discard fist half second of data
            h = reshape(D, NAI, []);
            h = mean(h(:,discardN+1:end),2)/obj.piezo.delta;
            %h = circshift(h, [-obj.piezo.hDelay 0]);
            h = h-mean(h);
            
            K = convmtx(h, NAI);
            K(1:length(h),:) = K(1:length(h),:) + K([NAI+1:end 1], :);
            K = K(1:NAI,:);
            DD = convmtx([1 -2 1], NAI);
            DD(:, 1:2) = DD(:,1:2) + DD(:, end-1:end);
            DD = DD(:,1:NAI);
            DD = DD'*DD;
            
            lambda = 1e6 * norm(K'*K, 'fro')/norm(DD, 'fro'); %1000000 1e6
            tic
            KK = pinv(K'*K + lambda*DD)*K';
            toc
            
            obj.piezo.h = h;
            obj.piezo.KK = KK;
            obj.piezo.starting = false;

            lsNAI = linspace(0,1,NAI+1);
            lsNAO = linspace(0,1,NAO+1);
            obj.piezo.target = obj.piezoDesiredX'-mean(obj.piezoDesiredX);
            obj.piezo.offset = mean(obj.piezoDesiredX);
            err_i = obj.piezo.KK*obj.piezo.target; err_i = err_i - mean(err_i);
            obj.piezoWaveform = (interp1(lsNAI, err_i([end 1:end]), lsNAO(1:end-1)) + obj.piezo.offset)';
            obj.hTaskPiezo.writeAnalogDataAsync(obj.piezoWaveform(1:end-1),[],[],[],[]);

            disp('Piezo response set');
%                         obj.abortAcq()
%                         obj.endAcq();
%                         figure, plot(h);
%                         keyboard
            
            %setup update event
            obj.hTaskPiezoFeedback.cfgSampClkTiming(obj.piezo.AIrate, 'DAQmx_Val_FiniteSamps', NAI);
            if obj.syncDaqsToFpga
                set(obj.hTaskPiezoFeedback,'sampClkTimebaseSrc', obj.refClockTermInt);
                set(obj.hTaskPiezoFeedback,'sampClkTimebaseRate', 5000000);
            end
            obj.hTaskPiezoFeedback.cfgInputBuffer(NAI);
            obj.hTaskPiezoFeedback.cfgDigEdgeStartTrig(obj.volClockFullTermInt);
            obj.hTaskPiezoFeedback.registerEveryNSamplesEvent(@obj.piezoCallbackUpdateAO,NAI,true);

            disp('Done setting H...'); toc;
            start(obj.hTimerPiezo);
        end
        function piezoCallbackUpdateAO(obj,task,evt)
            obj.hTaskPiezoFeedback.stop();
            if obj.piezo.starting
                obj.piezo.starting = false;
            else
                D = [evt.data];
                alpha = 0.5; %step size
                NAI = length(D);
                NAO = length(obj.piezoWaveform);
                lsNAI = linspace(0,1,NAI+1);
                lsNAO = linspace(0,1,NAO+1);
                
                err = D - obj.piezoDesiredX';
                err_m = mean(err); err = err - err_m;
                
                obj.piezo.offset = obj.piezo.offset - alpha*err_m;
                obj.piezo.target = obj.piezo.target - alpha*err;
                err_i = obj.piezo.KK*obj.piezo.target; %regularized deconvolution with pseudoinverse
                err_i = err_i - mean(err_i); %offset is solved separately because of difficulty measuring offset in h
                obj.piezoWaveform = max(0, min(10,interp1(lsNAI, err_i([end 1:end]), lsNAO(1:end-1))' + obj.piezo.offset));
                
%                 err_i =  obj.piezo.KK * err; 
%                 err_i = err_i - mean(err_i);
%                 err_i = interp1(lsNAI, err_i([end 1:end]), lsNAO(1:end-1));
%                 obj.piezoWaveform = max(0, min(10, obj.piezoWaveform - alpha.*(err_i' + err_m)));
                
                obj.hTaskPiezo.writeAnalogDataAsync(obj.piezoWaveform(1:end-1),[],[],[],[]);
                
                if obj.piezo.verbose && ishghandle(obj.hPiezoAx)
                    set(obj.hPiezoLines(1), 'Ydata', obj.piezoDesiredX([1:end 1]))
                    set(obj.hPiezoLines(2), 'Ydata',D([1:end 1]))
                    set(obj.hPiezoLines(3), 'Ydata',obj.piezoWaveform([1:end 1]))
                end
            end
            start(obj.hTimerPiezo);
        end
        
        
        % generate ao and prepare fpga
        function initAcq(obj)
            % make sure everything is clear
            obj.endAcq(true);
            
            %Prereqs
            if ~(obj.galvos.calFramePeriod==obj.framePeriod) || ~(all(obj.galvos.calPos == [obj.scanCenterPoint obj.tiling obj.aperture]))
                warning('Galvos not calibrated. You must accept a galvo calibration to initialize acquisition');
                return
            end
            
            %ensure the data directory exists
            if ~exist([obj.dataDir filesep obj.user], 'dir')
                mkdir([obj.dataDir filesep obj.user]);
            end
            
            % Clock routing
            obj.laserClockFullTermExt = sprintf('/%s/PFI%d', obj.piezoDaqName, obj.laserClockTerm);
            obj.laserClockFullTermInt = sprintf('/%s/PXI_Trig3', obj.piezoDaqName);
            obj.frameClockFullTermExt = sprintf('/%s/PFI%d', obj.scannerDaqName, obj.frameClockTerm);
            obj.frameClockFullTermInt = sprintf('/%s/PXI_Trig4', obj.scannerDaqName);
            obj.volClockFullTermInt = sprintf('/%s/PXI_Trig5', obj.scannerDaqName);
            obj.trgdVolClockFullTermInt = sprintf('/%s/PXI_Trig6', obj.beamDaqName);
            
            obj.hRouteRegistry = dabs.ni.daqmx.util.triggerRouteRegistry();
            %obj.hRouteRegistry.connectTerms(obj.frameClockFullTermInt, obj.frameClockFullTermExt); % frame clock out: DEV7/PFI0
            obj.hRouteRegistry.connectTerms(obj.trgdVolClockFullTermInt, obj.frameClockFullTermExt); %export triggered volume clock to DEV7/PFI0
            obj.hRouteRegistry.connectTerms(obj.laserClockFullTermExt, obj.laserClockFullTermInt); % laser clock in 
            
            % Create Galvo/beam/piezo tasks
            fprintf('Initializing DAQ tasks... ');
            
            obj.hTaskTriggerOutOnStart = most.util.safeCreateTask('SLAPMITriggerOut');
            obj.hTaskTriggerOutOnStart.createCOPulseChanTime(obj.beamDaqName, 0, '', 1e-3, 5e-3, 0); %low time, high time, initialdelay
            obj.hTaskTriggerOutOnStart.cfgImplicitTiming('DAQmx_Val_FiniteSamps',1);
            obj.hTaskTriggerOutOnStart.cfgDigEdgeStartTrig(obj.trgdVolClockFullTermInt);
            obj.hTaskTriggerOutOnStart.set('startTrigRetriggerable',false);
            if isempty(obj.acqStartPulseTerm)
                obj.hTaskTriggerOutOnStart.channels(1).set('pulseTerm',[]);
            else
                obj.hTaskTriggerOutOnStart.channels(1).set('pulseTerm',sprintf('PFI%d', obj.acqStartPulseTerm));
            end
            
            obj.hTaskPowerBeam = most.util.safeCreateTask('SLAPMIPowerBeam');
            obj.hTaskPowerBeam.createAOVoltageChan(obj.beamDaqName, 0);
            
            obj.hTaskPowerBeamOnDmd = most.util.safeCreateTask('SLAPMIPowerBeamOnDemand');
            obj.hTaskPowerBeamOnDmd.createAOVoltageChan(obj.beamDaqName, 0);
            
            obj.hTaskGalvoBeams = most.util.safeCreateTask('SLAPMIGalvosBeams');
            obj.hTaskGalvoBeams.createAOVoltageChan(obj.scannerDaqName, obj.scannerDaqChannels);
            
            obj.hTaskGalvo567OnDmd = most.util.safeCreateTask('SLAPMIGalvos567OnDmd');
            obj.hTaskGalvo567OnDmd.createAOVoltageChan(obj.scannerDaqName, obj.scannerDaqChannels(end-2:end));
            obj.hTaskGalvo567OnDmd.writeAnalogData([-2 -2 0]); %park galvos

            obj.hTaskPiezo = most.util.safeCreateTask('SLAPMIPiezo');
            obj.hTaskPiezo.createAOVoltageChan(obj.piezoDaqName, 0);
            
            obj.hTaskPiezoFeedback = most.util.safeCreateTask('SLAPMIPiezoFeedback');
            obj.hTaskPiezoFeedback.createAIVoltageChan(obj.piezoDaqName, 0);
            
            obj.hTaskPiezoLogger = most.util.safeCreateTask('SLAPMIPiezoLogger');
            obj.hTaskPiezoLogger.createAIVoltageChan(obj.beamDaqName, 0); %This is also the Piezo feedback for Scanimage
            
            if any(obj.logAIs)
                obj.hTaskAILogger = most.util.safeCreateTask('SLAPMIAILogger');
                obj.hTaskAILogger.createAIVoltageChan(obj.beamDaqName, find(obj.logAIs)); %#ok<FNDSB> %This is also the Piezo feedback for Scanimage
            end
            
            obj.hTaskFlipMirror = most.util.safeCreateTask('FlipMirror');
            obj.hTaskFlipMirror.createDOChan(obj.beamDaqName, 'port0/line0' , 'RasterFlip');
            obj.hTaskFlipMirror.writeDigitalData(false, 1, true); % FALSE for line scan
            
            obj.hTaskShutter = most.util.safeCreateTask('Shutter');
            obj.hTaskShutter.createDOChan(obj.piezoDaqName, 'port0/line7' , 'ShutterTTL');
            obj.hTaskShutter.writeDigitalData(false, 1, true); % FALSE for line scan
            
            obj.hTaskSLMTriggerOut = most.util.safeCreateTask('SLMTriggerOut');
            obj.hTaskSLMTriggerOut.createCOPulseChanTicks(obj.piezoDaqName, 1, '',obj.frameClockFullTermInt, 5, 5);

            obj.sampleRate = min([get(obj.hTaskPowerBeam, 'sampClkMaxRate'),...
                get(obj.hTaskGalvoBeams, 'sampClkMaxRate'),...
                get(obj.hTaskPiezo, 'sampClkMaxRate'),...
                500000]);
            
            disp('Done!');
            
            % AO generation
                %Galvo AO was already generated if prereqs are met
            
            %Prepare the live display GUI
            NLaserClcks = round(obj.framePeriod*obj.laser.repRate);
            obj.samplesExpected = NLaserClcks*obj.samplesPerLaserClk/2; %we pack pairs of PMT samples into a single FPGA sample
            obj.hGui.decimationMap = reshape(1:obj.samplesExpected, obj.samplesPerLaserClk/2, NLaserClcks);
            obj.hGui.decimationMap = obj.hGui.decimationMap(:,1:obj.hGui.dataDecimation:end);
            obj.hGui.decimationMap = obj.hGui.decimationMap(:);
            Ndecimated= ceil(NLaserClcks/obj.hGui.dataDecimation);
            obj.hGui.frameDataMap = reshape(1:NLaserClcks*obj.samplesPerLaserClk, obj.samplesPerLaserClk*obj.hGui.dataDecimation, Ndecimated);           
            
            % Write AO
            fprintf('Starting tasks... ');
            
            N_AO_galvo = size(obj.galvoWaveform,1)-1;
            
            obj.hTaskGalvoBeams.cfgSampClkTiming(obj.galvos.AOrate, 'DAQmx_Val_FiniteSamps', N_AO_galvo);
            if obj.syncDaqsToFpga
                set(obj.hTaskGalvoBeams,'sampClkTimebaseSrc', obj.refClockTermInt);
                set(obj.hTaskGalvoBeams,'sampClkTimebaseRate', 5000000);
            end
            obj.hTaskGalvoBeams.cfgOutputBuffer(N_AO_galvo);
            obj.hTaskGalvoBeams.cfgDigEdgeStartTrig(obj.frameClockFullTermInt);
            obj.hTaskGalvoBeams.set('startTrigRetriggerable',true);
            obj.hTaskGalvoBeams.writeAnalogData([obj.galvoWaveform(1:end-1,:) obj.beamWaveform(1:end-1,:)]);
            obj.hTaskGalvoBeams.start();
            
            %Piezo AO
            if obj.do3D
                %generate the piezo trace
                obj.piezo.starting = true;
                obj.piezo.AIrate = 10/(obj.framePeriod);
                obj.piezo.AOrate = 10*obj.piezo.AIrate;
                N_AI_piezo = round(obj.piezo.AIrate*obj.volumePeriod); %floor(obj.piezo.AIrate*obj.volumePeriod)-1;
                N_AO_piezo = round(obj.piezo.AOrate*obj.volumePeriod); %N_AI_piezo * round(obj.piezo.AOrate/obj.piezo.AIrate);

                if isempty(obj.piezoWaveform) %create a new piezo tuning
                    pOpts.Zmin = obj.Zmin/40;
                    pOpts.Zmax = obj.Zmax/40;
                    pOpts.duty = 0.85;
                    pOpts.bidi = obj.ZBiDi;
                    [obj.piezoWaveform, obj.powerWaveform] = piezo_trace_SLAPmi(N_AO_piezo,pOpts, obj.calib); %slow piezo (power waveform) is coupled to the piezo oscillation
                    obj.piezo.opts = pOpts;
                    obj.piezoDesiredX = interp1(linspace(0,1,N_AO_piezo+1), obj.piezoWaveform([1:end 1]), linspace(0,1,N_AI_piezo+1));
                    obj.piezoDesiredX = obj.calib.piezo.AO2AI(obj.piezoDesiredX(1:end-1));
                end
                
                %Piezo display
                if obj.piezo.verbose
                    if most.idioms.isValidObj(obj.hPiezoFig)
                        clf(obj.hPiezoFig)
                    else
                        obj.hPiezoFig = figure('name', 'Piezo Control', 'NumberTitle', 'off');
                    end
                   obj.hPiezoAx = axes('parent',obj.hPiezoFig); set(obj.hPiezoAx, 'xlim', [0 1], 'ylim', [0 10]);
                   obj.hPiezoLines(1) = plot(obj.hPiezoAx, linspace(0,1, N_AI_piezo+1), obj.piezoDesiredX([1:end 1]), 'b','parent',obj.hPiezoAx);
                   hold(obj.hPiezoAx,'on');
                   obj.hPiezoLines(2) = plot(obj.hPiezoAx, linspace(0,1, N_AI_piezo+1), nan(1,N_AI_piezo+1) , 'r','parent',obj.hPiezoAx);
                   obj.hPiezoLines(3) = plot(obj.hPiezoAx, linspace(0,1, length(obj.piezoWaveform)+1), obj.piezoWaveform([1:end 1]), 'k','parent',obj.hPiezoAx);
                   
                   legend(obj.hPiezoAx, {'Desired', 'Measured', 'AO'}); xlabel(obj.hPiezoAx,'Time'); ylabel(obj.hPiezoAx, 'Volts');
                end
                
                if length(obj.piezo.h)~=N_AI_piezo  %measure a new impulse response
                    obj.piezo.starting = true;
                    obj.piezo.h = [];
                    obj.piezo.KK = [];
                    %obj.piezo.hDelay = 0;
                    
                    %AO
                    obj.piezo.delta = 1; %amplitude of impulse for measuring impulse response
                    piezoAO = mean(obj.piezoWaveform)*ones(size(obj.piezoWaveform));
                    piezoAO(1:round(N_AO_piezo/N_AI_piezo)) = piezoAO(1) + obj.piezo.delta;
                    %AI
                    nAvg = max(7, ceil(4/obj.volumePeriod)); %go for 4 seconds but always do at least 7 cycles
                    obj.hTaskPiezoFeedback.cfgSampClkTiming(obj.piezo.AIrate, 'DAQmx_Val_FiniteSamps', nAvg*N_AI_piezo);
                    if obj.syncDaqsToFpga
                        set(obj.hTaskPiezoFeedback,'sampClkTimebaseSrc', obj.refClockTermInt);
                        set(obj.hTaskPiezoFeedback,'sampClkTimebaseRate', 5000000);
                    end
                    obj.hTaskPiezoFeedback.cfgInputBuffer(nAvg*N_AI_piezo);
                    obj.hTaskPiezoFeedback.cfgDigEdgeStartTrig(obj.volClockFullTermInt);
                    obj.hTaskPiezoFeedback.registerEveryNSamplesEvent(@obj.piezoCallbackSetH,nAvg*N_AI_piezo,true);
                    obj.hTaskPiezoFeedback.start();
                else %use old impulse response
                    %AO
                    piezoAO = obj.piezoWaveform;
                    %AI
                    obj.hTaskPiezoFeedback.cfgSampClkTiming(obj.piezo.AIrate, 'DAQmx_Val_FiniteSamps', N_AI_piezo);
                    if obj.syncDaqsToFpga
                        set(obj.hTaskPiezoFeedback,'sampClkTimebaseSrc', obj.refClockTermInt);
                        set(obj.hTaskPiezoFeedback,'sampClkTimebaseRate', 5000000);
                    end
                    obj.hTaskPiezoFeedback.cfgInputBuffer(N_AI_piezo);
                    obj.hTaskPiezoFeedback.cfgDigEdgeStartTrig(obj.volClockFullTermInt);
                    obj.hTaskPiezoFeedback.registerEveryNSamplesEvent(@obj.piezoCallbackUpdateAO,N_AI_piezo,true);
                    obj.hTaskPiezoFeedback.start();
                end
                
                obj.hTaskPiezo.cfgSampClkTiming(obj.piezo.AOrate, 'DAQmx_Val_FiniteSamps', N_AO_piezo-1);
                if obj.syncDaqsToFpga
                    set(obj.hTaskPiezo,'sampClkTimebaseSrc', obj.refClockTermInt);
                    set(obj.hTaskPiezo,'sampClkTimebaseRate', 5000000);
                end
                obj.hTaskPiezo.cfgOutputBuffer(N_AO_piezo-1); %N_AO_piezo-1 in order to ensure proper retriggering; not tested whether this is necessary
                obj.hTaskPiezo.cfgDigEdgeStartTrig(obj.volClockFullTermInt);
                obj.hTaskPiezo.set('startTrigRetriggerable',true);
                obj.hTaskPiezo.writeAnalogData(piezoAO(1:end-1));
                obj.hTaskPiezo.start();
                
                obj.hTaskPiezoLogger.cfgSampClkTiming(obj.piezo.AIrate, 'DAQmx_Val_FiniteSamps', N_AI_piezo-1);
                if obj.syncDaqsToFpga
                    set(obj.hTaskPiezoLogger,'sampClkTimebaseSrc', obj.refClockTermInt);
                    set(obj.hTaskPiezoLogger,'sampClkTimebaseRate', 5000000);
                end
                obj.hTaskPiezoLogger.cfgDigEdgeStartTrig(obj.trgdVolClockFullTermInt);
                obj.hTaskPiezoLogger.set('startTrigRetriggerable',true);
                obj.piezo.totalNSamples = min(1000, ceil(obj.nVolumes))*(N_AI_piezo-1); %there will be a buffer error for very long acquisitions using 'INF' frames
                obj.hTaskPiezoLogger.cfgInputBuffer(obj.piezo.totalNSamples);
                obj.hTaskPiezoLogger.start();
                
                obj.powerWaveform = repmat(obj.calib.B.B2V(obj.power.*obj.powerWaveformRaw), obj.nPlanes,1); %Look up the correct pockels cell voltage for the user's selected power level
                obj.powerWaveform = circshift(obj.powerWaveform, -[round(obj.sampleRate*obj.calib.pockels.wireDelaySlow) 0]);
                N_AO_beam = length(obj.powerWaveform);
            else %Any code specific to 2D mode only here
                N_AO_beam = size(obj.powerWaveformRaw,1);
                obj.powerWaveform = obj.calib.B.B2V(obj.power.*obj.powerWaveformRaw); %Look up the correct pockels cell voltage for the user's selected power level
                obj.powerWaveform = circshift(obj.powerWaveform, -[round(obj.sampleRate*obj.calib.pockels.wireDelaySlow) 0]);
            end
            

            % configure SLM trigger out
            %if obj.do3D 
                obj.SLMCounter.COrate=1/(10*obj.framePeriod);
                %N_CO_SLM=floor(obj.volumePeriod*obj.SLMCounter.COrate);
                N_CO_SLM=floor(obj.acqTime*obj.SLMCounter.COrate);

                obj.hTaskSLMTriggerOut.channels(1).set('pulseTerm',obj.SLMoutputTerminal);
                obj.hTaskSLMTriggerOut.cfgImplicitTiming('DAQmx_Val_FiniteSamps', N_CO_SLM)
                    if obj.syncDaqsToFpga
                        set(obj.hTaskSLMTriggerOut,'sampClkTimebaseSrc', obj.refClockTermInt);
                        set(obj.hTaskSLMTriggerOut,'sampClkTimebaseRate', 5000000);
                    end
                %obj.hTaskSLMTriggerOut.cfgDigEdgeStartTrig(obj.volClockFullTermInt);   
                obj.hTaskSLMTriggerOut.cfgDigEdgeStartTrig(obj.trgdVolClockFullTermInt);
                obj.hTaskSLMTriggerOut.set('startTrigRetriggerable',true);
                obj.hTaskSLMTriggerOut.start();
           % end
            
            %configure auxilliary AI logging 
            if any(obj.logAIs)
                obj.auxAI.nSamplesPerVolume = floor(obj.volumePeriod*obj.auxAI.AIrate); %number of AI samples per channel
                obj.hTaskAILogger.cfgSampClkTiming(obj.auxAI.AIrate, 'DAQmx_Val_FiniteSamps', obj.auxAI.nSamplesPerVolume-1);
                if obj.syncDaqsToFpga
                    set(obj.hTaskAILogger,'sampClkTimebaseSrc', obj.refClockTermInt);
                    set(obj.hTaskAILogger,'sampClkTimebaseRate', 5000000);
                end
                obj.hTaskAILogger.cfgDigEdgeStartTrig(obj.trgdVolClockFullTermInt);
                obj.hTaskAILogger.set('startTrigRetriggerable',true);
                obj.auxAI.totalNSamples = min(1e7, ceil(obj.nVolumes))*(obj.auxAI.nSamplesPerVolume-1); %there will be a buffer error for very long acquisitions using 'INF' frames
                obj.hTaskAILogger.cfgInputBuffer(obj.auxAI.totalNSamples);
                obj.hTaskAILogger.start();
            end
            
            %write power waveform
            obj.hTaskPowerBeam.stop();
            obj.hTaskPowerBeam.cfgSampClkTiming(obj.sampleRate, 'DAQmx_Val_FiniteSamps', N_AO_beam-1);
            if obj.syncDaqsToFpga
                set(obj.hTaskPowerBeam,'sampClkTimebaseSrc', obj.refClockTermInt);
                set(obj.hTaskPowerBeam,'sampClkTimebaseRate', 5000000);
            end
            obj.hTaskPowerBeam.cfgOutputBuffer(N_AO_beam-1);
            obj.hTaskPowerBeam.cfgDigEdgeStartTrig(obj.volClockFullTermInt);
            obj.hTaskPowerBeam.set('startTrigRetriggerable',true);
            obj.hTaskPowerBeam.writeAnalogData(obj.powerWaveform(1:end-1,:));
            obj.hTaskPowerBeam.start();
            
            disp('Done!');
            
            %% flush and configure fifos
            fprintf('Initializing FPGA... ');
            
            obj.takeFpga();
            
            obj.dispCorrupt = false;
            obj.streamCorrupt = false;
            
            obj.hFpga.Enable = false;
            obj.flushFifos();
            obj.hFpga.fifo_GalvoData.configure(10000000);
            obj.hFpga.fifo_PmtData.configure(125000000);
            obj.hFpga.fifo_ShortPmtData.configure(50000000);
            obj.hFpga.fifo_FrameTagData.configure(10000);
            
            obj.hFpga.fifo_GalvoData.start;
            obj.hFpga.fifo_PmtData.start;
            obj.hFpga.fifo_ShortPmtData.start;
            obj.hFpga.fifo_FrameTagData.start;
            obj.flushFifos();
            
            %% FPGA parameters
            obj.hFpga.EnableRawDataStream = obj.logData;
            obj.hFpga.EnableGalvoDataStream = obj.logData;
            obj.hFpga.SamplesPerLaserClk = ceil(obj.samplesPerLaserClk/2);
            obj.hFpga.SamplesPerLaserClk = ceil(obj.samplesPerLaserClk/2);
            obj.hFpga.FrameClockPeriod = 0;
            
            if isinf(obj.framesToCollect)
                obj.hFpga.FramesToCollect = 2^32-1;
            else
                obj.hFpga.FramesToCollect = obj.framesToCollect;
            end
            
            obj.hFpga.PXITrigLaserClk = 4;
            
            if isempty(obj.startTrigTerm)
                obj.hFpga.PXITrigExtStartTrig = 0;
            else
                obj.hRouteRegistry.connectTerms(['/' obj.scannerDaqName '/PFI' num2str(obj.startTrigTerm)], ['/' obj.scannerDaqName '/PXI_Trig3']);
                obj.hFpga.PXITrigExtStartTrig = 4;
            end
            
            obj.configOscope();
            
            disp('Done!');
            
            obj.status = 'Initialized';
        end
        
        % start actuators and arm FPGA for acquisition
        function armAcq(obj)
            % mex parameters
            if obj.logData
                obj.startNewFile();
            end
            
            % prepare acq start pulse
            obj.hTaskShutter.writeDigitalData(false, 1, true); %close shutter
            obj.hTaskTriggerOutOnStart.abort();
            obj.hTaskTriggerOutOnStart.start();
            
            % start actuators
            obj.hFpga.NumSlices = obj.nPlanes;
            obj.hFpga.FrameClockPeriod = floor(obj.framePeriod * 125e6);
            
            % reset timeout counters
            obj.hFpga.Reset = true;
            
            % start polling for FPGA data
            start(obj.hTimerPollData);
            
            % enable fpga for triggered acquisition
            obj.hFpga.Enable = true;
        end
        
        % trigger start of data acquisition
        function triggerAcq(obj)
            obj.hTaskShutter.writeDigitalData(true, 1, true); %open shutter
            if obj.logData %save public properties as metadata
                fields = fieldnames(obj);
                for f = fields'
                    savedata.(f{1}) = obj.(f{1});
                end
                if isfield(savedata.PSF, 'calib')
                    savedata.PSF = rmfield(savedata.PSF ,'calib');
                end
                %include the SLM properties in metadata
                if ~isempty(obj.hSI) 
                    if ~isempty(obj.hSI.hSlmScan.testPattern)
                        savedata.SLM.pattern = obj.hSI.hSlmScan.testPattern;
                        savedata.SLM.pixelToRefTransform = obj.hSI.hSlmScan.testPatternPixelToRefTransform; % #ok<STRNU>
                    end
                end
                %if ~isempty(obj.hSI) || (obj.hSLM.queueLength>0) || obj.do3D %check that 3D patterns are loaded
                    %savedata.SLM.frames=obj.hSLM.queuedFrames;
                    %savedata.SLM.pixelToRefTransform = obj.hSI.hSlmScan.testPatternPixelToRefTransform;
                %end
                %include time of acquisition
                savedata.timeNow = now;
                %save metadata
                save(obj.metaDataFileName, '-struct', 'savedata', '-v6') %-v6 flag makes the save much faster
            end
            obj.wasRunning = true;
            obj.status = 'Trigger Received'; %set here because it was creating a slowdown in the 'pollFPGAshortdata' loop
            obj.hFpga.AcqSoftStartTrig = true;
        end
        
        % abort data acquisition
        function abortAcq(obj)
            obj.hTaskShutter.writeDigitalData(false, 1, true); %close shutter
            
            %we now keep the beam on but the shutter closed until
            %acquisition ends
%             obj.hTaskPowerBeam.abort();
%             obj.hTaskPowerBeam.control('DAQmx_Val_Task_Unreserve');
%             obj.hTaskPowerBeamOnDmd.writeAnalogData(obj.calib.pockels.slow(1)); %turn off beam
%             obj.hTaskPowerBeam.writeAnalogData(obj.powerWaveform(1:end-1,:)); %write the buffer and restart the task
%             obj.hTaskPowerBeam.start();                            %in order to be ready for next acq
            
            obj.hTaskGalvoBeams.abort();
            obj.hTaskGalvoBeams.control('DAQmx_Val_Task_Unreserve');
            obj.hTaskGalvo567OnDmd.writeAnalogData([-2 -2 0]); %park galvos
            obj.hTaskGalvoBeams.writeAnalogData([obj.galvoWaveform(1:end-1,:) obj.beamWaveform(1:end-1,:)]);
            obj.hTaskGalvoBeams.start(); %write the buffer and restart the task in order to be ready for next acq
            
            % prepare acq start pulse
            obj.hTaskTriggerOutOnStart.abort();
            obj.hTaskTriggerOutOnStart.start();
            
            obj.hFpga.Abort = true;
            obj.hFpga.Reset = true;
            
            sampsread = 0;
            if obj.do3D
                sampsavail = obj.hTaskPiezoLogger.get('readAvailSampPerChan');
                if sampsavail
                    [obj.piezo.Zdata, sampsread] = obj.hTaskPiezoLogger.readAnalogData(obj.piezo.totalNSamples,[],0.5); %last arg is timeout in case aborted early
                end
            end
            
            if obj.logData
                %save piezo data and restart the piezo task
                if sampsread
                    if exist(obj.ZDataFileName, 'file')
                        disp('PIEZO DATA WRITE OVERWROTE A PREVIOUS FILE!!')
                    end
                    Zdata = obj.piezo; %#ok<NASGU>
                    save(obj.ZDataFileName, 'Zdata', '-v6');
                end
                
                AIsampsread = 0;
                if any(obj.logAIs)
                    sampsavail = obj.hTaskAILogger.get('readAvailSampPerChan');
                    AIdata.data = [];
                    AIdata.chans = find(obj.logAIs);
                    if sampsavail
                        [AIdata.data, AIsampsread] = obj.hTaskAILogger.readAnalogData(obj.auxAI.totalNSamples,[],0.5); %#ok<ASGLU> %last arg is timeout in case aborted early
                    end
                    if AIsampsread
                        if exist(obj.aiDataFileName, 'file')
                            warning('AI DATA WRITE OVERWROTE A PREVIOUS FILE!!')
                        end
                        save(obj.aiDataFileName, 'AIdata', '-v6');
                    end
                end
                
                
                if strcmp(obj.status, 'Finishing data logging...')
                    % stop logging of the current file immediately and
                    % prepare to start acq again
                    obj.finishFileCB(true);
                else
                    reenable = obj.hFpga.Enable;
                    obj.hFpga.Enable = false;
                    
                    if obj.stopGalvosWhileLogging
                        obj.hFpga.FrameClockPeriod = 0;
                    end
                    
                    obj.status = 'Finishing data logging...';
                    err = AcquisitionMex('finish', @()obj.finishFileCB(reenable));
                    if err
                        warning('Failed to send finish command to mex.');
                        obj.finishFileCB(reenable);
                    end
                end
            else
                if obj.retriggerable
                    obj.status = 'Armed';
                else
                    obj.hFpga.FrameClockPeriod = 0;
                    stop(obj.hTimerPollData);
                    obj.hFpga.Enable = false;
                    obj.status = 'Initialized';
                end
            end
        end
        
        function finishFileCB(obj,reenable)
            if obj.retriggerable
                obj.startNewFile();
                obj.hFpga.FrameClockPeriod = floor(obj.framePeriod * 125e6);
                obj.hFpga.Enable = reenable;
                obj.status = 'Armed';
            else
                stop(obj.hTimerPollData);
                obj.hFpga.Enable = false;
                obj.status = 'Initialized';
            end
        end
        
        % stop actuators and disarm FPGA
        function endAcq(obj,keepFpga)
            if nargin < 2
                keepFpga = false;
            end
            
            if most.idioms.isValidObj(obj.hTimerPollData)
                stop(obj.hTimerPollData);
            end
            obj.wasRunning = false;
            
            % stop frame clock
            obj.hFpga.FrameClockPeriod = 0;
            
            if most.idioms.isValidObj(obj.hTaskPowerBeam)
                obj.hTaskPowerBeam.abort();
                obj.hTaskPowerBeam.control('DAQmx_Val_Task_Unreserve');
                delete(obj.hTaskPowerBeam);
            end
            if most.idioms.isValidObj(obj.hTaskPowerBeamOnDmd)
                obj.hTaskPowerBeamOnDmd.writeAnalogData(obj.calib.pockels.slow(1)); %turn beam off
                delete(obj.hTaskPowerBeamOnDmd);
            end
            if most.idioms.isValidObj(obj.hTaskShutter)
                obj.hTaskShutter.writeDigitalData(false, 1, true); %close shutter
                delete(obj.hTaskShutter);
            end
            if most.idioms.isValidObj(obj.hTaskGalvoBeams)
                obj.hTaskGalvoBeams.abort();
                obj.hTaskGalvoBeams.control('DAQmx_Val_Task_Unreserve');
                delete(obj.hTaskGalvoBeams);
            end
            if most.idioms.isValidObj(obj.hTaskGalvo567OnDmd)
                obj.hTaskGalvo567OnDmd.writeAnalogData([-2 -2 0]); %park galvos
                delete(obj.hTaskGalvo567OnDmd);
            end
            if most.idioms.isValidObj(obj.hTaskPiezoFeedback)
                stop(obj.hTimerPiezo);
                obj.hTaskPiezoFeedback.abort();
                most.idioms.safeDeleteObj(obj.hTaskPiezoFeedback);
            end
            if most.idioms.isValidObj(obj.hTaskSLMTriggerOut)
                obj.hTaskSLMTriggerOut.abort();
                most.idioms.safeDeleteObj(obj.hTaskSLMTriggerOut);
            end
            
            % end SLM triggered series
            if obj.do3D
                %obj.hSlm.abortSlmQueue();
            end
            
            most.idioms.safeDeleteObj(obj.hTaskPiezo);
            most.idioms.safeDeleteObj(obj.hTaskPiezoLogger);
            most.idioms.safeDeleteObj(obj.hTaskAILogger);
            most.idioms.safeDeleteObj(obj.hRouteRegistry);
            most.idioms.safeDeleteObj(obj.hTaskFlipMirror);
            most.idioms.safeDeleteObj(obj.hTaskTriggerOutOnStart);
            
            err = AcquisitionMex('stop');
            
            if most.idioms.isValidObj(obj.hFpga)
                obj.hFpga.Enable = false;
                obj.hFpga.Reset = true;
                most.idioms.pauseTight(0.05);
            end
            
            if ~keepFpga
%                 obj.restoreSiFpga();
            end
            
            % status
            obj.status = 'Idle';
        end
        
        function pollFpgaShortData(obj,~,~)
            p = obj.hFpga.PmtDataFifoTimeouts;
            g = obj.hFpga.GalvoDataFifoTimeouts;
            if ~obj.streamCorrupt && ((p > 0) || (g > 0))
                obj.streamCorrupt = true;
                fprintf(2, 'Logged data is corrupted due to FIFO timeouts (pmt:%d, galvo:%d). Some data will be missing.\n',p,g);
            end
            
            [frameTags, frames] = obj.readAllAvailData();
            if ~isempty(frameTags)
                obj.latestFrameTag = frameTags(end);
                obj.latestFrame = frames{end};
                obj.wasRunning = true;
                if ~strcmp(obj.status, 'Acquiring Data')
                    obj.status = 'Acquiring Data';
                end
%                 fprintf('%d frames read. %.3fs\n', numel(frameTags), toc(t));
            else
                if strcmp(obj.hFpga.AcqStatus, 'Armed')
                    obj.wasRunning = true;
                    obj.status = 'Trigger Received';
                elseif strcmp(obj.hFpga.AcqStatus, 'Running')
                    obj.wasRunning = true;
                    obj.status = 'Acquiring Data';
                else
                    if obj.wasRunning
                        obj.wasRunning = false;
                        obj.abortAcq();
                    elseif obj.hFpga.Enable
                        obj.status = 'Armed';
                    end
                end
            end
        end
        
        function [frameTags, frames] = readAllAvailData(obj, readAll)
            frameTagData = obj.hFpga.fifo_FrameTagData.readAll();
            
            if isempty(frameTagData)
                frameTags =[];
                frames = [];
                return
            end
            frameTags = obj.decodeFrameTagData(frameTagData);
            
%             tic
%             frameTags = repmat(obj.decodeFrameTagData(0),numel(frameTagData),1);
%             for i = 1:numel(frameTagData)
%                 frameTags(i) = obj.decodeFrameTagData(frameTagData(i));
%             end
%             toc
            
            if (obj.hFpga.FrameTagFifoTimeouts > 0) || (obj.hFpga.ShortPmtDataFifoTimeouts > 0) || obj.dispCorrupt
                frames{1} = obj.hFpga.fifo_ShortPmtData.readAll();
                if ~obj.dispCorrupt
                    fprintf(2, 'Fifo timeouts detected. Live display data is corrupt. Datalogging will continue.\n');
                    obj.dispCorrupt = true;
                elseif (obj.hFpga.FrameTagFifoTimeouts == 0) && (obj.hFpga.ShortPmtDataFifoTimeouts == 0) && (numel(frameTags) == 0)
                    %recovered from corruption
                    obj.dispCorrupt = false;
                end
            else
                %discard all but last frame
                if nargin<2 || ~readAll
                    try
                        TS = [frameTags.totalSamples];
                        if numel(TS) > 1
                            obj.hFpga.fifo_ShortPmtData.read(sum(TS(1:end-1)),0);
                        end
                        frames = {obj.hFpga.fifo_ShortPmtData.read(TS(end))};
                    catch
                        disp('Frametags:')
                        disp(TS)
                        obj.abortAcq();
                        error('Error reading live acquisition data. Laser clock not connected?');
                    end
                else %read all frames
                    frames = cell(length(frameTags),1);
                    for frameIX = 1:length(frames)
                        frames{frameIX} = obj.hFpga.fifo_ShortPmtData.read(frameTags(frameIX).totalSamples,0);
                    end
                end
            end
        end
        
        function ft = decodeFrameTagData(obj, data)
            ft = struct('frameNumber', num2cell(bitand(data,uint64(1048575))),...
                'laserClcks', num2cell(bitshift(bitand(data,uint64(549754765312)),-20)),...
                'totalSamples', num2cell(bitshift(bitand(data,uint64(18446743523953737728)),-39)),...
                'sampsPerClk', obj.samplesPerLaserClk);
        end
        
        function flushFifos(obj)
            obj.hFpga.fifo_GalvoData.readAll();
            obj.hFpga.fifo_PmtData.readAll();
            obj.hFpga.fifo_ShortPmtData.readAll();
            obj.hFpga.fifo_FrameTagData.readAll();
        end
        
        function configOscope(obj)
            inputRg = repmat(obj.inputRange,1,4);
            inputRg(5:8) = 5;
            for ch = 0:7
                err = dabs.ni.oscope.configureChannel(ch, inputRg(ch+1), true, false);
                assert(err == 0, 'Error when attempting to configure NI 517x device. Code = %d', err);
            end
        end
        
        function incFileName(obj)
            obj.fileNameApp = obj.fileNameApp + 1;
        end
        
        function startNewFile(obj)
            while exist(obj.metaFileName,'file')
                obj.incFileName();
            end
            
            % put header string together
            obj.metaHeaderString = sprintf('headerVersion=0.1\n');
            obj.metaHeaderString = sprintf('%sinputRange=%s\n', obj.metaHeaderString, num2str(obj.inputRange));
            obj.metaHeaderString = sprintf('%sframesToCollect=%d\n', obj.metaHeaderString, obj.framesToCollect);
            obj.metaHeaderString = sprintf('%sframePeriod=%s\n', obj.metaHeaderString, num2str(obj.framePeriod));
            obj.metaHeaderString = sprintf('%spower=%s\n', obj.metaHeaderString, num2str(obj.power));
            
            err = AcquisitionMex('start', obj);
            assert(err == 0, 'Error initializing datalogger. Error %d.', err);
        end
        
        function [data, freq] = acqGalvoPosSamps(obj,N,freq)
            obj.takeFpga();
            obj.configOscope();
            obj.hFpga.GalvoRecordArm = true;
            if ~strcmp(obj.hFpga.AcqStatus, 'Galvo Record Armed');
                obj.hFpga.GalvoRecordArm = false;
                error('Could not enter galvo record mode. Ensure acquisition is inactive.');
            end
            
            decim = floor(125000000/freq);
            obj.hFpga.GalvoRecordDecimation = decim;
            freq = 125000000/decim;
            
            obj.hFpga.GalvoRecordSamples = N;
            
            %flush fifo
            obj.flushFifos();
            obj.hFpga.fifo_GalvoData.configure(10000000);
            obj.hFpga.fifo_GalvoData.start;
            
            obj.hFpga.AcqSoftStartTrig = true;
            most.idioms.pauseTight(N/freq*1.1);
            galvoData = obj.hFpga.fifo_GalvoData.readAll();
            obj.hFpga.GalvoRecordArm = false;
%             obj.restoreSiFpga();
            
            assert(length(galvoData) == N*2, 'Incorrect number of samples received.');
            
            %first element is four 16 bit samples, one from each galvo1-4
            dat1 = galvoData(1:2:end);
            posDat = double(typecast(dat1, 'int16'));
            
            %second element is two 14 bit samples, one from each galvo5-6, and a 36 bit clock that we are not interested in here
            dat2 = galvoData(2:2:end);
            g5 = typecast(uint16(bitand(bitshift(dat2,2),65532)),'int16');
            g6 = typecast(uint16(bitand(bitshift(dat2,-12),65532)),'int16');
            posDat2 = double([g5 g6]);
            
            data = [posDat(1:4:end) posDat(2:4:end) posDat(3:4:end) posDat(4:4:end) posDat2] * 2.5 / 2^15;
        end
        
        function freq = acqGalvoPosSampsAsync(obj,N,freq,extTrig)
            obj.takeFpga();
            obj.configOscope();
            obj.hFpga.PXITrigExtStartTrig = 0;
            obj.hFpga.GalvoRecordArm = true;
            if ~strcmp(obj.hFpga.AcqStatus, 'Galvo Record Armed');
                obj.hFpga.GalvoRecordArm = false;
                error('Could not enter galvo record mode. Ensure acquisition is inactive.');
            end
            
            decim = floor(125000000/freq);
            obj.hFpga.GalvoRecordDecimation = decim;
            freq = 125000000/decim;
            
            obj.hFpga.GalvoRecordSamples = N;
            
            %flush fifo
            obj.flushFifos();
            obj.hFpga.fifo_GalvoData.configure(125000000);
            obj.hFpga.fifo_GalvoData.start;
            
            if isempty(extTrig)
                obj.hFpga.AcqSoftStartTrig = true;
            else
                try
                    obj.hRouteRegistry = dabs.ni.daqmx.util.triggerRouteRegistry();
                    obj.hRouteRegistry.connectTerms(extTrig, ['/' obj.scannerDaqName '/PXI_Trig3']);
                    obj.hFpga.PXITrigExtStartTrig = 4;
                catch
                    fprintf(2, 'Failed to connect route for trigger.');
                    obj.hFpga.GalvoRecordArm = false;
                end
            end
        end
        
        function done = pollGalvoAsyncAcq(obj)
            done = obj.hFpga.GalvoSamplesRemaining == 0;
        end
        
        function data = endGalvoAsyncAcq(obj)
            galvoData = obj.hFpga.fifo_GalvoData.readAll();
            obj.hFpga.GalvoRecordArm = false;
            
            most.idioms.safeDeleteObj(obj.hRouteRegistry);
            nRS = obj.hFpga.GalvoRecordSamples;
%             obj.restoreSiFpga();
            
            assert(length(galvoData) == nRS*2, 'Incorrect number of samples received.');
            
            %first element is four 16 bit samples, one from each galvo1-4
            dat1 = galvoData(1:2:end);
            posDat = double(typecast(dat1, 'int16'));
            
            %second element is two 14 bit samples, one from each galvo5-6, and a 36 bit clock that we are not interested in here
            dat2 = galvoData(2:2:end);
            g5 = typecast(uint16(bitand(bitshift(dat2,2),65532)),'int16');
            g6 = typecast(uint16(bitand(bitshift(dat2,-12),65532)),'int16');
            posDat2 = double([g5 g6]);
            
            data = [posDat(1:4:end) posDat(2:4:end) posDat(3:4:end) posDat(4:4:end) posDat2] * 2.5 / 2^15;
        end
    end
    
    %% Prop access
    methods
        function v = get.pmtDataFileName(obj)
            v = [obj.dataDir filesep obj.user filesep obj.fileNameStem num2str(obj.fileNameApp) '.pdat'];
        end
        
        function v = get.galvoDataFileName(obj)
            v = [obj.dataDir filesep obj.user filesep obj.fileNameStem num2str(obj.fileNameApp) '.gdat'];
        end
        
        function v = get.ZDataFileName(obj)
            v = [obj.dataDir filesep obj.user filesep obj.fileNameStem num2str(obj.fileNameApp) '.zdat'];
        end
        
        function v = get.aiDataFileName(obj)
            v = [obj.dataDir filesep obj.user filesep obj.fileNameStem num2str(obj.fileNameApp) '.adat'];
        end
        
        function v = get.metaDataFileName(obj)
            v = [obj.dataDir filesep obj.user filesep obj.fileNameStem num2str(obj.fileNameApp) '.mdat'];
        end
        
        function v = get.metaFileName(obj)
            v = [obj.dataDir filesep obj.user filesep obj.fileNameStem num2str(obj.fileNameApp) '.txt'];
        end
        
        function v = get.fpgaSessionId(obj)
            v = obj.hFpga.session;
        end
        
        function v = get.pmtFifoNum(obj)
            v = obj.hFpga.fifo_PmtData.fifoNumber;
        end
        
        function v = get.galvoFifoNum(obj)
            v = obj.hFpga.fifo_GalvoData.fifoNumber;
        end
        
        function set.laserClkDelay(obj,v)
            n = floor(v / 8e-9);
            n(n<0) = 0;
            n(n>149) = 149;
            
            if n*8e-9 ~= v
                most.idioms.warn('Value coerced to acceptable value of %dns.', n*8);
            end
            
            obj.hFpga.LaserClkDelay = n;
            obj.laserClkDelay = n*8e-9;
        end
    end
end

