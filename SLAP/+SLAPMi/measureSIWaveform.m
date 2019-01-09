function [fsOut,cmdWvfm,fsIn,respWvfm,lineScanPeriod,lineAcquisitionPeriod] = measureSIWaveform(hSI,hSLAPMi)

    assert(~hSI.active, 'Cannot run test during active acquisition.');

    hWb = waitbar(0,'Preparing Waveform and DAQs...');
    hAO = [];
%     hDO = [];
    try
        %% prepare waveform
        zs = hSI.hStackManager.zs;
        sf = hSI.hRoiManager.currentRoiGroup.rois(1).get(zs(1));
        ss = hSI.hScan2D.scannerset;

        [lineScanPeriod,lineAcquisitionPeriod] = ss.linePeriod(sf);
        nx = ss.nsamples(ss.scanners{1},lineScanPeriod);           % total number of scan samples per line

        [ao_volts,~,~] = hSI.hRoiManager.currentRoiGroup.scanStackAO(ss,0,zs(1),'',0);
        xWvfm = ao_volts.G(:,1);

        cmdWvfm = xWvfm(nx*2+1:nx*4);
        testWvfm = repmat(cmdWvfm,20,1);     % 20 lines (10 cycles) are run
        fsOut = hSI.hScan2D.sampleRateCtl;
        T = (length(testWvfm)+1)/fsOut;
        
        
        %% prepare trigger task
%         hDO = most.util.safeCreateTask([hSI.hScan2D.name ' Soft Pulse2']);
%         hDO.createDOChan(hSI.hScan2D.mdfData.deviceNameGalvo,'port1/line0');
        
            
        %% prepare output task
        hAO = most.util.safeCreateTask([hSI.hScan2D.name '-ScannerTestOut2']);
        hAO.createAOVoltageChan(hSI.hScan2D.mdfData.deviceNameGalvo, hSI.hScan2D.mdfData.XMirrorChannelID, 'XMirrorChannel');
        hAO.control('DAQmx_Val_Task_Unreserve'); %Flush any previous data in the buffer
%         hAO.cfgDigEdgeStartTrig('PFI0', 'DAQmx_Val_Rising');
        hAO.cfgSampClkTiming(fsOut, 'DAQmx_Val_FiniteSamps', length(testWvfm));
        hAO.cfgOutputBuffer(length(testWvfm));
        hAO.set('startTrigRetriggerable',false);
        if ~hSI.hScan2D.simulated
            hAO.writeAnalogData(testWvfm);
        end
        hAO.control('DAQmx_Val_Task_Verify'); %%% Verify Task Configuration (mostly for trigger routing
        

        %% prepare and start input task
%         fsIn = hSLAPMi.acqGalvoPosSampsAsync(T*fsOut,fsOut,['/' hSLAPMi.scannerDaqName '/PFI0']);
        fsIn = fsOut * 10;
        fsIn = hSLAPMi.acqGalvoPosSampsAsync(ceil(T*fsIn),fsIn,sprintf('/%s/ao/StartTrigger',hSI.hScan2D.mdfData.deviceNameGalvo));


        %% start output task
        hAO.start();

        %% Trigger and wait
%         hDO.writeDigitalData([0;1;0],0.5,true);
        pt = T/10;
        for i = 0:11
            waitbar(i/10, hWb, 'Executing test...');
            t = tic;
            if hSLAPMi.pollGalvoAsyncAcq()
                break;
            end
            most.idioms.pauseTight(pt-toc(t));
        end
        assert(hSLAPMi.pollGalvoAsyncAcq(), 'Input did not complete in expected time.');

        %% read data and stop tasks
        waitbar(1, hWb, 'Parsing data...');
        data = hSLAPMi.endGalvoAsyncAcq();
        hAO.abort();

        %% parse and scale data
        sN = ceil(lineScanPeriod*fsIn);
        nCycles = 5;
        firstCycle = 5;
        respWvfm = data(1+sN*(2*(firstCycle-1)):sN*(2*(firstCycle-1)+2*nCycles),:);    % lines 17 and 18 (cycle 9) are extracted to plot
        
        respWvfm = respWvfm(:,3);
        respWvfm = mean(reshape(respWvfm, [2*sN nCycles]), 2);
        
        delete(hWb)
        most.idioms.safeDeleteObj(hAO);
%         most.idioms.safeDeleteObj(hDO);
    catch ME
        delete(hWb);
        most.idioms.safeDeleteObj(hAO);
%         most.idioms.safeDeleteObj(hDO);
        ME.rethrow
    end

    function cfg = daqMxTermCfgString(str)
        if length(str) > 4
            str = str(1:4);
        end
        cfg = ['DAQmx_Val_' str];
    end

end

