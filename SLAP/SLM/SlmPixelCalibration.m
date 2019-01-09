function SlmPixelCalibration(src,evt)
    persistent map frames frameGridRefSpaceXX frameGridRefSpaceYY chanIdx
    
    channel = 1;            % SI channel number used for acquisition
    initialDelay = 2;       % initial delay in frames to allow slm to reach thermal equilibrium
    slmResolution = 512;    % pixel resolution of SLM
    slmPixelLevels = 256;   % number of discrete steps the SLM pixels can be set to
    
    hSI = src.hSI;
    switch evt.EventName
        case 'focusStart'
            %%% setup
            fprintf('Initializing pixel Calibration\n');
            
            map = zeros(hSI.hRoiManager.linesPerFrame,hSI.hRoiManager.pixelsPerLine,slmPixelLevels,'int16');
            
            averageFactor = hSI.hDisplay.displayRollingAverageFactor;
            discardNFrames = 1;
            frames = repmat(averageFactor+discardNFrames,1,slmPixelLevels);
            frames = cumsum(frames);
            frames = frames + initialDelay;            
            [frameGridRefSpaceXX,frameGridRefSpaceYY] = hSI.hRoiManager.currentRoiGroup.rois(1).scanfields(1).meshgrid;
            
            [tf,chanIdx] = ismember(channel,hSI.hChannels.channelDisplay);
            
            if ~tf
                hSI.abort()
                error('Channel %d needs to be selected for display.',channel);
            end
            
            fprintf('Initialization done. Acquiring %d frames\n',frames(end));
            fprintf('Slm Pixel Value: %d\n',0);
            hSI.hSlmScan.writeTestPattern(zeros(slmResolution));
        case 'frameAcquired'
            [tf,idx] = ismember(hSI.hDisplay.lastFrameNumber,frames);
            
            if tf
                map(:,:,idx) = hSI.hDisplay.lastFrame{chanIdx};
                if idx < length(frames)
                    fprintf('Slm Pixel Value: %d\n',idx);
                    hSI.hSlmScan.writeTestPattern(ones(slmResolution)*idx);
                else
                    hSI.abort();
                    
                    [slmGridXX,slmGridYY] = meshgrid(1:slmResolution,1:slmResolution);
                    [slmGridRefSpaceXX,slmGridRefSpaceYY] = scanimage.mroi.util.xformMesh(slmGridXX,slmGridYY,hSI.hSlmScan.testPatternPixelToRefTransform);
                    
                    mapSlm = zeros(slmResolution,slmResolution,slmPixelLevels);
                    for idx = 1:size(map,3)
                        mapSlm(:,:,idx) = interp2(frameGridRefSpaceXX,frameGridRefSpaceYY,double(map(:,:,idx)),slmGridRefSpaceXX,slmGridRefSpaceYY,'linear',NaN);
                    end
                    
                    assignin('base','map',map);
                    assignin('base','mapSlm',mapSlm);
                    fprintf('SLM calibration: Acquisition completed. Check base workspace for resulting maps\n');
             
                    disp('fitting LUT...')
                    lut = SLM_lut(mapSlm);
                    assignin('base','lut',lut);
                end
            end
        otherwise
            %No-op
    end
end
