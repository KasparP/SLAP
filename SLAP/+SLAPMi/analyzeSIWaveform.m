function [xGalvoPosByPixel, measurements] = analyzeSIWaveform(hSI,hSLAPMi)
    [fsOut,cmdWvfm,fsIn,respWvfm,T,Ta] = SLAPMi.measureSIWaveform(hSI, hSLAPMi);
    measurements.T = T; measurements.Ta = Ta;measurements.respWvfm = respWvfm; 
    measurements.fsIn = fsIn; measurements.cmdWvfm = cmdWvfm; measurements.fsOut = fsOut;
    measurements.SIlinePhase = hSI.hScan2D.linePhase; measurements.N = hSI.hRoiManager.pixelsPerLine;
    
    bth = [cmdWvfm;respWvfm];
    yrg = max(bth) - min(bth);
    
    No = length(cmdWvfm);
    outTs = linspace(1/fsOut, No/fsOut, No)';
    Ni = length(respWvfm);
    inTs = linspace(1/fsIn, Ni/fsIn, Ni)';
    xrg = outTs(end);
    
    outTs = [outTs-(outTs(end)); outTs; outTs+(outTs(end))];
    inTs = [inTs-(inTs(end)); inTs; inTs+(inTs(end))];
    cmdWvfm = repmat(cmdWvfm,3,1);
    respWvfm = repmat(respWvfm,3,1);
    
    Td=(T-Ta)/2;
    vertlines = [Td Ta+Td];
    vertlines = [vertlines vertlines+T];
    
    %% plot not accounting phase adjust
%     figure;
%     plot(outTs, cmdWvfm);
%     hold on;
%     grid on;
%     plot(inTs, respWvfm);
%     plot(inTs, respWvfm-cmdWvfm);
%     title('Real time plot')
%     legend({'Command' 'Response' 'Error'});
%     xlabel('Time (s)');
%     ylabel('Position Command/Feedback Voltage (V)');
%     xlim([-xrg*.05 xrg*1.05]);
%     ylim([min(bth)-yrg/20 max(bth)+yrg/20]);
%     
%     for x = vertlines
%         plot([x x], [-15 15], '--k');
%     end
    
    
    
    %% plot accounting phase adjust
    sampShift = ceil(hSI.hScan2D.linePhase * fsIn);
    phsCorrRespWvfm = circshift(respWvfm,-sampShift);
    
    figure;
    plot(outTs, cmdWvfm);
    hold on;
    grid on;
    plot(inTs, phsCorrRespWvfm);
%     plot(inTs, phsCorrRespWvfm-cmdWvfm);
    title('ScanImage Galvo Command and Feedback Plot (Phase Adjusted)')
    legend({'Command' 'Response'});% 'Error'});
    xlabel('Time (s)');
    ylabel('Position Command/Feedback Voltage (V)');
    xlim([-xrg*.05 xrg*1.05]);
    ylim([min(bth)-yrg/20 max(bth)+yrg/20]);
    
    for x = vertlines
        plot([x x], [-15 15], '--k');
    end
    
    
    %% forward line pixel to position plot
    N = hSI.hRoiManager.pixelsPerLine;
    pixelTime = Ta / N;
    pixelTimes = linspace(pixelTime/2,pixelTime/2 + pixelTime*(N-1),N) + Td;
    pixelTimesR = linspace(pixelTime/2,pixelTime/2 + pixelTime*(N-1),N) + T + Td;
    xGalvoPosByPixel = [interp1(inTs,phsCorrRespWvfm,pixelTimes); fliplr(interp1(inTs,phsCorrRespWvfm,pixelTimesR))];
    
    figure;
    plot(xGalvoPosByPixel(1,:));
    hold on;
    grid on;
    plot(xGalvoPosByPixel(2,:))
    xlabel('Pixel Number');
    ylabel('Galvo Position Feedback Voltage (V)');
    title('ScanImage Pixel to Galvo Position Plot')
    legend({'Even Lines' 'Odd Lines'});
end

