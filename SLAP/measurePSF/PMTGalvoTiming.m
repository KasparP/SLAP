function PMTGalvoTiming (hSI, hSLAPMi, IM)
    [xGalvoPosByPixel, measurements] = SLAPMi.analyzeSIWaveform(hSI,hSLAPMi);
    galvoDelay(measurements, IM);
end
