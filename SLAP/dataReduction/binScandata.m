function scandata = binScandata(scandata, lxlbinfactor, framebinfactor)
%performs posthoc binning on a scandata file to make subsequent computations easier.

nCh = size(scandata.frames(1).pmtData,2);

avgIXs = [1:lxlbinfactor:length(scandata.line)] + [0:lxlbinfactor-1]';
scandata.Vx = mean(scandata.Vx(avgIXs),1)';
scandata.Vy = mean(scandata.Vy(avgIXs),1)';
scandata.line = scandata.line(1:lxlbinfactor:end);

for frame = 1:length(scandata.frames)
    for G = 6:-1:1
        tmpG = scandata.frames(frame).galvoData(:,G);
        galvoData(:,G) = mean(tmpG(avgIXs),1);
    end
    scandata.frames(frame).galvoData = galvoData;
    
    for Ch = nCh:-1:1
        tmpP = scandata.frames(frame).pmtData(:,Ch);
        pmtData(:,Ch) = mean(tmpP(avgIXs),1);
    end
    scandata.frames(frame).pmtData = pmtData;
    
    scandata.frames(frame).valid = all(scandata.frames(frame).valid(avgIXs),1);
end
    
    