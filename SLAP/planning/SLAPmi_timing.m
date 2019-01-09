function SLAPmi_timing (hSLAPmi)
slowrate = hSLAPmi.sampleRate; %the slowest clock we have to synchronize to
timePerLine = ceil((hSLAPmi.res/(hSLAPmi.laser.repRate/hSLAPmi.binFactor) + hSLAPmi.pockels.deadtime)*slowrate)/slowrate;
frameTime = 4*hSLAPmi.tiling*timePerLine;

%dummy check
if 1/frameTime>1400 %we are asking for too fast a galvo speed
    disp('Too high of a framerate was requested, binFactor was increased to compensate');
    hSLAPmi.binFactor = ceil(1/(1400*frameTime));
    timePerLine = ceil((hSLAPmi.res/(hSLAPmi.laser.repRate/hSLAPmi.binFactor) + hSLAPmi.pockels.deadtime)*slowrate)/slowrate;
    frameTime = 4*hSLAPmi.tiling*timePerLine;
end

hSLAPmi.framePeriod = frameTime;
hSLAPmi.acqTime = hSLAPmi.framesToCollect*hSLAPmi.framePeriod;

if hSLAPmi.do3D %volumePeriod is not a multiple of frameTime!
    hSLAPmi.nPlanes = round(hSLAPmi.volumePeriod/frameTime);
    hSLAPmi.volumePeriod = frameTime*hSLAPmi.nPlanes;
    if hSLAPmi.Zmax < hSLAPmi.Zmin
        hSLAPmi.Zmax = hSLAPmi.Zmin;
    end
    hSLAPmi.dZ = (hSLAPmi.Zmax - hSLAPmi.Zmin)/(0.5*hSLAPmi.nPlanes);     
else
    hSLAPmi.nPlanes = 1;
    hSLAPmi.volumePeriod = frameTime;
    hSLAPmi.dZ = [];
end
hSLAPmi.nVolumes = hSLAPmi.framesToCollect/hSLAPmi.nPlanes;


%deal with volume imaging?