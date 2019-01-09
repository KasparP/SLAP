function opts = calibrateTimingUsingAperture2(hSLAPMi)
%Oscillates the beam at high frequency at the edge of the aperture to
%figure out the delay between the galvo feedback and the PMT signal

%switch directory
hSLAPMi.user = 'Timing2'; TimingDir = [hSLAPMi.dataDir filesep hSLAPMi.user];
if ~exist(TimingDir, 'dir')
    mkdir(TimingDir);
end
cd(TimingDir);

%settings:
hSLAPMi.aperture = true;
hSLAPMi.framesToCollect = 50;
hSLAPMi.logData = true;
hSLAPMi.galvos.calFramePeriod = hSLAPMi.framePeriod;
hSLAPMi.galvos.calPos = [hSLAPMi.scanCenterPoint hSLAPMi.tiling];

freq =  1/hSLAPMi.framePeriod;
AOrate = hSLAPMi.galvos.AOrate;
NSampsOut = 4*round(AOrate/(4*freq));

XX = linspace(0, 4*pi, NSampsOut+1);
linelength = hSLAPMi.PSF.aperture.linelength;
waveform = linelength/6 * (sin(XX(1:end-1))-0.5);

calib = hSLAPMi.calib;
aperture = hSLAPMi.PSF.aperture;
offset = hSLAPMi.calib.galvos.offset;

for line = 1:4
    lname = ['line' int2str(line)];
    
    [Xv,Yv] = hSLAPMi.PSF.T2V(0.71*linelength, 0,aperture.(lname).offsetX, aperture.(lname).offsetY, hSLAPMi.PSF.(lname).g_angle);
    
    AO = zeros(NSampsOut, 6);
    switch line
        case 1
            AO(:,1) = Xv + waveform;
            AO(:,2) = Yv;
            AO(:,5) = offset.E(1);
            hSLAPMi.beamWaveform = calib.pockels.fast(2)*ones(size(AO,1),1); %shift the AO to compensate for wire transmission delay
        case 2
            AO(:,1) = Xv;
            AO(:,2) = Yv + waveform;
            AO(:,5) = offset.E(2);
            hSLAPMi.beamWaveform = calib.pockels.fast(2)*ones(size(AO,1),1); %shift the AO to compensate for wire transmission delay
        case 3
            AO(:,3) = Xv + waveform;
            AO(:,4) = Yv;
            AO(:,6) = offset.F(1);
            hSLAPMi.beamWaveform = calib.pockels.fast(1)*ones(size(AO,1),1); %shift the AO to compensate for wire transmission delay
        case 4
            AO(:,3) = Xv;
            AO(:,4) = Yv + waveform;
            AO(:,6) = offset.F(2);
            hSLAPMi.beamWaveform = calib.pockels.fast(1)*ones(size(AO,1),1); %shift the AO to compensate for wire transmission delay
    end
    
    hSLAPMi.lineIDs = line*ones(1,size(AO,1));
    hSLAPMi.galvoWaveform = AO;
    hSLAPMi.powerWaveformRaw = hSLAPMi.calib.B.scaleby(line)*ones(size(AO,1),1);
    
    
    fns{line} = hSLAPMi.galvoDataFileName;
    hSLAPMi.initAcq()
    hSLAPMi.armAcq();
    pause(0.1);
    hSLAPMi.triggerAcq();
    pause(0.2) %allow acquisition to complete
    try
        hSLAPMi.endAcq();
    catch ME
        hSLAPMi.abortAcq();
        hSLAPMi.endAcq();
        disp(getReport(ME))
        keyboard
    end
end

%reduce the data


for line = 1:4
    offset = 1;
    while ~(offset==0);
        opts.path = fns{line};
        [offset, Gn, opts] = SLAPMi_reduce_galvosRaw(opts);
        opts.AIdelay(Gn) = opts.AIdelay(Gn)+(offset*(2e-7));
        disp(opts.AIdelay);
    end
end
