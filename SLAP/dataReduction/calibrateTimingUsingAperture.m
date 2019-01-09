function calibrateTimingUsingAperture(hSLAPMi)
%Takes a SLAPMi image of a bright, uniform sample and compares to the known
%location of the aperture obtained with the camera in measurePSF.m
%Calculates the delay that lines up the edges of the aperture to the known locations on each axis 
%if using a fluorescent block, set the HPD detector to ~250V.

%switch directory
hSLAPMi.user = 'Timing'; TimingDir = [hSLAPMi.dataDir filesep hSLAPMi.user];
if ~exist(TimingDir, 'dir')
    mkdir(TimingDir);
end
cd(TimingDir);

%settings:
hSLAPMi.aperture = true;
hSLAPMi.framesToCollect = 50;
hSLAPMi.logData = true;

%set the scan length slightly longer to overlap the aperture more
linelength_old = hSLAPMi.PSF.aperture.linelength;
hSLAPMi.PSF.aperture.linelength = 1.2*linelength_old;

%calibrate the galvos at this new linelength
waitfor(SLAPmi_cal_Galvos(hSLAPMi))
%resp = msgbox('Hit OK once the galvos are finished calibrating');

hSLAPMi.initAcq()
hSLAPMi.armAcq();

try
    hSLAPMi.triggerAcq();
    pause(0.5) %allow acquisition to complete
    hSLAPMi.endAcq();
    hSLAPMi.PSF.aperture.linelength = linelength_old;
catch ME
    hSLAPMi.abortAcq();
    hSLAPMi.endAcq();
    hSLAPMi.PSF.aperture.linelength = linelength_old;
    
    disp('An error occurred:')
    disp(getReport(ME))
    
    keyboard
end

%reduce the data
scandata = SLAPMi_reduce;

%make a voltage space image of the aperture
data = hSLAPMi.PSF;

s = regionprops(data.aperture.bw(:,:,1), 'centroid');
apCenter_P = fliplr(s.Centroid);
apCenter_V = [data.p2vx(apCenter_P) data.p2vy(apCenter_P)];

coords.X = linspace(apCenter_V(1)-data.aperture.linelength, apCenter_V(1)+data.aperture.linelength, 500); 
coords.Y = linspace(apCenter_V(2)-data.aperture.linelength, apCenter_V(2)+data.aperture.linelength, 500);
coords.Z = 0;

[Xmesh, Ymesh] = ndgrid(coords.X, coords.Y);
Xgrid = reshape(data.v2px(Xmesh(:),Ymesh(:)), size(Xmesh));
Ygrid = reshape(data.v2py(Xmesh(:),Ymesh(:)), size(Xmesh));

for lineset = 1:2
    GI = griddedInterpolant(double(data.aperture.bw(:,:,lineset)));
    apIM{lineset} = reshape(GI(Xgrid(:), Ygrid(:)), size(Xgrid)); %#ok<AGROW>
end

scandata.refIMcoords = coords;
P = linePSF_full(scandata);

y = mean([scandata.frames.pmtData],2);

offsets = nan(1,4);
for line = 1:4
    line_ixs = scandata.line==line;
    %project the correct aperture image with the subset of the projection
    %for that line
    IM = apIM{ceil(line/2)};
    yE = P.P(line_ixs,:)*IM(:); yE = yE./max(yE);
    yL = y(line_ixs); yL = yL./max(yL); yL(isnan(yL)) = 0;
    
    yL = max(0,yL-0.05);
    yL(yL>0.5) = 0.5; yE(yE>0.5) = 0.5;
    %cross-correlate expected with actual
    yEcenter = sum(yE.*(1:length(yE))')./sum(yE);
    yCenter =  sum(yL.*(1:length(yL))')./sum(yL);
    %[xc, lags] = xcorr(yE, yL);
    %[maxcorr, lag] = max(xc);
    offsets(line) =  yEcenter- yCenter; %lags(lag);
    
    figure('Name', ['line' int2str(line)]), plot(yL), hold on, plot(yE), plot(interp1(yE, (1:length(yE))+offsets(line)))
end
offsets

keyboard


    
