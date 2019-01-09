function linePSF_testframes(data)
%Displays the testframe images and overlays reconstructed PSF for a PSF measurement (data)

%scanoptions contains the imaging metadata: 
    %the xyz location of the illumination line at each laser pulse
%calibration data obtained from measurePSF

%ensure that the scan and the psf measurements were obtained with the same
%calibration?

%load the measurement files if not provided
if nargin<1 || isempty(data)
    data = [];
    [fn, dr] = uigetfile('E:\SLAPmidata\PSF\*.psf', 'Select valid PSF files', 'multiselect', 'on');
    load([dr fn], '-mat')
end

calib = data.calib;
scandata.calib = calib;
scandata.calib.galvos.linelength = 4.3;
linelength = scandata.calib.galvos.linelength;
dZ = size(data.testframe,2);

Zcenter = ceil(dZ/2);
selframes=  false(size(data.testframe));
selframes(:,1,:) = true;
selframes = find(selframes(:)); 

scandata.line = [data.testframe(selframes).line];
scandata.Vx = [data.testframe(selframes).Vx];
scandata.Vy = [data.testframe(selframes).Vy]; 
scandata.frames(1).Z = [data.testframe(selframes).Z];
scandata.metadata.PSF = data;
scandata.metadata.calib = calib;
scandata.refIMcoords.X = linspace(calib.galvos.offset.raster.X-linelength, calib.galvos.offset.raster.X+linelength, 1024); %in V
scandata.refIMcoords.Y = linspace(calib.galvos.offset.raster.Y-linelength, calib.galvos.offset.raster.Y+linelength, 1024); %in V
scandata.refIMcoords.Z = unique([data.testframe.Z]); %reconstruct a volume with the right number of planes

opts.ignoreHash = true;
P = linePSF_full(scandata, opts);

[Xgrid,Ygrid] = ndgrid(P.coords{1}, P.coords{2});
Px = data.v2px(Xgrid(:), Ygrid(:));
Py = data.v2py(Xgrid(:), Ygrid(:));

for ind = 1:length(selframes)
    [sub1,sub2,sub3] = ind2sub(size(data.testframe), selframes(ind));
    disp(['Displaying testframe stack: ' int2str(ind)])
    R = nan([size(Xgrid) dZ]);
    %R = nan([size(Xmesh) dZ]);
    for Z = 1:dZ
        GI = griddedInterpolant(data.testframe(sub1,Z,sub3).im);
        imV = reshape(GI(Px,Py), size(Xgrid));
        R(:,:,Z) = imV;
    end
    R = R./prctile(R(:), 99.95);
    G = reshape(full(P.P(ind,:))./max(P.P(ind,:)), [length(P.coords{1}) length(P.coords{2}) length(P.coords{3})]);
    G = G./max(G(:));
    B = 0*R;
    
    figure('name', ['TestFrame ' int2str(selframes(ind)) '. Line:' int2str(scandata.line(ind)), ' X:' num2str(scandata.Vx(ind), 3), ' Y:' num2str(scandata.Vy(ind), 3)]),
    imshow3D(cat(4, R,G,B));
    drawnow
end
end