function calib = calculateZPSF_beads
%Measures the Z-X psf of the microscope for every line, and updates the
%calibration. Run 'measureZPSF_beads' first!

basedir = 'E:\SLAPMiData\';
%select the gdat files
[fns, dr] = uigetfile([basedir '\PSF\Zscans\*.gdat'], 'multiselect', 'on');
%sort the files
fns = sort_nat(fns);

%load the calibration file to save the PSF into
calib = [];
[calibfn, calibdr] = uigetfile([basedir '\Calibration\*.cal'], 'Select your calibration file');
load([calibdr filesep calibfn], '-mat');

%load the PSF file, this will be used to convert galvo voltages to R values
data = [];
[PSFfn, PSFdr] = uigetfile([basedir '\PSF\*.psf'], 'Select your PSF file');
load([PSFdr filesep PSFfn], '-mat');

opts.channels = 2;
opts= optionsGUI(opts);

opts.doPCA = false;
opts.doPlot = false;
%run SLAPMi_reduce on each of the files
for i = 1:length(fns)
    opts.path = [dr fns{i}(1:end-5)];
    [scandata, hFig] = SLAPMi_reduce(opts);
    frames = scandata.frames;
    
    pmtData = [frames.pmtData]; 
    pmtData = nanmean(pmtData,2);
    
    if i==1
        D = zeros(length(pmtData), length(fns));
        Vx = D; Vy = D;
    end
    ixs = 1:min(size(pmtData,1),  size(D,1));
    D(ixs,i) = pmtData(ixs, end);
    Vx(ixs,i) = scandata.Vx(ixs);
    Vy(ixs,i) = scandata.Vy(ixs);
    
    %drawnow  %show the figure from data reduction, then close it?
    close(hFig);
end

%measure pixel sizes
pixelsizeV = median(sqrt(diff(scandata.Vx).^2 + diff(scandata.Vy).^2)); %pixel size, in Volts
pixelsizeUM = pixelsizeV/calib.galvos.pixelsizeVperUM; %pixel size in microns; 1V=45.3um
Zstep = 0.2; %Zstep as set in measurePSF_beads

figure('name', 'Full XZ image'), imagesc(D');

%find local maxima in the smoothed image to get centroids
thresh = 0.9;
Dsmooth = D; Dsmooth(isnan(Dsmooth)) = 0;
Dsmooth = imgaussfilt(Dsmooth,1);
Dsmooth(Dsmooth<thresh) = 0;
maximaIM = imregionalmax(Dsmooth);

dZ = max(11, round(6/Zstep));
dR = max(11, round(1.7/pixelsizeUM));

[tmeshZ, tmeshR] = meshgrid((-dZ:dZ)./dZ, (-dR:dR)./dR);
thresh = (tmeshZ.^4 + tmeshR.^4)/2;
%for each centroid
%find the superresolution center, using COM
for line = 1:4
    lname = ['line' int2str(line)];
    line_ixs = scandata.line==line;
    
    %get the camera pixel size for this line (this will be negative if the beam is scanned in the negative R-axis)
    [T,Td] = data.V2T(Vx(line_ixs,:),Vy(line_ixs,:),calib.galvos.offset.(lname).X, calib.galvos.offset.(lname).Y, data.(lname).g_angle);
    R = data.(lname).interpRc_fromT(T,Td);
    pixelsizeP = median(diff(R(:, fix(end/2))));
    
    %Find the beads within the image
    D_line = D(line_ixs,:);
    [X,Z] = find(maximaIM(line_ixs,:));
    
    select = X>dR & X<size(D_line,1)-dR+1 & Z>5 & Z<size(D_line,2)-5+1;
    X = X(select); Z = Z(select);
    centered = zeros(2*dR+1,2*dZ+1, length(X));
    hCB = [];
    
    selected = true(1, length(X));
    hF = figure('name', lname, 'closeRequestFcn', @closereq,'units', 'normalized', 'pos', [0.0119    0.5933    0.9962    0.3150]);
    lineZ = nan(1, length(X));
    for c_ix = 1:length(X)
        x = X(c_ix); z = Z(c_ix);
        gi = griddedInterpolant(D_line);
        box =  gi({(x-dR:x+dR), (z-dZ:z+dZ)});
        imregion = false(size(box)); imregion(dR-5:dR+5, dZ-3:dZ+3) = true;
        [result] = regionprops(imregion, box, 'WeightedCentroid', 'MeanIntensity');
        x = result.WeightedCentroid(2); z = result.WeightedCentroid(1); I = result.MeanIntensity;
        
        lineZ(c_ix) = Z(c_ix) + z - dZ - 1;
        [Xq,Yq] = meshgrid(z-dZ:z+dZ, x-dR:x+dR);
        centered(:,:,c_ix) = interp2(box./I,Xq, Yq);
        
        hax = axes('parent', hF);
        imagesc(centered(:,:,c_ix), [0 3]); 
        set(hax, 'pos', [(c_ix-1)*(1/length(X)) 0.15 (1/length(X)) 0.8], 'xtick', [], 'ytick', [])
        %xlabel(int2str(c_ix))
        
        hCB(c_ix) = uicontrol(...
                'Units','Normalized',...
                'Parent',hF,...
                'Position', [(c_ix-0.7)*(1/length(X)) 0.05 (1/length(X)) 0.1],...
                'HorizontalAlignment','left',...
                'String',int2str(c_ix),...
                'Value', true,...
                'Style','checkbox'); %#ok<AGROW,NASGU>
    end

    waitfor(hF)

    PSF_im = trimmean(centered(:,:,selected),40,3);
    PSF_im = PSF_im - min(PSF_im(:));
    PSF_im(isnan(PSF_im)) = 0;
    %PSF_im = max(0, PSF_im - nanmedian(PSF_im(:)));
    %figure, imagesc(PSF_im)
    
    figure('Name', lname), 
    subplot(2,2,1),  imagesc(pixelsizeUM*[-dR dR], Zstep*[-dZ dZ], PSF_im'), ylabel('Z'), xlabel({'R'; 'PSF as measured'})
    
    %1D projections
    fwhmR = fwhm(pixelsizeUM.*(-dR:dR), PSF_im(:, ceil(end/2)));
    fwhmZ = fwhm(Zstep.*(-dZ:dZ),  PSF_im(ceil(end/2),:));
    
    subplot(2,2,3), plot(Zstep.*(-dZ:dZ), PSF_im(ceil(end/2),:)); xlabel(['Z (um) fwhm = ' num2str(fwhmZ, 4)]);
    
    
    subplot(2,2,4), plot(pixelsizeUM.*(-dR:dR), PSF_im(:, ceil(end/2))); xlabel(['R (um) fwhm = ' num2str(fwhmR, 4)]);
    
    
    
    %median filter the psf on each row, but don't reduce the maximum
    for z = 1:size(PSF_im,2)
       [maxval, maxel] = max(PSF_im(:,z));
       PSF_im(:,z) = medfilt2(PSF_im(:,z), [3 1]);
       PSF_im(maxel,z) = maxval;
    end
    PSF_im(PSF_im<thresh | isnan(PSF_im)) = 0; %PSF_im = imgaussfilt(PSF_im, [0.5 0.5]);
    subplot(2,2,2),  imagesc(pixelsizeUM*[-dR dR], Zstep*[-dZ dZ], PSF_im'), ylabel('Z'), xlabel({'R'; 'PSF cleaned and blurred'})    
    PSF_im = PSF_im./sum(PSF_im(:)); %normalize so that the PSF can be treated as a convolution kernel
    
   
    
    
    %create interpolants
    if pixelsizeP<0
        PSF_im = PSF_im(end:-1:1,:);
        pixelsizeP = -pixelsizeP;
    end
    calib.PSF.(['line' int2str(line)]).RZ2I = griddedInterpolant({pixelsizeP*(-dR:dR), Zstep*(-dZ:dZ)}, PSF_im);
    calib.PSF.(['line' int2str(line)]).Zcenter = nanmedian(lineZ);
end

done = false;
while ~done
    resp = input('Update calibration file? Y|N >>', 's');
    switch resp
        case 'Y'
            save([calibdr filesep calibfn], 'calib');
            disp('File saved.')
            done = true;
        case 'N'
            disp('File NOT saved.')
            done = true;
        otherwise
            disp('Bad response!');
    end
end


function closereq(src, evnt)
    for IX = 1:length(hCB)
        selected(IX) = get(hCB(IX), 'Value');
    end
    delete(src);
end
end