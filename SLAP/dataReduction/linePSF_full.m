function [P, scandata] = linePSF_full(datastruct, optsin)
%reconstructs a 3D PSF for a given X and Y scan
%this needs to be further processed to create a PSF appropriate for a Z
%scan
%%INPUT:
%datastruct is a scandata file OR a measurePSF data file

%For a 3D PSF, call linePSF_3D

opts.ignoreHash = false;
if nargin>1 %UPDATE WITH USER-SPECIFIED OPTIONS
    for field = fieldnames(optsin)'
         opts.(field{1}) = optsin.(field{1});
    end
end

%scandata contains the imaging metadata:
%the xyz location of the illumination line at each laser pulse
%savefile is the location to save the processed PSF data to
if exist('dataDirectory') == 2
    dataDir = dataDirectory;
else
    dataDir = pwd;
end
    
%load the measurement files if not provided
if ~nargin || isempty(datastruct)
    [fn, dr] = uigetfile(dataDir, 'Select a valid PSF or scandata file');
    datastruct = load([dr fn], '-mat');
    if isfield(datastruct, 'scandata') %scandata was passed
        scandata = datastruct.scandata;
        data = scandata.metadata.PSF;
        calib = scandata.metadata.calib;
    elseif isfield(datastruct, 'data') %a PSF file was passed
        scandata = [];
        data = datastruct.data;
        calib = data.calib;
    end
else
    if isfield(datastruct, 'frames') %scandata was passed
        scandata = datastruct;
        data = scandata.metadata.PSF;
        calib = scandata.metadata.calib;
    elseif isfield(datastruct, 'raster') %a PSF file was passed
        scandata = [];
        data = datastruct;
        calib = data.calib;
    end
end

if isempty(scandata) %if no scandata was passed, evaluate across the whole FOV
    %use ideal xyz positions
    scandata.res = 1000;
    scandata.tiling = 1;
    linelength = 4.3;
    
    scandata.refIMcoords.X = linspace(calib.galvos.offset.raster.X-linelength, calib.galvos.offset.raster.X+linelength, 1024); %in V
    scandata.refIMcoords.Y = linspace(calib.galvos.offset.raster.Y-linelength, calib.galvos.offset.raster.Y+linelength, 1024); %in V
    scandata.refIMcoords.Z = 0;
    
    AO = galvo_trace_SLAPmi(4*scandata.res*scandata.tiling, scandata.tiling, calib);
    
    lin_ixs{1} = repmat([true(1, scandata.res) false(1, scandata.res)], 1, 2*scandata.tiling); %for A and B (and E) galvos
    lin_ixs{2} = ~lin_ixs{1};  %for C and D (and F) galvos
    
    %indices for each line
    line = nan(1, 4*scandata.res*scandata.tiling);
    line(lin_ixs{1} & AO(:,5)'<(calib.galvos.offset.E(1)+0.1)) = 1;
    line(lin_ixs{1} & AO(:,5)'>(calib.galvos.offset.E(2)-0.1)) = 2;
    line(lin_ixs{2} & AO(:,6)'<(calib.galvos.offset.F(1)+0.1)) = 3;
    line(lin_ixs{2} & AO(:,6)'>(calib.galvos.offset.F(2)-0.1)) = 4;
    
    scandata.Vx = nan(size(line));
    scandata.Vx(line==1 | line==2) = AO(line==1 | line==2,1);
    scandata.Vx(line==3 | line==4) = AO(line==3 | line==4,3);
    
    scandata.Vy = nan(size(line));
    scandata.Vy(line==1 | line==2) = AO(line==1 | line==2, 2);
    scandata.Vy(line==3 | line==4) = AO(line==3 | line==4, 4);
    
    scandata.frames(1).Z = zeros(size(line));
    
    scandata.line = line;
else
    %extract scandata from galvo measurements
    if ~isfield(scandata, 'refIMcoords')
         %take only the first frame
        [fn, dr] = uigetfile([dataDir filesep '*.mat'], 'Select the reference image .mat file. Cancel to use default coordinates.');
        if ~fn
            if isfield(scandata, 'metadata') && scandata.metadata.aperture
                linelength = scandata.metadata.PSF.aperture.linelength;
            else
                linelength = calib.galvos.linelength;
            end
            scandata.refIMcoords.X = linspace(calib.galvos.offset.raster.X-linelength, calib.galvos.offset.raster.X+linelength, 1024); %in V
            scandata.refIMcoords.Y = linspace(calib.galvos.offset.raster.Y-linelength, calib.galvos.offset.raster.Y+linelength, 1024); %in V
            scandata.refIMcoords.Z = 0;
        else
            datastruct = load([dr fn], '-mat');
            scandata.refIMcoords = datastruct.refIM.M.coords;
            
            %test if X and Y are flipped?
%             Xtmp = scandata.refIMcoords.Y; Ytmp = scandata.refIMcoords.X;
%             scandata.refIMcoords.Ytmp = Xtmp; scandata.refIMcoords.X = Ytmp;
        end
    end
end

if scandata.refIMcoords.Z==0
    Zcoords = linspace(-6, 6, 13); %The Z coordinates at which to reconstruct the PSF
    disp(['Reconstructing PSF slices at: ' num2str(Zcoords)]);
    coords = input('If you prefer other planes, enter them now >>');
    if ~isempty(coords)
        Zcoords = coords;
    end
else
    %otherwise we reconstruct on a sampling grid the matches the reference image spacing
    Zcoords = suggestZcoords(scandata);
end

Pcoords = {scandata.refIMcoords.X, scandata.refIMcoords.Y, Zcoords};
if ~isfield(opts, 'Zlxl')
   opts.Zlxl = zeros(1, length(scandata.line));
end

hash = DataHash({Pcoords scandata.Vx scandata.Vy data.timeOfMeasurement opts.Zlxl});
savefile = [dataDir filesep 'PSF' filesep 'fullPSF' filesep hash '.mat'];
if ~opts.ignoreHash && exist(savefile, 'file')
    disp('Projection matrix has been previously generated, loading from file...');
    load(savefile);
    P.hash = hash;
else
    expected_nnz = length(scandata.refIMcoords.X)*30*min(6, length(Zcoords))*length(scandata.line);
    P = spalloc(length(scandata.refIMcoords.X)*length(scandata.refIMcoords.Y)*length(Zcoords), length(scandata.line),expected_nnz); %preallocate

    %The sampling grid
    %X is the row coordinate
    [Xmesh, Ymesh] = ndgrid(scandata.refIMcoords.X, scandata.refIMcoords.Y);
    Xgrid = repmat(reshape(data.v2px(Xmesh(:),Ymesh(:)), size(Xmesh)), 1, 1, length(Zcoords));
    Ygrid = repmat(reshape(data.v2py(Xmesh(:),Ymesh(:)), size(Xmesh)), 1, 1, length(Zcoords));
    Zgrid = repmat(reshape(Zcoords, [1 1 length(Zcoords)]), size(Xgrid, 1), size(Xgrid,2));
    Xmesh = repmat(Xmesh, [1 1 length(Zcoords)]);
    Ymesh = repmat(Ymesh, [1 1 length(Zcoords)]);
    
    %must convert X,Y voltages to pixels in the transformed space
    Rgrid = cell(1,4); Dgrid = cell(1,4);
    for line = 1:4
        theta = data.(['line' int2str(line)]).theta;
        Dgrid{line} = cosd(theta)*(Xgrid-600.5) - sind(theta)*(Ygrid-608.5) + 600.5;
        Rgrid{line} = sind(theta)*(Xgrid-600.5) + cosd(theta)*(Ygrid-608.5) + 608.5;
    end
    
    disp('Generating line PSF for lxl:')
    printstring = [];
    for lxl = 1:length(scandata.line) %for every line position
        if ~mod(lxl, 50)
            fprintf(repmat('\b', 1, length(printstring)));
            printstring = [int2str(lxl) ' of ' int2str(length(scandata.line))];
            fprintf(printstring);
        end
        
        Xv = scandata.Vx(lxl);
        Yv = scandata.Vy(lxl);
        Zpos = opts.Zlxl(lxl); 
        
        line = scandata.line(lxl);
        lname = ['line' int2str(line)];
        [T,Td] = V2T(Xv,Yv,calib.galvos.offset.(lname).X, calib.galvos.offset.(lname).Y, data.(lname).g_angle);
        
        %Get Dcenter, Rcenter
        Dc = data.(lname).interpDc_fromT(T,Td);
        Rc = data.(lname).interpRc_fromT(T,Td);
        Iscale =  data.(lname).interpId(Dc) .* data.(lname).interpIr(Rc); %brightness
        
        %select is the set of pixels in camera space that we want to reconstruct
        select = abs(Rgrid{line}-Rc)<11 & abs(Dgrid{line} - 600.5 - Dc)<520;
        Rq = Rgrid{line}(select);
        Dq = Dgrid{line}(select);
        Zq = Zgrid(select);
        
        %get differences for our specific R,D
        D = Dq - Dc; %where on the line are we, range: [0-1200] used to look up Z,I,R
        R = data.(lname).interpR_fromT(T*ones(size(Dq)),Td*ones(size(Dq)),Dq);
        Z = data.(lname).interpZ_fromT(T*ones(size(Dq)),Td*ones(size(Dq)),Dq);
        dR = Rq - R; %data.(lname).interpR(Xv*ones(size(Dq)), Yv*ones(size(Dq)), Dq); %How far from center we are, range [-1200:1200]
        
        %set dZ to 0, i.e. assume that the line lies in the same (curved) field as the raster image. This assumption is empirically valid!
        dZ = (Zq-Zpos); % - (Z-data.v2Z(Xmesh(select),Ymesh(select))); %first term accounts for nominal Z position, 2nd term is a correction for field curvature

        %increasing values of dZ indicate a more superficial part
        %of the PSF (i.e. when dZ>0, the PSF is deeper than the sample point)
        %increasing values of Zq indicate the PSF is deeper in the sample
        %increasing values of Z indicate the PSF is deeper in the sample
        
        nans = isnan(R) | isnan(D) | isnan(dZ);
        I = Iscale*data.(lname).Dprofile(round(D));
        I(~nans) = I(~nans).*calib.PSF.(['line' int2str(line)]).RZ2I(dR(~nans),dZ(~nans));
        I(nans | isnan(I) | I<0) = 0;
        P(select(:),lxl) = I; %#ok<SPRIX>
    end 
    P = P'; % The sparse matrix was most efficient to build by column so we actually made P' earlier
    
    fprintf('\nNormalizing P for equal line intensities')
    for line = 1:4
        lineixs = scandata.line==line;
        P(lineixs, :) = P(lineixs,:)./sum(sum(P(lineixs,:)));
    end
    
    P = struct('P', P);
    P.coords = Pcoords;
    P.hash = hash;
    
    disp(['Saving data to ' savefile ' ...'])
    save(savefile, 'P', '-v7.3');
end
disp('Done linePSF.')
end

function Z = suggestZcoords(scandata)
Z = scandata.refIMcoords.Z - mean(scandata.refIMcoords.Z);
Z = Z(abs(Z)<=6);
end