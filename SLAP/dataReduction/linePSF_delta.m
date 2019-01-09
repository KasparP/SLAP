function [P, scandata] = linePSF_delta(datastruct, optsin)
%reconstructs a delta PSF for a given set of 2D scan parameters
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
dataDir = 'E:\SLAPmidata\';
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
        end
    end
end

Zcoords = suggestZcoords(scandata);
Pcoords = {scandata.refIMcoords.X, scandata.refIMcoords.Y, Zcoords};
Z0 = ceil(length(Zcoords)/2);

hash = DataHash({Pcoords scandata.Vx scandata.Vy data.timeOfMeasurement 'delta'});
savefile = [dataDir filesep 'PSF\fullPSF\' hash '.mat'];
if ~opts.ignoreHash && exist(savefile, 'file')
    disp('Projection matrix has been previously generated, loading from file...');
    load(savefile);
    P.hash = hash;
else
    v2refX = griddedInterpolant(scandata.refIMcoords.X, 1:length(scandata.refIMcoords.X));
    v2refY = griddedInterpolant(scandata.refIMcoords.Y, 1:length(scandata.refIMcoords.X));
    
    disp('Generating delta PSF for lxl:')
    printstring = [];
    
    nD = 1200*3;
    Px = nan(1, 4*nD*length(scandata.line));
    Py = Px; Pj = Px; Pv = Px;
    P_ix = 0;
    for lxl = 1:length(scandata.line) %for every line position
        if ~mod(lxl, 50)
            fprintf(repmat('\b', 1, length(printstring)));
            printstring = [int2str(lxl) ' of ' int2str(length(scandata.line))];
            fprintf(printstring);
        end
        
        Xv = scandata.Vx(lxl);
        Yv = scandata.Vy(lxl);
        line = scandata.line(lxl);
        lname = ['line' int2str(line)];
        [T,Td] = data.V2T(Xv,Yv,calib.galvos.offset.(lname).X, calib.galvos.offset.(lname).Y, data.(lname).g_angle);
        
        %Get Dcenter, Rcenter
        Dc = data.(lname).interpDc_fromT(T,Td);
        Rc = data.(lname).interpRc_fromT(T,Td);
        Iscale =  data.(lname).interpId(Dc) .* data.(lname).interpIr(Rc); %brightness
                
        Dq = linspace(1,1200, nD); %(subtract Dc when looking up Dprofile!)
        Di = data.(lname).Dprofile(round(Dq));
        Rq = data.(lname).interpR_fromT(T*ones(size(Dq)),Td*ones(size(Dq)),Dq);
                
        theta = (data.(['line' int2str(line)]).theta);
        Rot = [cosd(theta) -sind(theta) ; sind(theta) cosd(theta)];
        pXY = Rot'*[(Dq-600.5) ; (Rq-608.5)];

        pX = pXY(1,:) + 600.5;
        pY = pXY(2,:) + 608.5;
        vX = data.p2vx(pX, pY);
        vY = data.p2vy(pX, pY);
        refX = v2refX(vX);
        refY = v2refY(vY);
        
        select = refX>1 & refX<length(scandata.refIMcoords.X) & refY>1 & refY<length(scandata.refIMcoords.Y);
        for lp = find(select)
            X = floor(refX(lp))+[0 1 0 1];
            Y = floor(refY(lp))+[0 0 1 1];
            Px(P_ix+1:P_ix+4) = X; Py(P_ix+1:P_ix+4) = Y; Pj(P_ix+1:P_ix+4) = lxl;
            It = Iscale.*Di(lp).*([1-mod(refX(lp),1) ; mod(refX(lp),1)] .* [1-mod(refY(lp),1) mod(refY(lp),1)]);
            Pv(P_ix+1:P_ix+4) = It(:);
            P_ix = P_ix+4;
        end
    end
    
    P = sparse(Pj(1:P_ix),sub2ind([length(scandata.refIMcoords.X), length(scandata.refIMcoords.Y), length(Zcoords)] , Px(1:P_ix), Py(1:P_ix), Z0*ones(1,P_ix)), Pv(1:P_ix), length(scandata.line), length(scandata.refIMcoords.X)*length(scandata.refIMcoords.Y)*length(Zcoords));

    fprintf('\nNormalizing P for equal line intensities')
    for line = 1:4
        lineixs = scandata.line==line;
        P(lineixs, :) = P(lineixs,:)./sum(sum(P(lineixs,:))); %#ok<SPRIX>
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