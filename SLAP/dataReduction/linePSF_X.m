function [P, scandata] = linePSF_X(datastruct, optsin)
%reconstructs a '3-plane' delta PSF for a given set of 2D scan parameters
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
dataDir = dataDirectory;
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

%see notes: http://www.evernote.com/l/ArR_rpb33hZIv7-qyeIbd5x1Fw3yJr2LbPc/
Popts.planes = [0 0.75 1.5]; %planes with nonzero weight for the PSF
Popts.weights = [1 20/80 12/80]; %The grand sum of the PSF in the given planes
Popts.dRs = [0 0 0.45]; %in um, the distance of the delta functions from the nominal center, for each plane

dRs = unique(Popts.dRs); nZ = length(dRs);

hash = DataHash({Pcoords scandata.Vx scandata.Vy data.timeOfMeasurement Popts});
try
    if ~exist([dataDir filesep 'PSF\fullPSF\'], 'dir')
        mkdir([dataDir filesep 'PSF'], 'fullPSF');
    end
end
savefile = [dataDir filesep 'PSF\fullPSF\' hash '.mat'];
if ~opts.ignoreHash && exist(savefile, 'file')
    disp('Projection matrix has been previously generated, loading from file...');
    load(savefile);
    P.hash = hash;
else
    v2refX = griddedInterpolant(scandata.refIMcoords.X, 1:length(scandata.refIMcoords.X));
    v2refY = griddedInterpolant(scandata.refIMcoords.Y, 1:length(scandata.refIMcoords.X));
    umPerPixel = sqrt(diff(data.p2vx(600 + [-1 1],600 + [0 0])).^2 + diff(data.p2vy(600 + [-1 1],600 + [0 0])).^2)./calib.galvos.pixelsizeVperUM; 

    disp('Generating X PSF for lxl:')
    printstring = [];
    
    nD = ceil(1200*3.3); %the 'stride' to sample the line at; should be >>1200, larger numbers are more accurate but slower
    for Z = nZ:-1:1 %preallocate
        Px{Z} = nan(1, 4*nD*length(scandata.line));
        Py{Z} = Px{Z}; Pj{Z} = Px{Z}; Pv{Z} = Px{Z};
        P_ix{Z} = 0;
    end
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
        
        for Z = 1:nZ %we will create a psf for the center plane and one for +Z
            dR = unique([dRs(Z) -dRs(Z)]);
            for dR_ix = 1:length(dR)
                pXY = Rot'*[(Dq-600.5) ; (Rq + (dR(dR_ix)./umPerPixel) -608.5)];
                
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
                    Px{Z}(P_ix{Z}+1:P_ix{Z}+4) = X; Py{Z}(P_ix{Z}+1:P_ix{Z}+4) = Y; Pj{Z}(P_ix{Z}+1:P_ix{Z}+4) = lxl;
                    It = Iscale.*Di(lp).*([1-mod(refX(lp),1) ; mod(refX(lp),1)] .* [1-mod(refY(lp),1) mod(refY(lp),1)]);
                    Pv{Z}(P_ix{Z}+1:P_ix{Z}+4) = It(:);
                    P_ix{Z} = P_ix{Z}+4;
                end
            end
        end
    end
    
    
    for Z = nZ:-1:1
        P{Z} = sparse(Pj{Z}(1:P_ix{Z}),sub2ind([length(scandata.refIMcoords.X), length(scandata.refIMcoords.Y)] , Px{Z}(1:P_ix{Z}), Py{Z}(1:P_ix{Z})), Pv{Z}(1:P_ix{Z}), length(scandata.line), length(scandata.refIMcoords.X)*length(scandata.refIMcoords.Y));
    end
    
    %assemble the full PSF by taking the appropriate Ps and multiplying by
    %the weights
    P = [P{2}.*(Popts.weights(3)./2)  P{1}.*Popts.weights(2) P{1}.*Popts.weights(1) P{1}.*Popts.weights(2) P{2}.*(Popts.weights(3)./2)];
    
    fprintf('\nNormalizing P for equal line intensities')
    for line = 1:4
        lineixs = scandata.line==line;
        P(lineixs, :) = P(lineixs,:)./sum(sum(P(lineixs,:))); %#ok<SPRIX>
    end
    
    P = struct('P', P);
    P.type = 'X'; %this is an X psf
    P.coords = Pcoords;
    P.hash = hash;
    P.opts = Popts;
    
    disp(['Saving data to ' savefile ' ...'])
    save(savefile, 'P', '-v7.3');
end
disp('Done linePSF.')
end

function Z = suggestZcoords(scandata)
Z = -1.5:0.75:1.5;
end