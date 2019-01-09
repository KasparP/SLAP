function SLAPMi_ROI
%A gui for quick backprojections of scandata, and for selection of ROIs for
%full reconstruction

hSI = []; hSLAPMi = [];
if evalin('base','exist(''hSI'',''var'')')
        hSI = evalin('base','hSI');
end
if evalin('base','exist(''hSLAPMi'',''var'')')
        hSLAPMi = evalin('base','hSLAPMi');
end

basedir = 'E:\SLAPMiData';

%state variables
alpha = 0.2;
Zpos = 1;
contrast = 100;
scanfield = [];
Xsz = 400; Ysz = 400; Zsz = 9;
imBuffer = 20*rand(Xsz,Ysz,Zsz);
doMax = false;
drawing = false;
curr_mask = 1;
refIM = [];
scandata = [];
P = [];
yAligned = [];
expected = [];
SLMmask = 1;
fig2locs = [];
fig2factor = [];
colors= {'r','g','b'};

opts.maxDFF = 2;
opts.fig2flag = false;
opts.channel = 1;
opts = optionsGUI(opts);

saveDir = 'E:\SLAPmidata\SLM\';
if ~exist(saveDir, 'dir')
    mkdir(saveDir)
end

%Initialize figure
handles.fig = figure('name', 'SLAPMi_ROI', 'pos', [50 50 1200 650], 'CloseRequestFcn',@exit);
set(handles.fig, 'resize', 'on', 'dockcontrols', 'off', 'numbertitle', 'off', 'visible', 'off')
set(handles.fig, 'WindowScrollWheelFcn', @scrollFcn);
set(handles.fig, 'colormap', [0 0 0.3 ; repmat(0:1/252:1, 3, 1)' ; 1 0 0]); %set the default colormap to Hi-Lo coloring

%figure Icon
warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
jframe=get(handles.fig,'javaframe');
icon = imread([fileparts(which('slapmi')) filesep 'roiIcon.bmp']);
javaImage = im2java(icon);
jIcon=javax.swing.ImageIcon(javaImage);
jframe.setFigureIcon(jIcon);

%initialize Axes
handles.ax1 = axes('units', 'pixels'); 
set(handles.ax1, 'pos', [8   55  400  400]);
handles.ax2 =  axes('units', 'pixels'); 
set(handles.ax2, 'pos', [450   55  400  200], 'visible', 'on');
handles.ax3 =  axes('units', 'pixels'); 
set(handles.ax3, 'pos', [450+30   255+20  400-30  200-20], 'visible', 'on');


%UICONTROLS

%Panels
handles.bottomPanel = uipanel('units', 'pixels',  'Position', [0 0 540 60]);

%CONTRAST SLIDERS
tmpX = 8;
tmpY = 12;
            uicontrol(...
                'Units','pixels',...
                'Parent',handles.fig,...
                'Style','text',...
                'Position',[tmpX tmpY+20 50 14],...
                'String','Contrast')
           handles.etContrast =  uicontrol(...
                'Units','pixels',...
                'Parent',handles.fig,...
                'Style','edit',...
                'Position',[tmpX+55 tmpY+20 35 18],...
                'callback', @(src,~)(setContrast(str2double(src.String))));
            handles.slContrast = uicontrol(...
                'Units','pixels',...
                'Parent',handles.fig,...
                'Style','slider',...
                'Position',[tmpX tmpY 90 16],...
                'Min',0,'Max',1,...
                'sliderstep', [1/1000 1/100],...
                'callback', @(src,~)(setContrast(src.Value))); 
          
            handles.pbSeg = uicontrol(...
                'Units','pixels',...
                'Parent',handles.fig,...
                'Callback',@(varargin)makeSeg,...
                'Position',[tmpX+158 tmpY+20 45 25],...
                'String','Segment',...
                'Tag','pbSeg');
            
            handles.pbScandata = uicontrol(...
                'Units','pixels',...
                'Parent',handles.fig,...
                'Callback',@(varargin)getScandata,...
                'Position',[tmpX+208 tmpY+20 107 25],...
                'String','get scandata',...
                'Tag','pbRefIM');
            handles.pbRefIM = uicontrol(...
                'Units','pixels',...
                'Parent',handles.fig,...
                'Callback',@(varargin)getReferenceImage,...
                'Position',[tmpX+208 tmpY-5 60 25],...
                'String','get refIM',...
                'Tag','pbRefIM');
            handles.puRefCh = uicontrol(...
                'Units','pixels',...
                'Style', 'popupmenu',...
                'Parent',handles.fig,...
                'Callback',@(varargin)setCh,...
                'Position',[tmpX+270 tmpY-5 45 25],...
                'String',{'Ch1'; 'Ch2'},...
                'Tag','puRefCh');
            
             handles.pbGrow = uicontrol(...
                'Units','pixels',...
                'Parent',handles.fig,...
                'Callback',@(varargin)(grow),...
                'Position',[tmpX+320 tmpY 50 20],...
                'String','Grow',...
                'Tag','pbGrow');
             handles.pbShrink = uicontrol(...
                'Units','pixels',...
                'Parent',handles.fig,...
                'Callback',@(varargin)(shrink),...
                'Position',[tmpX+375 tmpY 50 20],...
                'String','Shrink',...
                'Tag','pbShrink');
            
             handles.tbMaxProj = uicontrol(...
                'Units','pixels',...
                'Parent',handles.fig,...
                'style', 'togglebutton',...
                'Callback',@(src, evnt)(setMaxProj(src.Value)),...
                'Position',[tmpX+320 tmpY+23 105 20],...
                'String','MAX Projection',...
                'Tag','tbMaxProj');

            uicontrol(...
                'Units','pixels',...
                'Parent',handles.fig,...
                'Style','text',...
                'Position',[tmpX+443 tmpY+28 70 14],...
                'String','PATTERN')
            handles.tb1 =  uicontrol(...
                'Units','pixels',...
                'Parent',handles.fig,...
                'Style','toggle',...
                'Position',[tmpX+435 tmpY-3 30 30],...
                'String', 'R', ...
                'callback', @(varargin)(togglePattern(1)));
            handles.tb2 =  uicontrol(...
                'Units','pixels',...
                'Parent',handles.fig,...
                'Style','toggle',...
                'Position',[tmpX+465 tmpY-3 30 30],...
                'String', 'G', ...
                'callback', @(varargin)(togglePattern(2)));
            handles.tb3 =  uicontrol(...
                'Units','pixels',...
                'Parent',handles.fig,...
                'Style','toggle',...
                'Position',[tmpX+495 tmpY-3 30 30],...
                'String', 'B', ...
                'callback', @(varargin)(togglePattern(3)));
            
            
%Z control
            handles.slZ = uicontrol(...
                'Units','pixels',...
                'Parent',handles.fig,...
                'Style','slider',...
                'Position',[610 55 16 600],...
                'Min',1,'Max',Zsz,...
                'Value', round(Zsz/2),...
                'sliderstep', [1/Zsz 1/Zsz],...
                'callback', @(src, evnt)(updateZ(src.Value)));


             
           
            
set(handles.fig, 'ResizeFcn', @resizeFcn);
            
%initialize Screen
handles.screen = imshow(max(imBuffer,[],3), 'parent', handles.ax1, 'border', 'tight');
set(get(handles.screen, 'Parent'), 'Clim', [0 1]);
axis(handles.ax1, 'manual')
axis(handles.ax1, 'off');
hold(handles.ax1, 'on');

handles.Paint.selection = line(handles.ax1, [-1,-1,-1,-1,-1], [-1,-1,-1,-1,-1], 'color', 'c', 'linewidth', 2);
handles.mask = cell(1,3); 
updateCData;

%select scandata
[fn, basedir] = uigetfile('*.mat','Select reduced scandata',basedir);
if ~fn
    return
end
S = load([basedir filesep fn]);
scandata = S.scandata; clear S;

getReferenceImage;
getScandata;
updateCData;

handles.Paint.startpoint = [];
set(handles.fig, 'WindowButtonDownFcn', @manualpaint_BDF);
set(handles.fig, 'WindowButtonUpFcn', @manualpaint_BUF);
set(handles.fig, 'WindowButtonMotionFcn', @manualpaint_BMF);
set(handles.ax2, 'xtick', [], 'ytick', []);
resizeFcn;
set(handles.fig, 'visible', 'on')           

function manualpaint_BDF (src, event)
cp = get(handles.ax1,'CurrentPoint');
cp = cp(1, 1:2);

%If the point is within the axes
if all(cp>=1 & cp<=[Ysz Xsz])
    mod = get(src, 'SelectionType');
    switch mod
        case 'extend'
            %paint a circle on current point
            M = get(handles.mask{curr_mask}, 'Alphadata');
            [XX,YY] = meshgrid(1:Xsz,1:Ysz);
            M(sqrt((XX-cp(1)).^2 + (YY-cp(2)).^2)<5) = alpha;
            set(handles.mask{curr_mask}, 'Alphadata', M);
            drawing = true;
        otherwise
            handles.Paint.startpoint = cp;
    end
end
end
    function manualpaint_BUF (src, event)
        persistent X
        persistent Y
        drawing = false;
        mod = get(src, 'SelectionType');
        cp = get(handles.ax1,'CurrentPoint');
        cp = cp(1,1:2);
        if ~isempty(handles.Paint.startpoint)
            rect = [];
            X = double(repmat(1:Xsz, Ysz, 1));
            Y = double(repmat((1:Ysz)', 1, Xsz));
            %if we just clicked, not dragged
            if all(abs(cp-handles.Paint.startpoint)<2)
                rect = abs(Y-cp(1)) + abs(X-cp(2))<3;
                rect = rect';
            else %if all(cp>=1 & cp<=[F2P.acq.Ybound F2P.acq.Xbound])
                %define the rectangle between startpoint and endpoint
                sXY = [cp(1) handles.Paint.startpoint(1)];
                dXY = [cp(2) handles.Paint.startpoint(2)];
                rect = X>=min(sXY) & X<=max(sXY) & Y>=min(dXY) & Y<=max(dXY);
            end
            if ~isempty(rect)
                M = get(handles.mask{curr_mask}, 'Alphadata');
                mod = get(src, 'SelectionType');
                switch mod
                    case 'normal'
                        %add the pixels to the mask
                        M(rect) = alpha;
                    case 'alt'
                        %remove pixels from the mask
                        M(rect) = 0;
                end
                set(handles.mask{curr_mask}, 'alphadata', M);
            end
        end
        handles.Paint.startpoint = [];
        set(handles.Paint.selection, 'Xdata', -1, 'Ydata', -1);
        updateROI(curr_mask)
    end

    function manualpaint_BMF (src, event)
        mod = get(src, 'SelectionType');
        cp = get(handles.ax1,'CurrentPoint');
        cp = cp(1,1:2);
        switch mod
            case 'extend'
                if drawing
                    M = get(handles.mask{curr_mask}, 'Alphadata');
                    [XX,YY] = meshgrid(1:Xsz,1:Ysz);
                    M(sqrt((XX-cp(1)).^2 + (YY-cp(2)).^2)<20) = alpha;
                    set(handles.mask{curr_mask}, 'Alphadata', M);
                end
            otherwise
                drawing = false;
                if ~isempty(handles.Paint.startpoint);
                    x1=cp(1); y1=cp(2); x3=handles.Paint.startpoint(1); y3=handles.Paint.startpoint(2);
                    x2 = x1; y2 = y3;
                    x4 = x3; y4 = y1;
                    %color depends on click type
                    if strcmp(mod, 'alt')
                        C = 'r';
                    else
                        C = 'c';
                    end
                    set(handles.Paint.selection, 'Xdata', [x1 x2 x3 x4 x1], 'Ydata', [y1 y2 y3 y4 y1], 'color', C);
                else
                    set(handles.Paint.selection, 'Xdata', -1, 'Ydata', -1);
                end
        end
    end

    function getReferenceImage
        if ~isempty(hSI)
            basedir = hSI.hScan2D.logFilePath;
        end
        [fn, basedir] = uigetfile('*.mat','Select your aligned REF IM',basedir);
        if length(fn)>3 && strcmpi(fn(end-3:end), '.mat')%a single reference image was selected
            S = load([basedir filesep fn]);
            refIM = S.refIM;
        else
            disp('Using blank reference image...')
            pixelSize = 0.0044; %in volts
            refIM.data = ones(1280,1280,1);
            refIM.IM = refIM.data;
            refIM.M.coords.X = -639.5*pixelSize + scandata.metadata.calib.galvos.offset.raster.X + (1:1280) * pixelSize;
            refIM.M.coords.Y = -639.5*pixelSize + scandata.metadata.calib.galvos.offset.raster.Y + (1:1280) * pixelSize;
            refIM.M.coords.Z = 0;
            scanfield.pixelToRefTransform(1) = pixelSize;
            refIM.metadata.rois.RoiGroups.imagingRoiGroup.rois{1, 1}.scanfields(1) = scanfield;
        end
        
        %update the buffer and the scanfield, etc with the reference image
        channelSelect = min(size(refIM.data,4), opts.channel); 
        imBuffer = squeeze(refIM.data(:,:,:,channelSelect));
        imBuffer = double(permute(imBuffer, [2 1 3]));
        imBuffer = imBuffer./max(imBuffer(:));
        scanfield = refIM.metadata.rois.RoiGroups.imagingRoiGroup.rois{1, 1}.scanfields(1);
        scanfield.meshgrid = @(varargin)(meshgrid(refIM.M.coords.X, refIM.M.coords.Y));
        
        %update display
        Xsz = size(imBuffer,1); Ysz = size(imBuffer,2); Zsz = size(imBuffer,3);
        updateZ(Zsz/2);
        set(handles.ax1, 'xlim', 0.5+[0 Xsz], 'ylim', 0.5+[0 Ysz]);
        set(handles.slZ, 'Max', Zsz, 'Value', round(Zsz/2), 'sliderstep', [1/Zsz 1/Zsz]);
        updateCData;
    end


    function getScandata
        if ~isempty(hSLAPMi)
            basedir = [hSLAPMi.dataDir filesep hSLAPMi.user];
        end

        scandata.refIMcoords = refIM.M.coords;
        P = linePSF_delta(scandata);
        
        %compute the mask
        Vslm = scandata.metadata.SLM.pattern;
        SLMhigh  = abs(double(Vslm)-scandata.metadata.calib.SLM.lut(:,:,2))<abs(double(Vslm)-scandata.metadata.calib.SLM.lut(:,:,1));
        Tslm = scandata.metadata.calib.SLM.T(:,:,2).* SLMhigh + scandata.metadata.calib.SLM.T(:,:,1).*~SLMhigh;
        Tslm = max(Tslm, 5e-3); %Sometimes the SLM transmission seems to be underestimated
        scanfield = refIM.metadata.rois.RoiGroups.imagingRoiGroup.rois{1, 1}.scanfields(1);
        scanfield.meshgrid = @(varargin)(meshgrid(refIM.M.coords.X, refIM.M.coords.Y));
        pixelToRefTransform = scandata.metadata.SLM.pixelToRefTransform;
        Tim = mapSLMToImage(Tslm, pixelToRefTransform, scanfield);
        mask = imgaussfilt(Tim, 0.005/scanfield.pixelToRefTransform(1))';
        SLMmask = mask;
        
        %cut down P
        Psz = [length(P.coords{1}), length(P.coords{2}), length(P.coords{3})];
        sumP = sum(sum(reshape(full(sum(P.P,1)),Psz)));
        selectPlanes = sumP>max(sumP)/2;
        select = false(Psz); select(:,:,selectPlanes) = true;
        P.P = P.P(:, select(:));
        P.coords{3} = P.coords{3}(selectPlanes);
        expected=P.P*reshape(refIM.IM(:,:,Zpos).*SLMmask,[],1);
        
        y = [scandata.frames.pmtData];
        ymean = nanmean(y,1);
        [y0,yAligned] = yAlign(y, expected, scandata.line);
        
        cla(handles.ax3);cla(handles.ax2);
        plot(handles.ax3, ymean./nanmean(ymean), 'color', [0.3 0.3 0.3])
        hold(handles.ax3, 'on');
        legend(handles.ax3, 'Whole-frame Average')
        
        plot(handles.ax2, sqrt(expected./mean(expected)), 'color', [0.7 0.7 0.7], 'linewidth', 2)
        hold(handles.ax2, 'on');
        plot(handles.ax2, sqrt(y0./mean(y0)), 'color', [0.1 0.1 0.1])
        legend(handles.ax2, {'Measured', 'Expected'});
        
        for roi_ix = 1:3
            updateROI(roi_ix);
        end
        
        if opts.fig2flag
            showFig2(yAligned);
        end
    end

    function updateROI(N)
        if isfield(handles, ['plot2' int2str(N)])
            delete(handles.(['plot2' int2str(N)]));
            delete(handles.(['plot3' int2str(N)]));
        end
        
        select = get(handles.mask{N}, 'alphadata')>0;
        if any(select(:))
            P2 = P.P*reshape(refIM.IM(:,:,Zpos).*SLMmask.*select, [],1);
            handles.(['plot2' int2str(N)]) = plot(handles.ax2, sqrt(P2./mean(expected)), 'color', colors{N});
            %backprojection
            y2 = P2'*yAligned;
            handles.(['plot3' int2str(N)]) = plot(handles.ax3, y2./mean(y2), 'color', colors{N});
        end

    end

function updateCData
    if doMax
        I = sqrt(double(max(imBuffer,[],3)))'.*sqrt(SLMmask);
    else
        I = sqrt(double(imBuffer(:,:,Zpos))'.*sqrt(SLMmask));
    end
    set(handles.screen, 'CData', I);
    if isempty(handles.mask{1}) || ~all(size(get(handles.mask{1}, 'alphadata')) == [Xsz Ysz])
        for mask_ix = 1:length(handles.mask)
            delete(handles.mask{mask_ix});
            mask = zeros(Xsz,Ysz,3);
            mask(:,:,mask_ix) = 1;
            handles.mask{mask_ix} = imshow(mask, 'parent', handles.ax1);
            set(handles.mask{mask_ix}, 'alphadata', zeros(Xsz,Ysz), 'hittest', 'off');
        end
    end
end

function updateZ(val)
    Zpos = round(max(1,min(Zsz, val)));
    set(handles.slZ, 'Value', Zpos);
    updateCData;
end
    function setContrast(value)
        contrast = max(0, min(get(handles.slContrast, 'Max'), value));
        set(handles.etContrast, 'String', num2str(contrast,3));
        set(handles.slContrast, 'Value', contrast)
        set(handles.ax1, 'Clim', [0 max(1e-10, contrast)]);
    end
    function setMaxProj(value)
        doMax = logical(value);
        updateCData;
    end

    function scrollFcn(src, evnt)
        %updateZ(Zpos - evnt.VerticalScrollCount);
    end

    function grow
        S = strel('disk', 3);
        M = get(handles.mask{curr_mask}, 'alphadata')>0;
        M = imdilate(M,S);
        set(handles.mask{curr_mask},'alphadata', alpha*double(M));
    end
    function shrink
        S = strel('disk', 5);
        M = get(handles.mask{curr_mask}, 'alphadata')>0;
        M = imerode(M,S);
        set(handles.mask{curr_mask},'alphadata', alpha*double(M));
    end

    function togglePattern(n)
        curr_mask = n;
        for ix = 1:3
            set(handles.(['tb' int2str(ix)]), 'Value', ix==n)
        end
    end

    function setCh
       opts.channel = handles.puRefCh.Value; 
    end

    function exit(src, evnt)
        if strcmp(questdlg('Exit SLAPMi_ROI?', 'ROI GUI', 'Yes', 'No', 'No'), 'Yes')
            try
                delete(handles.fig);
                delete(handles.fig2);
            end
        end
    end

    function resizeFcn(src, evnt)
        p = get(handles.fig, 'pos');
        w = p(3); h = p(4);
        
        bpPos = get(handles.bottomPanel, 'pos');
        
        ax1H = h-bpPos(4);
        ax1W = min(ax1H, w/2-15);
        
        ax2W = w-(ax1W+40);
        
        set(handles.ax1, 'pos', [5 bpPos(4) ax1W ax1H]);
        set(handles.slZ, 'pos', [ax1W+5 bpPos(4) 15 ax1H])
        set(handles.ax2, 'pos', [ax1W+20 5 ax2W ceil((h-5)/2)], 'visible', 'on')
        set(handles.ax3, 'pos', [ax1W+40 5+ceil((h-5)/2) ax2W floor((h-5)/2)], 'visible', 'on')
    end

function [yRef,yAligned] = yAlign(yIn, expected, lineIDs)
    %get an estimate of the F0 line intesities
    nCut = round(size(yIn,2)/4);
    [~, sortorder] = sort(nansum(yIn,1));
    yOut = nanmean(yIn(:, sortorder(10:10+nCut)), 2);
    yOut(sum(~isnan(yIn),2)<30) = nan;
    
    %replace nans with prediction from surrounding measurements
    yRef = max(0, inpaint_nans(yOut));
    
    %align to expected
    yAligned = yIn;
    yAligned(isnan(yAligned)) = 0;
end

    function showFig2(yAligned)
           try
               delete(handles.fig2)
               delete(handles.fig2mask)
           end
           
           tau = 5;
           [~,DFF] = fastDFF(yAligned, tau);
           
           handles.fig2 = figure('Name', 'LXL SELECT (Shift-Click to add; Right Click to clear)', 'WindowButtonDownFcn', @fig2BDF, 'numbertitle', 'off');
           handles.fig2ax = axes(handles.fig2);
           fig2factor = size(yAligned,1)./size(DFF,1);
           handles.fig2im = imagesc( DFF, 'parent', handles.fig2ax);
           set(handles.fig2ax, 'pos', [0 0 1 1], 'xtick', [], 'ytick', [], 'clim', [0 min(max(DFF(:)),opts.maxDFF)])
           axis(handles.fig2ax, 'tight')
           hold(handles.fig2ax, 'on')
            mask = zeros(Xsz,Ysz,3);
            mask(:,:,1:2) = 1;
            handles.fig2mask = imshow(mask, 'parent', handles.ax1);
            set(handles.fig2mask, 'alphadata', zeros(Xsz,Ysz), 'hittest', 'off');
    end

    function fig2BDF (src, event)
        cp = get(handles.fig2ax,'CurrentPoint');
        cp = cp(1, 1:2);
        
        Xl = get(handles.fig2ax, 'xlim'); Xl = Xl(2);
        Yl = get(handles.fig2ax, 'ylim'); Yl = Yl(2);
        
        %If the point is within the axes
        if all(cp>=1 & cp<=[Xl Yl])
            mod = get(src, 'SelectionType');
            switch mod
                case 'extend' %add a line
                    fig2locs = [fig2locs round(cp(2))];
                case 'normal' %clear all and draw a line
                    fig2locs = round(cp(2));
                case 'alt' %clear all
                    fig2locs = [];
            end
            drawFig2locs;
        end
        
    end

    function drawFig2locs
        try
            delete(handles.fig2scatter)
        end
        handles.fig2scatter = scatter(30* ones(1,length(fig2locs)), fig2locs, 'parent', handles.fig2ax, 'markeredgecolor', 'y', 'linewidth', 2);
        set(handles.fig2mask, 'alphadata', 0.6*full(reshape(any(P.P(round(fig2locs.*fig2factor),:),1), Xsz,Ysz)))
    end


    function makeSeg
        ROIs = [];
        for ix = 1:3
            A = get(handles.mask{ix}, 'alphadata');
            if any(A(:))
                ROIs(:,:,end+1) = A>0;
            end
        end
        [fnsave drsave] = uiputfile(basedir);
        save([drsave filesep fnsave], 'ROIs');
    end
end




