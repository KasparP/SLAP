function SLMcontrol(varargin)
%A gui for controlling the SLM.

hSI = []; hSLAPMi = [];
if evalin('base','exist(''hSI'',''var'')')
        hSI = evalin('base','hSI');
end
if evalin('base','exist(''hSLAPMi'',''var'')')
        hSLAPMi = evalin('base','hSLAPMi');
end

%state variables
alpha = 0.2;
channel = 1;
Zpos = 1;
thresh = 500;
contrast = 100;
scanfield = [];
Xsz = 400; Ysz = 400; Zsz = 9;
imBuffer = 20*rand(Xsz,Ysz,Zsz);
doMax = false;
do3D = false;
drawing = false;
saveDir = 'E:\SLAPmidata\SLM\';
if ~exist(saveDir, 'dir')
    mkdir(saveDir)
end

%initialize ROIs
ROIs = struct('hScatter', {}, 'hLine', {}, 'center', {});
cell_radius = 10;

%Initialize figure
handles.fig = figure('name', 'SLM control', 'pos', [50 50 750 660], 'CloseRequestFcn',@exit);
set(handles.fig, 'resize', 'on', 'menubar', 'none', 'dockcontrols', 'off', 'toolbar', 'none', 'numbertitle', 'off')
set(handles.fig, 'WindowScrollWheelFcn', @scrollFcn);
set(handles.fig, 'colormap', [0 0 0.3 ; repmat(0:1/252:1, 3, 1)' ; 1 0 0]); %set the default colormap to Hi-Lo coloring


%figure Icon % for sandbox: comment out this section
%warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
%jframe=get(handles.fig,'javaframe');
%icon = imread([fileparts(which('slapmi')) filesep 'slmIcon.bmp']);
%javaImage = im2java(icon);
%jIcon=javax.swing.ImageIcon(javaImage);
%jframe.setFigureIcon(jIcon);


%initialize Axes
handles.ax1 = axes('units', 'pixels'); 
set(handles.ax1, 'pos', [8   55  600  600]);

%initialize Screen
handles.screen = imshow(max(imBuffer,[],3), 'parent', handles.ax1, 'border', 'tight');
set(get(handles.screen, 'Parent'), 'Clim', [0 50]);
axis manual
axis off
hold on

handles.Paint.selection = line([-1,-1,-1,-1,-1], [-1,-1,-1,-1,-1], 'color', 'c', 'linewidth', 2);
handles.mask = []; handles.masks=double.empty([0 0 0]);
updateCData;

%UICONTROLS
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
                'Min',-10,'Max',2000,...
                'sliderstep', [1/1000 1/100],...
                'callback', @(src,~)(setContrast(src.Value))); 
            
            uicontrol(...
                'Units','pixels',...
                'Parent',handles.fig,...
                'Style','text',...
                'Position',[tmpX+100 tmpY+22 100 14],...
                'String','Import SI Image')
            handles.pbImportCh1 = uicontrol(...
                'Units','pixels',...
                'Parent',handles.fig,...
                'Callback',@(varargin)getImageFromSI(hSI, 1),...
                'Position',[tmpX+110 tmpY 40 20],...
                'String','Ch1',...
                'Tag','pbImportCh1');
            handles.pbImportCh2 = uicontrol(...
                'Units','pixels',...
                'Parent',handles.fig,...
                'Callback',@(varargin)getImageFromSI(hSI, 2),...
                'Position',[tmpX+155 tmpY 40 20],...
                'String','Ch2',...
                'Tag','pbImportCh2');
            handles.pbRefIM = uicontrol(...
                'Units','pixels',...
                'Parent',handles.fig,...
                'Callback',@(varargin)getReferenceImage,...
                'Position',[tmpX+208 tmpY+5 60 25],...
                'String','get refIM',...
                'Tag','pbRefIM');
            handles.puRefCh = uicontrol(...
                'Units','pixels',...
                'Style', 'popupmenu',...
                'Parent',handles.fig,...
                'Callback',@(varargin)setCh,...
                'Position',[tmpX+270 tmpY+3 45 25],...
                'String',{'Ch1'; 'Ch2'},...
                'Tag','puRefCh');
            
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

%Right side controls
panelX = 635;tmpX = 2;

panelY = 200; tmpY = 62;
            handles.rightPanel = uipanel('units', 'pixels',  'Position', [panelX panelY 120 480]);
            
            handles.tbZStackPaint = uicontrol(...
                'Units','pixels',...
                'Parent',handles.rightPanel,...
                'style', 'togglebutton',...
                'Callback',@(src,~)(ZStackPaint(src.Value)),...
                'Position',[tmpX tmpY+370 110 20],...
                'String','Z Stack Paint',...
                'Tag','tbZStackPaint');
            
             handles.pbAutoPaint = uicontrol(...
                'Units','pixels',...
                'Parent',handles.rightPanel,...
                'Callback',@(varargin)(autoPaint),...
                'Position',[tmpX tmpY+350 75 20],...
                'String','AutoPaint',...
                'Tag','pbAutoPaint');
             handles.etThreshold =  uicontrol(...
                'Units','pixels',...
                'Parent',handles.rightPanel,...
                'Style','edit',...
                'String', num2str(thresh),...
                'Position',[tmpX+80 tmpY+351 30 18],...
                'callback', @(src, evnt)(updateThresh(str2double(src.String))));

             handles.pbGrow = uicontrol(...
                'Units','pixels',...
                'Parent',handles.rightPanel,...
                'Callback',@(varargin)(grow),...
                'Position',[tmpX tmpY+232 50 20],...
                'String','Grow',...
                'Tag','pbGrow');
             handles.pbShrink = uicontrol(...
                'Units','pixels',...
                'Parent',handles.rightPanel,...
                'Callback',@(varargin)(shrink),...
                'Position',[tmpX+60 tmpY+232 50 20],...
                'String','Shrink',...
                'Tag','pbShrink');
            handles.pbDespeckle = uicontrol(...
                'Units','pixels',...
                'Parent',handles.rightPanel,...
                'Callback',@(varargin)(despeckle),...
                'Position',[tmpX tmpY+210 70 20],...
                'String','Despeckle',...
                'Tag','pbDespeckle');
            handles.etDespeckle = uicontrol(...
                'Units','pixels',...
                'Parent',handles.rightPanel,...
                'Style','edit',...
                'Callback',@(varargin)(setDespeckle),...
                'Position',[tmpX+75 tmpY+210 35 20],...
                'String','30',...
                'Tag','etDespeckle');
            
             handles.tbManualPaint = uicontrol(...
                'Units','pixels',...
                'Parent',handles.rightPanel,...
                'style', 'togglebutton',...
                'Callback',@(src,~)(manualPaint(src.Value)),...
                'Position',[tmpX tmpY+328 65 20],...
                'String','ManualPaint',...
                'Tag','tbManualPaint');
            handles.tbMoveMask = uicontrol(...
                'Units','pixels',...
                'Parent',handles.rightPanel,...
                'style', 'togglebutton',...
                'Callback',@(src,~)(moveMask(src.Value)),...
                'Position',[tmpX+70 tmpY+328 40 20],...
                'String','Move',...
                'Tag','tbMoveMask');
            handles.tbSoma = uicontrol(...
                'Units','pixels',...
                'Parent',handles.rightPanel,...
                'style', 'togglebutton',...
                'Callback',@(src,~)(drawSoma(src.Value)),...
                'Position',[tmpX tmpY+305 65 20],...
                'String','Soma ROIs',...
                'Tag','tbSoma');
            handles.pbROI2Mask = uicontrol(...
                'Units','pixels',...
                'Parent',handles.rightPanel,...
                'Callback',@(varargin)(ROI2Mask),...
                'Position',[tmpX tmpY+283 70 20],...
                'String','ROIs>Mask',...
                'Tag','pbROI2Mask');
            
            uicontrol(...
                'Units','pixels',...
                'Parent',handles.rightPanel,...
                'Style','text',...
                'Position',[tmpX+65 tmpY+301 20 20],...
                'String','R=')
            handles.etCellRadius =uicontrol(...
                'Units','pixels',...
                'Parent',handles.rightPanel,...
                'Style','edit',...
                'Callback',@(varargin)(setRadius),...
                'Position',[tmpX+85 tmpY+305 25 20],...
                'String',int2str(cell_radius),...
                'Tag','etCellRadius');
            function setRadius; cell_radius = max(1,min(100, str2double(get(handles.etCellRadius, 'String')))); end
            
             handles.tbMaxProj = uicontrol(...
                'Units','pixels',...
                'Parent',handles.rightPanel,...
                'style', 'togglebutton',...
                'Callback',@(src, evnt)(setMaxProj(src.Value)),...
                'Position',[tmpX tmpY+255 110 20],...
                'String','MAX Projection',...
                'Tag','tbMaxProj');

            uicontrol(...
                'Units','pixels',...
                'Parent',handles.rightPanel,...
                'Style','text',...
                'Position',[tmpX+20 tmpY+190 70 14],...
                'String','SAVE PATTERN')
            handles.pbSave1 =  uicontrol(...
                'Units','pixels',...
                'Parent',handles.rightPanel,...
                'Position',[tmpX tmpY+160 30 30],...
                'String', '1', ...
                'callback', @(varargin)(savePattern(1)));
            handles.pbSave2 =  uicontrol(...
                'Units','pixels',...
                'Parent',handles.rightPanel,...
                'Position',[tmpX+32 tmpY+160 30 30],...
                'String', '2', ...
                'callback', @(varargin)(savePattern(2)));
            handles.pbSave3 =  uicontrol(...
                'Units','pixels',...
                'Parent',handles.rightPanel,...
                'Position',[tmpX+64 tmpY+160 30 30],...
                'String', '3', ...
                'callback', @(varargin)(savePattern(3)));
            uicontrol(...
                'Units','pixels',...
                'Parent',handles.rightPanel,...
                'Style','text',...
                'Position',[tmpX+20 tmpY+140 70 14],...
                'String','LOAD PATTERN')
            handles.pbLoad1 =  uicontrol(...
                'Units','pixels',...
                'Parent',handles.rightPanel,...
                'Position',[tmpX tmpY+110 30 30],...
                'String', '1', ...
                'callback', @(varargin)(loadPattern(1)));
            handles.pbLoad2 =  uicontrol(...
                'Units','pixels',...
                'Parent',handles.rightPanel,...
                'Position',[tmpX+32 tmpY+110 30 30],...
                'String', '2', ...
                'callback', @(varargin)(loadPattern(2)));
            handles.pbLoad3 =  uicontrol(...
                'Units','pixels',...
                'Parent',handles.rightPanel,...
                'Position',[tmpX+64 tmpY+110 30 30],...
                'String', '3', ...
                'callback', @(varargin)(loadPattern(3)));
            
            
            handles.pbAp =  uicontrol(...
                'Units','pixels',...
                'Parent',handles.rightPanel,...
                'Position',[tmpX tmpY+64 60 20],...
                'String', 'Aperture', ...
                'callback', @(src, evnt)(setAp));
            handles.etAp =  uicontrol(...
                'Units','pixels',...
                'Parent',handles.rightPanel,...
                'Style','edit',...
                'Position',[tmpX+75 tmpY+65 35 18]);
            handles.pbIntersect =  uicontrol(...
                'Units','pixels',...
                'Parent',handles.rightPanel,...
                'Position',[tmpX tmpY+40 110 20],...
                'String', 'Intersect Aperture', ...
                'callback', @(src, evnt)(intersectAp(str2double(get(handles.etAp, 'String')))));
            handles.pbWrite = uicontrol(...
                'Units','pixels',...
                'Parent',handles.rightPanel,...
                'Callback',@(varargin)(writeSLM),...
                'Position',[tmpX tmpY-10 110 40],...
                'String','Write SLM',...
                'Tag','pbWrite');
            
             handles.pbWrite3D = uicontrol(...
                'Units','pixels',...
                'Parent',handles.rightPanel,...
                'Callback',@(varargin)(writeSLM3D),...
                'Position',[tmpX tmpY-60 110 40],...
                'String','Write SLM 3D',...
                'Tag','pbWrite');
            
set(handles.fig, 'ResizeFcn', @resizeFcn);
            
%import image if available
if ~isempty(hSI)
    try
        getImageFromSI(hSI,2); 
    catch
        disp('SLMcontrol: No SI image to get')
    end
end

handles.etAp.String = '2'; setAp; %- for sandbox; comment out setAP

            
%Figure control
    function closeFcn(src, evnt)
        delete(src)
    end
function updateThresh(val)
    thresh = max(0,min(intmax('uint16'),val));
    set(handles.etThreshold, 'String', num2str(thresh));
    autoPaint;
end

%ZSTACKPAINT
    function ZStackPaint(val)
        if ~val %TURN OFF ZSTACKPAINT
            do3D = false;
            handles.masks = zeros(Xsz,Ysz,Zsz);
            delete(handles.mask);
            mask = cat(3, zeros(Xsz,Ysz), ones(Xsz,Ysz), ones(Xsz,Ysz));
            handles.mask = imshow(mask, 'parent', handles.ax1);
            set(handles.mask, 'alphadata', zeros(Xsz,Ysz), 'hittest', 'off');
            return
        end
        do3D = true;
        
        % read labels

        basedir='E:\SLAPMidata\SLM_segmentation\';
        [labelFn, labelDr] = uigetfile([basedir '*.h5'], 'Select an ilastik output');
        token=strsplit(labelFn,'_');
        bin=str2double(token{end-1}(end));
        dendrite=1;%label 1: dendrite, label 2: soma, label 3: bkg
        
        handles.masks=readlabel([labelDr labelFn],dendrite,bin); 
        CurrentMask=handles.masks(:,:,Zpos);
        set(handles.mask, 'alphadata', CurrentMask*0.3);
    end
    
    function im=readlabel(labelspath,targetlabel,bin)
        labels=h5read(labelspath, '/exported_data');
        im=squeeze(labels);
        im=(im==targetlabel);
        im=imresize(im,bin);
        se = strel('disk',5);
        for frame=1:size(im,3)
            im(:,:,frame)=imclose(bwareafilt(im(:,:,frame),[50 inf]),se);%despeckle and close
            %im(:,:,frame)=imclose(im(:,:,frame),se);
        end    
    end

%AUTOPAINT
    function autoPaint
        im = get(handles.screen, 'CData');
        newMask= autoPaintThresh(im, thresh);
        set(handles.mask, 'alphadata', newMask*alpha);
        updateMasks;

    end
    function imthresh = autoPaintThresh(im, thresh)
        sz = 5;
        SE = strel('disk', sz);
        
        for z = 1:size(im, 3)
            im(:,:,z) = ordfilt2(im(:,:,z),8,true(3));
        end
        immax = max(im,[],3);
        immax = imtophat(immax, ones(31));
        thresh = thresh + median(immax(:));
        
        imthresh = immax>thresh;
        imthresh = imdilate(imthresh, SE);
        imthresh = bwareafilt(imthresh, [(2*sz).^2 inf]);
        imthresh = imopen(imclose(imthresh,SE), SE);
    end
%MANUAL PAINT
%Manual Paint
    function manualPaint(val)
        if ~val %TURN OFF MANUAL PAINT MODE
            %set the cursor back
            set(handles.mask, 'hittest', 'off');
            set(handles.fig, 'WindowButtonDownFcn', '');
            set(handles.fig, 'WindowButtonUpFcn', '');
            set(handles.fig, 'WindowButtonMotionFcn', '');
            return
        end
        
        handles.tbMoveMask.Value = 0;
        handles.tbSoma.Value = 0;
        
        %make the mask clickable
        set(handles.mask, 'hittest', 'on');       
        %initialize variables
        handles.Paint.startpoint = [];
        %set the buttondownfunction
        set(handles.fig, 'WindowButtonDownFcn', @manualpaint_BDF);
        set(handles.fig, 'WindowButtonUpFcn', @manualpaint_BUF);
        set(handles.fig, 'WindowButtonMotionFcn', @manualpaint_BMF);
    end
function manualpaint_BDF (src, event)
cp = get(handles.ax1,'CurrentPoint');
cp = cp(1, 1:2);

%If the point is within the axes
if all(cp>=1 & cp<=[Ysz Xsz])
    mod = get(src, 'SelectionType');
    switch mod
        case 'extend'
            %paint a circle on current point
            M = get(handles.mask, 'Alphadata');
            [XX,YY] = meshgrid(1:Xsz,1:Ysz);
            M(sqrt((XX-cp(1)).^2 + (YY-cp(2)).^2)<5) = alpha;
            set(handles.mask, 'Alphadata', M);
            drawing = true;
            
            updateMasks;
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
                M = get(handles.mask, 'Alphadata');
                mod = get(src, 'SelectionType');
                switch mod
                    case 'normal'
                        %add the pixels to the mask
                        M(rect) = alpha;
                    case 'alt'
                        %remove pixels from the mask
                        M(rect) = 0;
                end
                set(handles.mask, 'alphadata', M);
                
                updateMasks;
            end
        end
        handles.Paint.startpoint = [];
        set(handles.Paint.selection, 'Xdata', -1, 'Ydata', -1);
    end

    function manualpaint_BMF (src, event)
        mod = get(src, 'SelectionType');
        cp = get(handles.ax1,'CurrentPoint');
        cp = cp(1,1:2);
        switch mod
            case 'extend'
                if drawing
                    M = get(handles.mask, 'Alphadata');
                    [XX,YY] = meshgrid(1:Xsz,1:Ysz);
                    M(sqrt((XX-cp(1)).^2 + (YY-cp(2)).^2)<5) = alpha;
                    set(handles.mask, 'Alphadata', M);
                    
                    updateMasks;
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
    function moveMask(val)
        if ~val %TURN OFF MOVE MODE
            %set the cursor back
            set(handles.mask, 'hittest', 'off');
            set(handles.fig, 'WindowButtonDownFcn', '');
            set(handles.fig, 'WindowButtonUpFcn', '');
            set(handles.fig, 'WindowButtonMotionFcn', '');
            return
        end
        
        %turn off competing toggles
        handles.tbManualPaint.Value = 0;
        handles.tbSoma.Value = 0;
        
        %make the mask clickable
        set(handles.mask, 'hittest', 'on');       
        %initialize variables
        handles.Move.startpoint = [];
        %set the buttondownfunction
        set(handles.fig, 'WindowButtonDownFcn', @moveMask_BDF);
        set(handles.fig, 'WindowButtonUpFcn', @moveMask_BUF);
        set(handles.fig, 'WindowButtonMotionFcn', @moveMask_BMF);
    end

function moveMask_BDF (src, event)
cp = get(handles.ax1,'CurrentPoint');
cp = cp(1, 1:2);
%If the point is within the axes
if all(cp>=1 & cp<=[Ysz Xsz])
    handles.Move.startpoint = cp;
    handles.Move.IM = get(handles.mask, 'alphadata');
end
end
    function moveMask_BUF (src, event)
        cp = get(handles.ax1,'CurrentPoint');
        cp = cp(1,1:2);
        if ~isempty(handles.Move.startpoint)
            %set the alphadata to the translation of the existing alphadata
            D = cp-handles.Move.startpoint;
            set(handles.mask, 'alphadata', imtranslate(handles.Move.IM, D));
        end
        handles.Move.startpoint = [];
    end

    function moveMask_BMF (src, event)
        cp = get(handles.ax1,'CurrentPoint');
        cp = cp(1,1:2);
        if ~isempty(handles.Move.startpoint);
            D = cp-handles.Move.startpoint;
            set(handles.mask, 'alphadata', imtranslate(handles.Move.IM, D));
        end
    end

    function drawSoma(val)
        set(handles.fig, 'WindowButtonDownFcn', '');
        set(handles.fig, 'WindowButtonUpFcn', '');
        set(handles.fig, 'WindowButtonMotionFcn', '');
        if ~val %TURN OFF SOMA MODE
            return
        end
        
        %turn off competing toggles
        handles.tbManualPaint.Value = 0;
        handles.tbMoveMask.Value = 0;
        
        %set the buttondownfunction
        set(handles.fig, 'WindowButtonDownFcn', @drawSoma_BDF);
    end
    function drawSoma_BDF(src, event)
        cp = get(handles.ax1,'CurrentPoint');
        cp = cp(1,1:2);
        if all(cp>=1 & cp<=[Ysz Xsz])
            mod = get(src, 'SelectionType');
            switch mod
                case 'normal' %left click to add
                    if isempty(ROIs) || min(sum((cell2mat({ROIs.center}') - repmat(cp, length(ROIs),1)).^2,2))>10
                        roiID = length(ROIs)+1;
                        ROIs(roiID,1).hScatter = scatter(cp(1),cp(2), '.r', 'linewidth', 2, 'parent', handles.ax1, 'userdata', 'ROI');
                        [ROIs(roiID).pixel_list,ROIs(roiID).nuc_list,ROIs(roiID).boundary_inner,ROIs(roiID).boundary_outer,ROIs(roiID).boundary_list]=donut_roi(get(handles.screen, 'CData'),cp([2 1]),cell_radius,[], false);
                        ROIs(roiID,1).hLine = plot(ROIs(roiID).boundary_outer(1,:), ROIs(roiID).boundary_outer(2,:), 'r.', 'parent', handles.ax1, 'userdata', 'ROI');
                        ROIs(roiID,1).center = cp;
                    end
                case 'alt'
                    d = sum((cell2mat({ROIs.center}') - repmat(cp, length(ROIs),1)).^2,2);
                    [minval, roiID] = min(d);
                    if minval<20 %delete nearest ROI
                        delete(ROIs(roiID).hScatter);
                        delete(ROIs(roiID).hLine);
                        ROIs(roiID) = [];
                    end

                otherwise
            end
        end
    end


    function redrawROIs
        objects = findobj(handles.ax1, 'userdata', 'ROI');
        try
            delete(objects);
        end
        for roiID = 1:length(ROIs)
            cp = ROIs(roiID).center;
            ROIs(roiID,1).hScatter = scatter(cp(1),cp(2), '.r', 'linewidth', 2, 'parent', handles.ax1, 'userdata', 'ROI');
            ROIs(roiID,1).hLine = plot(ROIs(roiID).boundary_outer(1,:), ROIs(roiID).boundary_outer(2,:), 'r.', 'parent', handles.ax1, 'userdata', 'ROI');
        end
    end


    function getImageFromSI(hSI, channel)
        if isempty(hSI)
            disp('No scanimage variable!');
            return
        end
        Zsz = length(hSI.hDisplay.rollingStripeDataBuffer);
        channelsAcquired = hSI.hDisplay.rollingStripeDataBuffer{1}{1}.roiData{1}.channels;
        [tf,chIdx] = ismember(channel,channelsAcquired);
        assert(tf,'Channel %d was not acquired',channel);
        sliceSz = size(hSI.hDisplay.rollingStripeDataBuffer{1}{1}.roiData{1}.imageData{chIdx}{1});
        Xsz = sliceSz(1); Ysz = sliceSz(2);
        im_SI = nan(Xsz,Ysz,Zsz);
        for sliceIdx = 1:Zsz
            % internally, ScanImage stores image data transposed (row-major order like in C)
            % transpose here to be consistent with Matlab's column-major order
            im_SI(:,:,sliceIdx) = hSI.hDisplay.rollingStripeDataBuffer{sliceIdx}{1}.roiData{1}.imageData{chIdx}{1}';
            %     zs(sliceIdx) = z;
            z = hSI.hDisplay.rollingStripeDataBuffer{sliceIdx}{1}.roiData{1}.zs(1);
            scanfield = hSI.hDisplay.rollingStripeDataBuffer{sliceIdx}{1}.roiData{1}.hRoi.get(z);
        end
        
        imBuffer = im_SI;
        Xsz = size(imBuffer,1); Ysz = size(imBuffer,2); Zsz = size(imBuffer,3);
        updateZ(Zsz/2);
        maxContrast =  max(imBuffer(:));
        set(handles.slContrast, 'Max' ,maxContrast);
        setContrast(maxContrast/2); 
        set(handles.ax1, 'xlim', 0.5+[0 Xsz], 'ylim', 0.5+[0 Ysz])
        setContrast(hSI.hChannels.channelLUT{channel}(2));      %prctile(imBuffer(:), 99));
        set(handles.slZ, 'Max', Zsz, 'Value', round(Zsz/2), 'sliderstep', [1/Zsz 1/Zsz]);
        updateCData;
    end

    function getReferenceImage
        basedir = 'E:\SLAPMiData';
        if ~isempty(hSI)
            basedir = hSI.hScan2D.logFilePath;
        end
        [fns, dr] = uigetfile({'*.mat;*.tif'},'Select the Tiff files in the stack, or your aligned REF IM',basedir, 'multiselect', 'on');
        if iscell(fns) %a stack of images was selected
            opts.fns = fns; opts.dr = dr;
            refIM = SLAPMi_refIM(opts);
        elseif length(fns)>4 && strcmpi(fns(end-3:end), '.mat')%a single reference image was selected
            load([dr filesep fns])
        else
            warning('Must select multiple tiff files or a single .mat file.') 
            return
        end
        
        %update the buffer and the scanfield, etc with the reference image
        channelSelect = min(size(refIM.data,4), channel); 
        imBuffer = squeeze(refIM.data(:,:,:,channelSelect));
        imBuffer = permute(imBuffer, [2 1 3]);
        scanfield = refIM.metadata.rois.RoiGroups.imagingRoiGroup.rois{1, 1}.scanfields(1);
        scanfield.meshgrid = @(varargin)(meshgrid(refIM.M.coords.X, refIM.M.coords.Y));
        
        %update display
        Xsz = size(imBuffer,1); Ysz = size(imBuffer,2); Zsz = size(imBuffer,3);
        updateZ(Zsz/2);
        set(handles.ax1, 'xlim', 0.5+[0 Xsz], 'ylim', 0.5+[0 Ysz]);
        maxContrast =  max(imBuffer(:));
        set(handles.slContrast, 'Max' ,maxContrast);
        setContrast(maxContrast/2); 
        set(handles.slZ, 'Max', Zsz, 'Value', round(Zsz/2), 'sliderstep', [1/Zsz 1/Zsz]);
        updateCData;
    end

function updateMasks
    if do3D
        M=get(handles.mask, 'alphadata')>0;
        handles.masks(:,:,Zpos)=M;
    end
end

function updateCData
    if doMax
        I = max(imBuffer,[],3);
    else
        I = imBuffer(:,:,Zpos);
    end
    set(handles.screen, 'CData', I);
    if isempty(handles.mask) || ~all(size(get(handles.mask, 'alphadata')) == [Xsz Ysz])
        delete(handles.mask);
        mask = cat(3, zeros(Xsz,Ysz), ones(Xsz,Ysz), ones(Xsz,Ysz));
        handles.mask = imshow(mask, 'parent', handles.ax1);
        set(handles.mask, 'alphadata', zeros(Xsz,Ysz), 'hittest', 'off');
    end
    
    % initialize masks (for 3D)
    if do3D
        CurrentMask=handles.masks(:,:,Zpos);
        set(handles.mask, 'alphadata', CurrentMask*0.3);
    else
        masks = zeros(Xsz,Ysz,Zsz);
        handles.masks=masks;
    end
end

function updateZ(val)
    Zpos = round(max(1,min(Zsz, val)));
    set(handles.slZ, 'Value', Zpos);
    updateCData;
end
    function setContrast(value)
        contrast = max(-5000, min(get(handles.slContrast, 'Max'), value));
        set(handles.etContrast, 'String', int2str(contrast));
        set(handles.slContrast, 'Value', contrast)
        set(handles.ax1, 'Clim', [min(0, contrast-2) contrast]);
    end
    function setMaxProj(value)
        doMax = logical(value);
        updateCData;
    end

    function setAp(val)
        if ~nargin
            val = str2double(handles.etAp.String);
        end
        pattern = SLM_drawPattern(val);
        if ~isempty(hSI) && ~isempty(scanfield)
            I = mapSLMToImage(pattern, hSI.hSlmScan.testPatternPixelToRefTransform, scanfield);
            set(handles.mask, 'alphadata', I*alpha)
        end
    end
    function intersectAp(val)
        if isnumeric(val) && val>0 && ~isempty(hSI)
            radius = val;
            [XX,YY] = meshgrid(linspace(-1,1,512));
            pattern = sqrt(XX.^2 + YY.^2)<radius;
            select = mapSLMToImage(pattern, hSI.hSlmScan.testPatternPixelToRefTransform, scanfield);
            M = get(handles.mask, 'alphadata')>0;
            M = M & select;
            set(handles.mask,'alphadata', alpha*double(M));
        end
    end
    function scrollFcn(src, evnt)
        updateZ(Zpos - evnt.VerticalScrollCount);
    end
    function writeSLM
        m_SI = get(handles.mask, 'alphadata');
        m_SI = double(m_SI>0);
        m_SLM = logical(mapImageToSLM(m_SI, hSI.hSlmScan.testPatternPixelToRefTransform, scanfield));
        SLM_drawPattern(m_SLM);
        if ~isempty(hSLAPMi)
            assignin('base', 'ROIs', ROIs);
            evalin('base', 'hSLAPMi.SLM.ROIs = ROIs');
        end

    end
    function writeSLM3D
            %m_SI=handles.masks;
            [XX,YY]=meshgrid(linspace(-1,1,512),linspace(-1,1,512));
            pattern=cat(3,sqrt(XX.^2 + YY.^2)<2,sqrt(XX.^2 + YY.^2)<0);
            %load('W:\JJ\SLM measurements\Masks\digitMasks\sequential\pattern.mat','pattern');
            m_SI=repmat(pattern,[1 1 50]);
            m_SLM=m_SI;
     %         for i=1:size(m_SI,3)
     %             m_SLM(:,:,i)=logical(mapImageToSLM(m_SI(:,:,i), hSI.hSlmScan.testPatternPixelToRefTransform, scanfield));
     %         end
            m_SLM=logical(m_SLM);

            SLM_drawPattern3D(m_SLM,true);
    end

    function grow
        S = strel('disk', 3);
        M = get(handles.mask, 'alphadata')>0;
        M = imdilate(M,S);
        set(handles.mask,'alphadata', alpha*double(M));
        updateMasks;


    end
    function shrink
        S = strel('disk', 5);
        M = get(handles.mask, 'alphadata')>0;
        M = imerode(M,S);
        set(handles.mask,'alphadata', alpha*double(M));

        updateMasks;

    end
    function despeckle
        M = get(handles.mask, 'alphadata')>0;
        M2 = bwareafilt(M,[str2double(get(handles.etDespeckle, 'String'))^2 inf]);
        set(handles.mask,'alphadata', alpha*double(M2));
        updateMasks;

    end
    function setDespeckle
        S = str2double(get(handles.etDespeckle, 'String')); 
        S = round(max(S, 1));
        set(handles.etDespeckle, 'String', int2str(S));
    end

    function savePattern(n)
       fn = [saveDir filesep 'mask' int2str(n)];
       m_SI = get(handles.mask, 'alphadata')>0;
       m_SI = double(m_SI>0);
       if ~isempty(hSI)
            m_SLM = logical(mapImageToSLM(m_SI, hSI.hSlmScan.testPatternPixelToRefTransform, scanfield)); % #ok<NASGU>
       else
            m_SLM =m_SI; % #ok<NASGU>
       end
       save(fn, 'm_SI', 'm_SLM', 'ROIs'); 
    end
    function loadPattern(n)
        fn = [saveDir filesep 'mask' int2str(n)];
        if exist([fn '.mat'], 'file')
            S = load(fn);
            m_SLM = S.m_SLM;
            m_SI = S.m_SI;
            ROIs = S.ROIs;
            redrawROIs;
            if size(m_SI) == [Xsz Ysz]
                set(handles.mask, 'alphadata', alpha*double(m_SI));
            elseif ~isempty(hSI)
                I = mapSLMToImage(m_SLM, hSI.hSlmScan.testPatternPixelToRefTransform, scanfield);
                I(isnan(I)) = 0;
                set(handles.mask, 'alphadata', I*alpha)
            else
                disp('There was an issue loading the SLM pattern');
            end
        end
    end

    function setCh
       channel = handles.puRefCh.Value; 
    end
    function exit(src, evnt)
        if strcmp(questdlg('Exit SLMcontrol?', 'SLMcontrol', 'Yes', 'No', 'No'), 'Yes')
            delete(handles.fig);
        end
    end

    function resizeFcn(src, evnt)
        p = get(handles.fig, 'pos');
        w = p(3); h = p(4);
        
        rPpos = get(handles.rightPanel, 'pos');
        set(handles.rightPanel, 'pos', [w-rPpos(3) (h-rPpos(4))/2 rPpos(3) rPpos(4)]);
        
        axpos = get(handles.ax1, 'pos');
        set(handles.ax1, 'pos', [axpos(1) axpos(2) w-axpos(1)-rPpos(3)-15 h-axpos(2)])
        
        set(handles.slZ, 'pos', [w-rPpos(3)-15 axpos(2) 15 h-axpos(2)])
    end

    function ROI2Mask
        M = get(handles.mask, 'alphadata')>eps;
        for roiID =1:length(ROIs)
            M(ROIs(roiID).pixel_list) = true;
        end
        set(handles.mask, 'alphadata', alpha*M);
    end
end




