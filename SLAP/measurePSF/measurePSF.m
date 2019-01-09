function measurePSF

%A GUI for aligning the SLAPMi scope and measuring the illumination PSF across the field
%of view

global calib
global hSLAPMi
global PSFalign
PSFalign = 1;

warning('off', 'imaq:pointgrey:roiChanged')
warning('off', 'MATLAB:timer:deleterunning')
if ~most.idioms.isValidObj(hSLAPMi)
    hSLAPMi = slapmi;
end

%PROPERTIES
basedir = 'E:\SLAPmidata\';
DZ = 5;
Zcenter = 200;
Z = 200;
hrect = [];
drawing = false;
scrolling = false;
generating = false; %true when generate_projections is running
FZ = [];
F1 = [];
h_im = [];
AX = [];
doScan = false;
galvoX = 0;
galvoY = 0;
calib = hSLAPMi.calib;
scantimer = [];
powerwarned = false;
Ap = 0.5;
rPOW = 0.04;
POW = 0.12; %line power. This cannot be very low! The line profile changes as power increases, it is stable above ~0.2


%Field of view plot
XYcenter = [950,600];
XYcorner = [950,500];
hFOVcenter = [];
hFOVcorner = [];
hFOV = [];
hFOVdots = [];
hFOVtext = [];
data = [];

%piezo parameters
steptime = 0.1;
AOrate = 10000;
N_AO = steptime*AOrate;

galvoRate = 1000000; %sample rate to galvos

%SETUP DAQ
hTaskPiezo = most.util.safeCreateTask('SLAPMIPiezo');
hTaskPiezo.createAOVoltageChan(hSLAPMi.piezoDaqName, 0);
hTaskPiezo.cfgSampClkTiming(10000, 'DAQmx_Val_FiniteSamps', N_AO);
hTaskPiezo.cfgOutputBuffer(N_AO);               

hTaskGalvoBeams = most.util.safeCreateTask('SLAPMIGalvosBeams');
hTaskGalvoBeams.createAOVoltageChan(hSLAPMi.scannerDaqName, [0 4 8 12 16 20 24]);

hTaskPowerBeam = most.util.safeCreateTask('SLAPMIPowerBeam');
hTaskPowerBeam.createAOVoltageChan(hSLAPMi.beamDaqName, 0); %Power for lines
hTaskPowerBeam.createAOVoltageChan(hSLAPMi.beamDaqName, 1); %Power for raster scan

hTaskFlipMirror = most.util.safeCreateTask('FlipMirror');
hTaskFlipMirror.createDOChan(hSLAPMi.beamDaqName, 'port0/line0' , 'RasterFlip');
hTaskFlipMirror.writeDigitalData(false, 1, true); % FALSE for line scan

hTaskShutter = most.util.safeCreateTask('Shutter');
hTaskShutter.createDOChan(hSLAPMi.piezoDaqName, 'port0/line7' , 'ShutterTTL');
hTaskShutter.writeDigitalData(true, 1, true); % OPEN


%SETUP VIDEO

vid = videoinput('pointgrey', 1, 'F7_Mono16_1920x1200_Mode7');
src = getselectedsource(vid);
src.Gamma = 1;
setSrcFor('preview');
triggerconfig(vid, 'manual');
vid.ROIPosition = [352 0 1216 1200];

F2 = figure('Toolbar','none','CloseRequestFcn', @F2close,...
    'Menubar', 'none',...
    'NumberTitle','Off',...
    'Name','Line Alignment Tool',...
    'pos',[1267 468 100  425], 'resize', 'off', 'windowbuttonupfcn', @(varargin)stopscroll);

%aperture controls
uicontrol(...
    'Units','Points',...
    'Parent',F2,...
    'Position', [5 275 40 20],...
    'HorizontalAlignment','left',...
    'String','Aperture',...
    'Style','text');
etAp = uicontrol(...
    'Units','Points',...
    'Parent',F2,...
    'Position', [40 280 45 20],...
    'HorizontalAlignment','center',...
    'String', num2str(Ap),...
    'Style','edit',...
    'Tag','etAp',...
    'callback',@(src, ~)setAperture(str2double(get(src,'string'))));

%generate projections
pbGenP =  uicontrol(...
    'Units','Points',...
    'Parent',F2,...
    'Callback',@(varargin)generate_projections,...
    'Position',[5 258 80 20],...
    'String','Generate P',...
    'Tag','pbGenP');

%ROI CONTROLS
pbROIset =  uicontrol(...
    'Units','Points',...
    'Parent',F2,...
    'Callback',@(varargin)drawROI,...
    'Position',[5 237 40 20],...
    'String','ROI Set',...
    'Tag','pbCopy1');
pbROIreset =  uicontrol(...
    'Units','Points',...
    'Parent',F2,...
    'Callback',@(varargin)resetROI,...
    'Position',[45 237 40 20],...
    'String','Reset',...
    'Tag','pbCopy1');
pbCalI1 =  uicontrol(...
    'Units','Points',...
    'Parent',F2,...
    'Callback',@(varargin)calibrateBcurve,...
    'Position',[5 215 25 20],...
    'String','EOM',...
    'Tag','pbCopy1');
pbCalI2 =  uicontrol(...
    'Units','Points',...
    'Parent',F2,...
    'Callback',@(varargin)calibrateIntensity,...
    'Position',[30 215 55 20],...
    'String','Line Powers',...
    'Tag','pbCopy1');

%Z CENTER CONTROLS
uicontrol(...
    'Units','Points',...
    'Parent',F2,...
    'Position', [12 183 10 20],...
    'HorizontalAlignment','left',...
    'String','Z',...
    'Style','text');
etZCenter = uicontrol(...
    'Units','Points',...
    'Parent',F2,...
    'Position', [30 188 30 20],...
    'HorizontalAlignment','center',...
    'String', num2str(Zcenter),...
    'Style','edit',...
    'Tag','etZCenter',...
    'callback',@(src, ~)moveZcenter(str2double(get(src,'string'))));
pbZup = uicontrol(...
    'Units','Points',...
    'Parent',F2,...
    'Position', [65 198 15 15],...
    'HorizontalAlignment','left',...
    'String','^',...
    'Style','pushbutton',...
    'Tag','pbZup', ...
    'callback',@(src, ~)moveZcenter([], +2));
pbZdn = uicontrol(...
    'Units','Points',...
    'Parent',F2,...
    'Position', [65 183 15 15],...
    'HorizontalAlignment','left',...
    'String','v',...
    'Style','pushbutton',...
    'Tag','pbZdn',...
    'callback',@(src, ~)moveZcenter([], -2));

%dZ CONTROLS
cbDoZ = uicontrol(...
    'Units','Points',...
    'Parent',F2,...
    'Position', [5 160 40 20],...
    'HorizontalAlignment','center',...
    'String','Z scan',...
    'Style','checkbox',...
    'Tag','cbDoZ',...
    'callback',@(varargin)toggleZScan);
etDZ = uicontrol(...
    'Units','Points',...
    'Parent',F2,...
    'Position', [50 160 30 20],...
    'HorizontalAlignment','center',...
    'String', num2str(DZ),...
    'Style','edit',...
    'Tag','etDZ',...
    'callback',@(src, ~)(updateDZ(str2double(get(src,'string')))));

%GALVO POSITION CONTROLS
uicontrol(...
    'Units','Points',...
    'Parent',F2,...
    'Position', [5 133 30 20],...
    'HorizontalAlignment','left',...
    'String','Beam',...
    'Style','text');
puBeam = uicontrol(...
    'Units','Points',...
    'Parent',F2,...
    'Position', [33 135 50 20],...
    'HorizontalAlignment','center',...
    'String', {'Line 1','Line 2','Line 3','Line 4','Raster Beam', 'Off'},...
    'Style','popup',...
    'Tag','etDZ',...
    'callback',@(varargin)updateBeam);
uicontrol('Units','Points','Parent',F2,'Position', [14 120 30 15],...
    'HorizontalAlignment','left','String','X','Style','text');
uicontrol('Units','Points','Parent',F2,'Position', [55 120 30 15],...
    'HorizontalAlignment','left','String','Y','Style','text');
etX = uicontrol('Units','Points','Parent',F2,'HorizontalAlignment','center',...
    'Position', [5 110 25 15],...
    'String', num2str(galvoX),...
    'Style','edit',...
    'Tag','etX',...
    'callback',@(varargin)updateGalvo);
pbXup = uicontrol(...
    'Units','Points',...
    'Parent',F2,...
    'Position', [30 120 15 15],...
    'HorizontalAlignment','left',...
    'String','^',...
    'Style','pushbutton',...
    'Tag','pbXup','Enable', 'off',...
    'buttondownfcn', @(varargin)bdXup);
pbXdn = uicontrol(...
    'Units','Points',...
    'Parent',F2,...
    'Position', [30 105  15 15],...
    'HorizontalAlignment','left',...
    'String','v',...
    'Style','pushbutton',...
    'Tag','pbXdn','Enable', 'off',...
    'ButtonDownFcn', @(varargin)bdXdn);

etY = uicontrol('Units','Points','Parent',F2,'HorizontalAlignment','center',...
    'Position', [47 110 25 15],...
    'String', num2str(galvoY),...
    'Style','edit',...
    'Tag','etX',...
    'callback',@(varargin)updateGalvo);
pbYup = uicontrol(...
    'Units','Points',...
    'Parent',F2,...
    'Position', [72 120 15 15],...
    'HorizontalAlignment','left',...
    'String','^',...
    'Style','pushbutton',...
    'Tag','pbYup','Enable', 'off',...
    'ButtonDownFcn', @(varargin)bdYup);
pbYdn = uicontrol(...
    'Units','Points',...
    'Parent',F2,...
    'Position', [72 105  15 15],...
    'HorizontalAlignment','left',...
    'String','v',...
    'Style','pushbutton',...
    'Tag','pbYdn','Enable', 'off',...
    'ButtonDownFcn', @(varargin)bdYdn);

%SHOW FOV?
cbShowFOV = uicontrol(...
    'Units','Points','Value', 1,...
    'Parent',F2,...
    'Position', [5 60 70 15],...
    'HorizontalAlignment','center',...
    'String','FOV Guide',...
    'Style','checkbox',...
    'Tag','cbDoZ',...
    'callback',@(src,~)FOVGuideVis(src.Value));

%SCAN TOGGLE
tbScan = uicontrol(...
    'Units','Points',...
    'Parent',F2,...
    'Position', [5 82 25 20],...
    'HorizontalAlignment','left',...
    'String','Scan',...
    'Style','togglebutton',...
    'Tag','pbZdn',...
    'callback',@(varargin)toggleScan);

%SET OFFSET BUTTON
pbSet0 = uicontrol(...
    'Units','Points',...
    'Parent',F2,...
    'Position', [35 82 25 20],...
    'HorizontalAlignment','left',...
    'String','Set 0',...
    'Style','pushbutton',...
    'Tag','pbSet0',...
    'callback',@(varargin)setOffset);

%SAVE and LOAD BUTTONS
pbSave = uicontrol(...
    'Units','Points',...
    'Parent',F2,...
    'Position', [65 82 25 20],...
    'HorizontalAlignment','left',...
    'String','Save',...
    'Style','pushbutton',...
    'Tag','pbZdn',...
    'callback',@(varargin)saveCalib);
pbLoad = uicontrol(...
    'Units','Points',...
    'Parent',F2,...
    'Position', [65 57 25 20],...
    'HorizontalAlignment','left',...
    'String','Load',...
    'Style','pushbutton',...
    'Tag','pbZdn',...
    'callback',@(varargin)loadCalib);

%POWER
uicontrol(...
    'Units','Points',...
    'Parent',F2,...
    'Position', [5 43 30 15],...
    'HorizontalAlignment','left',...
    'String','Power',...
    'Style','text');
etPower = uicontrol(...
    'Units','Points', ...
    'Parent',F2,...
    'Position', [7 33 20 15],...
    'String', '0.1',...
    'Style','edit',...
    'Tag','etPower',...
    'callback',@(src, ~)updatePower(max(0, min(1, str2double(get(src,'String'))))));
sPower = uicontrol(...
    'Units','Points', 'Min', 0, 'Max', 1, 'Value', 0.1,...
    'Parent',F2,...
    'Position', [35 35 50 15],...
    'Style','slider',...
    'Tag','sContrast',...
    'callback',@(src, ~)updatePower(get(src,'Value')));

%CONTRAST
uicontrol(...
    'Units','Points',...
    'Parent',F2,...
    'Position', [25 12 40 20],...
    'HorizontalAlignment','left',...
    'String','Contrast',...
    'Style','text');
sContrast = uicontrol(...
    'Units','Points', 'Min', 0, 'Max', 1, 'Value', 1,...
    'Parent',F2,...
    'Position', [5 5 80 15],...
    'Style','slider',...
    'Tag','sContrast',...
    'callback',@(src, ~)updateContrast(get(src,'Value')));


%initialize
updatePower(get(sPower, 'Value'))
ZScanOpen = true; toggleZScan(); %This switches us into Preview mode
moveZ(Z) %position the Piezo

%update the state to the raster beam
puBeam.Value = 5;
updateBeam;

    function kpf(event, keydata)
        switch keydata.Key
            case 'return'
                if drawing %check if drawing ROI
                    newpos = getPosition(hrect);
                    oldpos = vid.ROIPosition;
                    
                    %sanity checks
                    setpos = newpos + [oldpos(1:2) 0 0];
                    setpos = max([0 0 0 0], min([1900 1200 1900 1200], setpos));
                    setpos(3:4) = min(setpos(3:4), [1900 1200]-setpos(1:2));
                    
                    vid.ROIPosition = setpos;
                    set(AX, 'xlim', [0 setpos(3)])
                    set(AX, 'ylim', [0 setpos(4)])
                    
                    drawing = false;
                    delete(hrect)
                    drawFOV
                end
            case 'k'
                keyboard
            case 'a'
                imageAperture;
            case 'o'
                resetOffsets;
            case 'z'
                findZcenter;
        end
    end

    function drawROI
        if drawing
            drawing = false;
            delete(hrect)
        else
            drawing = true;
            hrect = imrect(AX);
        end
    end

    function resetROI
        newpos = [352 0 1216 1200];
        vid.ROIPosition = newpos;
        set(AX, 'xlim', [0 newpos(3)])
        set(AX, 'ylim', [0 newpos(4)])
        drawFOV
    end

    function toggleZScan
        if ~ZScanOpen
            %close the preview window
            delete(F1)
            F1 = [];
            
            ZScanOpen = true;
            FZ = figure('Toolbar','none',...
                'Menubar', 'none',...
                'NumberTitle','Off',...
                'Name','Z Scan');
            ROIpos = vid.ROIPosition;
            h_im = imshow(zeros(ROIpos(3), 3*ROIpos(4)), [0 65535], 'border', 'tight');
            AX = get(h_im, 'parent');
            colormap(AX, [0 0 0.3 ; repmat(0:1/252:1, 3, 1)' ; 1 0 0]);
            
            h_plot1 = figure; h_ax1 = axes('parent', h_plot1); %plot of intensity along the line at each Z
            
            %take a picture and Hough transform to get line axis
            Zpositions = Zcenter-DZ : Zcenter+DZ;
            set(h_ax1, 'ColorOrder', jet(length(Zpositions)), 'NextPlot', 'replacechildren');
            
            moveZ(Zcenter-DZ);
            
            setSrcFor('Zscan');
            vid.FramesPerTrigger = length(Zpositions);
            Ztrace = smooth([ones(1,100)*Zcenter-DZ ...
                        linspace(Zcenter-DZ, Zcenter+DZ, AOrate*(length(Zpositions)/src.FrameRate)) ...
                        MAA_curve(Zcenter+DZ,Zcenter-DZ, 0, 0,1, AOrate*steptime)],21);
            P = nan(0, length(Zpositions));
            
            stop(hTaskPiezo)
            hTaskPiezo.cfgOutputBuffer(length(Ztrace));
            hTaskPiezo.cfgSampClkTiming(10000, 'DAQmx_Val_FiniteSamps', length(Ztrace));
            hTaskPiezo.writeAnalogData(Ztrace(:)*10/400);
            
            theta = [];
            while true
                %take a picture at the top of the range
                %for f = 1:length(Zpositions)
                stop(hTaskPiezo)
                start(vid)
                start(hTaskPiezo)
                trigger(vid);
                
                try
                    frames = double(getdata(vid));
                catch
                    delete(FZ);
                    FZ = [];
                    return
                end
                if isempty(theta) %estimate angle
                    theta = getangle(sum(frames,4), puBeam.Value);
                end
                for f = 1:length(Zpositions)
                    I = imrotate(frames(:,:,1,f),theta, 'bilinear');
                    [~, maxrow] = max(sum(I,1));
                    
                    if maxrow<21 || maxrow>size(I,2)-21
                        P(:,f) = 0;
                    else
                        P(1:size(I,1),f) = smooth(sum(I(:, maxrow+[-20:20]),2));
                    end
                end
                if ZScanOpen
                    set(h_im, 'CData', [frames(:,:,1,1) frames(:,:,1,ceil(end/2)) frames(:,:,1,end)])
                    plot(h_ax1, P);
                else
                    delete(FZ);
                    FZ = [];
                    return
                end
            end
        else %close the Zscan
            ZScanOpen = false;
            F1 = figure('Toolbar','none',...
                'Menubar', 'none',...
                'NumberTitle','Off',...
                'Name','Line Alignment Tool', 'keypressfcn', @kpf);
            ROIpos = vid.ROIPosition;
            h_im = imshow(zeros(ROIpos(3),  ROIpos(4)), [0 255], 'border', 'tight');
            AX = get(h_im, 'parent');
            stop(vid)
            stop(hTaskPiezo)
            delete(FZ);
            FZ = [];
            moveZ(Zcenter)
            setSrcFor('preview');
            h_im = preview(vid, h_im);
            
            colormap(AX, [0 0 0.3 ; repmat(0:1/252:1, 3, 1)' ; 1 0 0]);
            hold(AX, 'on')
            drawFOV
        end
    end

function findZcenter
    %run this by hitting t he Z key while in the GUI
    doScan = true;
    setSrcFor('preview')
    vid.FramesPerTrigger = 5;
    
    Zs = Zcenter-DZ:Zcenter+DZ;
    I = nan(5,length(Zs));
        
    %raster beam
    set(puBeam, 'Value', 5)
    updateBeam();
    
    %get center position
    moveZ(Zcenter);
    updatePower(rPOW); start(vid); pause(0.01); trigger(vid); pause(vid.FramesPerTrigger/src.FrameRate + 0.02); updatePower(0);
    frame = getdata(vid);
    frame = median(double(frame),4);
    frame = frame - median(frame(:));
    fmed = imgaussfilt(medfilt2(frame, [3 3], 'symmetric'), [3 3]);
    [~,Xc] = max(max(fmed,[],2));
    [~,Yc] = max(max(fmed,[],1));

    for Zix = 1:length(Zs)
        moveZ(Zs(Zix));
        updatePower(rPOW); start(vid); pause(0.01); trigger(vid); pause(vid.FramesPerTrigger/src.FrameRate + 0.02); updatePower(0);
        frame = getdata(vid);
        frame = median(double(frame),4);
        frame = frame - median(frame(:));
        I(5,Zix) = sum(sum(frame(Xc-8:Xc+8, Yc-8:Yc+8)));
    end

    %other beams
    for beam = 1:4
        closeShutter;
        set(puBeam, 'Value', beam);
        updateBeam();
        updatePower(POW);
        for Zix = 1:length(Zs)
            moveZ(Zs(Zix));
            start(vid); openShutter; trigger(vid); pause(vid.FramesPerTrigger/src.FrameRate + 0.02);
            closeShutter;
            frame = getdata(vid);
            frame = median(double(frame),4);
            frame = frame - median(frame(:));
            I(beam,Zix) = sum(sum(frame(Xc-8:Xc+8, Yc-8:Yc+8)));
        end
    end
    figure('Name', 'Z intensity profiles'), plot(I'); legend({'1','2','3','4','Raster'});
    setSrcFor('preview')
    end


    function updateContrast(val)
        C = get(h_im, 'CData');
        maxC = intmax(class(C)); %max(C(:));
        %minC = min(C(:));
        set(AX, 'Clim', [0 val*maxC+1])
    end

    function F2close(obj, event)
        PSFalign = 0;
        delete(scantimer); scantimer = [];
        try
            closeShutter;
            parkLaser;
            stop(vid)
            stop(hTaskPiezo)
            delete(hTaskPiezo)
            stop(hTaskGalvoBeams)
            delete(hTaskGalvoBeams)
            stop(hTaskPowerBeam)
            delete(hTaskPowerBeam)
            stop(hTaskFlipMirror)
            delete(hTaskFlipMirror)
            stop(hTaskShutter)
            delete(hTaskShutter)
            
            most.idioms.safeDeleteObj(hSLAPMi);
            evalin('base','clear hSLAPMi');
        catch
            disp('An error occurred on closing')
        end
        if ~PSFalign
            if ishandle(F2)
                delete(F2)
            end
            if ishandle(F1)
                close(F1)
            end
        end
    end
        
    function updateDZ(val)
       DZ = max(0, min(50, val)); 
    end

    function moveZcenter(newZ, deltaZ)
        if isempty(newZ)
            newZ = Z;
        end
        if nargin<2
            deltaZ = 0;
        end
        moveZ(newZ+deltaZ);
        Zcenter = Z;
    end

    function moveZ(newZ)
            newZ = max(0,min(400, newZ)); %sanity check
            
            stop(hTaskPiezo)
            old = Z * 10/400;
            new = newZ * 10/400;
            
            piezoWaveform = MAA_curve(old,new,0,0,1,N_AO);
            hTaskPiezo.cfgOutputBuffer(N_AO);
            hTaskPiezo.cfgSampClkTiming(10000, 'DAQmx_Val_FiniteSamps', N_AO);
            hTaskPiezo.writeAnalogData(piezoWaveform(1:end-1)');
            start(hTaskPiezo)
            pause(steptime+(1/src.FrameRate))
            Z = newZ;
            set(etZCenter, 'String', num2str(Z));
    end

    function drawFOV
        ROIpos = vid.ROIPosition;
        
        [theta, r] = cart2pol(XYcorner(1)-XYcenter(1), XYcorner(2)-XYcenter(2));
        [X,Y] = pol2cart(theta -pi/8+ [0:pi/4:2*pi], r/(sqrt(2)*cos(pi/8)));
        
        delete([hFOVcenter hFOVcorner hFOV hFOVdots hFOVtext])
        hFOVcenter = scatter(AX, XYcenter(1)-ROIpos(1), XYcenter(2)-ROIpos(2), 'c');
        hFOVcorner = scatter(AX, XYcorner(1)-ROIpos(1), XYcorner(2)-ROIpos(2), 'c');
        
        hFOV = plot(AX,X+XYcenter(1)-ROIpos(1), Y+XYcenter(2)-ROIpos(2), 'c');
        set(hFOV, 'hittest', 'off')
        
        [X,Y] = pol2cart(theta + [pi/4:pi/4:(2*pi-pi/4)], r);
        hFOVdots = scatter(AX,X+XYcenter(1)-ROIpos(1), Y+XYcenter(2)-ROIpos(2), 'c.');
        hFOVtext = text(10,10, ['Line length: ' num2str(sqrt(2) * sqrt(sum((XYcorner-XYcenter).^2))/3.25, 4) ' um'], 'color', 'c', 'parent', AX);
        
        draggable(hFOVcenter, @(varargin)dragFOV, 'endfcn',@(varargin)endDragFOV)
        draggable(hFOVcorner, @(varargin)dragFOV, 'endfcn',@(varargin)endDragFOV)
    end

    function dragFOV
        center = [get(hFOVcenter, 'XData') get(hFOVcenter, 'YData')];
        corner = [get(hFOVcorner, 'XData') get(hFOVcorner, 'YData')];
        
        [theta, r] = cart2pol(corner(1)-center(1), corner(2)-center(2));
        [X,Y] = pol2cart(theta -pi/8 + [0:pi/4:2*pi], r/(sqrt(2)*cos(pi/8)));
   
        set(hFOV, 'Xdata', X+center(1), 'Ydata', Y+center(2))
        
        [X,Y] = pol2cart(theta + [pi/4:pi/4:(2*pi-pi/4)], r);
        set(hFOVdots, 'Xdata', X+center(1), 'Ydata', Y+center(2))
    end

    function FOVGuideVis(tf)
       if tf
           val = 'on';
       else
           val = 'off';
       end
        set(hFOVcenter, 'Visible', val)
        set(hFOVcorner, 'Visible', val)
        set(hFOV, 'Visible', val)
        set(hFOVdots, 'Visible', val)
    end

    function endDragFOV
        ROIpos = vid.ROIPosition;
        center = [get(hFOVcenter, 'XData') get(hFOVcenter, 'YData')];
        corner = [get(hFOVcorner, 'XData') get(hFOVcorner, 'YData')];
        XYcenter = center + ROIpos(1:2);
        XYcorner = corner + ROIpos(1:2);
        
        set(hFOVtext, 'String', ['Line length: ' num2str(sqrt(2) * sqrt(sum((XYcorner-XYcenter).^2))/3.25, 4) ' um']);
    end
    
    function updateBeam
        if doScan && get(puBeam, 'Value')<5
            doScan = false;
            toggleScan;
            hTaskShutter.writeDigitalData(true, 1, true); %OPEN SHUTTER
        elseif get(puBeam, 'Value')==6
            hTaskShutter.writeDigitalData(false, 1, true); %CLOSE SHUTTER
        else
            moveGalvo; 
            hTaskShutter.writeDigitalData(true, 1, true); %OPEN SHUTTER
        end
        %set the beam power
        updatePower(get(sPower, 'Value'));
    end

    function toggleScan        
        if doScan %we are already scanning
            delete(scantimer); scantimer = [];
            doScan = false;
            moveGalvo; %reset the galvo position, sets scan to False
        else
           resetScanTimer;
           doScan = true;
           ScanAmp = 6; 
           ScanRate = 400; %scanrate in Hz %400 is good
           nAO = galvoRate./ScanRate;
           T = linspace(-ScanAmp/2, ScanAmp/2, nAO/2+1);
           T = [T(1:end-1) T(end:-1:2)];
           
           %setup scan - change the parameters so only the target line shows up
           %AO = galvo_trace_SLAPmi(2000, 1, calib); %args: res, tiling, calib
           switch puBeam.Value %Turn the power down at the right time
               case 1
                   angle = 3*pi/8;
                   Va = T*sin(angle) + calib.galvos.offset.line1.X + galvoX;
                   Vb = T*cos(angle) + calib.galvos.offset.line1.Y + galvoY;
                   Vc = zeros(1, nAO); Vd = zeros(1, nAO);
                   Ve = ones(1,nAO)*calib.galvos.offset.E(1); Vf = zeros(1, nAO);
                   Vp = ones(1,nAO)*calib.pockels.fast(2); %THRU CUBE
               case 2
                   angle = pi/8;
                   Va = T*sin(angle) + calib.galvos.offset.line2.X + galvoX;
                   Vb = T*cos(angle) + calib.galvos.offset.line2.Y + galvoY;
                   Vc = zeros(1, nAO); Vd = zeros(1, nAO);
                   Ve = ones(1,nAO)*calib.galvos.offset.E(2); Vf = zeros(1, nAO);
                   Vp = ones(1,nAO)*calib.pockels.fast(2); %THRU CUBE
               case 3
                   angle = 3*pi/8;
                   Va = zeros(1, nAO); Vb = zeros(1, nAO);
                   Vc = T*sin(angle) + calib.galvos.offset.line3.X + galvoX;
                   Vd = T*cos(angle) + calib.galvos.offset.line3.Y + galvoY;
                   Ve = zeros(1, nAO); Vf = ones(1,nAO)*calib.galvos.offset.F(1); 
                   Vp = ones(1,nAO)*calib.pockels.fast(1); %THRU CUBE
               case 4
                   angle = pi/8;
                   Va = zeros(1, nAO); Vb = zeros(1, nAO);
                   Vc = T*sin(angle) + calib.galvos.offset.line4.X + galvoX;
                   Vd = T*cos(angle) + calib.galvos.offset.line4.Y + galvoY;
                   Ve = zeros(1, nAO); Vf = ones(1,nAO)*calib.galvos.offset.F(2); 
                   Vp = ones(1,nAO)*calib.pockels.fast(1); %THRU CUBE
               case 5 %raster scan; change the AOs completely
                   doScan = false;
                   return;
           end
           stop(hTaskGalvoBeams)
           hTaskGalvoBeams.cfgSampClkTiming(galvoRate, 'DAQmx_Val_ContSamps', nAO);
           hTaskGalvoBeams.cfgOutputBuffer(nAO);
           hTaskGalvoBeams.writeAnalogData([Va',Vb',Vc',Vd',Ve',Vf',Vp']);
           start(hTaskGalvoBeams)
        end
    end
    
    function updateGalvo
        if puBeam.Value<=4
            offsetX = calib.galvos.offset.(['line' int2str(puBeam.Value)]).X;
            offsetY = calib.galvos.offset.(['line' int2str(puBeam.Value)]).Y;
        elseif puBeam.Value==5
            offsetX = calib.galvos.offset.raster.X;
            offsetY = calib.galvos.offset.raster.Y;
        end
       galvoX = min(10,max(-10, str2double(get(etX, 'String'))+offsetX)) - offsetX; 
       galvoY = min(10,max(-10, str2double(get(etY, 'String'))+offsetY)) - offsetY;
       set(etX, 'String', num2str(galvoX));
       set(etY, 'String', num2str(galvoY));
       moveGalvo;
    end

    function bdXup
        set(etX, 'String', num2str(galvoX+0.01));
        updateGalvo;
        scrolling = true; pause(0.2)
        while scrolling
            set(etX, 'String', num2str(galvoX+0.1))
            pause(0.3)
            updateGalvo;
        end
    end
    function bdXdn
        set(etX, 'String', num2str(galvoX-0.01))
        updateGalvo;
        scrolling = true; pause(0.2)
        while scrolling
            set(etX, 'String', num2str(galvoX-0.1))
            pause(0.3)
            updateGalvo;
        end
    end
    function bdYup
        set(etY, 'String', num2str(galvoY+0.01))
        updateGalvo;
        scrolling = true; pause(0.2)
        while scrolling
            set(etY, 'String', num2str(galvoY+0.1))
            pause(0.3)
            updateGalvo;
        end
    end
    function bdYdn
        set(etY, 'String', num2str(galvoY-0.01))
        updateGalvo;
        scrolling = true; pause(0.2)
        while scrolling
            set(etY, 'String', num2str(galvoY-0.1))
            pause(0.3)
            updateGalvo;
        end
    end
    function stopscroll
        scrolling = false;
    end

    function updatePower(val)
        
        if ~powerwarned && ~doScan && ~generating && val>0.3
            powerwarned = true;
            button = questdlg('Not Scanning- do you really want such high power?','','Yes','No','No');
            if strcmp(button, 'No')
                val = 0.2;
                powerwarned = false;
            end
        elseif val<0.3
            powerwarned  = false;
        end
        val = round(val*1000)/1000;
        stop(hTaskPowerBeam)
        if puBeam.Value==5
            AOvals = min(2, max(0, [0 2*val]));
            hTaskFlipMirror.writeDigitalData(true, 1, true); %TRUE for raster scan
        elseif puBeam.Value==6 %OFF
            AOvals = [0 0]; %two values for the two lasers
        else
            B = calib.B.B2V(val*calib.B.scaleby(puBeam.Value));
            AOvals = min(2, max(0, [B 0]));
            hTaskFlipMirror.writeDigitalData(false, 1, true); %FALSE for line scan
        end
        hTaskPowerBeam.cfgSampClkTiming(10000, 'DAQmx_Val_FiniteSamps', 10);
        hTaskPowerBeam.cfgOutputBuffer(10);
        hTaskPowerBeam.writeAnalogData(repmat(AOvals, 10, 1));
        start(hTaskPowerBeam)
        set(etPower, 'String', num2str(val));
        set(sPower, 'Value', val);
    end

    function moveGalvo
        updatePower(get(sPower, 'Value'));
        stop(hTaskGalvoBeams)
        switch puBeam.Value
            %The AO channels are as follows:
            %AO 0,4,8,12 = ixs 1:4  = galvos A:D AB->Line1,2 CD->Line3,4,2D
            %AO 16 20 = ixs 5:6 = galvos E:F  E->AB    F->CD
            %AO 24 = ix 7 = Fast Pockels Cell  (ON for galvos A,B,E; OFF for galvos C,D,E)
            case 1
                Va = calib.galvos.offset.line1.X + galvoX;
                Vb = calib.galvos.offset.line1.Y + galvoY;
                Vc = 0; Vd = 0;
                Ve = calib.galvos.offset.E(1); Vf = 0;
                Vp = calib.pockels.fast(2); %THRU CUBE
            case 2
                Va = calib.galvos.offset.line2.X + galvoX;
                Vb = calib.galvos.offset.line2.Y + galvoY;
                Vc = 0; Vd = 0;
                Ve = calib.galvos.offset.E(2); Vf = 0;
                Vp = calib.pockels.fast(2); %THRU CUBE
            case 3
                Va = 0; Vb = 0;
                Vc = calib.galvos.offset.line3.X + galvoX;
                Vd = calib.galvos.offset.line3.Y + galvoY;
                Ve = 0; Vf = calib.galvos.offset.F(1);
                Vp = calib.pockels.fast(1); %REFLECT
            case 4
                Va = 0; Vb = 0;
                Vc = calib.galvos.offset.line4.X + galvoX;
                Vd = calib.galvos.offset.line4.Y + galvoY;
                Ve = 0; Vf = calib.galvos.offset.F(2);
                Vp = calib.pockels.fast(1); %REFLECT
            case 5
                Va = 0; Vb = 0;
                Vc = calib.galvos.offset.raster.X + galvoX;
                Vd = calib.galvos.offset.raster.Y + galvoY;
                Ve = 0; Vf = 4;
                Vp = 0; %OFF
        end
        AOvals = [Va,Vb,Vc,Vd,Ve,Vf,Vp];
        hTaskGalvoBeams.cfgSampClkTiming(10000, 'DAQmx_Val_FiniteSamps', 10);
        hTaskGalvoBeams.cfgOutputBuffer(10);
        hTaskGalvoBeams.writeAnalogData(repmat(AOvals, 10, 1));
        start(hTaskGalvoBeams)
        
    end
   
    function setOffset
        if puBeam.Value<=4
            calib.galvos.offset.(['line' int2str(puBeam.Value)]).X = calib.galvos.offset.(['line' int2str(puBeam.Value)]).X + galvoX;
            calib.galvos.offset.(['line' int2str(puBeam.Value)]).Y = calib.galvos.offset.(['line' int2str(puBeam.Value)]).Y + galvoY;
        elseif puBeam.Value==5
            calib.galvos.offset.raster.X = calib.galvos.offset.raster.X + galvoX;
            calib.galvos.offset.raster.Y = calib.galvos.offset.raster.Y + galvoY;
        end
        set(etX, 'String', num2str(0))
        set(etY, 'String', num2str(0))
        updateGalvo
    end

    function saveCalib
        [FileName,PathName] = uiputfile([basedir '\Calibration\*.cal'], 'Save your calibration file');
        if FileName
            save([PathName filesep FileName], 'calib');
        end
    end
    function loadCalib
        calib = SLAPMi.loadCalibrationData([]); %empty argument requests user selected filename
    end


function generate_projections
if generating
    disp('ABORTING GENERATION')
    generating = false;
    ZScanOpen = true;
    toggleZScan;
    src.Brightness = 2;
    return
end
delete(scantimer); scantimer = [];
resetROI;
stoppreview(vid);
resp = questdlg('Before proceeding, ensure that: Sample is in focus, Raster and Line lasers are both on, All apertures are removed. This will take ~10 minutes!','Ready?','Go','Cancel','Cancel');
if strcmpi(resp, 'Cancel')
    return
end
generating= true;
data = [];
delete(F1)
F1 = []; %kill the preview; it interferes with the camera?
set(tbScan, 'Value', 0)

%Set SLM to white
SLM_drawPattern(2);

Zstep = 2;
Zrange = [Z-4 Z+4];
Zs = linspace(Zrange(1), Zrange(2), ceil(diff(Zrange)/Zstep)+1); 
um2px = 3.256;  %pixels per micron

%pre drift
moveZ(Zs(ceil(length(Zs)/2))); 
[~, refcoords] = snapshot;

%ACQUIRE DATA
[data, T, Td] = measureLines(Zs, refcoords);  %LINE SCAN
data.raster = measureRaster(Zs, refcoords);  %RASTER SCAN % questdlg('Starting Raster Scan. Flip mirror if necessary. Hit OK to continue', 'Raster', 'OK', 'Cancel', 'OK')
moveZ(Zs(ceil(length(Zs)/2))); 

%post drift
[~,coords2] = snapshot;
drift = sqrt(sum((coords2-refcoords).^2));
disp(['The microscope drifted ' num2str(drift/um2px, 3) ' microns during acquisition'])

%POSTPROCESS RASTER SCAN
Icorrected = data.raster.I - repmat(min(data.raster.I,[],3), [1 1 size(data.raster.I,3)]);
data.raster.Zc = sum(repmat(reshape(Zs, 1, 1, []), [size(data.raster.Vx) 1]).*Icorrected, 3)./sum(Icorrected, 3);
data.raster.Zc = inpaint_nans(data.raster.Zc); %center of mass in Z
%data.raster.X = smoothn(data.raster.X);data.raster.Y = smoothn(data.raster.Y);
rasterX = nan(size(data.raster.Zc)); rasterY = rasterX;
for ix1 = 1:size(data.raster.X,1)
    for ix2 = 1:size(data.raster.X,2)
        rasterX(ix1,ix2) = interp1(Zs, squeeze(data.raster.X(ix1,ix2,:)), data.raster.Zc(ix1,ix2));
        rasterY(ix1,ix2) = interp1(Zs, squeeze(data.raster.Y(ix1,ix2,:)), data.raster.Zc(ix1,ix2));
    end
end
data.raster.X = rasterX; data.raster.Y = rasterY;
rasterValid = ~isnan(data.raster.X(:)) & ~isnan(data.raster.Y(:));
if ~all(rasterValid)
    warning('Some points of the raster scan were not found; this could cause errors in PSF estimation. Continuing by extrapolating missing values.')
    data.raster.X = inpaint_nans(data.raster.X); data.raster.Y = inpaint_nans(data.raster.Y);
end

nsamps = length(unique(T)); nsampsD = length(unique(Td));
data.p2vx = scatteredInterpolant(data.raster.X(:), data.raster.Y(:), data.raster.Vx(:));
data.p2vy = scatteredInterpolant(data.raster.X(:), data.raster.Y(:), data.raster.Vy(:));
data.v2px = griddedInterpolant(data.raster.Vx, data.raster.Vy, data.raster.X);
data.v2py = griddedInterpolant(data.raster.Vx, data.raster.Vy, data.raster.Y);
data.p2Z = scatteredInterpolant(data.raster.X(:), data.raster.Y(:), data.raster.Zc(:));
data.v2Z = griddedInterpolant(data.raster.Vx, data.raster.Vy, data.raster.Zc);

%POSTPROCESS LINES
for line = 1:4
    lname = ['line' int2str(line)];
    
    %Center Z plane
    A = repmat(reshape(1:length(Zs), 1, 1, []), [size(data.(lname).I, 1) size(data.(lname).I, 2) 1]);
    I = max(0,data.(lname).I); I = I - repmat(min(I,[],3), [1 1 size(I,3)]);
    Zc = sum(A.*I,3)./sum(I, 3);
    
    %now create an image of the intensity at zcenter
    data.(lname).I_Zc = max(0, nanmax((data.(lname).I(:,:,1:end-2) + data.(lname).I(:,:,2:end-1) + data.(lname).I(:,:,3:end))/3, [],3));
    Zc(data.(lname).I_Zc<max(max(data.(lname).I_Zc(:))/80,5e3)) = nan; %we don't trust the line in very dim regions
    refprof = nanmedian(data.(lname).I_Zc(:, floor(nsampsD/2)*nsamps + [1:nsamps]),2);
    
    %now for every t, get the Rcenter and Dcenter
    Dprofile = zeros(size(refprof,1), length(T));
    data.(lname).Dc = nan(length(T),1);
    xc_window = 400; 
    for t = 1:length(T)
        xc = xcovf(data.(lname).I_Zc(:,t), refprof, xc_window, 'biased');
        [~, maxix] = max(xc);
        data.(lname).Dc(t) = maxix-(xc_window+1); %offset in the center of the line (0 means exactly on the R-line)
        Dprofile(:,t) = circshift(data.(lname).I_Zc(:,t), -data.(lname).Dc(t));
    end
    
    %smooth the measured r values
    for t = 1:length(T)
        for z = 1:length(Zs)
            %nans= isnan(data.(lname).r(:,t,z));
            tmp = data.(lname).r(:,t,z);
            %tmp(nans) = nanmedian(tmp);
            tmp(~isnan(tmp)) = medfilt2(tmp(~isnan(tmp)), [120 1], 'symmetric');
            data.(lname).r(:,t,z) = tmp;
        end
    end
    
    %get the Rcenter at the middle of the line
    data.(lname).Rc = nan(length(T),1);
    for t = 1:length(T)
        d = round(size(Zc,1)/2 + data.(lname).Dc(t));
        data.(lname).Rc(t) = mean(data.(lname).r(d+(-100:100),t, round(Zc(d,t))));
    end

    data.(lname).Zc = interp1(1:length(Zs), Zs, Zc); %assign to be saved.
    data.(lname).Dprofile = max(0, smooth(trimmean(Dprofile,20,2), 0.02, 'loess')); %assign to be saved.
    data.(lname).linelength = sum(data.(lname).Dprofile> (0.8*min(data.(lname).Dprofile) + 0.2*max(data.(lname).Dprofile)));
    
    %postprocessing
    data.(lname).interpDc_fromT = griddedInterpolant(reshape(T, [nsamps nsampsD]), reshape(Td, [nsamps nsampsD]), reshape(data.(lname).Dc, [nsamps nsampsD]));
    data.(lname).interpRc_fromT = griddedInterpolant(reshape(T, [nsamps nsampsD]), reshape(Td, [nsamps nsampsD]), reshape(data.(lname).Rc, [nsamps nsampsD]));
    
    %1D interpolants (functions of Dc or Rc)
    I_r = nan(1,nsamps);
    Ravg = nan(1,nsamps);
    for ix = 1:nsamps
        select = sub2ind([nsamps nsampsD], ix*ones(1,nsampsD), 1:nsampsD);
        I_r(ix) = mean(nansum(data.(lname).I_Zc(:,select)));
        Ravg(ix) = nanmean(data.(lname).Rc(select));
    end
    I_d = nan(1,nsampsD);
    Davg = nan(1,nsampsD);
    for ix = 1:nsampsD
        select = sub2ind([nsamps nsampsD], 1:nsamps, ix*ones(1,nsamps));
        I_d(ix) = mean(nansum(data.(lname).I_Zc(:,select)));
        Davg(ix) = nanmean(data.(lname).Dc(select));
    end
    [~, sortorder] = sort(Davg);
    data.(lname).interpId = griddedInterpolant(Davg(sortorder),I_d(sortorder));
    [~, sortorder] = sort(Ravg);
    data.(lname).interpIr = griddedInterpolant(Ravg(sortorder),I_r(sortorder));
    
    figure('name', [lname ' Intensity']), plot(data.(lname).Dprofile)
    
    %3D interpolants
    %These are functions of (Vx,Vy,d)
    r_td = nan(length(T), size(Zc,1)); %T by D
    Z_td = nan(length(T), size(Zc,1)); %T by D
    for t = 1:length(T)
        for d = 1:size(Zc,1)
            if round(Zc(d,t))>0
                r_td(t,d) = data.(lname).r(d,t,round(Zc(d,t)));        
                Z_td(t,d) = data.(lname).Zc(d,t);
            end
        end
    end
    for i = 1:size(Z_td,1) %SMOOTH
        r_td(i,~isnan(r_td(i,:))) =  medfilt2(r_td(i,~isnan(r_td(i,:))), [1 15], 'symmetric');
        Z_td(i,~isnan(Z_td(i,:))) =  medfilt2(Z_td(i,~isnan(Z_td(i,:))), [1 15], 'symmetric');
    end
        
    data.(lname).interpR_fromT = griddedInterpolant({T(1:nsamps), Td(1:nsamps:end), 1:size(Zc,1)}, reshape(r_td, [nsamps nsampsD size(Zc,1)]));
    data.(lname).interpZ_fromT = griddedInterpolant({T(1:nsamps), Td(1:nsamps:end), 1:size(Zc,1)}, reshape(Z_td, [nsamps nsampsD size(Zc,1)]));
end

%line length
linelengthP = min([data.line1.linelength data.line2.linelength data.line3.linelength data.line4.linelength]) - 50;
data.linelength = sqrt(diff(data.p2vx(600+[-linelengthP/2 linelengthP/2],600+[0 0])).^2 + diff(data.p2vy(600+[-linelengthP/2 linelengthP/2],600+[0 0])).^2);

data.V2T = @V2T;
data.T2V = @T2V;
data.calib.galvos = calib.galvos; data.calib.B = calib.B; data.calib.PSF = calib.PSF; %save only some fields of the calibration
data.timeOfMeasurement = datestr(now, 30);
data.umPerPixel = sqrt(diff(data.p2vx(600 + [-1 1],600 + [0 0])).^2 + diff(data.p2vy(600 + [-1 1],600 + [0 0])).^2)./calib.galvos.pixelsizeVperUM; 

%calibrate aperture
SLM_drawPattern(0.4);
moveZ(Zs(ceil(length(Zs)/2)));
data.aperture = imageAperture(refcoords);

savefilename = [basedir filesep 'PSF' filesep 'PSFmeasured' datestr(now, 30) '.psf'];
save(savefilename, 'data', '-v7.3')
disp(['PSF file saved: ' savefilename])

%make sure that measured PSF is correct:
linePSF_testframes(data);
saveCalib; %save calibration file

setSrcFor('preview');
disp('Done PSF measurement!')
end

    function aperture = imageAperture(refpos)
        %take picture of aperture location; scan the beam
        %reference image:
        [driftIM, driftpos] = snapshot;
        setSrcFor('raster');
        
        updatePower(0)
        closeShutter;
        start(vid)
        pause(0.05)
        trigger(vid);
        ref = getdata(vid);
        ref = median(double(ref),4);
        ref = ref - median(ref(:));

        imAperture = cell(1,4);
        for linenum =1:4
            set(puBeam, 'Value', linenum)
            doScan = false;
            toggleScan; %ensure the beam is scanning
            start(vid)
            updatePower(0.2)
            openShutter;
            trigger(vid);
            
            IM = getdata(vid);
            updatePower(0)
            IM = median(double(IM),4);
            IM = IM - median(IM(:));
            IM = IM-ref;
            IM = imtranslate(IM, refpos-driftpos);
            
            imAperture{linenum} = IM;
        end
        
        Ap1IM = medfilt2(max(imAperture{1},imAperture{2}),[5 5]);
        thresh = (0.3*prctile(Ap1IM(:),99) + 0.7*prctile(Ap1IM(:),1));
        Ap12 = imfill(Ap1IM>thresh, 'holes');
        
        Ap2IM = medfilt2(max(imAperture{3},imAperture{4}),[5 5]);
        thresh = (0.3*prctile(Ap2IM(:),99) + 0.7*prctile(Ap2IM(:),1));
        Ap34 = imfill(Ap2IM>thresh, 'holes');
        
        aperture.raw = imAperture;
        aperture.bw = cat(3, Ap12, Ap34);

        radius = sqrt(max(sum(Ap12(:)), sum(Ap34(:)))/pi);
        aperture.linelength = 1.04*sqrt(diff(data.p2vx(600 + [-radius radius],600 + [0 0])).^2 + diff(data.p2vy(600 + [-radius radius],600 + [0 0])).^2);

        %figure out the range of the aperture
        for line = 1:4
            lname = ['line' int2str(line)];
            bw = aperture.bw(:,:, ceil(line/2));
            frot = imrotate(bw, data.(['line' int2str(line)]).theta, 'nearest','crop');
            Rmin = find(any(frot,1), 1, 'first'); Rmax =  find(any(frot,1), 1, 'last');
            Dmin = find(any(frot,2), 1, 'first'); Dmax =  find(any(frot,2), 1, 'last');
            Rcenter = (Rmin+Rmax)/2; Dcenter = (Dmin+Dmax-size(bw,2))/2;
            [T,Td] = data.V2T(data.(lname).V(:,1),data.(lname).V(:,2),calib.galvos.offset.(lname).X, calib.galvos.offset.(lname).Y, data.(lname).g_angle);
            interpT = scatteredInterpolant(data.(lname).Rc, data.(lname).Dc, T);
            interpTd = scatteredInterpolant(data.(lname).Rc, data.(lname).Dc, Td);
            centerT = interpT(Rcenter, Dcenter);
            centerTd = interpTd(Rcenter, Dcenter);
            [aperture.(lname).offsetX, aperture.(lname).offsetY] = data.T2V(centerT, centerTd,calib.galvos.offset.(lname).X, calib.galvos.offset.(lname).Y, data.(lname).g_angle);
        end
        %show a picture
        figure, imshow(double(cat(3, aperture.bw, Ap34)))
        disp(['Aperture size in Volts: ' num2str(aperture.linelength)])
    end

    function calibrateBcurve
        calib.B.scaleby = ones(1,4); %reset calibration
        calib.B.B2V = griddedInterpolant([0 1], [0 2]);
        calib.B.V2B = griddedInterpolant([0 2], [0 1]);
        calib.B.curve = [];
        
        vid.FramesPerTrigger = 3;
        src.Shutter = 30;
        src.Gain = 4;
        
        %reference
        set(puBeam, 'Value', 1)
        updatePower(0)
        closeShutter;
        start(vid)
        pause(0.05)
        trigger(vid);
        ref = getdata(vid);
        ref = median(double(ref),4);
        
        %get a power curve for one line
        set(puBeam, 'Value', 1)
        updatePower(0)
        powers = 0:0.02:1;
        powcurve = nan(size(powers));
        doScan = false;
        toggleScan;
        for powerix = 1:length(powers)
            start(vid)
            updatePower(powers(powerix))
            pause(0.05)
            trigger(vid);
            IMdata = getdata(vid);
            updatePower(0)
            
            powcurve(powerix) = sum(sum(median(double(IMdata),4) - ref));
            if any(IMdata(:)>(0.9*intmax('uint16')))
                keyboard
               break 
            end
            pause(0.05);
        end
        
        %fit (sin^2(x))^2 to the power curve
        powcurve = powcurve(1:powerix)./nanmax(powcurve(1:powerix));
        
        [maxpow, maxix] = max(smooth(powcurve,5));
        calib.pockels.slow = [0 2*powers(maxix)];
        
        figure, plot(powcurve)
        
        F = @(x,xdata)(powcurve(1) + x(1)*sin((2*xdata)*(pi/2)/diff(calib.pockels.slow)).^x(2));
        x0 = [2 3]; lb = [0 2]; ub = [4 4];
        [x,resnorm,~,exitflag,output] = lsqcurvefit(F,x0,powers(1:powerix)', powcurve', lb, ub);
        
        powcurve = F(x,powers(1:powerix));
        [maxpow, maxix] = max(powcurve);
        powcurve = powcurve./(maxpow);
        
        hold on, plot(powcurve)
        
        calib.B.V2B = griddedInterpolant(2*powers(1:powerix), powcurve);
        calib.B.B2V = griddedInterpolant([powcurve(1:maxix) 100], 2*[powers(1:maxix) powers(maxix)]);
        calib.B.curve = [2*powers(1:maxix)' powcurve(1:maxix)']; 
    end
    function setAperture (val)
        %sanity check
        val = max(min(val, 2),0);
        set(etAp, 'String', num2str(val));
        Ap = val;
        SLM_drawPattern(val);
    end

    function calibrateIntensity %update the calib.B variable
        %calib.B.scaleby = ones(1,4); %reset calibration
        
        vid.FramesPerTrigger = 5;
        src.Gain = 4;
        src.Brightness = 2;
        src.Shutter = 30;
        powIn = get(sPower, 'Value');
        powNow = max(powIn, 0.2);
        
        iterate = true;
        iter = 0;
        while iterate
            iter = iter+1;
            done = false;
            while ~done
                %reference
                set(puBeam, 'Value', 1)
                updatePower(0)
                closeShutter;
                start(vid)
                pause(0.05)
                trigger(vid);
                ref = getdata(vid);
                ref = median(double(ref),4);
                done = true;
                
                doScan = false;
                toggleScan;
                intLine = nan(1,4);
                for line = 1:4
                    set(puBeam, 'Value', line)
                    updateBeam;
                    pause(0.1);
                    
                    start(vid)
                    openShutter;
                    updatePower(powNow)
                    pause(0.05)
                    trigger(vid);
                    IMdata = getdata(vid);
                    closeShutter; 
                    if sum(IMdata(:)>65000)>50
                        disp('Too much power, restart!');
                        powNow = powNow*0.8;
                        done = false;
                        break
                    elseif sum(IMdata(:)>1000)<100
                        disp('Too little power, restart!');
                        powNow = powNow*1.4;
                        done = false;
                        break
                    end
                    pause(0.05)
                    intLine(line) = sum(sum(median(double(IMdata),4) - ref));
                end
            end
            
            updatePower(0);
            scaleby = min(intLine)./intLine;
            scaleby
            calib.B.scaleby = calib.B.scaleby.*scaleby;
            calib.B.scaleby = calib.B.scaleby./min(calib.B.scaleby);
            if min(scaleby)>0.99 || iter>4
                iterate = false; %done
            end
        end
    end

    function closeShutter
        hTaskShutter.writeDigitalData(false, 1, true); % CLOSE SHUTTER
        pause(0.008);
    end

    function openShutter
        hTaskShutter.writeDigitalData(true, 1, true); % OPEN SHUTTER
        pause(0.008);
    end

function parkLaser(startvals)
       stop(hTaskGalvoBeams)
       endvals = [8,8,8,8,0,2,0];
       if nargin<1 || isempty(startvals)
           startvals = endvals;
       elseif get(puBeam, 'Value')~=5 %if startvals are given and we are imaging the lines, don't move the main galvos
           endvals(1:4) = startvals(1:4);
       end
       
       for ix = 7:-1:1
           AOvals(:, ix) = linspace(startvals(ix), endvals(ix), 10);
       end
       hTaskGalvoBeams.cfgSampClkTiming(10000, 'DAQmx_Val_FiniteSamps', 10);
       hTaskGalvoBeams.cfgOutputBuffer(10);
       hTaskGalvoBeams.writeAnalogData(repmat(AOvals, 10, 1));
       start(hTaskGalvoBeams)
end

    function resetScanTimer %a timer to turn off scanning if user accidentally leaves it on
        delete(scantimer);
        scantimer = timer;
        scantimer.TimerFcn = @(myTimerObj, thisEvent)(parkLaser);
        scantimer.startDelay = 1000;
        start(scantimer)
    end

    function setSrcFor(StringIn)
        stop(vid);
        switch StringIn
            case 'lines'
                vid.FramesPerTrigger = 3;
                src.Brightness = 2;
                src.Exposure = 0;
                src.Gain = 4;
                src.Shutter = 0.5;
            case 'raster'
                vid.FramesPerTrigger = 3;
                src.Brightness = 5;
                src.Exposure = 0;
                src.Gain = 5;
                src.Shutter = 30;
            case 'snapshot'
                vid.FramesPerTrigger = 3;
                src.Brightness = 2;
                src.Exposure = 0;
                src.Gain = 4;
                src.Shutter = 30;
            case 'preview'
                vid.FramesPerTrigger = 3;
                src.Brightness = 2;
                src.Exposure = 0;
                src.Gain = 4;
                src.Shutter = 30;
            case 'Zscan'
                src.Brightness = 2;
                src.Exposure = 0;
                src.Gain = 4;
                src.Shutter = 30;
        end
    end

    function [fmed, coords] = snapshot %takes a picture of the raster beam at center of field, to align for drift during psf measurement 
        stop(vid)
        setSrcFor('snapshot');
        if get(puBeam, 'Value')~=5
            set(puBeam, 'Value', 5); updateBeam;
            pause(0.8);
        end
        galvoX = 0; galvoY = 0; moveGalvo;
        updatePower(rPOW)
        start(vid)
        trigger(vid);
        pause(0.1);
        IM = median(double(getdata(vid)),4);
        updatePower(0); parkLaser;
        IM = IM - median(IM(:));
        
        fmed = imgaussfilt(medfilt2(IM, [3 3], 'symmetric'), [3 3]);
        [~,Xc] = max(max(fmed,[],2));
        [~,Yc] = max(max(fmed,[],1));
        box = fmed(Xc + (-10:10), Yc+(-10:10));
        sX = sum(box,1); sY = sum(box,2);
        coords(1) = Xc + sum((-10:10)'.*sY)/sum(sY);
        coords(2) = Yc + sum((-10:10).*sX)/sum(sX);
    end

    function raster = measureRaster(Zs, refpos)
        rasterpos = linspace(-3.75, 3.75, 8);
        wh = waitbar(0, 'Initializing Raster Scan...');
        
        stop(vid)
        setSrcFor('raster');
        set(puBeam, 'Value', 5)
        updatePower(0)
        closeShutter;
        start(vid)
        trigger(vid);
        ref = median(double(getdata(vid)),4);
        ref = ref - median(ref(:));
        ref = ref/src.Shutter;
        
        raster.X = nan(length(rasterpos),length(rasterpos), length(Zs));
        raster.Y = nan(length(rasterpos),length(rasterpos), length(Zs));
        raster.Vx = repmat((rasterpos + calib.galvos.offset.raster.X)', 1, length(rasterpos));
        raster.Vy = repmat(rasterpos + calib.galvos.offset.raster.Y, length(rasterpos), 1);
        
        for z = 1:length(Zs)
            %drift
            moveZ(Zs(ceil(length(Zs)/2)));
            [~, driftpos] = snapshot; disp(['Drift:' num2str(driftpos-refpos)])
            setSrcFor('raster');
            
            moveZ(Zs(z));
            pause(0.1);
            for X = 1:length(rasterpos)
                waitbar((z-1)/length(Zs) + (X/length(rasterpos))/length(Zs), wh, ['Performing raster scan, z=' int2str(z) ', X=' int2str(X)]);
                for Y = 1:length(rasterpos)
                    Vx = raster.Vx(X,Y);
                    Vy = raster.Vy(X,Y);
                    
                    repeat = 1;
                    while repeat<5
                        stop(hTaskGalvoBeams)
                        AOvals = [0,0,Vx,Vy,-1,-1,0];
                        hTaskGalvoBeams.cfgSampClkTiming(10000, 'DAQmx_Val_FiniteSamps', 10);
                        hTaskGalvoBeams.cfgOutputBuffer(10);
                        hTaskGalvoBeams.writeAnalogData(repmat(AOvals, 10, 1));
                        start(hTaskGalvoBeams)
                        updatePower(rPOW)
                        
                        start(vid)
                        pause(0.01)
                        trigger(vid);
                        pause(vid.FramesPerTrigger/src.FrameRate + 0.02)
                        updatePower(0)
                        
                        frame = getdata(vid);
                        saturating = sum(frame(:)> intmax('uint16')*0.99);
                        good = sum(frame(:)>intmax('uint16')/20)-saturating;
                        if saturating>1
                            if src.Shutter<0.05
                                disp([int2str(saturating) ' pixels saturated']);  repeat = inf;
                            else
                                pause(0.01)
                                src.Shutter = src.Shutter * (2/3);
                                repeat = repeat+1; disp(['Shutter set to: ' num2str(src.Shutter)]);
                            end
                        elseif good<5
                            if src.Shutter>30
                                disp(['Only ' int2str(good) ' pixels good']); repeat = inf;
                            else
                                pause(0.01)
                                src.Shutter = src.Shutter * (3/2);
                                repeat = repeat+1; disp(['Shutter set to: ' num2str(src.Shutter)]);
                            end
                        else
                            repeat = inf; %done
                        end
                    end
                    
                    frame = median(double(frame),4);
                    frame = frame - median(frame(:));
                    frame = (frame/src.Shutter)-ref;
                    fmed = imgaussfilt(medfilt2(frame, [3 3], 'symmetric'), [3 3]);
                    fmed = imtranslate(fmed, refpos-driftpos);
                    
                    [~,Xc] = max(max(fmed,[],2));
                    [~,Yc] = max(max(fmed,[],1));
                    if Xc>10 && Xc<=(size(frame,1)-10) && Yc>10 && Yc<=(size(frame,2)-10)
                        box = frame(Xc + (-10:10), Yc+(-10:10));
                        sX = sum(box,1); sY = sum(box,2);
                        raster.X(X,Y,z) = Xc + sum((-10:10)'.*sY)/sum(sY);
                        raster.Y(X,Y,z) = Yc + sum((-10:10).*sX)/sum(sX);
                        raster.I(X,Y,z) = sum(sX); %sum intensity? max intensity?
                    else
                        disp(['Couldn''t find raster beam at Z=' int2str(Z) ' X=' int2str(X) ' Y=' int2str(Y)])
                    end
                end
            end
        end
        delete(wh)  %delete the waitbar
    end

    function [data, T, Td] = measureLines(Zs, refpos)
        wh = waitbar(0, 'Initializing...');
        tilingfactor = 1;
        %POW = 0.06; %line power
        fprintf('Line Power: %d. Volts: %d \n', [POW calib.B.B2V(POW*calib.B.scaleby(1))]);
        um2v = calib.galvos.pixelsizeVperUM; %1/45.3;% scanner volts per micron
        scansize = (250 + 15)*um2v; %line length in volts
        scansizeD = (1-1/tilingfactor)*scansize + 15*um2v;
        nsamps = 7; %10;
        nsampsD = 2+ tilingfactor; %number of samples in the perpendicular direction
        [T,Td] = ndgrid(linspace(-scansize/2, scansize/2, nsamps), linspace(-scansizeD/2, scansizeD/2, nsampsD));
        T = T(:); Td=Td(:);
        
        %%%%%%%%% GET REFERENCE IMAGE
        stop(vid)
        setSrcFor('lines')
        set(puBeam, 'Value', 1)
        updatePower(0)
        closeShutter;
        start(vid)
        pause(0.05)
        trigger(vid);
        ref = getdata(vid);
        ref = median(double(ref),4);
        ref = ref - median(ref(:));
        ref = ref/src.Shutter;
        
        Ithresh = min(1000, max(100, 6*std(ref(:)))); %std of reference image?
        data.refIMsize = size(ref);
        
        
        data.testframe = struct([]); %We save testframes to reconstruct to ensure PSF measurement is accurate
        testframe_ixs = [2 length(T)-1];
        data.testframe(4, length(Zs), length(testframe_ixs)).t = nan; %extize
        
        %%%%%%%% GET LINE IMAGES
        lines = [1 2 3 4]; %randperm(4);
        for line = lines
            %drift
            moveZ(Zs(ceil(length(Zs)/2)));
            [~, driftpos] = snapshot; disp(['Drift: ' num2str(driftpos-refpos)])
            setSrcFor('lines');
            
            set(puBeam, 'Value', line); updateBeam; pause(1); %allow the mirror to move back to position
            switch line
                case 1
                    angle = -pi/8;
                    [Va, Vb] = T2V(T,Td, calib.galvos.offset.line1.X, calib.galvos.offset.line1.Y, angle);
                    Vc = T*0; Vd = T*0;
                    Ve = calib.galvos.offset.E(1); Vf = 0;
                    Vp = calib.pockels.fast(2);
                    data.line1.V = [Va Vb];
                case 2
                    angle = -3*pi/8;
                    [Va, Vb] = T2V(T,Td, calib.galvos.offset.line2.X, calib.galvos.offset.line2.Y, angle);
                    Vc = T*0; Vd = T*0;
                    Ve = calib.galvos.offset.E(2); Vf = 0;
                    Vp = calib.pockels.fast(2);
                    data.line2.V = [Va Vb];
                case 3
                    angle = -pi/8;
                    Va = T*0; Vb = T*0;
                    [Vc, Vd] = T2V(T,Td, calib.galvos.offset.line3.X, calib.galvos.offset.line3.Y, angle);
                    Ve = 0; Vf = calib.galvos.offset.F(1);
                    Vp = calib.pockels.fast(1);
                    data.line3.V = [Vc Vd];
                case 4
                    angle = -3*pi/8;
                    Va = T*0; Vb = T*0;
                    [Vc, Vd] = T2V(T,Td, calib.galvos.offset.line4.X, calib.galvos.offset.line4.Y, angle);
                    Ve = 0; Vf = calib.galvos.offset.F(2);
                    Vp = calib.pockels.fast(1);
                    data.line4.V = [Vc Vd];
            end
            %Sanity check
            if any(Va<-10 | Vb<-10 | Vc<-10 | Vd<-10 | Va>10 | Vb>10 | Vc>10 | Vd>10)
                keyboard
            end
            data.(['line' int2str(line)]).g_angle = angle;
            data.(['line' int2str(line)]).r = nan(size(ref,1),length(T), length(Zs));
            data.(['line' int2str(line)]).I = nan(size(ref,1),length(T), length(Zs));
            data.(['line' int2str(line)]).sigma = nan(size(ref,1),length(T), length(Zs));
            data.(['line' int2str(line)]).theta = nan;
            
            for z = 1:length(Zs)
                stop(vid)
                moveZ(Zs(z));
                pause(0.1);
                updatePower(POW)
                for t = 1:length(T)
                    %move the line, take a picture
                    waitbar((z-1)/length(Zs) + (t/length(T))/length(Zs), wh, ['Performing line scan line=' int2str(line) ', z=' int2str(z) ', t=' int2str(t)]);
                    repeat = 1;
                    while repeat<6
                        stop(hTaskGalvoBeams)
                        AOvals = [Va(t),Vb(t),Vc(t),Vd(t),Ve,Vf,Vp];
                        hTaskGalvoBeams.cfgSampClkTiming(10000, 'DAQmx_Val_FiniteSamps', 10);
                        hTaskGalvoBeams.cfgOutputBuffer(10);
                        hTaskGalvoBeams.writeAnalogData(repmat(AOvals, 10, 1));
                        start(hTaskGalvoBeams)
                        openShutter;
                        start(vid)
                        pause(0.02)
                        trigger(vid);
                        pause(vid.FramesPerTrigger/src.FrameRate+0.01)
                        closeShutter;
                        frame = getdata(vid);
                        
                        saturating = sum(frame(:)> intmax('uint16')*0.99);
                        good = sum(frame(:)>intmax('uint16')/30)-saturating;
                        if saturating>800
                            if src.Shutter<0.05
                                disp([int2str(saturating) ' pixels saturated']);  repeat = inf;
                            else
                                pause(0.01)
                                src.Shutter = src.Shutter * (2/3);
                                disp(['Frame sum: ' int2str(sum(frame(:)))])
                                repeat = repeat+1; disp(['Shutter set to: ' num2str(src.Shutter)]);
                            end
                        elseif good<2000
                            if src.Shutter>30
                                disp(['Only ' int2str(good) ' pixels good']); repeat = inf;
                            else
                                pause(0.01)
                                src.Shutter = src.Shutter * (3/2);
                                repeat = repeat+1; disp(['Shutter set to: ' num2str(src.Shutter)]);
                            end
                        else
                            repeat = inf; %done
                        end
                    end
                    frame = median(double(frame),4);
                    frame = frame - median(frame(:));
                    frame = (frame/src.Shutter) -ref;
                    frame = imtranslate(frame, refpos-driftpos);
                    
                    %Save TEST FRAMES we'll use these to verify the quality of the psf measurement
                    if any(t==testframe_ixs)
                        testframeix = find(t==testframe_ixs, 1);
                        data.testframe(line, z, testframeix).im = frame;
                        data.testframe(line, z, testframeix).line = line;
                        data.testframe(line, z, testframeix).Brightness = src.Shutter;
                        data.testframe(line, z, testframeix).Z = Zs(z);
                        data.testframe(line, z, testframeix).z = z;
                        data.testframe(line, z, testframeix).T= T(t);
                        data.testframe(line, z, testframeix).Td= Td(t);
                        data.testframe(line, z, testframeix).t = t;
                        if line<=2
                            data.testframe(line, z, testframeix).Vx = AOvals(1);
                            data.testframe(line, z, testframeix).Vy = AOvals(2);
                        else
                            data.testframe(line, z, testframeix).Vx = AOvals(3);
                            data.testframe(line, z, testframeix).Vy = AOvals(4);
                        end
                    end
                    
                    %identify the image angle
                    if isnan(data.(['line' int2str(line)]).theta)
                        data.(['line' int2str(line)]).theta = getangle(frame, line);
                    end
                    %rotate the image
                    frot = imrotate(frame, data.(['line' int2str(line)]).theta, 'bilinear','crop');
                    
                    %get position and intensity in rotated coordinate frame
                    [~, bestcol] = max(sum(frot,1));
                    if bestcol<10 || bestcol>size(frot,2)-10
                        keyboard
                    end
                    dsupport = mean(frot(:,bestcol-10:bestcol+10),2)>Ithresh;
                    dsupport = max(1, find(dsupport,1,'first')-5):min(length(dsupport), find(dsupport,1,'last')+5);
                    for d = dsupport
                        ldata = max(0, frot(d, bestcol+(-24:24)));
                        CoM = sum((-24:24) .* (ldata-min(ldata)))/sum(ldata-min(ldata));
                        if abs(CoM)<=13 %center of mass should be in the center of the line
                            data.(['line' int2str(line)]).r(d,t,z) = CoM+bestcol;
                            data.(['line' int2str(line)]).I(d,t,z) = sum(ldata(25+round(CoM)+(-10:10)));
                        else
                            [];
                            %bad fit
                        end
                    end
                end
            end
        end
        delete(wh)  %delete the waitbar
    end
end

    
    function [T,Td] = V2T(Vx,Vy,Xoff,Yoff,angle)
        T = (Vx - Xoff)*cos(angle) - (Vy-Yoff)*sin(angle);
        Td = (Vx - Xoff)*sin(angle) + (Vy-Yoff)*cos(angle);
    end
    function [Vx,Vy] = T2V(T,Td,Xoff,Yoff,angle)
        Vx = T*cos(angle) + Td*sin(angle) + Xoff;
        Vy = -T*sin(angle) + Td*cos(angle) + Yoff;
    end
    
    function theta = getangle(image, line)
    %returns the angle in degrees (as passed to imrotate) that image should be rotated to make the line vertical
    image = max(0, image-median(image(:)));
    if line>0 && line<5
        im_angles = [21  65   -24   -67]; %initial guesses for angles
        center =  im_angles(line); %initial guess
        scale = 5; %we will do a pyramid search
        nsamps = 7;
        bestangle = center;
        bestcontrast = -inf;
        for iter = 1:5
            angles = bestangle + linspace(-scale,scale, nsamps);
            for angle = angles
               frot = imrotate(image, angle, 'bilinear','crop');
               contrast = max(sum(frot,1));
               if contrast>bestcontrast
                   bestangle = angle;
                   bestcontrast = contrast;
               end
            end
            scale = 2*scale/(nsamps-1);
        end
        theta = bestangle;
    else
        disp('getangle: bad linenumber argument')
        theta = 0;
    end
    disp(['Line ' int2str(line) ' Angle: ' num2str(theta)]);
    end
    
    function resetOffsets
    global calib
    if strcmpi(questdlg('Reset Offsets?','reset', 'Yes','No', 'No'), 'Yes')
        X12r = (calib.galvos.offset.line1.X - calib.galvos.offset.line2.X);
        calib.galvos.offset.line1.X = X12r/2;
        calib.galvos.offset.line2.X =  -X12r/2;
        
        Y12r = (calib.galvos.offset.line1.Y - calib.galvos.offset.line2.Y);
        calib.galvos.offset.line1.Y = Y12r/2;
        calib.galvos.offset.line2.Y =  -Y12r/2;
        
        X34r = (calib.galvos.offset.line3.X - calib.galvos.offset.line4.X);
        calib.galvos.offset.line3.X = X34r/2;
        calib.galvos.offset.line4.X =  -X34r/2;
        
        Y34r = (calib.galvos.offset.line3.Y - calib.galvos.offset.line4.Y);
        calib.galvos.offset.line3.Y = Y34r/2;
        calib.galvos.offset.line4.Y =  -Y34r/2;
    end
    end