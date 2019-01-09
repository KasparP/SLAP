function hFig = scandataPlot3D(scandata, optsin)
F = scandata.frames;
Nchans = size(F(1).pmtData,2);
nVol = floor(scandata.metadata.nVolumes);

P = cell2mat(reshape({F(:,1:nVol).pmtData}, 1, size(F,1), nVol));
P = permute(P, [2 1 3]);

Z = cell2mat(reshape({F(:,1:nVol).Z}, 1, size(F,1), nVol));
Z = permute(Z, [2 1 3]);

%properties
S = 1; %slice number
maxC = max(P(:));
C = maxC/2; %contrast level
cmap = gray(256); cmap(end,:) = [1 0 0];

M = cell(1,4); %Movie for upstroke and downstroke, and each line

if scandata.metadata.ZBiDi
    nVol = 2*nVol;
    volTimes = (scandata.metadata.volumePeriod/2)*(0:nVol-1);
    Zdes = linspace(scandata.metadata.Zmin,scandata.metadata.Zmax,floor(size(P,1)/2))/40;
    uplines = 1:floor(size(P,1)/2);
    downlines = floor(size(P,1)/2)+1:size(P,1);
else %unidirectional Zscan
    volTimes = scandata.metadata.volumePeriod*(0:nVol-1);
    duty = 0.85;  %duty = scandata.metadata.piezo.opts.duty;
    Zdes = linspace(scandata.metadata.Zmin,scandata.metadata.Zmax,ceil(duty*size(P,1)))/40;
    uplines = [(-1:0)+size(P,1) 1:floor(duty*size(P,1))];
    downlines = [];
end
%movie of projections in each direction on each axis
hFig = figure('name', [scandata.metadata.galvoDataFileName ' movie'], 'pos',  [25  900   1550   200], 'WindowScrollWheelFcn', @scrollfun);
hSliderFrame = uicontrol('parent', hFig, 'Style', 'slider','Min',1,'Max',nVol,'Value',S,'SliderStep',[1/(nVol-1) 10/(nVol-1)],'Position', [0 0 200 15], 'Callback', @framefun);
hSliderContrast = uicontrol('parent', hFig, 'Style', 'slider','Min',0,'Max',maxC,'Value',C,'SliderStep',min(1, max(1/256,[1/(maxC-1) 10/(maxC-1)])),'Position', [200 0 100 15], 'Callback', @contrastfun); %#ok<NASGU>
hPlayTimer = timer('ExecutionMode', 'fixedSpacing', 'Period', 0.1, 'TimerFcn', @timerFunc);
hPlayButton = uicontrol('parent', hFig, 'Style', 'toggle','Value',0,'Position', [300 0 80 15], 'String', 'Play', 'Callback', @playButton); %#ok<NASGU>
hExportButton = uicontrol('parent', hFig, 'Style', 'pushbutton','Position', [380 0 50 15], 'String', 'Export', 'Callback', @exportButton); %#ok<NASGU>

for lineID = 1:4
    select = scandata.line==lineID;

    Pl = P(:,select,:);
    Zl = Z(:,select,:);
    
    if scandata.metadata.ZBiDi
        M{lineID} = nan(length(Zdes), size(Pl,2),2*size(Pl,3));
        for lxl = 1:size(Pl,2)
            for vol = 1:size(Pl,3)
                if sum(~isnan(Pl(uplines,lxl,vol)))>2 && sum(~isnan(Pl(downlines,lxl,vol)))>2
                    M{lineID}(:,lxl,2*vol-1) = interp1(Zl(uplines,lxl,vol), Pl(uplines,lxl,vol), Zdes);
                    M{lineID}(:,lxl,2*vol) = interp1(Zl(downlines,lxl,vol), Pl(downlines,lxl,vol), Zdes);
                end
            end
        end
    else %unidirectional Zscan
        M{lineID} = nan(length(Zdes), size(Pl,2),size(Pl,3));
        for lxl = 1:size(Pl,2)
            for vol = 1:size(Pl,3)
                M{lineID}(:,lxl,vol) = interp1(Zl(uplines,lxl,vol), Pl(uplines,lxl,vol), Zdes);
            end
        end
    end
    
    M{lineID}(isnan(M{lineID})) = 0;
    %M(:,line) = alignLateral(M(:,line));

    hAx(lineID) = axes('pos', [0 (lineID-1)*0.25 1 0.25]);
    hIm(lineID) = imagesc(M{lineID}(:,:,S));
    if lineID<4
        line([0 1e5], [0.5 0.5], 'color', 'r', 'parent', hAx(lineID));
    end
    colormap(cmap)
    set(hAx(lineID), 'visible', 'off');
end

annotation('textbox', [0.9 0.15 0.1 0.1], 'String', 'Line 1', 'color', 'c', 'linestyle', 'none', 'horizontalalignment', 'right')  
annotation('textbox', [0.9 0.4 0.1 0.1], 'String', 'Line 2', 'color', 'c','linestyle', 'none','horizontalalignment', 'right')    
annotation('textbox', [0.9 0.65 0.1 0.1], 'String', 'Line 3', 'color', 'c','linestyle', 'none','horizontalalignment', 'right')    
annotation('textbox', [0.9 0.9 0.1 0.1], 'String', 'Line 4', 'color', 'c','linestyle', 'none', 'horizontalalignment', 'right')   

%Frame time
axes(hAx(4));
hTime = text(25,5,'t = 0 sec', 'color', 'c');
%add scalebar
axes(hAx(1));
pxPer5umR = 5/scandata.pixelSizeUM;
pxPer5umZ = 5/scandata.metadata.dZ;
hScale = scalebar('border', 'UL', 'Xlen', pxPer5umR, 'Ylen', pxPer5umZ);
hScale.Position =[size(M{1},2)*0.03 1];
hScale.hTextX_Pos = [-2,12]; hScale.hTextY_Pos = [-1.4 12];
set(hScale.hLineX, 'color', 'c', 'linewidth', 2); set(hScale.hTextX, 'color', 'c', 'String', '5 um','verticalAlignment', 'top', 'horizontalAlignment', 'left');
set(hScale.hLineY, 'color', 'c', 'linewidth', 2); set(hScale.hTextY, 'color', 'c', 'String', '5 um', 'verticalAlignment', 'bottom', 'horizontalAlignment', 'left');

% for I =1:4  %for figure generation
%     figure('name', int2str(I), 'pos', [0 700 1600 1.6*750/4])
%     hAxI = axes('pos', [0 0 1 1]);
%     imagesc(M{I}(:,:,1), [0 20]);
%     colormap gray
%     axis off
% end

    function scrollfun(hObj, eventdata)
        UPDN = eventdata.VerticalScrollCount;
        S = max(1, min(nVol, S+UPDN));
        set(hSliderFrame, 'Value', S);
        updateImage();
    end

    function framefun (hObj,event)
        S = round(get(hObj,'Value'));
        updateImage();
    end

    function updateImage
        for l = 1:4
            set(hIm(l),'cdata',M{l}(:,:,S))
        end
        %add frametime
        set(hTime, 'String', ['t= ' num2str(volTimes(S), 2) 'sec'])

    end
    function contrastfun(hObj,event)
        C = max(hObj.Value, 1e-4);
        for l = 1:4
            set(hAx(l),'clim',[0 C])
        end
    end
    function timerFunc(obj, event)
        S = mod(S, nVol)+1;
        updateImage();
    end
    function playButton(obj, event)
        if obj.Value
           start(hPlayTimer);
           set(obj, 'String', 'Pause');
        else
            stop(hPlayTimer);
            set(obj, 'String', 'Play');
        end
    end
    function exportButton(obj, event)
        [fn, dr] = uiputfile([scandata.metadata.fileNameStem '.avi'], 'Select a name for the movie');
        set(hPlayButton, 'value', false); set(hPlayButton, 'String', 'Play'); stop(hPlayTimer);
        
        set([hSliderFrame hSliderContrast hPlayButton hExportButton], 'visible', 'off');
        
        w = VideoWriter([dr filesep fn]);
        w.FrameRate = 5;
        
        open(w);
        rect = get(hFig, 'pos');
        for v = 1:nVol
            S = v; updateImage();
            frame = getframe(hFig);
            writeVideo(w,frame);
        end
        close(w)
        
        set([hSliderFrame hSliderContrast hPlayButton hExportButton], 'visible', 'on');
    end


    function mov = alignLateral(mov)
        maxlag = 10;
        nV = size(mov{2,1},3);
        for v = 1:nV
            v1filt = imgaussfilt(mov{1}(:,:,v),1); v3filt = imgaussfilt(mov{1}(:,:,min(v+1,end)),1);
            v2filt = imgaussfilt(mov{2}(:,:,v),1);
            
            lags = nan(1,size(v2filt,1));
            for strip = 1:size(v2filt,1)
               [xc1, lags1] = xcorr(v2filt(strip,:), v1filt(strip,:), maxlag);
               [xc2, lags2] = xcorr(v2filt(strip,:), v3filt(strip,:), maxlag);
               [b1, l1] = max(xc1); [b2, l2] = max(xc2); 
               lags(strip) = (lags1(l1)*b1 + lags2(l2)*b2)/(b1+b2);
               mov{2}(strip,:,v) = interp1(mov{2}(strip,:,v), (1:size(v2filt,2))-lags(strip));
            end
            
            
%             CC = normxcorr2_general(mov{1,line}(:,:,v),mov{2,line}(:,:,v), 0.8*numel(mov{1,line}(:,:,v)));
%             [maxval, maxind] = max(CC(:));
%             [row,col] = ind2sub(size(CC),maxind);
%             rowshift = row-ceil(size(CC,1)/2); colshift = col-ceil(size(CC,2)/2);
%             mov{2,line}(:,:,v) = imtranslate(mov{2,line}(:,:,v), [-colshift -rowshift]);
        end
    end
end