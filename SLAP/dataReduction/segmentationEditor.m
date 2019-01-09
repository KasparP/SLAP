function refIM = segmentationEditor(refIM)

opts = refIM.segopts;

%data prep
bw = refIM.bw;
core = refIM.labels == 2;
shaft = refIM.labels==4;
bw = bw|core|shaft;

Inorm = imtophat(refIM.IM, strel('disk',ceil(1/opts.XYscale)));
i1 = refIM.pts(:,1);
i2 = refIM.pts(:,2);
i3 = refIM.pts(:,3);
newPt = false(size(i3));
pointTypes = refIM.labels(sub2ind(size(refIM.labels), i1,i2,i3));
distThresh = round(opts.distThreshUM/opts.XYscale);
ysz = [size(bw,1) size(bw,2)];

%state variables
imBlock = ones([size(refIM.IM) 3]);
Z = 1;
nZ = size(refIM.IM,3);
nS = length(i1);
addedTF = [];

%handles
hF = figure('Name', 'Segmentation Editor', 'numbertitle', 'off');
hS = []; 
hIM = imshow(squeeze(imBlock(:,:,1,:)));
hold on, hT = text(size(imBlock,1)/2, size(imBlock,2)/2, 'Loading...', 'color', 'k', 'fontsize', 20, 'HorizontalAlignment', 'center');
hAx = get(hIM, 'parent');
hButtonRecalc = uicontrol('Style', 'pushbutton','Position', [5 5 80 25],'String','Recalculate', 'Callback' , @(varargin)(calculateSupport)); %#ok<NASGU>
hButtonSave = uicontrol('Style', 'pushbutton','Position', [95 5 80 25],'String','Save', 'Callback' , @(varargin)(saveSeg)); %#ok<NASGU>

set (hF, 'WindowScrollWheelFcn', @scrollFunc);
set (hF, 'WindowButtonDownFcn', @bdf);
set (hF, 'WindowKeyPressFcn', @kpf);
set (hF, 'CloseRequestFcn', @exit);

drawnow;
calculateSupport;
delete(hT); hT = [];

if nargout
    waitfor(hF)
end
    

    function exit(src, evnt)
        if strcmp(questdlg('Exit segmentationEditor?', 'SEGEDIT', 'Yes', 'No', 'No'), 'Yes')
            delete(hF);
        end
    end
    function saveSeg
        [fn dr] = uiputfile([dataDirectory '*.mat']);
        if fn
            save([dr filesep fn], 'refIM');
            disp(['SEG file saved to: ' dr filesep fn]);
        end
    end

    function scrollFunc(obj, evnt)
        UPDN = evnt.VerticalScrollCount;
        Z = max(1,min(nZ, Z+UPDN));
        redraw;
    end

    function kpf(src, event)
       switch event.Key
           case {'d', 'rightarrow'}
               xl = get(hAx, 'xlim');
               maxmove = ysz(2) - xl(2);
               set(hAx, 'xlim', xl+min(20,maxmove));
           case {'a', 'leftarrow'}
               xl = get(hAx, 'xlim');
               maxmove = xl(1);
               set(hAx, 'xlim', xl-min(20,maxmove));
           case {'w', 'uparrow'}
               yl = get(hAx, 'ylim');
               maxmove = yl(1);
               set(hAx, 'ylim', yl-min(20,maxmove));
           case {'s', 'downarrow'}
               yl = get(hAx, 'ylim');
               maxmove = ysz(1)-yl(2);
               set(hAx, 'ylim', yl+min(20,maxmove));
           case 'hyphen'
               zoom(0.8);
           case 'equal'
               zoom(1/0.8);
           case {'2', 'rightbracket'}
               evnt.VerticalScrollCount = 1;
               scrollFunc([], evnt);
           case {'1', 'leftbracket'}
               evnt.VerticalScrollCount = -1;
               scrollFunc([], evnt);
           case 'k'
               keyboard
       end
    end
    function bdf(src, event)
        cp = get(hAx,'CurrentPoint');
        cp = cp(1, 1:2);
        %If the point is within the axes
        if all(cp>=1 & cp<=ysz)
            mod = get(src, 'SelectionType');
            switch mod
                case 'normal'
                    %add a point
                    i1(end+1) = round(cp(2));
                    i2(end+1) = round(cp(1));
                    i3(end+1) = Z; 
                    addedTF(end+1,:) = false;
                    pointTypes(end+1) = 5;
                    newPt(end+1) = true;
                case 'alt'
                    %delete nearest point
                    select = find(i3==Z);
                    d = sqrt((i1(select)-cp(2)).^2 + (i2(select)-cp(1)).^2);
                    [minval,minix] = min(d); minix = select(minix);
                    if minval<6
                        i1(minix) = [];
                        i2(minix) = [];
                        i3(minix) = [];
                        addedTF(minix,:) = [];
                        pointTypes(minix) = [];
                        newPt(minix) = [];
                    end
                otherwise
                    keyboard
            end
            redraw;
        end
    end

    function redraw
        set(hIM, 'CData', squeeze(imBlock(:,:,Z,:)));
        
        hold(hAx, 'on');
        delete(hS); %clear the scatterpoints
        %delete(hT);
        hS = []; %hT = [];
        
        markercolors = {'r', 'c', 'g','b','m'};
        for pt = 1:max(pointTypes)
            select = i3==Z & pointTypes==pt;
            hS(end+1) = scatter(hAx,i2(select & newPt),i1(select & newPt), 'marker', 'o', 'SizeData', 10, 'markeredgecolor','none', 'markerfacecolor', 'r');
            hS(end+1) = scatter(hAx,i2(select),i1(select), 'marker', 'o', 'SizeData', 10, 'markeredgecolor', markercolors{pt});
            
            select = i3<Z & addedTF(:,Z) & pointTypes == pt;
            hS(end+1) = scatter(hAx,i2(select & newPt),i1(select & newPt), 'marker', '^', 'SizeData', 16, 'markeredgecolor','none', 'markerfacecolor', 'r');
            hS(end+1) = scatter(hAx,i2(select),i1(select), 'marker', '^', 'SizeData', 16, 'markeredgecolor', markercolors{pt});
            
            
            select = i3>Z & addedTF(:,Z) & pointTypes == pt;
            hS(end+1) = scatter(hAx,i2(select & newPt),i1(select & newPt), 'marker', 'v', 'SizeData', 10, 'markeredgecolor','none', 'markerfacecolor', 'r');
            hS(end+1) = scatter(hAx,i2(select),i1(select), 'marker', 'v', 'SizeData', 10, 'markeredgecolor', markercolors{pt});
        end
    end

    function calculateSupport
        %we are going to calculate influence of seeds over nearby pixels with a
        %'flood fill' operation
        added = cell(length(i1),1);
        addedTF = false(length(i1),size(refIM.IM,3)); 
        
        %recalculate pointTypes
        select = pointTypes==5;
        pointTypes(select) = refIM.labels(sub2ind(size(refIM.labels), i1(select),i2(select),i3(select)));
        
        %recalculate all the 'added' points
        intensities = Inorm(sub2ind(size(Inorm), i1,i2,i3));
        intensities(pointTypes==2 | pointTypes==4 | pointTypes==5) = 2*intensities(pointTypes==2 | pointTypes==4 | pointTypes==5);
        [sorted,sortorder] = sort(intensities, 'descend');
        
        %calculate added points
        for ix = sortorder'
            %if this is a shaft, shafts veto over a long distance, everything else vetos over a normal distance
            %if this is a core point, everything vetos over a short distance
            %otherwise everything vetos over a normal distance
            distThreshShaft = distThresh; distThreshOther = distThresh;
            if pointTypes(ix)==4
                distThreshShaft = 1.5*distThresh;
            elseif pointTypes(ix)==2 || newPt(ix)
                distThreshOther = distThresh/1.5; distThreshShaft = distThresh/1.5;
            end
            
            for z2 = i3(ix)+1:min(i3(ix)+3, size(refIM.IM,3))
                %all potential neighbours
                n4 = (i3==z2) | addedTF(:,z2);
                select = n4 & pointTypes==4;
                n4(select) = sqrt((i2(ix)-i2(select)).^2 + (i1(ix)-i1(select)).^2)<distThreshShaft;
                select = n4 & pointTypes~=4;
                n4(select) = sqrt((i2(ix)-i2(select)).^2 + (i1(ix)-i1(select)).^2)<distThreshOther;
                if ~any(n4) && bw(i1(ix),i2(ix),z2) %check if the point is in BW
                    added{ix} = [added{ix} z2];
                    addedTF(ix,z2) = true;
                else
                    break
                end
            end
            for z2 = i3(ix)-1:-1:max(1, i3(ix)-3)
                %all potential neighbours
                n4 = (i3==z2) | addedTF(:,z2);
                select = n4 & pointTypes==4;
                n4(select) = sqrt((i2(ix)-i2(select)).^2 + (i1(ix)-i1(select)).^2)<distThreshShaft;
                select = n4 & pointTypes~=4;
                n4(select) = sqrt((i2(ix)-i2(select)).^2 + (i1(ix)-i1(select)).^2)<distThreshOther;
                if ~any(n4) && bw(i1(ix),i2(ix),z2) %check if the point is in BW
                    added{ix} = [added{ix} z2];
                    addedTF(ix,z2) = true;
                else
                    break
                end
            end
        end
        
        nh_size = 3*distThresh; %maximum radius of seed influence, in pixels. This should be >2*dist_thresh
        I = cell(4*length(i1),1); J = I; V = I;
        cellIX = 1;
        for zz = 1:size(bw,3)
            i_ixs = find(i3==zz  |  cellfun(@(x)any(x==zz), added));
            bw_padded = padarray(bw(:,:,zz), [nh_size nh_size]);
            shaft_padded = padarray(shaft(:,:,zz), [nh_size nh_size]);
            nonshaft_padded = bw_padded & ~shaft_padded;
            IMpadded = padarray(refIM.IM(:,:,zz), [nh_size nh_size]);
            
            nhood = false(2*nh_size+1,2*nh_size+1,length(i_ixs));
            nhoodbw = false(2*nh_size+1,2*nh_size+1,length(i_ixs));
            nhoodIM = zeros(2*nh_size+1,2*nh_size+1,length(i_ixs));
            for seed = 1:length(i_ixs)
                if pointTypes(i_ixs(seed))==4 %i3, not zz here
                    nhood(:,:,seed) = shaft_padded((0:2*nh_size)+i1(i_ixs(seed)), (0:2*nh_size)+i2(i_ixs(seed)));
                else
                    nhood(:,:,seed) = nonshaft_padded((0:2*nh_size)+i1(i_ixs(seed)), (0:2*nh_size)+i2(i_ixs(seed)));
                end
                nhoodbw(:,:,seed) = bw_padded((0:2*nh_size)+i1(i_ixs(seed)), (0:2*nh_size)+i2(i_ixs(seed)));
                nhoodIM(:,:,seed) = IMpadded((0:2*nh_size)+i1(i_ixs(seed)), (0:2*nh_size)+i2(i_ixs(seed)));
            end
            nhoodIM = nhoodIM./repmat(max(max(nhoodIM,[],1),[],2), size(nhood,1),size(nhood,2),1);
            influence = zeros(2*nh_size+1,2*nh_size+1,length(i_ixs));
            influence(nh_size+1,nh_size+1,:) = nhoodIM(nh_size+1,nh_size+1,:);
            
            SE = ones(3)/9;
            %Flood fill:
            for i = 1:2*distThresh
                influence = influence + convn(influence.*nhoodIM, SE, 'same').*nhood;
            end
            for i = 1:distThresh
                influence = influence + convn(influence.*nhoodIM, SE, 'same').*nhoodbw;
            end
            influence = influence./sum(sum(influence,1),2); %so seeds at edges don't get swamped
            influence(:,:,pointTypes(i_ixs)==4) = 2*influence(:,:,pointTypes(i_ixs)==4); %the shaft points are bigger; compensate
            influence(isnan(influence)) = 0;
            
            %accumulate indices for sparse matrix, #voxels x #seeds
            for seed = 1:length(i_ixs)
                [row, col, val] = find(influence(:,:,seed));
                seg_ixs = sub2ind(size(bw), row+i1(i_ixs(seed))-(nh_size+1),col+i2(i_ixs(seed))-(nh_size+1),zz*ones(size(col)));
                %accumulators for sparse matrix:
                I{cellIX} = seg_ixs; J{cellIX} = i_ixs(seed)*ones(size(seg_ixs)); V{cellIX} = val; cellIX = cellIX+1;
            end
        end
        
        %error checking: make sure every pixel is influenced by a seed
        no_seed = bw;
        no_seed(unique(cell2mat(I))) = false;
        if any(no_seed(:))
            no_seed = bwareaopen(no_seed, 4, 4);%discard isolated pixels
            [L,num] = bwlabeln(no_seed);
            i1B = nan(num,1); i2B = nan(num,1); i3B = nan(num,1); pTB = nan(num,1);
            for Lix = 1:num
                L_inds = find(L(:)==Lix);
                [~,maxind] = max(refIM.IM(L_inds));
                [i1B(Lix), i2B(Lix), i3B(Lix)] = ind2sub(size(refIM.IM), L_inds(maxind));
                pTB(Lix) = refIM.labels(L_inds(maxind));
                
                I{cellIX} = L_inds; J{cellIX} = (length(i1)+Lix)*ones(size(I{cellIX})); V{cellIX} = ones(size(I{cellIX}));
                cellIX = cellIX+1;
            end
            i1 = [i1 ; i1B]; i2 = [i2 ; i2B]; i3 = [i3 ; i3B]; pointTypes = [pointTypes ; pTB];
            added{length(i1)} = []; addedTF(length(i1),1) = false; newPt(length(i1)) = false;
        end
        
        i = cell2mat(I); j = cell2mat(J); v = cell2mat(V);
        %apply boundary sharpening
        v = v.^opts.sharpness;
        sumIM = sum(sparse(i,j,v, numel(bw), length(i1)),2);
        select = v>(0.05.*sumIM(i));
        i = i(select); j = j(select); v=v(select);
        v = v./sumIM(i);
        refIM.seg = sparse(i,j,v.*refIM.IM(i), numel(bw), length(i1));
        
        %force contiguous segments
        for ix = 1:size(refIM.seg,2)
             lininds = find(refIM.seg(:,ix));
             inds = sub2ind(size(refIM.IM), i1(ix),i2(ix),i3(ix));
             l2 = 1;
             growing = true;
             while growing
                inds = unique([inds; inds+1; inds-1; inds+size(refIM.IM,1); inds-size(refIM.IM,1); inds+size(refIM.IM,1)*size(refIM.IM,2); inds-size(refIM.IM,1)*size(refIM.IM,2)]);
                inds = intersect(lininds,inds);
                growing = length(inds)-l2;
                l2 = length(inds);
             end
             reject = setdiff(lininds, inds);
             refIM.seg(reject,ix) = 0;
        end
        
        %remove '0' segments
        select = any(refIM.seg,1);
        refIM.seg = refIM.seg(:,select);
        i1 = i1(select);
        i2 = i2(select);
        i3 = i3(select);
        pointTypes = pointTypes(select);
        newPt = newPt(select);
        
        %ensure that the segments add up to the reference image within the support
        valid = full(any(refIM.seg,2));
        refIM.seg(valid,:) = refIM.seg(valid,:).*(refIM.IM(valid)./(sum(refIM.seg(valid,:),2))); %, 1, size(refIM.seg,2));
        refIM.pts = [i1 i2 i3];
        refIM.ptType = pointTypes;

        calcBlock;
        redraw;
    end

    function calcBlock
        nS = size(refIM.seg,2);
        RGB = rand(3,nS); RGB = RGB./repmat(sum(RGB,1), 3,1);
        S_RGB = sqrt([refIM.seg*RGB(1,:)' refIM.seg*RGB(2,:)' refIM.seg*RGB(3,:)']);
        S_RGB = 1.5* S_RGB./max(S_RGB(:));
        imBlock = reshape(full(S_RGB), [size(refIM.IM,1) size(refIM.IM,2) size(refIM.IM,3) 3]);
        if isfield(refIM, 'mask') && any(refIM.mask(:)<1)
            imBlock = imBlock + repmat(0.2*(1-refIM.mask), 1,1, size(refIM.IM,3), 3);
        end
    end

end