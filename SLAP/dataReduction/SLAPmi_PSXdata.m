function SLAPMi_PSXdata(optsin)
%assembles data for activity reconstruction in a convenient format

%OPTIONS
opts.doFastMotion = true;   tooltips.doFastMotion = 'Correct for fast motion with cross-correlation';
opts.distThreshSeg = 1.4;   tooltips.distThreshSeg = 'Spatial scale of reference segmentation';
opts.maxSeedsRef = Inf;     tooltips.maxSeedsRef = 'Maximum number of seeds in the reference segmentation; usually Inf';
opts.manualROIs = false;    tooltips.manualROIs = 'Set to true to load a SLAPMI_ROI file to merge the autosegmentation';

if nargin %UPDATE WITH USER-SPECIFIED OPTIONS
    for field = fieldnames(optsin)'
         opts.(field{1}) = optsin.(field{1});
    end
else
      opts = optionsGUI(opts, tooltips);
end

basedir = 'E:\SLAPmidata\';

%ask for refIM/segmentation data
[fnRef, drRef] = uigetfile([basedir '*.mat'] , 'Select your refIM');
refIM = []; load([drRef fnRef]);

%select scandata
[fns, dr] = uigetfile([drRef filesep '..' filesep '*.mat'], 'Select your reduced scandata', 'multiselect', 'on');
if ~iscell(fns)
    fns = {fns};
end

%if ~isfield(refIM, 'labels')
    opts.fn = [drRef fnRef];
    refIM = SLAPMi_ClassifyPixels(opts);
    if opts.manualROIs
        [ROIfn, ROIdr] = uigetfile([drRef filesep '..' filesep '*.mat'], 'Select your ROI file');
        ROIs = [];
        load([ROIdr filesep ROIfn], 'ROIs');
        refIM = SLAPMi_ROI_seg(ROIs, refIM);
        disp(['Saving file to: ' drRef fnRef(1:end-8) '_SEG_ROIs.mat'])
        save([drRef fnRef(1:end-8) '_SEG_ROIs'], 'refIM')
    end
%end
if ~isfield(refIM, 'seg')
    disp('Performing reference segmentation...');
    segopts.maxseeds = opts.maxSeedsRef;
    segopts.visualize = true;
    segopts.distThreshUM = opts.distThreshSeg; %The scale of the segmentation, in microns
    segopts.savefile = [drRef fnRef(1:end-4)];
    if max(refIM.labels(:))<4
        input('The state of the art segmentation is Core-Shell, your segmentation is out of date. Continue? >>')
        refIM = SLAPMi_segmentOnce(refIM, segopts);
    else
        refIM = SLAPMi_segmentCoreShell(refIM, segopts);
    end
    
    %get the mask
    scandata = [];
    load([dr fns{end}]);
    refIM.mask = getSLMmask(scandata, refIM);
    
    disp(['Saving file to: ' drRef fnRef(1:end-8) '_SEG.mat'])
    save([drRef fnRef(1:end-8) '_SEG'], 'refIM')
    
    %Segmentation editor?
    if strcmp(questdlg('Run segmentationEditor?', 'SEGEDIT', 'Yes', 'No', 'No'), 'Yes')
            segmentationEditor(refIM);
            disp('To continue, rerun SLAPMi_PSXdata with your saved segmentation.')
            return
    end
end

% set alignment type
alignData = struct;
alignfun = @SLAPMi_alignDTW_2D_SLM;
alignOffset= nan(length(fns),3);
alignError = nan(length(fns),1);
alignCensor = nan(length(fns),1);
alignFast = nan(length(fns),1);
scanTime = nan(length(fns),1);

hFig = []; hFigM = [];
for fnum = 1:length(fns)
    processFile(fnum,[]);
end

%look for motion smoothness violations
if length(fns)>1
motionspline = csaps(scanTime, alignOffset',1-1e-11);
vals= fnval(motionspline, scanTime);
E = sum((alignOffset.*[1 1 2] - vals'.*[1 1 2]).^2, 2);
violAlign = E>3;
if any(violAlign) %figure showing motion violations
    %make a figure that reports alignment error
    hAlignFig = figure('Name', 'Alignment Violations');
    
    subplot(2,1,1)
    for dim = 1:3
        hold on; scatter(scanTime, alignOffset(:,dim));
    end
    legend({'X','Y','Z'})
    for dim = 1:3
        hold on; scatter(scanTime(violAlign), alignOffset(violAlign,dim), 'markeredgecolor', 'none', 'markerfacecolor', 'r');
        plot(linspace(min(scanTime),max(scanTime), 1000), fnval(motionspline, linspace(min(scanTime),max(scanTime), 1000)));
    end
    xlabel('Time'); ylabel('Position (pixels)'); 
    
    subplot(2,1,2)
    scatter(scanTime, alignError, 'markerfacecolor', 'b')
    set(gca, 'ylim', [0 Inf]); xlabel('Time'); ylabel('Alignment Error (a.u.)')
    
    for fnum = find(violAlign)' %REPROCESS VIOLATIONS WITH FORCING
        processFile(fnum, round(vals(:,fnum)));
    end
    
    figure(hAlignFig)
    subplot(2,1,1);
     for dim = 1:3
        hold on; scatter(scanTime(violAlign), alignOffset(violAlign,dim), 'markeredgecolor', 'none', 'markerfacecolor', 'g');
     end
     subplot(2,1,2);
     hold on; scatter(scanTime(violAlign), alignError(violAlign), 'markeredgecolor', 'none', 'markerfacecolor', 'g');
end

%Save alignment dataset
alignData.scanTime = scanTime; alignData.alignError = alignError; alignData.alignOffset = alignOffset; alignData.fastAlignCensor = alignCensor; alignData.fastAlignError = alignFast;
alignData.fileNames = fns;
save([dr filesep fns{1}(1:end-11) '-' int2str(length(fns)) 'files_ALIGNDATA.mat'], 'alignData', '-v7.3');

hCensor = figure('Name', 'Censored Motion'); subplot(1,2,1);  imagesc(alignData.fastAlignCensor); subplot(1,2,2); imagesc(alignData.fastAlignError);
saveas(hCensor, [dr filesep fns{1}(1:end-11) '-' int2str(length(fns)) 'files_CENSOR.fig']);

ETH = imtophat(alignData.alignError, ones(17,1));
cutoff = 5*std(ETH(ETH<prctile(ETH, 90)));
if any(ETH>cutoff)
    figure('name', 'PROBLEMDATA ALIGNMENT ERRORS'); 
    plot(ETH); hold on, plot(alignData.alignError, ':');
    hold on, plot([0 length(ETH)], [cutoff cutoff], 'r', 'linewidth', 2);
    xlabel('File number'); ylabel('Error');
    msgbox(['Some of the generated problemdata files failed to align within tolerance:' ; reshape(fns(ETH>cutoff), [],1)]);
    disp(reshape(fns(ETH>cutoff), [],1))
end
end
    function processFile(fnum, forceAlign)
        fn = fns{fnum};
        load([dr fn]);
        disp(['Processing: ' fn])
        try %#ok<TRYNC> %delete the figures from the previous file
            delete(hFig);
            delete(hFigM);
        end
        
        %ALIGN
        scandata.refIMcoords = refIM.M.coords;
        if opts.doFastMotion
            [scandata,hFigM] = SLAPMi_Motion_noREF (scandata, opts);
            M = imtophat(mean(scandata.motion.noREF.error,1), ones(1,500));
            M = M-nanmean(M);
            thresh =  4*std(M(M<prctile(M,90)));
            M(isnan(M)) = Inf;
            scandata.motion.noREF.censor = smooth(M,5)>thresh;
            scandata.motion.noREF.error = M;
            if any(scandata.motion.noREF.censor)
                figure('name', ['Uncorrected Motion above threshold: ' fns{fnum}]), plot(M);
                hold on, plot(smooth(M,4));
                hold on, plot([1 length(M)], [thresh thresh], 'r');
                drawnow;
            end
            alignFast(fnum,1:length(scandata.frames)) = scandata.motion.noREF.error;
            alignCensor(fnum,1:length(scandata.frames)) = scandata.motion.noREF.censor;
        end
        
        alignopts.dense = false; %whether the reference image is dense
        alignopts.PSFtype = 'X'; 
        alignopts.Zrange = 3; 
        alignopts.Zres=  1; 
        alignopts.prior = 0.0005; 
        [scandata, P, S, hFig]= alignfun(scandata, refIM, forceAlign, alignopts);
        scanTime(fnum) = scandata.metadata.timeNow;
        alignOffset(fnum,:) = scandata.motion.dXYZ;
        alignError(fnum,:) = scandata.motion.error;

        try %#ok<TRYNC>
            saveas(hFig, [dr filesep fn(1:end-11) 'ALIGNMENT.fig']);
            disp(['Alignment Figure saved to: ' dr filesep fn(1:end-11) 'ALIGNMENT.fig']);
            visopts.doPCA = false;
        end
        
        for ch = 1 %:size(scandata.frames(1).pmtData,2) %TODO: add channel 2 processing
            if ch==2 %TODO
               %change S.IM to channel 2 
               
               %change S.seg to match S.IM
               S.seg(valid,:) = S.seg(valid,:).*(S.IM(valid)./sum(S.seg(valid,:),2));
               
               %adjust P appropriately
            end
            
            %Fuse dim or overlapping seeds, considering mask
            S = fuseSeeds(S);
            
            %discard seeds far from mask
            PS = P.P*(repmat(S.mask(:),length(P.coords{3}),1).*S.seg);
            
            sel = nansum(PS,1)>nansum(PS(:))./1e4;  %seeds must meet a minimum brightness
            includeBrightOutsideMask = true;
            if ~includeBrightOutsideMask
                inMask = any(repmat(S.mask(:)>0.9, size(S.IM,3),1).*S.seg); %seeds that are within the mask
                disp([int2str(sum(sel & inMask)) ' bright objects outside the mask were discarded'])
                sel = sel & inMask;
            end
            includeBorder = true;
            if includeBorder
                maskmask = S.mask(:)>0.8;
                inMask = any(repmat(maskmask, size(S.IM,3),1).*S.seg); %seeds that are within the mask
                outMask = any(repmat(~maskmask, size(S.IM,3),1).*S.seg);
                disp([int2str(sum(~sel & inMask & outMask)) ' dim objects at the border of the mask were included'])
                sel = sel | (inMask & outMask);
            end
            
            %figure('Name', [int2str(sum(select)) '/' int2str(sum(~select)) ' Accepted/Rejected seeds for: ' fn]),
            %imshow(imfuse(sqrt(reshape(full(sum(S.seg(:,select),2)), 1280,1280)), sqrt(reshape(full(sum(S.seg(:,~select),2)), 1280,1280))))
            S.seg(:,~sel) = 0;
            
            %ensure that the segments add up to the reference image within the support
            valid = any(S.seg,2);
            S.seg(valid,:) = S.seg(valid,:).*(S.IM(valid)./sum(S.seg(valid,:),2));
            
            disp(['There were ' int2str(sum(any(S.seg,1))) ' seeds in the reconstruction volume.'])
            
            %SAVE
            toSave.scandata = scandata;
            for f = 1:length(toSave.scandata.frames)
                toSave.scandata.frames(f).pmtData = toSave.scandata.frames(f).pmtData(:,ch); %discard other channels for solver
            end
            toSave.S = S;
            toSave.P = P;
            toSave.fnP = ['P-' P.hash]; %deprecate this if we continue to store the P matrix in the problem data
            fnP = [dr filesep toSave.fnP '.mat'];
            if ~exist(fnP, 'file')
                disp(['Saving P matrix to: ' fnP]);
                save(fnP, 'P', '-v7.3');
            end
            disp(['Saving ProblemData to: ' dr filesep fn(1:end-11) 'PROBLEMDATA_CH' int2str(ch) '.mat']);
            save([dr filesep fn(1:end-11) 'PROBLEMDATA_CH' int2str(ch) '.mat'], '-struct', 'toSave', '-v7.3');
        end
        
        disp(['Done file: ' fn])
        drawnow;
    end
end

function mask = getSLMmask(scandata, refIM)
%convert SLM voltages to transmission
Vslm = scandata.metadata.SLM.pattern;
SLMhigh  = abs(double(Vslm)-scandata.metadata.calib.SLM.lut(:,:,2))<abs(double(Vslm)-scandata.metadata.calib.SLM.lut(:,:,1));
Tslm = scandata.metadata.calib.SLM.T(:,:,2).* SLMhigh + scandata.metadata.calib.SLM.T(:,:,1).*~SLMhigh;
Tslm = max(Tslm, 1e-2); %Sometimes the SLM transmission seems to be underestimated, let's assume it's kind of large

%interpolate onto image space
scanfield = refIM.metadata.rois.RoiGroups.imagingRoiGroup.rois{1, 1}.scanfields(1);
scanfield.meshgrid = @(varargin)(meshgrid(refIM.M.coords.X, refIM.M.coords.Y));
pixelToRefTransform = scandata.metadata.SLM.pixelToRefTransform;
Tim = mapSLMToImage(Tslm, pixelToRefTransform, scanfield);

%smooth the mask
mask = imgaussfilt(Tim, 0.005/scanfield.pixelToRefTransform(1))';
end
