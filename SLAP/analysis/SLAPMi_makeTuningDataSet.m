function dataset = SLAPMi_makeTuningDataSet

%select files
disp('Select your saved RECON files for analysis');
[fns, dr] = uigetfile('*.mat', 'Select your saved RECON files for analysis', 'multiselect', 'on');
if ~iscell(fns)
    fns = {fns};
end
fns = sort_nat(fns);

% %ALIGNMNET data
% alignFN = dir([dr filesep '/../*ALIGNDATA*.mat']);
% if length(alignFN)==1
%     alignDR =  alignFN.folder;
%     alignFN = alignFN.name;
% else
%     disp('There should be exactly one ALIGNDATA file in the problemdata directory. Select your aligndata manually, or hit cancel to use default settings')
%     [alignFN, alignDR] = uigetfile([dr filesep '/../*ALIGNDATA*.mat'], 'Select Aligndata or cancel to use all files and default plane');
% end     
        
%SEGMENTATION data
disp('select your reference Segmentation _SEG file');
[fnRef, drRef] = uigetfile([dr filesep '*.mat'], 'select your reference Segmentation _SEG file');

%LOAD STIMULUS
stimfiles = dir([dr  filesep '..' filesep '*Timings*']);
if ~isempty(stimfiles)
    stimList = [];
    for sf = 1:length(stimfiles)
        SF = load([dr filesep '..' filesep stimfiles(sf).name]);
        stimList = [stimList ; SF.time_stamp]; %#ok<AGROW>
    end
else
    disp('Could not find Stimulus Timing files, please select manually');
    [stimfileNames, stimfileDR] = uigetfile([dr filesep '*Timing*.mat'], 'Could not find Stimulus Timing files, please select manually', 'multiselect', 'on');
    if ~iscell(stimfileNames)
        stimfileNames = {stimfileNames};
    end
    stimList = [];
    for sf = 1:length(stimfileNames)
        SF = load([stimfileDR filesep stimfileNames{sf}]);
        stimList = [stimList ; SF.time_stamp]; %#ok<AGROW>
    end
end

%load each file, identify stimulus time
warningshown = false;

for fnum = 1:length(fns)
    disp(['Loading file: ' int2str(fnum) ' of ' int2str(length(fns)) '...'])
    try
        load([dr fns{fnum}]);
    catch
        keyboard
        continue
    end
    if fnum==1
        if isfield(sys_recon.output,'SLM_mask')
            sys_recon.opts.SLM_threshold = prctile(sys_recon.output.SLM_mask(:),5);
        end
        if isfield(sys_recon.output,'SLM_mask') && ~isfield(sys_recon.output,'OutsideSLMSeeds')
        sys_recon.opts.SLM_threshold = prctile(sys_recon.output.SLM_mask(:),5);
        InsideSLM = sys_recon.output.SLM_mask > sys_recon.opts.SLM_threshold;
        nPlanes = size(sys_recon.input.ref_image,3);
        SLM_extended = repmat(InsideSLM,1,1,nPlanes);
        S_temp = sys_recon.input.S(:,[1:end-2,end]).*SLM_extended(:);
        sys_recon.output.OutsideSLMSeeds = ~any(S_temp);
        else
            OutsideSLMSeeds = sys_recon.output.OutsideSLMSeeds;
        end
        dataset.P = sys_recon.input.P;
        dataset.X = nan([size(sys_recon.output.F) length(fns)]);
        dataset.Ysum = nan([size(sys_recon.input.y,2) length(fns)]);
        dataset.motion.error = nan([length(sys_recon.input.scandata.motion.noREF.error) length(fns)]);
        dataset.motion.shifts = nan([3 size(sys_recon.input.y,2) length(fns)]);
        dataset.dPhotons = dataset.X;
        dataset.fusedInto = nan(size(dataset.X,1)-2,length(fns));
        dataset.F0_photons = nan(size(dataset.X,1)-2,length(fns));
        dataset.Z = zeros(1, length(fns));
        dataset.spikes = nan([size(sys_recon.output.spikes) length(fns)]);
        dataset.stimTime = nan(1,length(fns));
        dataset.alignError = nan(1,length(fns));
    end
    indices = false(size(sys_recon.input.ref_image));
    indices(:,:,(size(indices,3)+1)/2) = 1;
    centerPlaneSeeds = any(sys_recon.input.S(indices(:),:));
    valid = sys_recon.input.ref.fusedInto>0;
    dataset.fusedInto(1:length(sys_recon.input.ref.fusedInto),fnum) = sys_recon.input.ref.fusedInto';
    dataset.X(valid,:,fnum) = sys_recon.output.F(sys_recon.input.ref.fusedInto(valid), :);
    dataset.Ysum(:,fnum) = sum(sys_recon.input.y,1);
    dataset.motion.error(:,fnum) = sum(sys_recon.input.scandata.motion.noREF.error,1);
    shiftfacs = pca(sys_recon.input.scandata.motion.noREF.shift);
    dataset.motion.shifts(:,:,fnum) = shiftfacs';
    PS = full(sum(sys_recon.output.PS)');
    dataset.F0_photons(valid,fnum) = PS(valid);
    dataset.dPhotons(valid,:,fnum) = PS(valid).*sys_recon.output.F(sys_recon.input.ref.fusedInto(valid), :);
    if isfield(sys_recon.input.scandata.motion, 'dXYZ')
        dataset.Z(fnum) = sys_recon.input.scandata.motion.dXYZ(3);
    elseif ~warningshown
        warningshown = true;
        warning('Solve was performed with defunct motion correction, you should rerun PSXdata! Assuming 0 offset.');
    end
    dataset.X(~centerPlaneSeeds,:,fnum) = nan;
    dataset.F0_photons(~centerPlaneSeeds,fnum) = nan;
    dataset.dPhotons(~centerPlaneSeeds,:,fnum) = nan;
    dataset.spikes(valid,:,fnum) = sys_recon.output.spikes(sys_recon.input.ref.fusedInto(valid), :);
    dataset.spikes(~centerPlaneSeeds,:,fnum) = nan;
    %stimulus timing
    dataset.stimTime(fnum) = sys_recon.input.scandata.metadata.timeNow;
    try
        dataset.alignError(fnum) = sys_recon.input.scandata.motion.error;
    catch
        disp('This solve was performed on a defunct problemdata file'); keyboard;
    end
end
dataset.X(OutsideSLMSeeds,:,:) = nan;
dataset.spikes(OutsideSLMSeeds,:,:) = nan;
dataset.dPhotons(OutsideSLMSeeds,:,fnum) = nan;

%get stimulus IDs
dist = repmat(stimList(:,1), 1, length(dataset.stimTime)) - repmat(dataset.stimTime, size(stimList,1),1);
distP = dist; distM = dist;
distP(distP<0) = nan; distM(distM>0) = nan;
mindistP = nanmin(distP, [], 1);
mindistM = nanmax(distM,[],1);
if nanvar(mindistM)<nanvar(mindistP)
    dataset.stimulus.stimDelay = nanmedian(mindistM);
    dataset.stimulus.delayVariance = nanvar(mindistM);
else
    dataset.stimulus.stimDelay = nanmedian(mindistP);
    dataset.stimulus.delayVariance = nanvar(mindistP);
end
dataset.stimulus.stim = nan(length(fns),1);
dataset.stimulus.timeError = nan(length(fns),1);
dataset.stimulus.stimTime = nan(length(fns),1);
for fnum = 1:length(fns)
    [dataset.stimulus.timeError(fnum), minIx] =  min(abs(dataset.stimTime(fnum) - stimList(:,1) + dataset.stimulus.stimDelay));
    dataset.stimulus.timeError(fnum) = 86400*dataset.stimulus.timeError(fnum); %convert to seconds
    dataset.stimulus.stim(fnum) = stimList(minIx,2);
    dataset.stimulus.stimTime(fnum) = stimList(minIx,1);
end
figure,  scatter(86400*(dataset.stimulus.stimTime-min(dataset.stimulus.stimTime)), dataset.stimulus.stim);
xlabel('Stimulus time (s)');
ylabel('Stimulus ID');
title(['mean absolute stimulus timing error = ' num2str(mean(abs(dataset.stimulus.timeError))) ' seconds']);
refIM = [];
load([drRef filesep fnRef]);
dataset.filenames = fns;
dataset.dr = dr;
dataset.refIM = refIM;
dataset.solverOpts = sys_recon.opts;
dataset.SLM = sys_recon.output.SLM_mask;
dataset.Zcenter = mode(dataset.Z) + ceil(size(dataset.refIM.IM,3)/2);  %should this be -mode(...)??

if ~all(isnan(dataset.alignError))
    ETH = imtophat(dataset.alignError, ones(17,1));
    cutoff = 5*std(ETH(ETH<prctile(ETH, 90)));
    if any(ETH>cutoff)
        figure('name', 'PROBLEMDATA ALIGNMENT ERRORS');
        plot(ETH); hold on, plot(dataset.alignError, ':');
        hold on, plot([0 length(ETH)], [cutoff cutoff], 'r', 'linewidth', 2);
        xlabel('File number'); ylabel('Error');
        
        [dataset] = QCSLAPMiDataset(dataset);
    end
end
    
try
    mkdir([dr filesep 'SLAPMi Tuning']);
end
save([dr filesep 'SLAPMi Tuning' filesep 'TUNING_' int2str(length(fns)) 'FILES_' fns{1}(1:end-4)], 'dataset','-v7.3')

% %User Quality Control
% [dataset, selection] = QCSLAPMiDataset(dataset);
% if sum(selection<length(selection))
%     save([dr filesep 'SLAPMi Tuning' filesep 'TUNING_' int2str(sum(selection)) 'FILES_' fns{1}(1:end-4)], 'dataset','-v7.3')
% end
end