function getStimDelay
    %get filenames
    [fns, dr]= uigetfile('*.mat','Select a few contiguous scandata files', 'MultiSelect', 'On');
    
    %load stimtimes
    disp('Determining delay... (this goes faster if you select fewer scandata files)')
    times1 = [];
    for i = 1:length(fns)
        load([dr fns{i}])
        times1 = [times1 scandata.metadata.timeNow];
    end
    
    %read in stimtimes
    stimfiles = dir([dr filesep '*Timings*']);
    if ~isempty(stimfiles)
        stimList = [];
        for sf = 1:length(stimfiles)
            SF = load([dr filesep stimfiles(sf).name]);
            stimList = [stimList ; SF.time_stamp]; %#ok<AGROW>
        end
    else
        error('no timing fiiles in folder');
    end
    
    dist = repmat(stimList(:,1), 1, length(times1)) - repmat(times1, size(stimList,1),1);
    distP = dist; distM = dist;
    distP(distP<0) = nan; distM(distM>0) = nan;
    mindistP = nanmin(distP, [], 1);
    mindistM = nanmax(distM,[],1);
    if var(mindistM)<var(mindistP)
        StimDelay = nanmedian(mindistM);
        delayvariance = nanvar(mindistM);
    else
        StimDelay = nanmedian(mindistP);
        delayvariance = nanvar(mindistP);
    end
    StimDelay
    if delayvariance>1e-12
        warning('Trouble identifying stimulus delay')
        keyboard
    end

    %resave the timings files
    disp('adding stimulus identities to Timing File...')
     for sf = 1:length(stimfiles)
        load([dr filesep stimfiles(sf).name])
        save([dr filesep stimfiles(sf).name], 'reftime', 'stim_count', 'time_stamp', 'StimDelay', 'delayvariance')
     end
end