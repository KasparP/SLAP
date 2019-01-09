function SLAPMi_dataReductionServer(basedir)
%repeatedly poll a network location for SLAPMi data files to reduce

%on startup, ensure latest SLAPmi is pulled from repo
startdr = pwd;
cd(fileparts(which('slapmi')));
git('pull');
cd(startdr);

%options
opts.align = false;         tooltips.align = 'align to remove motion, good for in vivo data';
opts.channels = 1;         tooltips.channels = 'which channels to process; if empty do all channels';
opts.ignoreSaved = false;   tooltips.ignoreSaved = 'if data has been reduced before, reduce it again anyways';
opts.amp_setup = 'MPPC';    tooltips.amp_setup = 'supported: ''MPPC''';
opts.doRegression = false;  tooltips.doRegression = 'reject dark photons using regression';
opts.doNND = false;         tooltips.doNND = 'nonnegative deconvolution of PMT traces';
opts.doPlot = false;         tooltips.doPlot = 'show plots';
opts.doPCA  = false;         tooltips.doPCA = 'show a plot with NMF';
opts = optionsGUI(opts, tooltips);


if ~nargin || isempty(basedir)
    basedir = 'Z:\dataReduction';
end
disp(['SLAPMi_dataReductionServer is running...'])
while true
    try
        checkFolder(basedir, opts)
    catch ME
        warning(['SLAPMi_dataReductionServer ancounterted an error: ' datestr(now)])
        disp(ME)
    end
    pause(15);
    %disp(['SLAPMi_dataReductionServer finished a cycle: ' datestr(now)])
end

end

function checkFolder(D, opts)
    S = dir(D);
    fns = {S(3:end).name};
    folders = [S(3:end).isdir];
    for f_ix = find(folders)
        checkFolder ([D filesep fns{f_ix}]);
    end
    
    %for every gdat file, check if the reduced version exists
    S = dir([D filesep '*.gdat']);
    for f_ix = 1:length(S)
       if ~exist([D filesep S(f_ix).name(1:end-5) '_REDUCED.mat'], 'file')
           opts.path = [D filesep S(f_ix).name];
           SLAPMi_reduce(opts);
       end
    end
end
