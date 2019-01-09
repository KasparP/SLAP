function hFig = scandataPlot(scandata, opts)
if nargin<2
    opts = [];
end
if nargin %scandata supplied
    if scandata.metadata.do3D
        hFig = scandataPlot3D(scandata, opts);
    else
        hFig = scandataPlot2D(scandata, opts);
    end
else %load files
    basedir = 'E:\SLAPMiData';
    [fns, dr] = uigetfile([basedir filesep '*.mat'], 'Select your scandata files', 'multiselect', 'on');
    if ~iscell(fns)
        fns = {fns};
    end
    hFig = [];
    for ix = 1:length(fns)
        load([dr filesep fns{ix}]);
        if scandata.metadata.do3D
            hFig = [hFig scandataPlot3D(scandata)]; %#ok<AGROW>
        else
            hFig = [hFig scandataPlot2D(scandata)]; %#ok<AGROW>
        end
    end
end
end

