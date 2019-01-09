function refIM = SLAPMi_classifyPixels(opts)
%Label a SLAPMi reference image using Ilastik in batch mode

%DEFAULT OPTIONS
if ~nargin || ~isfield(opts, 'maxseeds')
    opts.maxseeds = 250; %maximum number of seed regions per XY plane; the sparsity of the activity
end

%REFERENCE IMAGE
if ~isfield(opts, 'fn')
    [fn, dr] = uigetfile('*.mat', 'Select imagestack data');
    opts.fn = [dr fn];
end
[dr, fn] = fileparts(opts.fn);
dr = [dr filesep];
opts.savefile = [dr fn '_ILSTK.mat'];

refIM = [];
load([dr fn]);

%ilastikPath = 'C:\Program Files\ilastik-1.2.0rc10\';
ilastikPath = 'C:\Program Files\ilastik-1.2.2rc10';
projectDr = 'E:\SLAPmidata\Ilastik\';
[projectFn, projectDr] = uigetfile([projectDr '*.ilp'], 'Select an ilastik project or classifier output');
if strcmpi(projectFn(end-3:end), '.ilp')
    projectPath = [projectDr projectFn];
    IMpath = [dr fn '.tif'];
    if ~exist(IMpath, 'file')
        [IMfn, IMdr] = uigetfile([dr filesep '*.tif'], 'TIF file not found, please select your image stack');
        IMpath = [IMdr IMfn];
    end
    if strfind(projectFn, 'autocontext')
        labelspath = [IMpath(1:end-4) '_Simple Segmentation Stage 2.h5'];
        if ~exist(labelspath, 'file')
            status = dos(['"' ilastikPath '\run-ilastik.bat" --headless --export_source="simple segmentation stage 2" --project="' projectPath '" ' IMpath]);
        end
    else
        labelspath = [IMpath(1:end-4) '_Simple Segmentation.h5'];
        if ~exist(labelspath, 'file')
            status = dos(['"' ilastikPath '\run-ilastik.bat" --headless --export_source="Simple Segmentation" --project="' projectPath '" ' IMpath]);
        end
    end
    labels = h5read(labelspath, '/exported_data');
elseif strcmpi(projectFn(end-2:end), '.h5') %LOAD LABELS
    labels = h5read([projectDr filesep projectFn], '/exported_data');
end

if size(labels,4)>1
    labels = squeeze(labels(1,:,:,:));
end
labels = permute(labels,[2 1 3]); %Ilastik transposes the first two dimensions
figure('name', ['Ilastik labels for ' fn]), imshow3D(double(labels))

refIM.labels = labels;