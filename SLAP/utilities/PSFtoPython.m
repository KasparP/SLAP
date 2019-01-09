function PSFtoPython (filename)
%converts a psf file to a python-readable format
if ~nargin
    [fns, dr] = uigetfile('E:/Slapmidata/PSF/*.psf', 'multiselect', 'on');
    if ~iscell(fns)
        fns = {fns};
    end
else
    [dr, fns{1}] =fileparts(filename);
end

for fnum = 1:length(fns)
    data = [];
    load([dr filesep fns{fnum}], '-mat')
    
    %base of struct
    fieldnames = {'p2vx', 'p2vy','v2px','v2py','p2Z','v2Z'};
    for i =1:length(fieldnames)
        data.(fieldnames{i}) = struct(data.(fieldnames{i}));
    end
    
    
    %interpolants for each line
    for linenum = 1:4
        linename = ['line' int2str(linenum)];
        fieldnames = {'interpDc_fromT',   'interpRc_fromT',    'interpId',    'interpIr',    'interpR_fromT',    'interpZ_fromT'};
        for i =1:length(fieldnames)
            data.(linename).(fieldnames{i}) = struct(data.(linename).(fieldnames{i}));
        end
        
        data.calib.PSF.(linename).RZ2I = struct(data.calib.PSF.(linename).RZ2I);
    end
    
    save([dr filesep fns{fnum}(1:end-4) '.h5'], 'data', '-v7.3');
end
   