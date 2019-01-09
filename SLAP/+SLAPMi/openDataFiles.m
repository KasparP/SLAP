function [metaData, pmtFile, galvoData, saveData, Zdata] = openDataFiles(fileName)
%% file names
metaFile = [fileName '.txt'];
pmtFile = [fileName '.pdat'];
galvoFile = [fileName '.gdat'];
saveFile = [fileName '.mdat'];
Zfile =  [fileName '.zdat'];

%% read savedata
saveData = load(saveFile, '-mat');

%% read meta data file
fid = fopen(metaFile);
rows = textscan(fid,'%s','Delimiter','\n');
fclose(fid);
rows = rows{1};

%% process meta data into a struct
befvars = []; %#ok
befvars = who();
rows = strcat(rows,';');
cellfun(@eval,rows(2:end));
vars = who();
vars = setdiff(vars,befvars);
cellfun(@(x)evalin('caller', ['metaData.' x '=' x ';']),vars);

%% read pmt data file
% fid = fopen(pmtFile);
% dat = fread(fid, '*uint64');
% fclose(fid);
% 
% %% process pmt data file
% % each element is a u64 packed with 4 i16 values:
% % [ch0 sample 1, ch0 sample 2, ch1 sample 1, ch1 sample 2]
% dat = typecast(dat, 'uint32');
% switch ch
%     case 0 %all channels
%         pmtData.data = double([typecast(dat(1:2:end), 'int16') typecast(dat(2:2:end), 'int16')]) * (metaData.inputRange/2) / 2^15;
%     case 2
%         pmtData.data = double(typecast(dat(2:2:end), 'int16')) * (metaData.inputRange/2) / 2^15;
%     case 1
%         pmtData.data = zeros(length(dat), 2);
%         pmtData.data = double(typecast(dat(1:2:end), 'int16')) * (metaData.inputRange/2) / 2^15;
% end
%         
% pmtData.dt = 1/250000000;
% clear cdat;

%% read galvo data file
fid = fopen(galvoFile);
dat = fread(fid, '*uint64');
fclose(fid);

%% process galvo data file
% data consists of pairs of corresponding values. The first element is a
% u64 packed with 4 i16 values giving a sample of each galvo pos channel.
% the second element is packed with 2 14bit samples and a 36 bit clock.
posDat1 = double(typecast(dat(1:2:end), 'int16'));
dat = dat(2:2:end);
galvoData.t = double(bitshift(dat,-28)) * 8e-9;
dat = double([typecast(uint16(bitand(bitshift(dat,2),65532)),'int16')...
    typecast(uint16(bitand(bitshift(dat,-12),65532)),'int16')]);

galvoData.data = [posDat1(1:4:end) posDat1(2:4:end) posDat1(3:4:end) posDat1(4:4:end) dat] * 2.5 / 2^15;

%% read Z data file
Zdata = [];
if exist(Zfile, 'file')
    load(Zfile, '-mat');    
end
end

