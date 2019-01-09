function rigidImage = processRasterStack (doMotion)

basedir = 'C:\Users\podgorskik\Dropbox (HHMI)\Kaspar-Emilia\2016-08-11';
[fns, dr] = uigetfile([basedir filesep '*.tif'], 'multiselect', 'on');
[savefn savedr] = uiputfile('*.tif', 'Select a name for the merged stack');

if ~nargin
    doMotion = ~strcmpi(input('Correct for sample motion? [Y]|N >>', 's'), 'N');
end


for fnum = 1:length(fns)
    disp(['Processing imageset: ' int2str(fnum) ' of ' int2str(length(fns))])
    fileName=[dr filesep fns{fnum}];
    
    vol=double(ScanImageTiffReader(fileName).data());
    
    Lx = size(vol,1);
    Ly = size(vol,2)/2;
    N = size(vol,3);
    
    IMfastY = vol(:,1:Ly,:);
    IMfastX = rot90(vol(:, Ly + (1:Ly),:));
    
    if fnum==1
        rigidImage = zeros(Lx,Ly, length(fns));
    end
   
    
    if doMotion
        concensus = makeRigidImage(IMfastX, IMfastY); %rigidimages alignment
    else
       %only correct for slow intensity noise, etc 
       
       %adjust for bidi scan phase errors
       [IMfastX, shiftX(fnum)] = fixphase(IMfastX, 1);
       [IMfastY, shiftY(fnum)] = fixphase(permute(IMfastY, [2 1 3]),0); IMfastY = permute(IMfastY, [2 1 3]);
       
       %first, get the offset between the X and Y blocks
       medFX = double(median(IMfastX,3)).^3; %cube to emphasize signal over low-amplitude noise
       medFY = double(median(IMfastY,3)).^3; %cube to emphasize signal over low-amplitude noise
       

       
       output = dftregistration_wprior(fft2(imgaussfilt(medFX,2)), fft2(imgaussfilt(medFY,2)), 1);
       
       rowshift = output(3);
       colshift = output(4);
       
       %overlap images
       buffer = 1+max(abs([rowshift colshift]));
       concensusX = nan(Lx+2*buffer,Ly+2*buffer);
       concensusY = nan(Lx+2*buffer,Ly+2*buffer);
       concensusX(buffer+(1:Lx), buffer+(1:Ly)) = double(median(IMfastX,3));
       concensusY(buffer+rowshift+(1:Lx), buffer+colshift+(1:Ly)) = double(median(IMfastY,3));
       
       figure, imshow(concensusX,[])
       for row = 1:size(concensusX,1)
            valid = ~isnan(concensusX(row,:)+concensusY(row,:));
            concensusX(row,:) = concensusX(row,:) - (mean(concensusX(row,valid)) - mean(concensusY(row,valid)));
       end
       %figure, imshow(concensusX,[])
       
       for col = 1:size(concensusY,2)
           valid = ~isnan(concensusX(:,col)+concensusY(:,col));
           concensusY(:,col) = concensusY(:,col) + (mean(concensusX(valid,col)) - mean(concensusY(valid,col)));
       end
       %figure, imshow(concensusY,[])
       
       concensus = (concensusX+concensusY)/2; %arithmetic mean
       %concensus = sqrt(max(0,concensusX).*max(0,concensusY)); %geometric mean
       concensus = concensus(buffer+(1:Lx), buffer+(1:Ly));
       
       %figure, imshow(concensus,[]);
    end
    rigidImage(:,:,fnum) = concensus;
end

rigidImage = rigidImage - (-1);
rigidImage = uint16(rigidImage.*floor(double(intmax('uint16'))./max(rigidImage(:))));

options.overwrite = true;
errorcode = saveastiff(rigidImage, [savedr savefn], options);
if errorcode
    keyboard
end


end


function [IMset, rowshift] = fixphase(IMset, rowshift)
if nargin<2 || isempty(rowshift)    
    H2 = floor((size(IMset,1)-1)/2);
    IM1 = median(IMset(1:2:2*H2-1,:,:),3);
    IM2 = median(IMset(2:2:2*H2,:,:),3);
    IM3 =  median(IMset(3:2:2*H2+1,:,:),3);
    
    IM1fft = fft2(IM1.^3);
    IM2fft = fft2(IM2.^3);
    IM3fft = fft2(IM3.^3);
    
    p1 = dftregistration(IM1fft(1,:), IM2fft(1,:), 4);
    p3 = dftregistration(IM3fft(1,:), IM2fft(1,:), 4);
    
    rowshift = ((p1(4) + p3(4)) /2);
end
if abs(rowshift>0) && abs(rowshift)<2.5
    IMset(1:2:end,:,:) = imtranslate(IMset(1:2:end,:,:), [-rowshift 0]);
end
end