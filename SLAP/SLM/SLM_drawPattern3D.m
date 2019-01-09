function patternOut = SLM_drawPattern3D(pattern,overdrive)
    %loads a pattern onto the SLM for registration
    %if PATTERN is a scalar, draw a circular aperture of diameter PATTERN
    %if PATTERN is a 512x512 binary, draw it as on/off
    %if PATTERN is the string 'test', draw the test pattern
    
    if any(isnan(pattern(:)))
        error('Nans in SLM pattern');
    end
    
    %check is Scanimage is running
    if evalin('base','exist(''hSLM'',''var'') && most.idioms.isValidObj(hSLM)')
        handle = evalin('base', 'hSLM');
    elseif evalin('base','exist(''hSI'',''var'') && most.idioms.isValidObj(hSI)')
        handle = evalin('base','hSI.hSlmScan.hSlm');
    else
        handle = dabs.meadowlark.MeadowlarkSlmOdp();
        assignin('base', 'hSLM', handle);
    end
        
    if evalin('base','exist(''hSLAPMi'',''var'')')
        hSLAPMi = evalin('base','hSLAPMi');
        calib = hSLAPMi.calib;
    else
        dataDir = 'E:\SLAPmidata';
        calib = SLAPMi.loadCalibrationData([dataDir filesep 'Calibration\calibration.cal']);
    end
    lut = calib.SLM.lut;
    lut=SLM_lut_recalib(lut,dlmread('C:\slm4411_at1064_P8.lut'));
    handle.overdriveEnable=true;
    handle.trueFrames=5;
    handle.lutFile=[];
    %handle.lutFile='C:\slm4411_at1064_P8.lut';

    %LED locations for 2 spot uncaging
    LED1 = [286, 214]; LED2 = [214,284];

    %Parse Pattern
    if numel(pattern)==1 %DRAW AN APERTURE
        radius = pattern;
        [XX,YY] = meshgrid(linspace(-1,1,512));
        select = sqrt(XX.^2 + YY.^2)<radius;
        pattern = zeros(size(select), 'uint8');
        pattern(select) = lut(cat(3, false(size(select)), select));
        pattern(~select) = lut(cat(3, ~select, false(size(select))));
    elseif all([size(pattern,1) size(pattern,2)] == [512 512])
        stack=size(pattern,3);
        pattern = repmat(lut(:,:,1),[1 1 stack]).*(1-pattern) + repmat(lut(:,:,2),[1 1 stack]).*pattern;
    elseif strcmpi(pattern, 'test')
        pattern = imread('C:\BLINK_PCIe\Image_Files\test.bmp');
    elseif strcmpi(pattern, 'spots')
        [XX,YY] = meshgrid(1:512);
        select= (mod(XX,10)<5 & mod(YY,10)<5);
        pattern = zeros(size(select), 'uint8');
        pattern(select) = lut(cat(3, false(size(select)), select));
        pattern(~select) = lut(cat(3, ~select, false(size(select))));
    elseif strcmpi(pattern, 'leds+')
        select = true(512);
        select(LED1(1)+(-5:5), LED1(2) +(-5:5)) = false;
        select(LED2(1)+(-5:5), LED2(2) +(-5:5)) = false;
        pattern = zeros(size(select), 'uint8');
        pattern(select) = lut(cat(3, false(size(select)), select));
        pattern(~select) = lut(cat(3, ~select, false(size(select))));
    elseif strcmpi(pattern, 'leds-')
        select = false(512);
        select(LED1(1)+(-5:5), LED1(2) +(-5:5)) = true;
        select(LED2(1)+(-5:5), LED2(2) +(-5:5)) = true;
        pattern = zeros(size(select), 'uint8');
        pattern(select) = lut(cat(3, false(size(select)), select));
        pattern(~select) = lut(cat(3, ~select, false(size(select))));
    else
        error('Pattern argument didn''t parse!') 
    end
    
    pattern = uint8(pattern); %just in case
    
    % write SLM queue
    handle.abortQueue;
    handle.writeQueue(pattern);
    handle.startQueue;
    
    % initialize overdrive
    
    t1 = tic;
    fprintf(1,'Writing a series of SLM patterns... ');

    % writeImage(obj,im,overdrive,wait_for_trigger,external_pulse)
    %handle.writeImage(pattern,overdrive,true,false)
    
    disp(['Done! Duration = ' num2str(toc(t1)) 's.'])
 
    %invert the LUT to get a pattern with 'expected' brightness
    lut = round(lut);
    patternOut = (double(mean(pattern,3)) - lut(:,:,1))./(lut(:,:,2)-lut(:,:,1));
    
   