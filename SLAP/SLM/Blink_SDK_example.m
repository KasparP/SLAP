% Example usage of Blink_SDK_C.dll
% Meadowlark Optics Spatial Light Modulators
% March 28 2015

% Load the DLL
% Blink_SDK_C.dll, Blink_SDK.dll, FreeImage.dll and wdapi1021.dll
% should all be located in the same directory as the program referencing the
% library
% Matlab only supports C-style headers so this is a 'sanitized' version of
% the normal header file
cd('C:\Users\scanimage\Desktop\tmp\Blink SDK');
loadlibrary('Blink_SDK_C.dll',str2func('Blink_SDK_C'));


% Basic parameters for calling Create_SDK
bit_depth = 8;
slm_resolution = 512;
num_boards_found = libpointer('uint32Ptr', 0);
constructed_okay = libpointer('int32Ptr', 0);
is_nematic_type = 1;
RAM_write_enable = 1;
use_GPU = 1;
max_transients = 20;

% OverDrive Plus Parameters
% Matlab automatically escapes backslashes (unlike most languages)
lut_file = 'C:\Program Files\Meadowlark Optics\Blink OverDrive Plus\LUT Files\slm4411_at1064_regional_encrypt.txt';

% Basic SLM parameters
true_frames = 3;

% Blank calibration image
cal_image = imread('C:\Program Files\Meadowlark Optics\Blink OverDrive Plus\SDK\512white.bmp');

% Arrays for image data
ramp_0 = imread('C:\Program Files\Meadowlark Optics\Blink OverDrive Plus\SDK\ramp_0_512.bmp');
ramp_1 = imread('C:\Program Files\Meadowlark Optics\Blink OverDrive Plus\SDK\ramp_1_512.bmp');
on=imread('C:\Program Files\Meadowlark Optics\Blink OverDrive Plus\Other Image Files\0.bmp');
on=on./255; % normalize to 0 - 1.0
off=imread('C:\Program Files\Meadowlark Optics\Blink OverDrive Plus\Other Image Files\1.bmp');

sdk = calllib('Blink_SDK_C', 'Create_SDK', bit_depth, slm_resolution, num_boards_found, constructed_okay, is_nematic_type, RAM_write_enable, use_GPU, max_transients, lut_file);

if constructed_okay.value == 0
    disp('Blink SDK was not successfully constructed');
    disp(calllib('Blink_SDK_C', 'Get_last_error_message', sdk));
    calllib('Blink_SDK_C', 'Delete_SDK', sdk);
else
    disp('Blink SDK was successfully constructed');
    fprintf('Found %u SLM controller(s)', num_boards_found.value);
    % Set the basic SLM parameters
    calllib('Blink_SDK_C', 'Set_true_frames', sdk, true_frames);
    % A blank calibration file must be loaded to the SLM controller
    calllib('Blink_SDK_C', 'Write_cal_buffer', sdk, 1, cal_image);
    % A linear LUT must be loaded to the controller for OverDrive Plus
    %calllib('Blink_SDK_C', 'Load_linear_LUT', sdk, 1);
    calllib('Blink_SDK_C', 'Load_LUT_file', sdk, 1,'C:\slm4411_at1064_P8.lut');
    
    % Turn the SLM power on
    calllib('Blink_SDK_C', 'SLM_power', sdk, 1);
    
    % Loop between our ramp images
    success=calllib('Blink_SDK_C', 'Write_overdrive_image', sdk, 1, on, 1, 0);
    %sucess=calllib('Blink_SDK_C', 'Write_image', sdk, 1, cal_image,size(cal_image,1), 1, 0);
    
    % Always call Delete_SDK before exiting
    %calllib('Blink_SDK_C', 'Delete_SDK', sdk);
end

%unloadlibrary('Blink_SDK_C');