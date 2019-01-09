clearvars;

[~,~,Z] = peaks(100);
Z = single(Z);
Z_index = uint8((Z - min(Z(:))) * (255 / (max(Z(:)) - min(Z(:)))));
Z_color = uint8(ind2rgb(Z_index, hsv(256)*256));
Z_alpha = Z_color; Z_alpha(:,:,4) = repmat(0:1/100:1-1/100, 100, 1)*255;
Z_color_multiframe = reshape([Z_color(:)*0.2 Z_color(:)*0.6 Z_color(:)], 100, 100, 3, 3);

% 8-bit, grayscale image
clear options;
saveastiff(uint8(Z_index), 'Z_uint8.tif');

% Lossless LZW compression
clear options;
options.comp = 'adobe';
saveastiff(uint8(Z_index), 'Z_uint8_compress.tif', options);

% Overwrite to an existing file
clear options;
options.overwrite = true;
options.compress = 'adobe';
saveastiff(uint8(Z_index), 'Z_uint8_compress.tif', options);

% Allow message printing.
clear options;
options.message = true;
saveastiff(uint8(Z_index), 'Z_uint8_Message.tif', options);

% 16 bit, grayscale image
clear options;
saveastiff(uint16(single(Z_index)/255*(2^16-1)), 'Z_uint16.tif');

% 32 bit, grayscale image
clear options;
saveastiff(uint32(single(Z_index)/255*(2^32-1)), 'Z_uint32.tif');

% 32 bit single, grayscale image
clear options;
saveastiff(Z, 'Z_single.tif');

% Color image
clear options;
options.color = true;
saveastiff(Z_color, 'Z_rgb.tif', options);

% Color image with alpha channel
clear options;
options.color = true;
saveastiff(Z_alpha, 'Z_alpha.tif', options);

% Save each R, G and B chanels of the color image, separately.
clear options;
saveastiff(Z_color, 'Z_rgb_channel.tif');

% Save the multi-frame RGB color image
clear options;
options.color = true;
saveastiff(Z_color_multiframe, 'Z_rgb_multiframe.tif', options);

% 32 bit single, 50x50x50 volume data
clear options;
saveastiff(single(rand(50, 50, 50)), 'volume_50x50x50.tif');

% Append option is ignored if path dose not exist.
clear options;
options.append = true;
saveastiff(Z_index, 'Z_uint8_append.tif', options);

% You can append any type of image to an existing tiff file.
clear options;
options.append = true;
saveastiff(single(rand(10, 10, 3)), 'Z_uint8_append.tif', options);
options.color = true;
saveastiff(Z_color_multiframe, 'Z_uint8_append.tif', options);

% Save image to a sub directory
clear options;
saveastiff(uint8(Z_index), 'sub1/sub2/Z_uint8.tif');
options.append = true;
saveastiff(uint8(Z_index), 'sub1/sub2/Z_uint8.tif', options);

% Save complex number images
saveastiff(fft2(Z_index), 'Z_fft2_gray.tif');
saveastiff(fft2(Z_color), 'Z_fft2_color.tif', options);
clear options;
options.color = true;
options.append = true;
saveastiff(fft2(Z_index), 'Z_fft2_append.tif');
saveastiff(fft2(Z_color), 'Z_fft2_append.tif', options);
saveastiff(fft2(Z_color), 'Z_fft2_append.tif', options);
saveastiff(fft2(Z_alpha), 'Z_fft2_append.tif', options);

% Load multiframe tiff
multiframe = loadtiff('volume_50x50x50.tif');
complexframe = loadtiff('Z_uint8_append.tif');
subdir = loadtiff('sub1/sub2/Z_uint8.tif');
fftgray = loadtiff('Z_fft2_gray.tif');
fftcolor = loadtiff('Z_fft2_color.tif');
fftappend = loadtiff('Z_fft2_append.tif');

% Big Tiff File (64 bit)
[~,maxArraySize]=computer; 
is64bitMatlab = maxArraySize> 2^31;
if is64bitMatlab
    % Use 'big' option to save larget than 4GB files .
    clear options;
    options.big = true; % Use BigTIFF format
    saveastiff(single(zeros(24000, 24000)), 'BigTiff(2GB+2GB).btf', options);
    options.append = true;
    saveastiff(single(zeros(24000, 24000)), 'BigTiff(2GB+2GB).btf', options); % 4GB Big TIFF file
    bigtiff = loadtiff('BigTiff(2GB+2GB).btf');
end
