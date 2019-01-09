[XX,YY] = meshgrid(linspace(-1,1, 512));
pattern = mod(sqrt(XX.^2 + YY.^2)*512, 256);

pattern(10:10:end, 10:10:end) = 128;

pattern(254:257, 254:401) = 20;
pattern(254:321, 254:257) = 20;
pattern(255:256, 255:400) = 100;
pattern(255:320, 255:256) = 100;

pattern = uint8(pattern);
figure, imshow(pattern,[])
imwrite(pattern, 'c:\blink_pcie\image_files\test.bmp')