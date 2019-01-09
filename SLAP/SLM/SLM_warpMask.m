function SLM_warpMask
%draws a mask used for phase alignment and dewarping in scanimage

res = [512 512];

im = rand(res);
im = imgaussfilt(im, 5);
im = mod(round(800*im), 2);

im = im>0.5;
%im = medfilt2(im, [3 3])>0.5;
%im = medfilt2(im, [3 3])>0.5;

%im = 256*double(im);
SLM_drawPattern(im);


    
    