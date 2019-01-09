

function ddisplay = SevenSegment(digit,mode,matrixHeight,poligonWidth,distans,matrixBand,stat)
%Seven segment image builder
%Format:
%
%SevenSegmentImage=SevenSegment(digit)
%SevenSegmentImage=SevenSegment(digit,mode,matrixHeight,poligonWidth,distans,matrixBand,stat)
%
%Input:
%   %all dimensions in pixels
%   digit:        the digit to display
%   mode:         'free'- the output image length and decimal point, according to the input digit.
%                 'x.y' - the output image length fixed to 'x' digits include the sign before
%                 decimal point and 'y' digits after it; segments that not in use, will fill by zero.
%   matrixHeight: is a height of output image, length is calculated
%                 from other input parameters.
%   poligonWidth: is a width of single seven segment poligon
%   distans:      is a distans between seven segment digits
%   matrixBand:   is a white area between edge of image and digits
%   stat          1 for display dimantions, 0 for not
%   
%
%Output:
%                 the image with seven segment representation of input digit
%                 if stat=1: list of output dimantions in Command Window
%
%Usage example:
%               SevenSegmentImage=sevenSegment(-20.56,'4.3',100,10,20,10,0)
%               imshow(SevenSegmentImage)
%               SevenSegmentImage=sevenSegment(-20.56,'free',100,10,20,10,0)
%               imshow(SevenSegmentImage)
%               
%
% Autor: Alex Zeleznyak       Version 1.0
% Copyright (c) 2010 Alex Zeleznyak. 
% Please email me (Alex74zel@gmail.com) if you find bugs, or have 
% suggestions or questions!
%
%If you use this code, you can donate PayPal Alex74zel@gmail.com
%

%default input parameters
if  nargin == 1
    matrixHeight=100;
    poligonWidth=10;
    distans=20;
    matrixBand=10;
    stat=0;
    mode='free';
end
if not(ischar(mode))
    error('ErrorTests:legalityTest','mode must be a string')
end
if length(digit)>1 || isempty(digit)
    error('ErrorTests:legalityTest','digit must be a real number' )
end
if matrixHeight<1 || poligonWidth<1 || distans<1 || matrixBand<0 || stat>1 || stat<0  
    error('ErrorTests:legalityTest','Illegal input parameters' )
end
%chack input parameters proportion
if (matrixHeight-2*matrixBand)< poligonWidth*3+1
    error('ErrorTests:resolutionTest','Input parameters under posible resolution \n matrixHeight-2*matrixBand must be big that poligonWidth*3+1')
end

numstring=num2str(digit);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(mode,'free')
    %fixed format calculations
else
    %calculation of number of symbols (include sign) that different from 
    %zero befor decimal point
    digit_before=length(num2str(fix(digit)));
    if digit<0 && fix(digit)==0
        digit_before=digit_before+1;
    end
    %calculation of number of symbols that different from 
    %zero after decimal point
    digit_after=length(num2str(digit))-digit_before-1;
    %convert mode to number
    nummode=str2double(mode);
    %calculation number of symbols before decimal point according to fixed
    %mode restriction
    before=int8(fix(nummode));
    %calculation number of symbols after decimal point according to fixed
    %mode restriction
    after=int8((nummode-double(before))*10);
        
    %Add zeros before decimal point to fill format restriction
    if before>digit_before
        dif=before-digit_before;
        for i=1:dif
            numstring=strcat('0',numstring);
        end
    elseif before<digit_before
        error('ErrorTests:legalityTest','lenth of input number is out of format restrictions' )                             
    end 
    %Add zeros after decimal point to fill format restriction
    if after>digit_after 
        if digit_after == -1 %if the input digit is integer
            numstring=strcat(numstring,'.');
            digit_after = 0;
        end
        dif=after-digit_after;
        for i=1:dif
            numstring=strcat(numstring,'0');
        end
    %in case conflict between symbols afer decimal point of input digit and
    %format restriction, digits after decimal point displayed according to
    %input number
    elseif after<digit_after
        warning('digits after decimal point displayed according to input number' )                             
    end 
    %sign forfvarding
    if digit<0
        sign=strfind(numstring,'-');
        numstring=strcat(numstring(sign),numstring(1:sign-1),numstring(sign+1:end));  
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%flag of point exist
afterPoint=0;


digLength=length(numstring);
%chack for point in input digit
point=numstring=='.';
if max(point)>0
    %pointIndex=find(point==1);
    ispoint=1;
else
    ispoint=0;
end
%calculation of seven segment digit width
digW=ceil((matrixHeight-2*matrixBand+poligonWidth)/2);
%calculation of seven segment poligon length
poliL=digW;  
%calculation of seven segment poligon width
poliW=poligonWidth;

%calculation of seven segment digit height
digH=matrixHeight-2*matrixBand;
%calculation of output image band
band=matrixBand;
%memory allocation for output image
ddisplay=ones(matrixHeight,digW*(digLength-ispoint)+band*2+distans*(digLength-1)+poliW*ispoint)*255;
%horizontal poligon definition
 horpol=ones(poliW,poliL)*255; 
%verticalp poligon definition
 verpol=ones(poliL,poliW)*255;
%point poligon definition
 pointpol=ones(poliW,poliW)*255;
%loop for avery symbol in the digit
for i=1:digLength
    %a,b,c,d,e,f,g;p seven segment elements
    a=horpol;%    -a-
    b=verpol;%  f|   |b
    c=verpol;%    -g-
    d=horpol;%  e|   |c
    e=verpol;%    -d-    p
    f=verpol;
    g=horpol;
    p=pointpol;
    %color segment difeniton
    switch numstring(i)
        case '0'
            a=0;b=0;c=0;d=0;e=0;f=0;
        case '1'
            b=0;c=0;
        case '2'
            a=0;b=0;d=0;e=0;g=0;
        case '3'
            a=0;b=0;d=0;c=0;g=0;
        case '4'
            b=0;c=0;g=0;f=0;
        case '5'
            a=0;c=0;d=0;g=0;f=0;
        case '6'
            a=0;c=0;d=0;e=0;f=0;g=0;
        case '7'
            a=0;b=0;c=0;
        case '8'
            a=0;b=0;c=0;d=0;e=0;f=0;g=0;
        case '9'
            a=0;b=0;c=0;d=0;f=0;g=0;
        case '-'
            g=0;
        case '.'
            p=0;
            afterPoint=1;
    end
    %calculation position of carrent symbol of input digit in reference to output image X coordinates (start from band)
    digitPlace=(i-1-afterPoint)*(poliL+distans)+afterPoint*(poliW+distans)+band;
    if a==0
        %display segment a as black
        ddisplay(1+band:band+poliW,digitPlace+1:poliL+digitPlace)=a;
    end
    if g==0
        %display segment g as black
        ddisplay(band+poliL-poliW+1:band+poliL,digitPlace+1:poliL+digitPlace)=g;
    end
    if d==0
        %display segment d as black
        ddisplay(band+2*poliL-2*poliW+1:band+2*poliL-poliW,digitPlace+1:poliL+digitPlace)=d;
    end
    if f==0
        %display segment f as black
        ddisplay(band+1:band+poliL,digitPlace+1:poliW+digitPlace)=f;
    end
    if e==0
        %display segment e as black
        ddisplay(band+poliL-poliW+1:band+2*poliL-poliW,digitPlace+1:poliW+digitPlace)=e;
    end
    if b==0
        %display segment b as black
        ddisplay(band+1:band+poliL,poliL-poliW+1+digitPlace:poliL+digitPlace)=b;
    end
    if c==0
        %display segment c as black
        ddisplay(band+poliL-poliW+1:band+2*poliL-poliW,poliL-poliW+1+digitPlace:poliL+digitPlace)=c;
    end
    if p==0
        %display point
        ddisplay(band+2*poliL-2*poliW+1:band+2*poliL-poliW,digitPlace+1-afterPoint*(poliW+distans)+(poliL+distans):poliW+digitPlace-afterPoint*(poliW+distans)+(poliL+distans))=p;
    end
    
end
%dimantions of output 
if stat==1
    disp(['Digits width       = ',num2str(digW)]);
    disp(['Digit height       = ',num2str(digH)]);
    disp(['Poligon width      = ',num2str(poligonWidth)]);
    disp(['Poligon length     = ',num2str(digW)]);
    disp(['Distans            = ',num2str(distans)]);
    disp(['Band               = ',num2str(matrixBand)]);
    disp('Output image size : Rows   Columns');
    disp(['                    ',num2str(size(ddisplay,1)),'      ',num2str(size(ddisplay,2))]);
    
end
end

