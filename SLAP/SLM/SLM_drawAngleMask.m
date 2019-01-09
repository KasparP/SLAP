function SLM_drawAngleMask(freq, angles, aperture)
if nargin<3
    aperture = 0.55;
end
res = [512 512];
t = 0.9;
[xx,yy] = meshgrid(linspace(0,2*pi,res(1))*freq, linspace(0,2*pi,res(2))*freq);

mask = false(res);

for theta = (angles+[0 45 90 135])*pi/180
    x = xx*sin(theta);
    y = yy*cos(theta);
    mask = (mask | (sin(x+y)>t));
end

%figure, imshow(mask)
SLM_drawPattern(mask);
end
