function CM = cmap_KP(n)
if ~nargin
    n=255;
end
cmap = [0 0 0; 28 16 68 ; 79,18,123 ; 129,37,129; 181,54,122;  229,89,100; 251,135,97 ; 254,194,135 ; 251,253,191];
CM = interp1(linspace(0,1,size(cmap,1)),cmap./255,linspace(0,1,n));
end