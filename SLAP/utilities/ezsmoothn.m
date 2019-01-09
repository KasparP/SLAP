function z = ezsmoothn(y)

%EZSMOOTHN Robust smoothing for 1-D to N-D data (easy version of SMOOTHN)
%   EZSMOOTHN provides a fast, automatized and robust discretized spline
%   smoothing for real or complex data of arbitrary dimension. EZSMOOTHN
%   can deal with missing data.
%
%   Z = EZSMOOTHN(Y) automatically smoothes the uniformly-sampled array Y.
%   Y can be any real or complex N-D array. Non finite data (NaN or Inf)
%   are treated as missing values.
%
%   EZSMOOTHN is an easy and simplified version of the SMOOTHN function. It
%   contains around 40 lines of code to make the understanding of SMOOTHN
%   much easier.
%
%   References (please refer to the two following papers)
%   ---------- 
%   1) Garcia D, Robust smoothing of gridded data in one and higher
%   dimensions with missing values. Computational Statistics & Data
%   Analysis, 2010;54:1167-1178. 
%   <a
%   href="matlab:web('http://www.biomecardio.com/pageshtm/publi/csda10.pdf')">PDF download</a>
%   2) Garcia D, A fast all-in-one method for automated post-processing of
%   PIV data. Exp Fluids, 2011;50:1247-1259.
%   <a
%   href="matlab:web('http://www.biomecardio.com/pageshtm/publi/media11.pdf')">PDF download</a>
%
%   Examples:
%   --------
%   %--- Example #1: smooth a curve ---
%   x = linspace(0,100,2^8);
%   y = cos(x/10)+(x/50).^2 + randn(size(x))/4;
%   y([70 75 80]) = [5.5 5 6];
%   z = ezsmoothn(y); % Regular smoothing
%   plot(x,y,'r'), hold on, plot(x,z,'k','LineWidth',2)
%
%   %--- Example #2: smooth a surface ---
%   xp = 0:.02:1;
%   [x,y] = meshgrid(xp);
%   f = exp(x+y) + sin((x-2*y)*3);
%   fn = f + randn(size(f))*0.5;
%   fs = ezsmoothn(fn);
%   subplot(121), surf(xp,xp,fn), zlim([0 8]), axis square
%   subplot(122), surf(xp,xp,fs), zlim([0 8]), axis square
%
%   %--- Example #3: smooth an image with missing data ---
%   n = 256;
%   y0 = peaks(n);
%   y = y0 + randn(size(y0))*2;
%   I = randperm(n^2);
%   y(I(1:n^2*0.5)) = NaN; % lose 1/2 of data
%   y(40:90,140:190) = NaN; % create a hole
%   z = ezsmoothn(y); % smooth data
%   subplot(2,2,1:2), imagesc(y), axis equal off
%   title('Noisy corrupt data')
%   subplot(223), imagesc(z), axis equal off
%   title('Recovered data ...')
%   subplot(224), imagesc(y0), axis equal off
%   title('... compared with the actual data')
%
%   %--- Example #4: smooth volumetric data ---
%   [x,y,z] = meshgrid(-2:.2:2);
%   xslice = [-0.8,1]; yslice = 2; zslice = [-2,0];
%   vn = x.*exp(-x.^2-y.^2-z.^2) + randn(size(x))*0.06;
%   subplot(121), slice(x,y,z,vn,xslice,yslice,zslice,'cubic')
%   title('Noisy data')
%   v = ezsmoothn(vn);
%   subplot(122), slice(x,y,z,v,xslice,yslice,zslice,'cubic')
%   title('Smoothed data')
%
%   %--- Example #5: smooth a cardioid ---
%   t = linspace(0,2*pi,1000);
%   x = 2*cos(t).*(1-cos(t)) + randn(size(t))*0.1;
%   y = 2*sin(t).*(1-cos(t)) + randn(size(t))*0.1;
%   z = ezsmoothn(complex(x,y));
%   plot(x,y,'r.',real(z),imag(z),'k','linewidth',2)
%   axis equal tight
%
%   %--- Example #6: smooth a 2-D vector field ---
%   [x,y] = meshgrid(linspace(0,1,24));
%   Vx = cos(2*pi*x+pi/2).*cos(2*pi*y);
%   Vy = sin(2*pi*x+pi/2).*sin(2*pi*y);
%   Vx = Vx + sqrt(0.05)*randn(24,24); % adding Gaussian noise
%   Vy = Vy + sqrt(0.05)*randn(24,24); % adding Gaussian noise
%   I = randperm(numel(Vx));
%   Vx(I(1:30)) = (rand(30,1)-0.5)*5; % adding outliers
%   Vy(I(1:30)) = (rand(30,1)-0.5)*5; % adding outliers
%   Vx(I(31:60)) = NaN; % missing values
%   Vy(I(31:60)) = NaN; % missing values
%   Vs = ezsmoothn(complex(Vx,Vy)); % automatic smoothing
%   subplot(121), quiver(x,y,Vx,Vy,2.5), axis square
%   title('Noisy velocity field')
%   subplot(122), quiver(x,y,real(Vs),imag(Vs)), axis square
%   title('Smoothed velocity field')
%
%   See also SMOOTHN, SMOOTH1Q.
%
%   -- Damien Garcia -- 2014/02, revised 2014/02/26
%   website: <a
%   href="matlab:web('http://www.biomecardio.com')">www.BiomeCardio.com</a>

sizy = size(y);
n = prod(sizy); % total number of elements in y
N = sum(sizy~=1); % rank tensor

Lambda = zeros(sizy);
d = ndims(y);
for i = 1:d
    siz0 = ones(1,d);
    siz0(i) = sizy(i);
    Lambda = bsxfun(@plus,Lambda,...
        2-2*cos(pi*(reshape(1:sizy(i),siz0)-1)/sizy(i)));
end

%-- Weights
W0 = ones(sizy);
I = isfinite(y); % missing data (NaN or Inf values)
W0(~I) = 0; % weights for missing data are 0
y(~I) = mean(y(I));
W = W0;

zz = y;
for k = 1:4
    tol = Inf;
    while tol>1e-4
        DCTy = dctn(W.*(y-zz)+zz);
        fminbnd(@GCVscore,-10,30,optimset('TolX',.1));
        tol = norm(zz(:)-z(:))/norm(z(:));
        zz = z;
    end
    h = (sqrt(1+sqrt(1+16*s))/sqrt(2)/sqrt(1+16*s))^N;
    W = bisquare(y,z,I,h).*W0; % weights
end
    function GCVs = GCVscore(p) % GCV score
        s = 10^p;
        Gamma = 1./(1+s*Lambda.^2);
        z = idctn(Gamma.*DCTy);
        RSS = norm(sqrt(W(:)).*(y(:)-z(:)))^2;
        TrH = sum(Gamma(:));
        GCVs = RSS/n/(1-TrH/n)^2;
    end
end

function W = bisquare(y,z,I,h)
r = y-z; % residuals
MAD = median(abs(r(I)-median(r(I)))); % median absolute deviation
u = abs(r/(1.4826*MAD)/sqrt(1-h)); % studentized residuals
W = (1-(u/4.685).^2).^2.*((u/4.685)<1); % bisquare weights
end



%% DCTN
function y = dctn(y)

%DCTN N-D discrete cosine transform.
%   Y = DCTN(X) returns the discrete cosine transform of X. The array Y is
%   the same size as X and contains the discrete cosine transform
%   coefficients. This transform can be inverted using IDCTN.
%
%   Reference
%   ---------
%   Narasimha M. et al, On the computation of the discrete cosine
%   transform, IEEE Trans Comm, 26, 6, 1978, pp 934-936.
%
%   Example
%   -------
%       RGB = imread('autumn.tif');
%       I = rgb2gray(RGB);
%       J = dctn(I);
%       imshow(log(abs(J)),[]), colormap(jet), colorbar
%
%   The commands below set values less than magnitude 10 in the DCT matrix
%   to zero, then reconstruct the image using the inverse DCT.
%
%       J(abs(J)<10) = 0;
%       K = idctn(J);
%       figure, imshow(I)
%       figure, imshow(K,[0 255])
%
%   -- Damien Garcia -- 2008/06, revised 2011/11
%   -- www.BiomeCardio.com --

y = double(y);
sizy = size(y);
y = squeeze(y);
dimy = ndims(y);

% Some modifications are required if Y is a vector
if isvector(y)
    dimy = 1;
    if size(y,1)==1, y = y.'; end
end

% Weighting vectors
w = cell(1,dimy);
for dim = 1:dimy
    n = (dimy==1)*numel(y) + (dimy>1)*sizy(dim);
    w{dim} = exp(1i*(0:n-1)'*pi/2/n);
end

% --- DCT algorithm ---
if ~isreal(y)
    y = complex(dctn(real(y)),dctn(imag(y)));
else
    for dim = 1:dimy
        siz = size(y);
        n = siz(1);
        y = y([1:2:n 2*floor(n/2):-2:2],:);
        y = reshape(y,n,[]);
        y = y*sqrt(2*n);
        y = ifft(y,[],1);
        y = bsxfun(@times,y,w{dim});
        y = real(y);
        y(1,:) = y(1,:)/sqrt(2);
        y = reshape(y,siz);
        y = shiftdim(y,1);
    end
end
        
y = reshape(y,sizy);

end

%% IDCTN
function y = idctn(y)

%IDCTN N-D inverse discrete cosine transform.
%   X = IDCTN(Y) inverts the N-D DCT transform, returning the original
%   array if Y was obtained using Y = DCTN(X).
%
%   Reference
%   ---------
%   Narasimha M. et al, On the computation of the discrete cosine
%   transform, IEEE Trans Comm, 26, 6, 1978, pp 934-936.
%
%   Example
%   -------
%       RGB = imread('autumn.tif');
%       I = rgb2gray(RGB);
%       J = dctn(I);
%       imshow(log(abs(J)),[]), colormap(jet), colorbar
%
%   The commands below set values less than magnitude 10 in the DCT matrix
%   to zero, then reconstruct the image using the inverse DCT.
%
%       J(abs(J)<10) = 0;
%       K = idctn(J);
%       figure, imshow(I)
%       figure, imshow(K,[0 255])
%
%   See also DCTN, IDSTN, IDCT, IDCT2, IDCT3.
%
%   -- Damien Garcia -- 2009/04, revised 2011/11
%   -- www.BiomeCardio.com --

y = double(y);
sizy = size(y);
y = squeeze(y);
dimy = ndims(y);

% Some modifications are required if Y is a vector
if isvector(y)
    dimy = 1;
    if size(y,1)==1
        y = y.';
    end
end

% Weighing vectors
w = cell(1,dimy);
for dim = 1:dimy
    n = (dimy==1)*numel(y) + (dimy>1)*sizy(dim);
    w{dim} = exp(1i*(0:n-1)'*pi/2/n);
end

% --- IDCT algorithm ---
if ~isreal(y)
    y = complex(idctn(real(y)),idctn(imag(y)));
else
    for dim = 1:dimy
        siz = size(y);
        n = siz(1);
        y = reshape(y,n,[]);
        y = bsxfun(@times,y,w{dim});
        y(1,:) = y(1,:)/sqrt(2);
        y = ifft(y,[],1);
        y = real(y*sqrt(2*n));
        I = (1:n)*0.5+0.5;
        I(2:2:end) = n-I(1:2:end-1)+1;
        y = y(I,:);
        y = reshape(y,siz);
        y = shiftdim(y,1);            
    end
end
        
y = reshape(y,sizy);

end