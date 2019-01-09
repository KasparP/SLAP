function [C,numberOfOverlapPixels] = normxcorr2_FFT(T, A, FFT_T, FFT_A, requiredNumberOfOverlapPixels)
%Like NORMXCORR2_GENERAL, but takes the FFT of the inputs in and always does FFT-based registration
if nargin<5
    requiredNumberOfOverlapPixels = 1;
end

sizeA = size(A);
sizeT = size(T);

% Find the number of pixels used for the calculation as the two images are
% correlated.  The size of this image will be the same as the correlation
% image.  
numberOfOverlapPixels = local_sum(ones(sizeA),sizeT(1),sizeT(2));

local_sum_A = local_sum(A,sizeT(1),sizeT(2));
local_sum_A2 = local_sum(A.*A,sizeT(1),sizeT(2));

% Note: diff_local_sums should be nonnegative, but it may have negative
% values due to round off errors. Below, we use max to ensure the radicand
% is nonnegative.
diff_local_sums_A = ( local_sum_A2 - (local_sum_A.^2)./ numberOfOverlapPixels );
clear local_sum_A2;
denom_A = max(diff_local_sums_A,0); 
clear diff_local_sums_A;

% Flip T in both dimensions so that its correlation can be more easily
% handled.
rotatedT = rot90(T,2);
local_sum_T = local_sum(rotatedT,sizeA(1),sizeA(2));
local_sum_T2 = local_sum(rotatedT.*rotatedT,sizeA(1),sizeA(2));
clear rotatedT;

diff_local_sums_T = ( local_sum_T2 - (local_sum_T.^2)./ numberOfOverlapPixels );
clear local_sum_T2;
denom_T = max(diff_local_sums_T,0); 
clear diff_local_sums_T;

denom = sqrt(denom_T .* denom_A);
clear denom_T denom_A;

xcorr_TA = xcorr2_fast(T,A,FFT_T,FFT_A);
clear A T;
numerator = xcorr_TA - local_sum_A .* local_sum_T ./ numberOfOverlapPixels;
clear xcorr_TA local_sum_A local_sum_T;

% denom is the sqrt of the product of positive numbers so it must be
% positive or zero.  Therefore, the only danger in dividing the numerator
% by the denominator is when dividing by zero. We know denom_T~=0 from
% input parsing; so denom is only zero where denom_A is zero, and in these
% locations, C is also zero.
C = zeros(size(numerator));
tol = 1000*eps( max(abs(denom(:))) );
i_nonzero = find(denom > tol);
C(i_nonzero) = numerator(i_nonzero) ./ denom(i_nonzero);
clear numerator denom;

% Remove the border values since they result from calculations using very
% few pixels and are thus statistically unstable.
% By default, requiredNumberOfOverlapPixels = 0, so C is not modified.
if( requiredNumberOfOverlapPixels > max(numberOfOverlapPixels(:)) )
    error(['ERROR: requiredNumberOfOverlapPixels ' num2str(requiredNumberOfOverlapPixels) ...
    ' must not be greater than the maximum number of overlap pixels ' ...
    num2str(max(numberOfOverlapPixels(:))) '.']);
end
C(numberOfOverlapPixels < requiredNumberOfOverlapPixels) = 0;

end

%-------------------------------
% Function  local_sum
%
function local_sum_A = local_sum(A,m,n)

% This algorithm depends on precomputing running sums.

% If m,n are equal to the size of A, a faster method can be used for
% calculating the local sum.  Otherwise, the slower but more general method
% can be used.  The faster method is more than twice as fast and is also
% less memory intensive. 
if( m == size(A,1) && n == size(A,2) )
    s = cumsum(A,1);
    c = [s; repmat(s(end,:),m-1,1) - s(1:end-1,:)];
    s = cumsum(c,2);
    clear c;
    local_sum_A = [s, repmat(s(:,end),1,n-1) - s(:,1:end-1)];
else
    % Break the padding into parts to save on memory.
    B = zeros(size(A,1)+2*m,size(A,2));
    B(m+1:m+size(A,1),:) = A;
    s = cumsum(B,1);
    c = s(1+m:end-1,:)-s(1:end-m-1,:);
    d = zeros(size(c,1),size(c,2)+2*n);
    d(:,n+1:n+size(c,2)) = c;
    s = cumsum(d,2);
    local_sum_A = s(:,1+n:end-1)-s(:,1:end-n-1);
end

end

%-------------------------------
% Function  xcorr2_fast
%
function cross_corr = xcorr2_fast(FFT_T,FFT_A)

Fa = fft2(rot90(a,2),optimalSize(1),optimalSize(2));
Fb = fft2(b,optimalSize(1),optimalSize(2));
xcorr_ab = real(ifft2(Fa .* Fb));

xcorr_ab = xcorr_ab(1:outsize(1),1:outsize(2));


T_size = size(FFT_T);
A_size = size(FFT_A);
outsize = A_size + T_size - 1;
cross_corr = real(ifft2(FFT_T .* FFT_A));
cross_corr = cross_corr(1:outsize(1),1:outsize(2));
end
