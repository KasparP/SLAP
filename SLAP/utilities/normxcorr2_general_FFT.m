function [C,local_sum_A, local_sum_A2] = normxcorr2_general_FFT(T, A, fftArot, local_sum_A, local_sum_A2, requiredNumberOfOverlapPixels)
%this function has been stripped down for speed; ensure that neither T nor
%A are flat and that T and A are strictly positive

%NORMXCORR2_GENERAL Normalized two-dimensional cross-correlation.
%   [C,numberOfOverlapPixels] = NORMXCORR2_GENERAL(TEMPLATE,A) computes the
%   normalized cross-correlation of matrices TEMPLATE and A. The resulting
%   matrix C contains correlation coefficients and its values may range
%   from -1.0 to 1.0.
%
%   [C,numberOfOverlapPixels] =
%   NORMXCORR2_GENERAL(TEMPLATE,A,requiredNumberOfOverlapPixels) sets to 0
%   all locations in C computed from positions where A and T overlap less
%   than requiredNumberOfOverlapPixels.
%   Larger values of requiredNumberOfOverlapPixels zero-out pixels on a
%   larger border around C. 
%   Thus, larger values remove less stable computations but also limit the
%   capture range. 
%   If the template is smaller than the image and it is desired that the
%   computation only be carried out when the template is fully overlapping
%   the image, requiredNumberOfOverlapPixels should be set to
%   numel(template).
%   The default is set to 0, meaning no modifications to C.
%
%   Limitations of normxcorr2:
%   The documentation of normxcorr2 states that, "The matrix A must be
%   larger than the matrix TEMPLATE for the normalization to be
%   meaningful." It is implemented following the details of the paper "Fast
%   Normalized Cross-Correlation", by J. P. Lewis, Industrial Light &
%   Magic. This approach assumes the template is small relative to the
%   image and proceeds to calculate the normalization across the entire
%   template. This leads to correct computations wherever the template is
%   wholly overlapping with the image, but the computation is incorrect in
%   the borders of the output (the border size is proportional to the
%   template size).  This problem is therefore worse for larger templates
%   to the point that, when the template is the same size as the image, the
%   only correct value is at the center pixel (where the images are fully
%   overlapping). Thus, if normxcorr2 is used for such things as
%   registering images of the same size, the result will be incorrect.
%   
%   The new normxcorr2_general:
%   normxcorr2_general is more general than normxcorr2 in that it gives
%   correct results everywhere regardless of the relative size of A and
%   TEMPLATE. It accomplishes this by computing the normalized correlation
%   only in the overlap regions between the two matrices. Thus, the result
%   is correct for all locations of correlation.  The result is the same as
%   if the NCC were carried out in the spatial domain (which would take a
%   long time to compute for large matrices).
% 
%   Class Support
%   -------------
%   The input matrices can be numeric. The output matrix C is double.
%
%   Example
%   -------
%   This example correlates an input with itself using normxcorr2 (the
%   built-in Matlab version) and normxcorr2_general (the general version).
%   Because the template is not small compared with the input image (they
%   are the same size in this case), the output of normxcorr2.m is
%   incorrect for most pixels.  On the other hand, the general version is
%   correct at all locations, which can be easily verified analytically or
%   visually.
%   
%   Note that the image processing toolbox (IPT) is needed to run this
%   example since normxcorr2 is part of that toolbox.  However,
%   normxcorr2_general does not require the IPT.
% 
%   input = repmat([1:6 5:-1:1],11,1);
%   normxcorr2_output = normxcorr2(input,input);
%   normxcorr2_general_output = normxcorr2_general(input,input);
%   figure;
%   subplot(2,2,1), imagesc(input); title('Input pattern');
%   subplot(2,2,3), imagesc(normxcorr2_output); title('Output of Matlab built-in normxcorr2');
%   subplot(2,2,4), imagesc(normxcorr2_general_output); title('Output of normxcorr2\_general');
%
%   See also NORMXCORR2.
%
%   References: Dirk Padfield. "Masked FFT registration". In Proc. Computer
%   Vision and Pattern Recognition, 2010.
%
%   Author: Dirk Padfield, GE Global Research, padfield@research.ge.com
%

%   Input-output specs
%   ------------------
%   T:    2-D, real, full matrix
%         logical, uint8, uint16, or double
%         no NaNs, no Infs
%         prod(size(T)) >= 2
%
%   A:    2-D, real, full matrix
%         logical, uint8, uint16, or double
%         no NaNs, no Infs
%         prod(size(A)) >= 2
%
%   C:    double

sizeA = size(A);
sizeT = size(T);

% Find the number of pixels used for the calculation as the two images are
% correlated.  The size of this image will be the same as the correlation
% image.  
numberOfOverlapPixels = local_sum(ones(sizeA),sizeT(1),sizeT(2));

if isempty(local_sum_A)
    local_sum_A = local_sum(A,sizeT(1),sizeT(2));
    local_sum_A2 = local_sum(A.*A,sizeT(1),sizeT(2));
end
% Note: diff_local_sums should be nonnegative, but it may have negative
% values due to round off errors. Below, we use max to ensure the radicand
% is nonnegative.
diff_local_sums_A = ( local_sum_A2 - (local_sum_A.^2)./ numberOfOverlapPixels );

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

outsize = sizeA + sizeT - 1;
xcorr_TA = freqxcorr(T,fftArot,outsize);

clear A T;
numerator = xcorr_TA - local_sum_A .* local_sum_T ./ numberOfOverlapPixels;
clear xcorr_TA local_sum_T;

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


%-------------------------------
% Function  freqxcorr
%
function xcorr_ab = freqxcorr(T,fftArot,outsize)
% Calculate correlation in frequency domain
Ft = fft2(T,size(fftArot,1),size(fftArot,2));
xcorr_ab = real(ifft2(fftArot .* Ft));
xcorr_ab = rot90(xcorr_ab(1:outsize(1),1:outsize(2)),2);
