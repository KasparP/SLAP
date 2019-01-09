function B = ordfilt2(varargin)
%ORDFILT2 2-D order-statistic filtering.
%   B=ORDFILT2(A,ORDER,DOMAIN) replaces each element in A by the
%   ORDER-th element in the sorted set of neighbors specified by
%   the nonzero elements in DOMAIN.  
%
%   B = ORDFILT2(A,ORDER,DOMAIN,S), where S is the same size as
%   DOMAIN, uses the values of S corresponding to the nonzero
%   values of DOMAIN as additive offsets.
%
%   B = ORDFILT2(...,PADOPT) controls how the matrix boundaries
%   are padded.  PADOPT may be 'zeros' (the default) or
%   'symmetric'.  If PADOPT is 'zeros', A is padded with zeros at
%   the boundaries.  If PADOPT is 'symmetric', A is symmetrically
%   extended at the boundaries.
%
%   Class Support
%   -------------
%   The class of A may be numeric or logical.  The class of B is 
%   the same as the class of A, unless the additive offset form of 
%   ORDFILT2 is used, in which case the class of B is double.
%
%   Example
%   -------
%   Use a maximum filter on snowflakes.png with a [5 5] neighborhood.  This is
%   equivalent to imdilate(image,strel('square',5)).
%  
%       A = imread('snowflakes.png');
%       B = ordfilt2(A,25,true(5));
%       figure, imshow(A), figure, imshow(B)
%  
%   Remarks
%   -------
%   DOMAIN is equivalent to the structuring element used for
%   binary image operations. It is a matrix containing only 1's
%   and 0's; the 1's define the neighborhood for the filtering
%   operation.
%
%   For example, B=ORDFILT2(A,5,ONES(3,3)) implements a 3-by-3
%   median filter; B=ORDFILT2(A,1,ONES(3,3)) implements a 3-by-3
%   minimum filter; and B=ORDFILT2(A,9,ONES(3,3)) implements a
%   3-by-3 maximum filter.  B=ORDFILT2(A,4,[0 1 0; 1 0 1; 0 1 0])
%   replaces each element in A by the maximum of its north, east,
%   south, and west neighbors. 
%
%   See also MEDFILT2.

%   Copyright 1993-2014 The MathWorks, Inc.

[A,order,domain,s,padopt,~] = ParseInputs(varargin{:});

if(isempty(A))
    B = A;
    return;
end

domainSize = size(domain);
center = floor((domainSize + 1) / 2);
[r,c] = find(domain);
r = r - center(1);
c = c - center(2);
padSize = max(max(abs(r)), max(abs(c)));
originalSize = size(A);
if (strcmp(padopt, 'zeros'))
    A = padarray(A, padSize * [1 1], 0, 'both');
elseif (strcmp(padopt, 'ones'))
    % padopt of 'ones' is for support of medfilt2; it is
    % undocumented
    A = padarray(A, padSize * [1 1], 1, 'both');
elseif (strcmp(padopt, 'replicate'))
    A = padarray(A, padSize * [1 1], 'replicate', 'both');
else
    A = padarray(A, padSize * [1 1], 'symmetric', 'both');
end
Ma = size(A,1);
offsets = c*Ma + r;

% make sure that offsets are valid
if (~isreal(offsets) || any(floor(offsets) ~= offsets) || any(~isfinite(offsets)))
    %should never get here
    error(message('images:ordfilt2:badOffsets'))
end

if isempty(s)
  %ORDFILT2(A,ORDER,DOMAIN)
  B = ordfilt2mex(A, order, offsets, [padSize padSize] + 1, ...
             originalSize, domainSize);
else
  %ORDFILT2(A,ORDER,DOMAIN,S,PADOPT)
  B = ordfilt2mex(A, order, offsets, [padSize padSize] + 1, ...
           originalSize, domainSize, s);
end


%%%
%%% ParseInputs
%%%
function [A,order,domain,s,padopt,msg] = ParseInputs(varargin)

s = [];
padopt = 'zeros';
msg = '';

narginchk(3,5);

A = varargin{1};
order = varargin{2};
domain = varargin{3};
options = {'zeros', 'ones', 'symmetric', 'replicate'};
% padopt of 'ones' is for supporting medfilt2; it is undocumented.

if (nargin == 4)
  if (ischar(varargin{4}))
    padopt = validatestring(varargin{4},options,mfilename,'PADOPT',4);
  else
    s = varargin{4};
  end
    
elseif (nargin == 5)
  s = varargin{4};
  padopt = validatestring(varargin{5},options,mfilename,'PADOPT',5);  
end


% make sure that arguments are valid
validateattributes(order,{'double'},{'real','scalar','integer'},mfilename, ...
              'ORDER',2);

% check domain
validateattributes(domain,{'numeric','logical'}, {'2d'}, mfilename, ...
              'DOMAIN',3);
domain = logical(domain);
          
% check order
if ((order<1) || order>sum(domain(:)))
    error(message('images:ordfilt2:orderNotValid'));
end


if ~isempty(s)
  validateattributes(A, ...
      {'uint8','uint16','uint32','int8','int16','int32','single','double','logical'},...
      {'2d','real','nonsparse'}, mfilename, 'A', 1);  
  A = double(A);
  % check offsets
  if(~isequal(size(s), size(domain)))
      error(message('images:ordfilt2:sizeMismatch'));
  end
  s = s(domain ~= 0);
  validateattributes(s, {'double'}, {'real'}, mfilename, 'S', 4);
else
  validateattributes(A, ...
      {'uint8','uint16','uint32','int8','int16','int32','single','double','logical'},...
      {'2d','real','nonsparse'}, mfilename, 'A', 1);
end




