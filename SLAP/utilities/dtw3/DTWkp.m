function [Ta, D] = DTWkp (A,B,maxlag, doplot)
%B is a fixed signal we are registering to
%A is the signal to warp

%Ta is the warping vector; i.e. call interp1(A,Ta) to get a signal that matches B. 
%D is the overall similarity of A to B (min cost), bounded below at 0

L = length(A); 
if L~=length(B)
    error('bad inputs')
end

%stepcost = [1 1 1.0; 0 1 1.2 ; 1 0 1.2];  
stepcost = [1 1 1; 1 2 2; 2 1 2; 1 3 5; 3 1 5];

repeat = true;
lambda = max(max(A),max(B))/(10*maxlag);
while repeat
    N = lambda*(1:L)';
    
    %distance matrix
    %should be log probability of the moving signal being a sample from the fixed signal
    
    %scale B to the magnitude of A
    %B = B * mean(A)/mean(B);
    %S = 1-log(bsxfun(@poisspdf, A', B))
    S = pdist2([A(:) N ],[B(:) N]).^2;
    
    if nargin>2 && ~isempty(maxlag)
        for i = 1:length(S)
            S(i, [1:i-maxlag-1 i+maxlag+1:end]) = inf;
        end
    end
    
    [p,q,C] = dpfast(S, stepcost); %Call Dan Ellis Mex file
    
    if any(abs(p-q)>=maxlag)
        lambda = lambda*1.1;
    else
        repeat = false;
    end
end

D = C(end,end); %similarity of A and B

Ta = interp1(q,p, 1:length(A));
Ta = smooth(Ta,5);

%PLOTTING
if nargin>3 && doplot
    C(C>1e100) = nan; figure, imshow(C,[]), hold on, plot(p,q,'r')
    figure, plot(B); hold on, plot(interp1(A,Ta));
    figure, plot(Ta-(1:length(Ta))')
end