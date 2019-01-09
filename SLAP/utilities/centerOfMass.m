function index = centerOfMass (input, dim)
    if any(input<0)
        warning('Center of mass encountered negative data!');
    end
    if nargin<2
        dim=2;
    end
    if dim==1
        input = input';
    end
        
    index  = input*(1:size(input,2))' ./ sum(input, 2);
    
    if dim==1
        index = index';
    end