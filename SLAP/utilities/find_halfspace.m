 function middle = find_halfspace(x,target)
    % returns the number of the last sample
    % in the ASCENDING array data that is
    % < x using a simple half space search
    first = 1;
    last = length(x);
    while first <= last
        middle = floor((first + last)/2);
        if target < x(middle)
            last = middle - 1;
        elseif target > x(middle)
            first = middle + 1;
        else
            break
        end
    end
end