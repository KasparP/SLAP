function calibrateTimingUsingBeads(scandata, P)
%fine calibrate the slapmi timing. scandata is a scan of some small number
%of beads


if nargin<2
    P = linePSF_full(scandata);
end

select = false(length(P.coords{1}), length(P.coords{2}), length(P.coords{3}));
select(:,:, ceil(end/2)) = true;
P.P = P.P(:,select); P.coords{3} = 0;

Y = mean([scandata.frames.pmtData],2);
shifts = -20:20;

for line = 4:-1:1
    line_ixs = scandata.line==line;
    Pline{line} = P.P(line_ixs,:)'*Y(line_ixs);
end

for line = 1:4
    sel3  = (1:4)~=line;
    P3Y=prod(cat(3,Pline{sel3}),3);
    score = nan(1, length(shifts));
    for shiftix = 1:length(shifts)
        P4Y = P.P(line_ixs,:)'*interp1(Y(line_ixs), (1:sum(line_ixs))+shifts(shiftix))';
        
        score(shiftix) = nansum(P3Y.*P4Y);
        
%         figure, imshow(reshape(P4Y, length(P.coords{1}), length(P.coords{2})),[])
    end
    keyboard
end
%get 

keyboard