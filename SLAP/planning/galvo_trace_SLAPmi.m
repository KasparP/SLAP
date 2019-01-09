function [AO, lineID] = galvo_trace_SLAPmi (res, tiling, calib, noplot)

%generates the desired galvo traces
    %res: number of samples for the trace
    %tiling:    the number of strips in which to acquire the image, larger
    %           numbers produce a bigger FOV
    %calib: calbration data specifying the X and Y offsets, etc.
    %noplot: if true, suppresses plot of galvo path
if nargin<2
    tiling = 1;
end
if nargin<4
    noplot = true;
end

offset = calib.galvos.offset; %a structure containing the X and Y offsets for each line
linelength = calib.galvos.linelength;
scanlength = tiling*linelength; %the scan length in volts

theta = pi/8; %the workign axis of the galvos (45degree) + tilt of the line (22.5deg)
R1 = [cos(theta) -sin(theta) ; sin(theta) cos(theta)]'; % rotation matrix
theta = 3*pi/8;
R2 = [cos(theta) -sin(theta) ; sin(theta) cos(theta)]'; % rotation matrix

%Make traces for lines 1 and 2
[Atrace, Btrace, L1_starts, L1_ends, L2_starts, L2_ends] = makeTraces('line1','line2');
L1_corners = repmat([offset.line1.X offset.line1.Y], 5, 1) + ([[1 1 -1 -1 1]' [1 -1 -1 1 1]']*scanlength/2)*R1;
L2_corners = repmat([offset.line2.X offset.line2.Y], 5, 1) + ([[1 1 -1 -1 1]' [1 -1 -1 1 1]']*scanlength/2)*R2;

%make traces for lines 3 and 4
[Ctrace, Dtrace, L3_starts, L3_ends, L4_starts, L4_ends] = makeTraces('line3','line4');
L3_corners = repmat([offset.line3.X offset.line3.Y], 5, 1) + ([[1 1 -1 -1 1]' [1 -1 -1 1 1]']*scanlength/2)*R1;
L4_corners = repmat([offset.line4.X offset.line4.Y], 5, 1) + ([[1 1 -1 -1 1]' [1 -1 -1 1 1]']*scanlength/2)*R2;

%correct the phase of the CD lines
Ctrace = circshift(Ctrace, [0 res/(4*tiling)]);
Dtrace = circshift(Dtrace, [0 res/(4*tiling)]);

%generate traces for the E and F galvos, which select between the cylinder assemblies entering the A/B and
%C/D galvos, respectively
buffer = round(res/50);
E = [ones(1,(2*tiling-1)*res-1)*offset.E(1) MAA_curve(offset.E(1), offset.E(2), 0, 0, 1, res-buffer) ...
    ones(1,(2*tiling-1)*res+buffer-2)*offset.E(2) MAA_curve(offset.E(2), offset.E(1), 0, 0, 1, res-buffer) ones(1,buffer-2)*offset.E(1)];
E = E(1:(4*tiling):end);
F = [ones(1,(2*tiling-1)*res-1)*offset.F(1) MAA_curve(offset.F(1), offset.F(2), 0, 0, 1, res-buffer) ...
    ones(1,(2*tiling-1)*res+buffer-2)*offset.F(2) MAA_curve(offset.F(2), offset.F(1), 0, 0, 1, res-buffer) ones(1,buffer-2)*offset.F(1)];
F = F(1:(4*tiling):end);
F = circshift(F, [0 res/(4*tiling)]);

%Combine the traces to return
AO = [Atrace ; Btrace ; Ctrace ; Dtrace ; E ; F]';

lineID = nan(size(AO,2), 1);
lin_ixs = repmat([true(1, res/(4*tiling)) false(1, res/(4*tiling))], 1, 2*tiling);
lineID(lin_ixs & abs(E-offset.E(1))<1e-10) = 1;
lineID(lin_ixs & abs(E-offset.E(2))<1e-10) = 2;
lineID(~lin_ixs & abs(F-offset.F(1))<1e-10) = 3;
lineID(~lin_ixs & abs(F-offset.F(2))<1e-10) = 4;

%Plot line1 and line2 scans
if ~noplot
    figure,
    subplot(1,2,1)
    scatter(L1_starts(:,1), L1_starts(:,2)), hold on, scatter(L1_ends(:,1), L1_ends(:,2), 'r');
    hold on, scatter(L2_starts(:,1), L2_starts(:,2)), hold on, scatter(L2_ends(:,1), L2_ends(:,2), 'r');
    plot(Atrace, Btrace);
    plot(L1_corners(:,1), L1_corners(:,2), 'k'), plot(L2_corners(:,1), L2_corners(:,2), 'k')
    axis image
    subplot(1,2,2)
    scatter(L3_starts(:,1), L3_starts(:,2)), hold on, scatter(L3_ends(:,1), L3_ends(:,2), 'r');
    hold on, scatter(L4_starts(:,1), L4_starts(:,2)), hold on, scatter(L4_ends(:,1), L4_ends(:,2), 'r');
    plot(Ctrace, Dtrace);
    plot(L3_corners(:,1), L3_corners(:,2), 'k'), plot(L4_corners(:,1), L4_corners(:,2), 'k')
    axis image
end

    function [X_ab, Y_ab, A_starts, A_ends, B_starts, B_ends] = makeTraces(ID1, ID2)
        A_starts = nan(tiling, 2); A_ends = nan(tiling, 2);
        A_tops = repmat([offset.(ID1).X offset.(ID1).Y], tiling,1)+[(scanlength/2)*ones(tiling,1)  linelength*(-(tiling-1)/2:(tiling-1)/2)']*R1 ;  %[Xpos  Ypos]
        A_bottoms = repmat([offset.(ID1).X offset.(ID1).Y], tiling, 1) + [-(scanlength/2)*ones(tiling,1)  linelength*(-(tiling-1)/2:(tiling-1)/2)']*R1 ;
        A_starts(1:2:end,:) = A_tops(1:2:end,:);    A_ends(1:2:end,:) = A_bottoms(1:2:end,:);
        A_starts(2:2:end,:) = A_bottoms(2:2:end,:); A_ends(2:2:end,:) = A_tops(2:2:end,:);
        
        %we can change which direction we come into the boxes:
        [A_starts, A_ends] = deal(A_ends, A_starts); %switch these variables
        
        B_starts = nan(tiling, 2); B_ends = nan(tiling, 2);
        B_tops = repmat([offset.(ID2).X offset.(ID2).Y], tiling, 1) + [(scanlength/2)*ones(tiling,1)  linelength*((tiling-1)/2:-1:-(tiling-1)/2)']*R2 ;  %[Xpos  Ypos]
        B_bottoms = repmat([offset.(ID2).X offset.(ID2).Y], tiling, 1) + [-(scanlength/2)*ones(tiling,1)  linelength*((tiling-1)/2:-1:-(tiling-1)/2)']*R2 ;
        B_starts(1:2:end,:) = B_tops(1:2:end,:);    B_ends(1:2:end,:) = B_bottoms(1:2:end,:);
        B_starts(2:2:end,:) = B_bottoms(2:2:end,:); B_ends(2:2:end,:) = B_tops(2:2:end,:);
        if ~mod(tiling,2)
            [B_starts, B_ends] = deal(B_ends, B_starts);
        end
        
        X_ab = makeTrace(A_starts(:,1), A_ends(:,1), B_starts(:,1), B_ends(:,1));
        Y_ab = makeTrace(A_starts(:,2), A_ends(:,2), B_starts(:,2), B_ends(:,2));
    end

    function trace = makeTrace(starts1,ends1,starts2,ends2)
        trace = [];
        for ix = 1:length(starts1)
            L = linspace(starts1(ix), ends1(ix), res+1);
            if ix<length(starts1)
                C = MAA_curve(ends1(ix), starts1(ix+1), ends1(ix)-starts1(ix), ends1(ix+1)-starts1(ix+1), 1, res);
            else
                C = MAA_curve(ends1(ix), starts2(1), ends1(ix)-starts1(ix), ends2(1)-starts2(1), 1, res);
            end
            trace = [trace L(1:end-1) C(1:end-1)]; %#ok<AGROW>
        end
        for ix = 1:length(starts2)
            L = linspace(starts2(ix), ends2(ix), res+1);
            if ix<length(starts2)
                C = MAA_curve(ends2(ix), starts2(ix+1), ends2(ix)-starts2(ix), ends2(ix+1)-starts2(ix+1), 1, res);
            else
                C = MAA_curve(ends2(ix), starts1(1), ends2(ix)-starts2(ix), ends1(1)-starts1(1), 1, res);
            end
            trace = [trace L(1:end-1) C(1:end-1)]; %#ok<AGROW>
        end
        trace = trace(1:4*tiling:end);
    end
end
