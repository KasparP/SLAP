function [AOfast, AOslow] = pockels_trace_SLAPmi(NSampsOut, hSLAPMi)

calib = hSLAPMi.calib;
tiling = hSLAPMi.tiling;

deadsamps = calib.pockels.deadsamps;
fast = calib.pockels.fast;
slow = calib.pockels.slow;

% make fast pockels trace
AOfast = nan(NSampsOut, 1);
step = step_settle(1,NSampsOut/(4*tiling),(4*tiling*deadsamps)/NSampsOut);
AOfast(:,1) = fast(1) + (fast(2)-fast(1))*repmat([step -step+1 step -step+1], 1, tiling);

%make slow pockels trace
AOslow = hSLAPMi.calib.B.scaleby(hSLAPMi.lineIDs);
NSampsSlow = round((hSLAPMi.sampleRate/hSLAPMi.galvos.AOrate)*NSampsOut);
x1 = linspace(0,1, length(AOslow)+1); x2 = linspace(0,1,NSampsSlow+1); 
AOslow = interp1(x1(1:end-1), AOslow, x2(1:end-1))';

end