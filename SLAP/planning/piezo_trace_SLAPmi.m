function [trace, powertrace] = piezo_trace_SLAPmi(res, optsin, calib)

opts.Zmin = -1;
opts.Zmax = +1;
opts.duty = 0.8;
opts.bidi = true;
if nargin>1 %UPDATE WITH USER-SPECIFIED OPTIONS
    for field = fieldnames(optsin)'
         opts.(field{1}) = optsin.(field{1});
    end
end
%sanity checks
if opts.duty>=1 || opts.duty<=0
    error('Invalid duty cycle')
end

opts.duty = round(opts.duty*res)/res;


if ~opts.bidi
    line = linspace(-0.5,0.5, opts.duty*res+1);
    connector = MAA_curve(0.5, -0.5, 1/opts.duty, 1/opts.duty, (1-opts.duty), (1-opts.duty)*res);
    trace = [line(1:end-1) connector(1:end-1)];
    trace = (trace - min(trace))/(max(trace)-min(trace));
    trace = trace'*(opts.Zmax-opts.Zmin) + opts.Zmin;
    
    %power trace can be zeroed at the turnaround points but currently isn't
    powertrace = ones(size(trace)) * calib.pockels.slow(2); %always on
else
    %produce a near-triangle wave with minimum absolute acceleration
    line = linspace(-0.5,0.5, opts.duty*res+1);
    connector = MAA_curve(0.5, 0.5, 1/opts.duty, -1/opts.duty, (1-opts.duty), (1-opts.duty)*res);
    trace = [line(1:end-1) connector(1:end-1) -line(1:end-1) -connector(1:end-1)];
    
    trace = (trace - min(trace))/(max(trace)-min(trace));
    trace = trace*(opts.Zmax-opts.Zmin) + opts.Zmin;
    trace = circshift(trace, [0 floor(length(connector)/2)-1]);
    trace = trace(1:2:end)';

    %power trace can be zeroed at the turnaround points but currently isn't
    powertrace = ones(size(trace)) * calib.pockels.slow(2); %always on
end

end