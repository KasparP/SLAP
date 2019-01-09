function deleteTimers(varargin)
    ts = [timerfind('Userdata', 8); timerfind('Userdata', 1)];
    disp(['Stopping ' int2str(length(ts)) ' timers'])
    for t =1:length(ts)
        stop(ts(t));
        delete(ts(t));
    end
    
end