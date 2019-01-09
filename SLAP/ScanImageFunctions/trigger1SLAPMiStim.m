function trigger1SLAPMiStim(hSLAPMi, delay)
if nargin<2
    delay = 1;
end
    ts = timerfind('Userdata', '8stimTimer');
    delete(ts);
    ts = timerfind('Userdata', 1);
    delete(ts);
    
    hTimer = timer( 'ExecutionMode', 'fixedRate', 'StartDelay', delay, 'TasksToExecute', 1, 'TimerFcn', @tFunc, 'StopFcn', @stopFunc, 'Userdata', 1);
    start(hTimer)
    
function stopFunc(x,y)
    delete(hTimer);
    disp('Done stimuli');
end

function tFunc(x,y)
   triggerSLAPMiStim(hSLAPMi);
end
end