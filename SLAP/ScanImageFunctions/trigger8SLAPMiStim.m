function trigger8SLAPMiStim(hSLAPMi, hSI)
    ts = timerfind('Userdata', 8);
    length(ts)
    delete(ts);
    
    frameRate = hSI.hRoiManager.scanVolumeRate;
    disp(frameRate)
    
    frameInterval = 1/frameRate;
    period = round(5/frameInterval)*frameInterval; %make ISI close to 5 seconds
    
    hTaskTriggerOutOnDmd = most.util.safeCreateTask('TriggerOutOnDemand');
    hTaskTriggerOutOnDmd.createDOChan(hSLAPMi.beamDaqName, 'PFI0' , 'TriggerOutOnDmd');
    hTimer = timer( 'ExecutionMode', 'fixedRate', 'StartDelay', period-1, 'Period', period, 'TasksToExecute', 8, 'TimerFcn', @tFunc, 'StopFcn', @stopFunc, 'Userdata', 8);
    start(hTimer)
    
function stopFunc(x,y)
    delete(hTaskTriggerOutOnDmd);
    delete(hTimer);
    disp('Done stimuli');
end

function tFunc(x,y)
    hTaskTriggerOutOnDmd.writeDigitalData(false, 1, true);
    hTaskTriggerOutOnDmd.writeDigitalData(true, 1, true);
    disp(['Stimulus triggered at ' datestr(now)])
    hTaskTriggerOutOnDmd.writeDigitalData(false, 1, true); %turn off beam
end
end