function SoftTrigOverdrive
%read patterns
readdir='W:\JJ\mask_SLM\';
writedir='W:\JJ\mask_SLM_lut\';
pattern=zeros(512,512,100,'uint8');
for i=1:100
   pattern(:,:,i)=imread([readdir num2str(i) '.bmp']);
end

% load lut
dataDir = 'E:\SLAPmidata';
calib = SLAPMi.loadCalibrationData([dataDir filesep 'Calibration\calibration.cal']);
lut = calib.SLM.lut;
lut=SLM_lut_recalib(lut,dlmread('C:\slm4411_at1064_P8.lut'));

% refit pattern
stack=size(pattern,3);
pattern2=double(pattern);
pattern2 = repmat(lut(:,:,1),[1 1 stack]).*(1-pattern2) + repmat(lut(:,:,2),[1 1 stack]).*pattern2;
pattern2=uint8(pattern2);

% calculate transient frames, but don't start the frame queue
hODP = hSI.hSlmScan.hSlm();
transientFrames = hODP.calculateTransientFrames(pattern2);

% creat timer
timing=zeros(1,100);
t0=0;
t = timer;
t.StartFcn = @startFrames;
t.TimerFcn = @updateFrames;
t.StopFcn = @logFrames;
t.Period = 0.01;
t.TasksToExecute = 100;
t.ExecutionMode = 'fixedDelay';
start(t)

    function startFrames(obj,event)
        disp([event.Type ' executed '...
            datestr(event.Data.time,'dd-mmm-yyyy HH:MM:SS.FFF')]);
        t0=event.Data.time(end);
    end

    function updateFrames(obj,event)
        frame=get(obj,'TasksExecuted');
        hODP.writeTransientFrames(transientFrames{frame});
        timing(frame)=event.Data.time(end)-t0;
    end

    function logFrames(obj,event)
        disp([event.Type ' executed '...
            datestr(event.Data.time,'dd-mmm-yyyy HH:MM:SS.FFF')]);
       save('timing','timing');
       save('pattern','pattern');
       delete(obj);
    end
end