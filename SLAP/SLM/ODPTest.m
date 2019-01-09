%%
dirpath='W:\JJ\masks3\';
m_SI=false(1280,1280,100);
for i=1:100
    m_SI(:,:,i)=logical(imread([dirpath num2str(i) '.bmp']));
end


%%
writedir='W:\JJ\mask_SLM\';
for i=1:size(m_SI,3)
    m_SLM(:,:,i)=logical(mapImageToSLM(m_SI(:,:,i), hSI.hSlmScan.testPatternPixelToRefTransform, scanfield));
    imwrite(m_SLM(:,:,i),[writedir num2str(i) '.bmp']);
end

%%
dataDir = 'E:\SLAPmidata';
calib = SLAPMi.loadCalibrationData([dataDir filesep 'Calibration\calibration.cal']);
lut = calib.SLM.lut;
lut=SLM_lut_recalib(lut,dlmread('C:\slm4411_at1064_P8.lut'));
pattern=m_SLM;
%%
pattern=rand(512,512,100);
stack=size(pattern,3);
%pattern = repmat(lut(:,:,1),[1 1 stack]).*(1-pattern) + repmat(lut(:,:,2),[1 1 stack]).*pattern;
pattern=uint8(pattern);
%%
%hODP = hSI.hSlmScan.hSlm();
% calculate transient frames, but don't start the frame queue
transientFrames = hODP.calculateTransientFrames(pattern);
% write the transient frames to the SLM manually
%%
for j=1:5
 for idx = 1:100
     hODP.writeTransientFrames(transientFrames{idx});
     pause(0.01);
 end
end
%%
% pattern
pattern=rand(512,512,100);
figure();
timing=zeros(100,1);
t = timer;
t.StartFcn = @(~,thisEvent)disp([thisEvent.Type ' executed '...
    datestr(thisEvent.Data.time,'dd-mmm-yyyy HH:MM:SS.FFF')]);
t.TimerFcn = @(~,thisEvent)disp(['Event ' num2str(get(t,'TasksExecuted')) ' executed '...
      datestr(thisEvent.Data.time,'dd-mmm-yyyy HH:MM:SS.FFF')]);
%t.TimerFcn = @(~,thisEvent) hODP.writeTransientFrames(transientFrames{get(t,'TasksExecuted')});
%t.TimerFcn = @(~,thisEvent) imshow(pattern(:,:,get(t,'TasksExecuted')),[]);
t.StopFcn = @(~,thisEvent)disp([thisEvent.Type ' executed '...
    datestr(thisEvent.Data.time,'dd-mmm-yyyy HH:MM:SS.FFF')]);
t.Period = 0.01;
t.TasksToExecute = 100;
t.ExecutionMode = 'fixedRate';
start(t)