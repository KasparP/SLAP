function numericMicroscope(hSI)

chans = hSI.hDisplay.lastFrameChannels;
nChan = length(chans);
ydata = nan(30, nChan);

hF = figure; hAx = axes; hL = plot(ydata);
updateY;

set(hF, 'closeRequest', @closehF)

t = timer;
t.TimerFcn = @(~,thisEvent)(updateY);
t.Period = 0.3;
t.ExecutionMode = 'fixedRate';
start(t)


function closehF(varargin)
    delete(hF);
    stop(t); 
    delete(t)
end

    function updateY(varargin)
        for chan = chans(:)'
            ydata(:, chan) = [ydata(2:end,chan) ; sum(sum(hSI.hDisplay.lastFrame{chan}))];
            set(hL(chan), 'ydata', ydata(:,chan));
        end
        
    end
end