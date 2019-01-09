function Y = step_settle(speed, nsamples, dutyrange)

%produces a trace that steps from 0 to 1 over nsamples, achieving the final
%level within dutyrange x 100% of the total samples
if dutyrange<0 || dutyrange>1 ||speed<0
    error('Bad argument')
end

x = linspace(0,1,floor(nsamples*dutyrange));
curve = x.^(2^speed);

Y = [fliplr(-curve+1) ones(1, nsamples-length(curve))];
