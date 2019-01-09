function counts = NND (data, tau)
%performs fast nonnegative deconvolution on pmt signal to solve for minimum MSE
%photon rate
% counts:   estimated photon rate
% data  :   The data to be deconvolved
% tau   :   The time constant of the PMT, in data samples

T = length(data);
counts = zeros(1, T);
counts(end) = max(0,data(end));
cutoff = ceil(8*tau);       %how long the effect of a timepoint can be
k = exp(-(0:cutoff)/tau);   %the convolution kernel
recent = [T nan(1,T-1)];    %stored locations where we assigned counts
recent_ix = 1;

%the points that could potentially be assigned counts:
points = data(1:end-1)>(k(2).*[0; data(1:end-2)]) & data(1:end-1)>0;
%dividing these points up into runs, for speed
runstarts = find(points & ~[false ; points(1:end-1)]);
runends = find(points & ~[points(2:end) ; false]);
run_id = length(runends);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%UNCOMMENT FOR DISPLAY
%h = figure; plot(data);
%display_deconv;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while run_id>0
    %%%%%%%%%%%%%%%%%%%%
    %FOR DISPLAY ONLY  % 
    %delete([h1 h2]);  %
    %display_deconv;   %
    %pause(0.5)        %
    %drawnow;          %
    %%%%%%%%%%%%%%%%%%%%
    
    %initialization
    converged = false;
    oldtop = 0;
    oldbottom = 0;
    t = runends(run_id);
    t1 = t;
    accum = 0;
    
    while ~converged
        if recent_ix && recent(recent_ix)<(t+cutoff)
            t2 = recent(recent_ix);
            C_max = counts(t2)./k(t2-t+1);
        else
            t2 = min(t+cutoff, T+1);
            C_max = inf;
        end
        
        %b: kernel
        b = k(t1-t+1:t2-t);
        top = b*data(t1:t2-1)+oldtop; %this is the numerator of the least squares fit for an exponential
        bottom = sum(b.^2)+oldbottom; %this is the denominator of the fit
        done = false;
        
        while ~done
            %the error function is (data-kernel.*C)^2
            best_C = max(top/bottom, 0);  %C=top/bottom sets the derivative of the error to 0
            
            if best_C>(C_max+accum) %does not meet nonnegative constraint. Continue to adjust previous solutions.
                if counts(t2)
                    accum = accum + counts(t2)/k(t2-t+1);
                    counts(t2) = 0;
                end
                t1 = t2;
                oldtop = top;
                oldbottom = bottom;
                recent_ix = recent_ix-1;
                
                done = true;
            else                    %converged!
                %now that we have found the MSE counts for times t<end, check if
                %this will be swamped by the next timepoint in the run
                if  t==runstarts(run_id) || data(t-1)<(best_C/k(2)) %C_max won't necessarily get swamped
                    if recent_ix && t2<=t+cutoff
                            counts(t2) = counts(t2) - (best_C-accum)*k(t2-t+1);
                    end
                    run_start = runstarts(run_id);
                    init_ix = recent_ix + 1;
                    recent_ix = recent_ix + 1 + t - run_start;
                    recent(init_ix:recent_ix) = t:-1:run_start;
                    
                    counts(runstarts(run_id):t) = [data(run_start:t-1) ; best_C] - [0 ; k(2)*data(run_start:t-1)];
                    
                    done = true;
                    converged = true;
                else %C_max will get swamped
                    %in this situation, we know that this point will be removed
                    %as we continue to process the run. To save time:
                    t = t-1;
                    runends(run_id) = t;
                    accum = accum/k(2);
                    top = top*k(2)+data(t); % %this is the correct adjustment to the derivative term above
                    bottom = bottom*(k(2)^2)+1; % %this is the correct adjustment to the derivative term above
                end
            end
        end
    end   
    run_id = run_id-1;
end

% function display_deconv
%     hold on, 
%     h1 = plot(counts, 'g');
%     fit = conv(counts, k);
%     hold on
%     h2 = plot(fit(1:T), 'r');
% end
end