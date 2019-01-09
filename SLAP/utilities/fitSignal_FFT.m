% *******************************************************
% delay-matching between two signals (complex/real-valued)
% M. Nentwig
%
% * matches the continuous-time equivalent waveforms
%   of the signal vectors (reconstruction at Nyquist limit =>
%   ideal lowpass filter)
% * Signals are considered cyclic. Use arbitrary-length 
%   zero-padding to turn a one-shot signal into a cyclic one.
%
% * output:
%   => coeff: complex scaling factor that scales 'ref' into 'signal'
%   => delay 'deltaN' in units of samples (subsample resolution)
%      apply both to minimize the least-square residual      
%   => 'shiftedRef': a shifted and scaled version of 'ref' that 
%      matches 'signal' 
%   => (signal - shiftedRef) gives the residual (vector error)
%
% *******************************************************
function [coeff, shiftedRef, deltaN] = fitSignal_FFT(signal, ref)
    n=length(signal);
    % xyz_FD: Frequency Domain
    % xyz_TD: Time Domain
    % all references to 'time' and 'frequency' are for illustration only

    forceReal = isreal(signal) && isreal(ref);
    
    % *******************************************************
    % Calculate the frequency that corresponds to each FFT bin
    % [-0.5..0.5[
    % *******************************************************
    binFreq=(mod(((0:n-1)+floor(n/2)), n)-floor(n/2))/n;

    % *******************************************************
    % Delay calculation starts:
    % Convert to frequency domain...
    % *******************************************************
    sig_FD = fft(signal);
    ref_FD = fft(ref, n);

    % *******************************************************
    % ... calculate crosscorrelation between 
    % signal and reference...
    % *******************************************************
    u=sig_FD .* conj(ref_FD);
    if mod(n, 2) == 0
        % for an even sized FFT the center bin represents a signal
        % [-1 1 -1 1 ...] (subject to interpretation). It cannot be delayed. 
        % The frequency component is therefore excluded from the calculation.
        u(length(u)/2+1)=0;
    end
    Xcor=abs(ifft(u));

    %  figure(); plot(abs(Xcor));
    
    % *******************************************************
    % Each bin in Xcor corresponds to a given delay in samples.
    % The bin with the highest absolute value corresponds to
    % the delay where maximum correlation occurs.
    % *******************************************************
    integerDelay = find(Xcor==max(Xcor));
    
    % (1): in case there are several bitwise identical peaks, use the first one
    % Minus one: Delay 0 appears in bin 1
    integerDelay=integerDelay(1)-1;

    % Fourier transform of a pulse shifted by one sample
    rotN = exp(2i*pi*integerDelay .* binFreq);

    uDelayPhase = -2*pi*binFreq;
    
    % *******************************************************
    % Since the signal was multiplied with the conjugate of the
    % reference, the phase is rotated back to 0 degrees in case
    % of no delay. Delay appears as linear increase in phase, but
    % it has discontinuities.
    % Use the known phase (with +/- 1/2 sample accuracy) to 
    % rotate back the phase. This removes the discontinuities.
    % *******************************************************
    %  figure(); plot(angle(u)); title('phase before rotation');
    u=u .* rotN;
    
    % figure(); plot(angle(u)); title('phase after rotation');
    
    % *******************************************************
    % Obtain the delay using linear least mean squares fit
    % The phase is weighted according to the amplitude.
    % This suppresses the error caused by frequencies with
    % little power, that may have radically different phase.
    % *******************************************************
    weight = abs(u); 
    constRotPhase = 1 .* weight;
    uDelayPhase = uDelayPhase .* weight;
    ang = angle(u) .* weight;
    r = [constRotPhase; uDelayPhase] .' \ ang.'; %linear mean square
    
    %rotPhase=r(1); % constant phase rotation, not used.
    % the same will be obtained via the phase of 'coeff' further down
    fractionalDelay=r(2);
    
    % *******************************************************
    % Finally, the total delay is the sum of integer part and
    % fractional part.
    % *******************************************************
    deltaN = integerDelay + fractionalDelay;

    % *******************************************************
    % provide shifted and scaled 'ref' signal
    % *******************************************************
    % this is effectively time-convolution with a unit pulse shifted by deltaN
    rotN = exp(-2i*pi*deltaN .* binFreq);
    ref_FD = ref_FD .* rotN;
    shiftedRef = ifft(ref_FD);
    
    % *******************************************************
    % Again, crosscorrelation with the now time-aligned signal
    % *******************************************************
    coeff=sum(signal .* conj(shiftedRef)) / sum(shiftedRef .* conj(shiftedRef));
    shiftedRef=shiftedRef * coeff;

    if forceReal
        shiftedRef = real(shiftedRef);
    end
end