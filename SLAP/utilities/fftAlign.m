function [y_align,lag]=fftAlign(x, y)

fty = fft(y);
ftx = fft(x);

xc = ifft(conj(fty).*ftx);
[~,i] = max(abs(xc));
y_align = circshift(y,i);
lag = i;

end