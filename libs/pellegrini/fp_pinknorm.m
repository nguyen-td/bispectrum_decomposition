function signal_t = fp_pinknorm(signal)
% Normalize signal to 1/f shape 
%
% Copyright (c) 2023 Franziska Pellegrini and Stefan Haufe

time_pnts = size(signal,1);

%get 1/f shape 
n1=2*ceil((time_pnts-1)/2);
scal=sqrt(1./[1:(n1-2)/2]);
ff=zeros(n1,1);
ff(2:(n1)/2,:)=scal;

%fft of signal
signal_fft = fft(signal);

%1/f transform of fft of signal and get back to time domain 
signal_t=2*real(ifft(signal_fft.*ff));
signal_t=signal_t(1:time_pnts,:);