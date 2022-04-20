clear
close all

% Declare the input for a cos(2*pi*f*t)
n = 1000;
vecTime = linspace(0,1,n);
freq = 1/0.03;  % Hz
vecInputV = sin(2*pi*freq*vecTime);

Fs = n;
df = Fs/length(vecInputV);
vecFreq = -Fs/2:df:Fs/2-df;
fftVin = fftshift(fft(vecInputV))/n;
plot(vecFreq, abs(fftVin), "-b.")  % Plot the input fft
hold on 
title("FFT of the step input")
xlabel("Omega (rad/s)")
ylabel("V (dB)")
grid on


% vecOmega = wspace(vecTime);  % Create the omega vector
% fftVin = abs(fft(vecInputV))/n; % Input fft in dB
% plot(vecOmega, fftVin, "-b.")  % Plot the input fft
% hold on 
% title("FFT of the step input")
% xlabel("Omega (rad/s)")
% ylabel("V (dB)")
% grid on



% This function constructs a linearly-spaced vector of angular
% frequencies that correspond to the points in an FFT spectrum.
function w = wspace(t,nt)
% This function constructs a linearly-spaced vector of angular
% frequencies that correspond to the points in an FFT spectrum.
% The second half of the vector is aliased to negative
% frequencies.
%
% USAGE
%
% w = wspace(tv);
% w = wspace(t,nt);
%
% INPUT
%
% tv - vector of linearly-spaced time values
% t - scalar representing the periodicity of the time sequence
% nt - Number of points in time sequence
%      (should only be provided if first argument is scalar)
%
% OUTPUT
%
% w - vector of angular frequencies
%
% EXAMPLE
%
%   t = linspace(-10,10,2048)';   % set up time vector
%   x = exp(-t.^2);               % construct time sequence
%   w = wspace(t);                % construct w vector
%   Xhat = fft(x);                % calculate spectrum
%   plot(w,abs(Xhat))             % plot spectrum
%
% AUTHOR:  Thomas E. Murphy (tem@umd.edu)

if (nargin<2)
    nt = length(t);
    dt = t(2) - t(1);
    t = t(nt) - t(1) + dt;
end

if (nargin == 2)
    dt = t/nt;
end

w = 2*pi*(0:nt-1)'/t;
kv = find(w >= pi/dt);
w(kv) = w(kv) - 2*pi/dt;
end