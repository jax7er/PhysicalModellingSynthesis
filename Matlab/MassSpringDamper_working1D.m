% A Single Mass-Spring-Damper Simple Harmonic Oscillator as originally
% demonstrated in PhyMod
% DTM 21/01/2016

function Copy_of_MassSpringDamper()
clc;
close all; 

num_masses = 1000; % mass 1 and num_masses are fixed

% System Definitions
Fs = 44100;         % Sample Rate for audio output
T = 1/Fs;           % Sample period
% seconds = input('How many seconds would you like me to last? ');
seconds = 2;
numSamples = 44100*seconds;          % Number of samples in output
frequency = 100;
excitation_index = round(num_masses/2, 0);
excitation_amplitude = 1;
      
xc = zeros(1, num_masses);
% xc(25:29) = hanning(5)';
xc(excitation_index) = excitation_amplitude;
xp = xc;
xn = zeros(1, num_masses);

% m = 100 .* ones(1, num_masses) ./ 1000000000; % mass given in kg*10e-9
% 
% k =  m .* (2*pi*frequency)^2; % omega0 = sqrt(k/m) = 2*pi*frequency
%                               % k = m*omega0^2 = m*(2*pi*frequency)^2
% z = 100 .* ones(1, num_masses) ./ 1000000000; % z given in N*s/(kg*10e-9)
% 
% k_ = (k .* T^2) ./ m;
% z_ = (z .* T) ./ m;

% set all masses, spring constants and damping factors to be equal
m1 = 100 / 1000000000;
% k1 = (m1 * (2*pi*frequency)^2) * (T^2 / m1);
k1 = 10 * (T^2 / m1);
z1 = (1000 / 1000000000) * (T / m1);

out = zeros(1, numSamples);

% Main Loop
for n = (1 : numSamples) 
%     for i = (2 : num_masses - 1)
%         xn(i) = xc(i) * (k_(i) + z_(i) + k_(i - 1) + z_(i - 1) + 2) - xc(i + 1) * (k_(i) + z_(i)) - xc(i - 1) * (k_(i - 1) + z_(i - 1)) ...
%               - xp(i) * (z_(i) + z_(i - 1) + 1)                     + xp(i + 1) * z_(i)           + xp(i - 1) * z_(i - 1);
%     end
    
    for i = (2 : num_masses - 1)
        xn(i) = xc(i) * 2 * (1 - k1 - z1) + (xc(i + 1) + xc(i - 1)) * (k1 + z1) ...
              - xp(i) * (1 - 2 * z1)      - (xp(i + 1) + xp(i - 1)) * z1;
    end
    
%     out(n) = xn(excitation_index);
    out(n) = xn(2);
    
%     Resamples the output at a rate of Fs/200 to plot the mass displacement in real-time.
%     if mod(n, 100) == 0   
%         plot(xn); 
%         axis([0 num_masses -0.1 0.1]);
%         pause(0.0001);
%     end
    
    xp = xc;
    xc = xn;
end
    
% Plot time domain response of output
figure(2);
plot(out);
xlabel('Time (samples)');
ylabel('Amplitude');

fftSize = 8192;
f = (0:fftSize-1)*(Fs/fftSize);

% Plot frequency domain response of output
figure(3);
semilogx(f,20*log10(abs(fft(out,fftSize))/max(abs(fft(out,fftSize)))));
axis([0 8000 -40 0]);
xlabel('Frequency (Hz)');
ylabel('Magnitude Response (dB)');

sound(out ./ max(out), Fs);
