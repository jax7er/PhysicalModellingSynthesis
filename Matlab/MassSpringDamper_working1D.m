
% prepare Matlab for simulation
clc;       % clear console
close all; % close all open figures

% define audio parameters
Fs = 44100;                   % sample rate for audio output
Ts = 1 / Fs;                  % sample period
seconds = 4;                  % number of seconds of output to be calculated
Ns = 44100 * seconds;         % number of samples in output
all_samples_range = (1 : Ns); % range containing all masses

% define mass-spring-damper system parameters
Nm = 100;                           % total number of masses
mobile_masses_range = (2 : Nm-1);   % range containing all masses with a mobile mass on each side
f0 = 1000;                          % resonant freqency of system in Hz
m = 0.001 * ones(1, Nm);            % mass given in kg
m_ = [m m(Nm)];                     % length Nm+1 with last element equal to the right-most mass
k_ = (2*pi*f0)^2 * m_;              % k given in N/m or kg/m*s^2 -> omega0 = 2*pi*f0 = sqrt(k*l/m))
z_crit = (4 * m_ .* k_).^(-1/2);    % z^2 - 4*m*k = 0 -> http://hyperphysics.phy-astr.gsu.edu/hbase/oscda.html#c1
z_ = z_crit;                        % z given in N*s/m or kg/m*s

% scale spring constant and damping factor by appropriate factors to simplify the loop calculation
k = (k_ ./ m_) * Ts^2;
z = (z_ ./ m_) * Ts; 

% define excitation of the system
Ei = round(Nm/2, 0); % index of center of excitation
Ea = Nm / 10;         % amplitude of excitation

% initialise the displacements of all masses
x_next = zeros(1, Nm); % set all next displacements to zero
x_curr = x_next;       % set all current displacements to zero
x_curr(Ei) = Ea;       % insert excitation into current displacements
x_prev = x_curr;       % set previous displacements to be equal to current

% define index of mass whose displacement will be taken as output
Oi = Nm/5;          % index of mass to use as output
O = zeros(1, Ns); % output buffer for samples

% initialise loop parameters
plot_enable = 1;          % enables the real-time plot of the masses' displacement
plot_update_delay = 0.01; % seconds between updating real-time graph
dispstat('', 'init');     % initialise dispstat so that % progress is displayed properly
f0_scale = 2;

% main loop
for n = all_samples_range
    if mod(n, Ns/4) == 0
        x_next = zeros(1, Nm); % set all next displacements to zero
        x_curr = x_next;       % set all current displacements to zero
        x_curr(Ei) = Ea;       % insert excitation into current displacements
        x_prev = x_curr;       % set previous displacements to be equal to current
        f0 = f0 * f0_scale;
        k = k * f0_scale;
    end
    
    % calculate next displacement of first mass simplified with x(i-1) terms disappearing
    x_next(1) = x_curr(1) * (2 - k(2) - z(2) - k(1) - z(1)) + x_curr(2) * (k(2) + z(2)) ...
              - x_prev(1) * (1 -        z(2) -        z(1)) - x_prev(2) *         z(2);
          
    % calculate next displacement of middle masses
    for i = mobile_masses_range
        x_next(i) = x_curr(i) * (2 - k(i+1) - z(i+1) - k(i) - z(i)) + x_curr(i+1) * (k(i+1) + z(i+1)) + x_curr(i-1) * (k(i) + z(i)) ...
                  - x_prev(i) * (1 -          z(i+1) -        z(i)) - x_prev(i+1) *           z(i+1)  - x_prev(i-1) *         z(i);
    end
    
    % calculate next displacement of last mass simplified with x(i+1) terms disappearing
    x_next(Nm) = x_curr(Nm) * (2 - k(Nm+1) - z(Nm+1) - k(Nm) - z(Nm)) + x_curr(Nm-1) * (k(Nm) + z(Nm)) ...
               - x_prev(Nm) * (1 - z(Nm+1) -                   z(Nm)) - x_prev(Nm-1) *          z(Nm);
    
    % set output sample to new displacement of the output mass
    O(n) = x_next(Oi);
    
    % plot the mass displacement in real-time
    if (plot_enable == 1) && (mod(n, Fs*plot_update_delay) == 0)
        subplot(2,1,1), plot(x_next); 
        axis([1 Nm -Ea/4 Ea/4]), xlabel('Mass'), ylabel('Amplitude');
        title(sprintf('m=%fkg, k=%fN/m, z=%fN*s/m, f0=%dHz', m(1), k(1), z(1), f0));
        subplot(2,1,2), plot(O); 
        axis([1 Ns -Ea/4 Ea/4]), xlabel('Sample'), ylabel(sprintf('Amplitude of mass %d', Oi));
        pause(plot_update_delay);
    end
    
    % display progress of simulation in percent
    if mod(n, Ns/100) == 0
        dispstat(sprintf('%d%%', 100*n/Ns));
    end
    
    % set up arrays for next iteration
    x_prev = x_curr;
    x_curr = x_next;
end

% check whether there are any output samples
if abs(max(O)) > 0 
    Onorm = O / max(O);
    
    % Plot time domain response of output
    figure(2);
    plot(Onorm);
    xlabel('Time (samples)');
    ylabel('Amplitude');

    fftSize = 8192;
    f = (0:fftSize-1)*(Fs/fftSize);

    % Plot frequency domain response of output
    figure(3);
    semilogx(f, 20 * log10(abs(fft(Onorm, fftSize)) / max(abs(fft(Onorm, fftSize)))));
    axis([0 8000 -40 0]);
    xlabel('Frequency (Hz)');
    ylabel('Magnitude Response (dB)');

    sound(Onorm, Fs);
else
    disp('All output samples are 0');
end
