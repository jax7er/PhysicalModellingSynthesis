% prepare Matlab for simulation
clc;       % clear console
close all; % close all open figures

% define audio parameters
Fs = 44100;                   % sample rate for audio output
Ts = 1 / Fs;                  % sample period
seconds = 1;                  % number of seconds of output to be calculated
Ns = 44100 * seconds;         % number of samples in output
all_samples_range = (1 : Ns); % range containing all masses

% define mass-spring-damper system parameters
Nm = 100;                           % total number of masses
mobile_masses_range = (2 : Nm-1);   % range containing all masses with a mobile mass on each side
m = 0.001 * ones(1, Nm);            % mass given in kg
m_ = [m m(Nm)];                     % length Nm+1 with last element equal to the right-most mass
k_ = -10000 * ones(1, Nm+1);              % k given in N/m or kg/m*s^2 -> f0 = sqrt(k*l/m)/(2*pi)
f0_1 = sqrt(2 * k_ / m(1)) / (2*pi); % first resonant freqency of system in Hz -> f0 = sqrt(k*l/m)/(2*pi)
% z_crit = (4 * m_ .* k_).^(-1/2);    % z^2 - 4*m*k = 0 -> http://hyperphysics.phy-astr.gsu.edu/hbase/oscda.html#c1
z_ = 0.01 * ones(1, Nm+1);%*z_crit;                    % z given in N*s/m or kg/m*s

% scale spring constant and damping factor by appropriate factors to simplify the loop calculation
K = (k_ ./ m_) * Ts^2;
Z = (z_ ./ m_) * Ts; 

% define excitation of the system
Ei = round(Nm/2, 0); % index of center of excitation
Ea = Nm / 10;         % amplitude of excitation

% initialise the displacements of all masses
x_next = zeros(1, Nm); % set all next displacements to zero
x_curr = x_next;       % set all current displacements to zero
x_curr(Ei) = Ea;       % insert excitation into current displacements
x_prev = x_curr;       % set previous displacements to be equal to current

% define index of mass whose displacement will be taken as output
Oi = round(Nm/5, 0);          % index of mass to use as output
O = zeros(1, Ns); % output buffer for samples

% initialise loop parameters
k_scale_enable = 0;          % enables multiple pitches
k_scale_factor = 4;
plot_enable = 0;           % enables the real-time plot of the masses' displacement
plot_update_delay = 0.01;  % seconds between updating real-time graph
delay_enable = 0;
delay_interval = 0.25;
delay_gain = 0.5;
reverb_enable = 0;            % enables reverb in the output
reverb_density = 0.9;      % \
reverb_feedback = 1;       % |- reverb parameters
reverb_feedthrough = 0.75; % /
figure_enable = 1;
write_audio = 1;
dispstat('', 'init');      % initialise dispstat so that % progress is displayed properly

dispstat('Simulating... ', 'keepthis');

% main loop
for n = all_samples_range    
    % scale k by a fixed factor after a number of samples
    if (k_scale_enable == 1) && (mod(n, Ns/4) == 0) && (n < Ns)
        x_next = zeros(1, Nm); % set all next displacements to zero
        x_curr = x_next;       % set all current displacements to zero
        x_curr(Ei) = Ea;       % insert excitation into current displacements
        x_prev = x_curr;       % set previous displacements to be equal to current
        k_ = k_ * k_scale_factor;
        K = (k_ ./ m_) * Ts^2;
    end
    
    % calculate next displacement of first mass simplified with x(i-1) terms disappearing
    x_next(1) = (x_curr(1) * (2 - K(2) - Z(2) - K(1) - Z(1)) + x_curr(2) * (K(2) + Z(2)) ...
               - x_prev(1) * (1 -        Z(2) -        Z(1)) - x_prev(2) *         Z(2));
          
    % calculate next displacement of middle masses
    for i = mobile_masses_range
        x_next(i) = (x_curr(i) * (2 - K(i+1) - Z(i+1) - K(i) - Z(i)) + x_curr(i+1) * (K(i+1) + Z(i+1)) + x_curr(i-1) * (K(i) + Z(i)) ...
                   - x_prev(i) * (1 -          Z(i+1) -        Z(i)) - x_prev(i+1) *           Z(i+1)  - x_prev(i-1) *         Z(i));
    end
    
    % calculate next displacement of last mass simplified with x(i+1) terms disappearing
    x_next(Nm) = (x_curr(Nm) * (2 - K(Nm+1) - Z(Nm+1) - K(Nm) - Z(Nm)) + x_curr(Nm-1) * (K(Nm) + Z(Nm)) ...
                - x_prev(Nm) * (1 - Z(Nm+1) -                   Z(Nm)) - x_prev(Nm-1) *          Z(Nm));
    
    % set output sample to new displacement of the output mass
    O(n) = x_next(Oi);
    
    % plot the mass displacement in real-time
    if (plot_enable == 1) && (mod(n, Fs*plot_update_delay) == 0)
        subplot(2,1,1), plot(x_next); 
        axis([1 Nm -Ea/4 Ea/4]), xlabel('Mass'), ylabel('Amplitude');
        title(sprintf('m=%fkg, k=%fN/m, z=%fN*s/m, f0_1=%dHz', m(1), K(1), Z(1), f0_1));
        subplot(2,1,2), plot(O); 
        axis([1 Ns -Ea/4 Ea/4]), xlabel('Sample'), ylabel(sprintf('Amplitude of mass %d', Oi));
        pause(plot_update_delay);
    end
    
    % set up arrays for next iteration
    x_prev = x_curr;
    x_curr = x_next;    
    
    % display progress of simulation in percent
    if mod(n, Ns/100) == 0
        dispstat(sprintf('%d%%', 100*n/Ns));
    end
end

% check whether there are any output samples
if abs(max(O)) > 0         
    % add reverb to the output signal -> https://uk.mathworks.com/help/releases/R2016b/coder/examples/reverberation-using-matlab-classes.html - 'Reverberation Using MATLAB Classes' - Matlab documentation
    if reverb_enable == 1
        dispstat('Reverberating... ', 'keepthis');
        
        % set up reverb parameters
        reverb = Reverb(Fs);
        reverb.Density = reverb_density;
        reverb.Feedback = reverb_feedback;
        reverb.Feedthrough = reverb_feedthrough;
        
        % perform reverberation
        for i = (1 : Ns)
            reverb.update(O(i));
            O(i) = reverb.output();
            
            % display progress of reverbaration in percent
            if mod(i, Ns/100) == 0
                dispstat(sprintf('%d%%', 100*i/Ns));
            end
        end
    end
    
    % add delay to the output signal -> https://www.elec.york.ac.uk/internal_web/bsc/yr3/modules/Audio_Algorithms/AADI_SP_lec6-2015-V1-1_ONLINE.pdf - 'AADI: Lecture 6: MATLAB simulation of synthetic reverberation' - Tony Tew
    if delay_enable == 1
        dispstat('Delaying... ', 'keepthis');
        
        % set up comb filter parameters
        comb_Tc = delay_interval;
        comb_M = round(comb_Tc / Ts, 0);
        comb_g = delay_gain;
        comb_zeros = [zeros(1, comb_M) 1]; 
        comb_poles = [1 zeros(1, comb_M-1) -comb_g];
        
        % create new vector with (comb_Tc * Fs) elements more than Onorm and place filtered result inside
        Odelay = [O zeros(1, comb_Tc * Fs)];
        Odelay = filter(comb_zeros, comb_poles, Odelay);
        
        % take samples from index (comb_Tc * Fs) to remove inital delay
        O = Odelay(comb_Tc * Fs : end);
    end
    
    % normalise output signal
    Onorm = O / max(O);
    
    if figure_enable == 1
        % Plot frequency domain response of output
        figure(3);
        fftSize = Ns;%8192;
        fft_ = fft(Onorm, fftSize);
        f = (0:fftSize-1)*(Fs/fftSize);
        semilogx(f, 20 * log10(abs(fft_) / max(abs(fft_))));
        axis([1 20000 -40 0]);
        xlabel('Frequency (Hz)');
        ylabel('Magnitude Response (dB)');
        
        % Plot time domain response of output
        figure(2);
        plot(Onorm);
        axis([0 Ns -1 1]);
        xlabel('Time (samples)');
        ylabel('Amplitude');
    end

    % play audio output
    sound(Onorm, Fs);
    
    % save audio output as a .wav file
    if write_audio == 1
        filename = sprintf('audio/Nm=%d_m=%f_k=%f_z=%f_delay=%d_reverb=%d.wav', Nm, m_(1), k_(1), z_(1), delay_enable, reverb_enable);
        audiowrite(filename, Onorm, Fs);
    end
    
    dispstat('Done', 'keepthis');
else
    dispstat('All output samples are 0', 'keepthis');
end
